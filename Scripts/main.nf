nextflow.enable.dsl = 2

params.config_file = params.config_file ?: "/Volumes/Bibhus_SSD/RNA_seq/Pipeline/config.yaml"
def cfg = new groovy.yaml.YamlSlurper().parse(new File(params.config_file))

assert cfg.input_dir           : "Missing input_dir in ${params.config_file}"
assert cfg.output_dir          : "Missing output_dir in ${params.config_file}"
assert cfg.files?.pattern      : "Missing files.pattern in ${params.config_file}"
assert cfg.reference?.star_index : "Missing reference.star_index in ${params.config_file}"

params.input_dir    = cfg.input_dir
params.output_dir   = cfg.output_dir
params.threads      = (cfg.threads ?: 4) as int
params.read_pattern = cfg.files.pattern
params.star_index   = cfg.reference.star_index

process FASTP {
    tag { sample_id }
    cpus params.threads
    publishDir "${params.output_dir}/fastp", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trim.fq.gz"),
          path("${sample_id}_R2.trim.fq.gz"),
          path("${sample_id}.fastp.json"),
          path("${sample_id}.fastp.html")

    script:
    """
    fastp \
      -i ${reads[0]} \
      -I ${reads[1]} \
      -o ${sample_id}_R1.trim.fq.gz \
      -O ${sample_id}_R2.trim.fq.gz \
      -j ${sample_id}.fastp.json \
      -h ${sample_id}.fastp.html \
      --detect_adapter_for_pe \
      -w ${task.cpus}
    """
}

process STAR_ALIGN {
    tag { sample_id }
    cpus params.threads
    publishDir "${params.output_dir}/star", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}.Aligned.sortedByCoord.out.bam")

    script:
    """
    TMPDIR=\$(mktemp -d "${PWD}/${sample_id}.STAR.XXXXXX")
    trap 'rm -rf "\$TMPDIR"' EXIT

    STAR \
    --runThreadN ${task.cpus} \
    --genomeDir "${params.star_index}" \
    --readFilesIn "${r1}" "${r2}" \
    --readFilesCommand zcat \
    --outTmpDir "\$TMPDIR" \
    --outFileNamePrefix "${sample_id}." \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM MD
    """
}

workflow {
    read_pairs_ch = Channel
        .fromFilePairs("${params.input_dir}/${params.read_pattern}", size: 2, checkIfExists: true)

    trimmed_ch = FASTP(read_pairs_ch)

    star_input_ch = trimmed_ch.map { it ->
        tuple(it[0], it[1], it[2])
    }

    STAR_ALIGN(star_input_ch)
}