nextflow.enable.dsl = 2

params.config_file = params.config_file ?: "/Users/bibhusimkhada/Desktop/RNA_seq/Pipeline/config.yaml"
def cfg = new groovy.yaml.YamlSlurper().parse(new File(params.config_file))

assert cfg.input_dir              : "Missing input_dir in ${params.config_file}"
assert cfg.output_dir             : "Missing output_dir in ${params.config_file}"
assert cfg.files?.pattern         : "Missing files.pattern in ${params.config_file}"
assert cfg.reference?.star_index  : "Missing reference.star_index in ${params.config_file}"
assert cfg.reference?.rRNA        : "Missing reference.rRNA in ${params.config_file}"
assert cfg.tools?.bbduk           : "Missing tools.bbduk in ${params.config_file}"

params.input_dir    = cfg.input_dir
params.output_dir   = cfg.output_dir
params.threads      = (cfg.threads ?: 4) as int
params.read_pattern = cfg.files.pattern
params.star_index   = cfg.reference.star_index
params.rRNA         = cfg.reference.rRNA
params.bbduk        = cfg.tools.bbduk

process FASTP {
    tag { sample_id }
    cpus params.threads
    publishDir "${params.output_dir}/1_fastp", mode: 'copy', overwrite: true

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

process BBDUK_RRNA {
    tag { sample_id }
    cpus params.threads
    publishDir "${params.output_dir}/2_bbduk", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.clean.fq.gz"),
          path("${sample_id}_R2.clean.fq.gz")

    script:
    """
    ${params.bbduk} \
      in1=${r1} \
      in2=${r2} \
      out1=${sample_id}_R1.clean.fq.gz \
      out2=${sample_id}_R2.clean.fq.gz \
      ref=${params.rRNA} \
      k=31 \
      hdist=1 \
      stats=${sample_id}.bbduk.stats.txt \
      t=${task.cpus}
    """
}

process STAR_ALIGN {
    tag { sample_id }
    cpus params.threads
    publishDir "${params.output_dir}/3_star", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("${sample_id}.Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR \
      --runThreadN ${task.cpus} \
      --genomeDir "${params.star_index}" \
      --readFilesIn "${r1}" "${r2}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${sample_id}." \
      --outSAMtype BAM SortedByCoordinate
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/${params.read_pattern}", size: 2, checkIfExists: true)

    trimmed_ch = FASTP(read_pairs_ch)

    rrna_input_ch = trimmed_ch.map { sample_id, r1, r2, json, html ->
        tuple(sample_id, r1, r2)
    }

    rrna_filtered_ch = BBDUK_RRNA(rrna_input_ch)

    star_input_ch = rrna_filtered_ch.map { sample_id, r1, r2 ->
        tuple(sample_id, r1, r2)
    }

    STAR_ALIGN(star_input_ch)
}