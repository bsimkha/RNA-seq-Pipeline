nextflow.enable.dsl = 2

params.config_file = params.config_file ?: "/Volumes/Bibhu's SSD/RNA_seq/Pipeline/config.yaml"
def cfg = new groovy.yaml.YamlSlurper().parse(new File(params.config_file))

assert cfg.input_dir      : "Missing input_dir in ${params.config_file}"
assert cfg.output_dir     : "Missing output_dir in ${params.config_file}"
assert cfg.files?.pattern : "Missing files.pattern in ${params.config_file}"

params.input_dir    = cfg.input_dir
params.output_dir   = cfg.output_dir
params.threads      = (cfg.run?.threads ?: 4) as int
params.read_pattern = cfg.files.pattern

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

workflow {
    read_pairs_ch = Channel
        .fromFilePairs("${params.input_dir}/${params.read_pattern}", size: 2, checkIfExists: true)

    FASTP(read_pairs_ch)
}