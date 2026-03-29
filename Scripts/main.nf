nextflow.enable.dsl = 2

// Load config.yaml
params.config = "/Volumes/Bibhu's SSD/RNA_seq/Pipeline/config.yaml"
config = new groovy.yaml.YamlSlurper().parse(new File(params.config))

params.input_dir  = config.input_dir
params.output_dir = config.output_dir
params.threads    = config.run?.threads ?: 4

params.read_pattern = "${params.input_dir}/*_{1,2}.fq.gz"