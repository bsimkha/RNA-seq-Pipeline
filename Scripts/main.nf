nextflow.enable.dsl = 2

params.config_file = params.config_file ?: "path/to/config.yaml"
def cfg = new groovy.yaml.YamlSlurper().parse(new File(params.config_file))

assert cfg.input_dir              : "Missing input_dir in ${params.config_file}"
assert cfg.output_dir             : "Missing output_dir in ${params.config_file}"
assert cfg.files?.pattern         : "Missing files.pattern in ${params.config_file}"
assert cfg.reference?.star_index  : "Missing reference.star_index in ${params.config_file}"
assert cfg.reference?.rRNA        : "Missing reference.rRNA in ${params.config_file}"
assert cfg.tools?.bbduk           : "Missing tools.bbduk in ${params.config_file}"
assert cfg.tools?.featurecounts   : "Missing tools.featurecounts in ${params.config_file}"
assert cfg.reference?.annotation  : "Missing reference.annotation in ${params.config_file}"

params.input_dir    = cfg.input_dir
params.output_dir   = cfg.output_dir
params.threads      = (cfg.threads ?: 4) as int
params.read_pattern = cfg.files.pattern
params.star_index   = cfg.reference.star_index
params.rRNA         = cfg.reference.rRNA
params.bbduk        = cfg.tools.bbduk
params.featurecounts = cfg.tools.featurecounts
params.gtf          = cfg.reference.annotation

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

process FASTQC {
    tag { sample_id }
    cpus 2
    publishDir "${params.output_dir}/3_fastqc", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*.html"), path("*.zip")

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    """
}

process STAR_ALIGN {
    tag { sample_id }
    cpus params.threads
    publishDir "${params.output_dir}/4_star", mode: 'copy', overwrite: true

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

process FEATURECOUNTS {
    tag { sample_id }
    cpus params.threads
    publishDir "${params.output_dir}/5_featurecounts", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id),
          path("${sample_id}.counts.txt"),
          path("${sample_id}.counts.txt.summary")

    script:
    """
    ${params.featurecounts} \
      -T ${task.cpus} \
      -a "${params.gtf}" \
      -o ${sample_id}.counts.txt \
      ${bam}
    """
}

process MERGE_COUNTS {
    tag "merge_counts"
    cpus params.threads
    publishDir "${params.output_dir}/6_counts_matrix", mode: 'copy', overwrite: true

    input:
    path(count_files)

    output:
    path("combined_counts.txt")

    script:
    """
    cat > merge_counts.R << 'RS'
    files <- list.files()
    files <- files[endsWith(files, ".counts.txt")]

    dfs <- lapply(files, function(f) {
      read.table(f, sep = "\\t", header = TRUE, comment.char = "#", check.names = FALSE)
    })

    counts <- do.call(cbind, lapply(dfs, function(df) df[, 7, drop = FALSE]))
    ref <- dfs[[1]]

    result <- data.frame(
      Geneid = ref[, 1],
      length = ref[, 6],
      counts,
      check.names = FALSE
    )

    colnames(result) <- c(
      "Geneid",
      "length",
      sub("_counts.txt", "", basename(files), fixed = TRUE)
    )

    write.table(result, "combined_counts.txt", sep = "\\t", quote = FALSE, row.names = FALSE)
    RS

    Rscript merge_counts.R
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/${params.read_pattern}", size: 2, checkIfExists: true)

    trimmed_ch = FASTP(read_pairs_ch)

    rrna_input_ch = trimmed_ch.map { sample_id, r1, r2, json, html ->
        tuple(sample_id, r1, r2)
    }

    rrna_filtered_ch = BBDUK_RRNA(rrna_input_ch)

    FASTQC(rrna_filtered_ch)

    star_input_ch = rrna_filtered_ch.map { sample_id, r1, r2 ->
        tuple(sample_id, r1, r2)
    }

    star_out = STAR_ALIGN(star_input_ch)

    per_sample_counts = FEATURECOUNTS(star_out)

    count_files = per_sample_counts.map { sample_id, counts, summary -> counts }.collect()

    MERGE_COUNTS(count_files)
}
