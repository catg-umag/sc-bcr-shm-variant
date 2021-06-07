nextflow.enable.dsl = 2


workflow Preprocessing {
  take:
    reads           // channel [subject, FASTA files]
    cell_index
    cell_whitelist

  main:
    // correct and extract barcodes
    extractBarcodes(reads, cell_whitelist)
    
    // sequence cleanup
    sequenceCleanup(extractBarcodes.out.reads)

    splitFastqBySample(sequenceCleanup.out.reads, cell_index)
      | transpose
      | map { [it[0].simpleName.replaceAll(/_L\d+$/, '')] + it }
      | groupTuple
      | mergeLanes

    fastQC(mergeLanes.out.reads)
      | collect
      | multiQC

  emit:
    mergeLanes.out.reads
}


/*
 * Extracts barcodes from FASTQ applying correction
 */
process extractBarcodes {
  tag "$name"
  publishDir "output/bc_corrections/", pattern: '*.json', mode: 'copy'
  label 'julia'
  cpus 2
  memory 16.GB

  input:
  tuple val(name), path(reads)
  path(cell_whitelist)

  output:
  tuple val(name), path("barcode_corrections_${name}.json"), emit: corrections
  tuple val(name), path("${name}.r1.fastq"), path("${name}.r2.fastq"), emit: reads

  script:
  (r1, r2) = reads.sort { it.name }
  """
  export JULIA_NUM_THREADS=${task.cpus}
  extract_barcodes.jl \
    -o ${name}.r1.fastq -O ${name}.r2.fastq \
    -i $r1 -I $r2 \
    -w $cell_whitelist \
    -c barcode_corrections_${name}.json
  """
}


/*
 * Trims adapters (automatic) and clean remaining FASTQ files
 */
process sequenceCleanup {
  tag "$name"
  label 'fastp'
  publishDir 'output/qc/fastp/', pattern: "${name}_report*", mode: 'copy'
  cpus 8

  input:
  tuple val(name), path(fastq_r1), path(fastq_r2)

  output:
  tuple val(name), path("*_clean.*.fastq.gz"), emit: reads
  tuple path("${name}_report.html"), path("${name}_report.json"), emit: report
  
  script:
  """
  fastp \
    -i ${fastq_r1} -I ${fastq_r2} \
    -o ${name}_clean.R1.fastq.gz -O ${name}_clean.R2.fastq.gz \
    --adapter_sequence=TTTCTTATATGGG \
    --qualified_quality_phred 30 --unqualified_percent_limit 25 \
    --trim_poly_g --trim_poly_x --poly_x_min_len 6 \
    --trim_front1 12 --trim_tail1 2 --trim_front2 2 --trim_tail2 2 \
    --cut_front --cut_tail --cut_mean_quality 30 \
    --thread ${task.cpus} \
    -j ${name}_report.json -h ${name}_report.html
  """
}


/*
 * Splits FASTQ files by subject
 */
process splitFastqBySample {
  tag "$name"
  label 'julia'
  memory 2.GB

  input:
  tuple val(name), path(reads)
  path cell_index

  output:
  tuple path("*.r1.fastq"), path("*.r2.fastq"), path("*_info.csv")

  script:
  (r1, r2) = reads.sort { it.name }
  exp = name.replaceAll(/_L\d+$/, '')
  """
  export JULIA_NUM_THREADS=${task.cpus}
  split_fastq_by_sample.jl --r1 $r1 --r2 $r2 -i $cell_index -n $exp
  """
}


/*
 * Merges files from same experiment/subject but different lane
 */
process mergeLanes {
  tag "$name"
  publishDir "output/fastq_by_sample/$exp/$subject/", mode: 'copy'
  cpus 4

  input:
  tuple \
    val(name), \
    path(r1, stageAs: 'r1_*.fastq'), \
    path(r2, stageAs: 'r2_*.fastq'), \
    path(info, stageAs: 'info_*.csv')

  output:
  tuple val(name), path("${name}.r?.fastq.gz"), emit: reads
  tuple val(name), path("${name}_info.csv"), emit: info

  script:
  (exp, subject) = name.split('_')
  if (r1.size() > 1)
    """
    cat $r1 | pigz -c -p ${task.cpus} > ${name}.r1.fastq.gz
    cat $r2 | pigz -c -p ${task.cpus} > ${name}.r2.fastq.gz
    cat $info > ${name}_info.csv
    """
  else
    """
    pigz -c -p ${task.cpus} $r1 > ${name}.r1.fastq.gz
    pigz -c -p ${task.cpus} $r2 > ${name}.r2.fastq.gz
    cp $info ${name}_info.csv
    """
}


process fastQC {
  tag "$name"
  label 'fastqc'
  publishDir "output/qc/fastqc", mode: 'copy'
  cpus 2

  input:
  tuple val(name), path(fastq_files)

  output:
  path("fastqc_${name}")

  script:
  """
  mkdir fastqc_${name}
  fastqc $fastq_files -o fastqc_${name} - t ${task.cpus}
  """
}


process multiQC {
  label 'multiqc'
  publishDir "output/qc/multiqc", mode: 'copy'
  
  input:
  path fastqc_reports

  output:
  tuple path('multiqc_data'), path('multiqc_report.html')

  script:
  """
  cp ${workflow.projectDir}/conf/multiqc_config.yaml .
  multiqc .
  """
}