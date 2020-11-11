#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { build_cutadapt_args } from './lib/nextflow/utils'


/*
 * Prepare inputs
 */
cell_index = file(params.cell_index)
cell_whitelist = file(params.cell_whitelist)
region_positions = file(params.region_positions)

channel
  .fromPath(params.references)
  .splitFasta(record: [id: true, sequence: true])
  .collectFile(storeDir: 'output/references') {
    ["${it.id.replaceFirst(/_[HL]C/, '') }.fasta", ">${it.id}\n${it.sequence}"]
  }
  .map { [it.baseName, it] }
  .set { references }

channel
  .fromList(params.data)
  .map { ["${it.name}_L${it.lane}", file("$HOME" + it.fastq.replace('$HOME', ''))]}
  .set { reads }

channel
  .fromPath(params.external_consensus_path)
  .map { [it.baseName, it] }
  .set { external_consensus_with_name }


/*
 * Workflow
 */
workflow {
  getBarcodeCorrections(reads, cell_whitelist)

  extractBarcodes(reads.join(getBarcodeCorrections.out))
  
  // trimAdaptersCutAdapt(extractBarcodes.out.reads)
  // splitFastqBySample(trimAdapters.out.reads.join(extractBarcodes.out.barcodes), cell_index)

  // no barcode trimming (temporary)
  extractBarcodes.out.reads
    .map { [it[0], [it[1], it[2]]] }
    .join(extractBarcodes.out.barcodes)
    .set { reads_with_barcodes }

  splitFastqBySample(reads_with_barcodes, cell_index)

  splitFastqBySample.out
    .transpose()
    .map { [it[0].simpleName.replaceAll(/_L\d+$/, '')] + it }
    .groupTuple()
    .set { grouped_lanes }

  mergeLanes(grouped_lanes)

  mergeLanes.out.reads
    .filter { it[0] =~ /K?B_S\d+/ }
    .join(references)
    .set { reads_with_ref }

  sortAndConvert(mapping(reads_with_ref))
  
  sortAndConvert.out.join(references)
    .join(mergeLanes.out.info)
    .set { consensus_pre_info }

  makeConsensus(consensus_pre_info)

  filterConsensus(makeConsensus.out.flatten(), region_positions)

  filterConsensus.out
    .map { [it.baseName.replaceAll(/_.*/, ""), it] }
    .combine(external_consensus_with_name, by: 0)
    .map { it[1..2] }
    .set { cons_with_ext }

  searchSequences(cons_with_ext)
}


/*
 * Corrects barcodes and saves them in a JSON file.
 */
process getBarcodeCorrections {
  tag "$name"
  publishDir "output/bc_corrections/", pattern: '*.json', mode: 'copy'
  label 'julia'
  cpus 1
  memory 16.GB

  input:  
  tuple val(name), path(reads)
  path cell_whitelist

  output:
  tuple val(name), path("corrections_${name}.json")

  script:
  (r1, r2) = reads.sort { it.simpleName }
  """
  export JULIA_NUM_THREADS=${task.cpus}
  get_barcode_corrections.jl $r1 $cell_whitelist -o corrections_${name}.json
  """
}


/*
 * Extracts barcodes from FASTQ applying correction
 */
process extractBarcodes {
  tag "$name"
  label 'julia'
  cpus 2
  memory 16.GB

  input:
  tuple val(name), path(reads), path(corrections)

  output:
  tuple val(name), path("${name}.r1.fastq.gz"), path(r2), emit: reads
  tuple val(name), path("${name}_barcodes.csv"), emit: barcodes

  script:
  (r1, r2) = reads.sort { it.simpleName }
  """
  export JULIA_NUM_THREADS=${task.cpus}
  extract_barcodes.jl -r ${name}.r1.fastq.gz -b ${name}_barcodes.csv -c $corrections $r1
  """
}


/*
 * Clean the sequence of primer sequences
 */
process trimAdaptersCutAdapt {
  tag "$name"
  publishDir 'output/cutadapt/', pattern: '*.txt', mode: 'copy'
  cpus 8

  input:
  tuple val(name), path(fastq_r1), path(fastq_r2)

  output:
  tuple val(name), path("*_clean.*.fastq.gz"), emit: reads
  path "${name}_report.txt", emit: report
  
  script:
  cutadapt_primers = build_cutadapt_args(params.primers)
  """
  cutadapt \
    -e 0.12 \
    --times 3 \
    --overlap 5 \
    -m 20 \
    -j ${task.cpus / 2} \
    -f fastq \
    -o ${name}_clean.R1.fastq.gz \
    -p ${name}_clean.R2.fastq.gz \
    ${cutadapt_primers} \
    ${fastq_r1} ${fastq_r2} > ${name}_report.txt
  """
}


/*
 * Trims adapters (automatic) and clean remaining FASTQ files
 */
process trimAdaptersFastp {
  tag "$name"
  publishDir 'output/fastp/', pattern: "${name}_report*", mode: 'copy'
  cpus 4

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
    -c -g \
    -w ${task.cpus} \
    -j ${name}_report.json -h ${name}_report.html
  """
}


/*
 * Splits FASTQ files by subject
 */
process splitFastqBySample {
  tag "$name"
  label 'julia'
  cpus 4
  memory 16.GB

  input:
  tuple val(name), path(reads), path(barcodes)
  path cell_index

  output:
  tuple path("*.r1.fastq.gz*"), path("*.r2.fastq.gz"), path("*_info.csv")

  script:
  (r1, r2) = reads.sort { it.simpleName }
  exp = name.replaceAll(/_L\d+$/, '')
  """
  export JULIA_NUM_THREADS=${task.cpus}
  split_fastq_by_sample.jl --r1 $r1 --r2 $r2 -b $barcodes -i $cell_index -n $exp
  pigz -p ${task.cpus} *.fastq
  """
}


/*
 * Merges files from same experiment/subject but different lane
 */
process mergeLanes {
  tag "$name"
  publishDir "output/fastq_by_sample/$exp/$subject/", mode: 'copy'

  input:
  tuple \
    val(name), \
    path(r1, stageAs: 'r1_*.fastq.gz'), \
    path(r2, stageAs: 'r2_*.fastq.gz'), \
    path(info, stageAs: 'info_*.csv')

  output:
  tuple val(name), path("${name}.r?.fastq.gz"), emit: reads
  tuple val(name), path("${name}_info.csv"), emit: info

  script:
  (exp, subject) = name.split('_')
  if (r1.size() > 1)
    """
    cat $r1 > ${name}.r1.fastq.gz
    cat $r2 > ${name}.r2.fastq.gz
    cat $info > ${name}_info.csv
    """
  else
    """
    cp $r1 ${name}.r1.fastq.gz
    cp $r2 ${name}.r2.fastq.gz
    cp $info ${name}_info.csv
    """
}


/*
 * Maps reads against reference
 */
process mapping {
  tag "$subject"
  cpus 8

  input:
  tuple val(subject), path(reads), path(reference)

  output:
  tuple val(subject), file("${subject}.sam")

  script:
  """
  minimap2 -a -xsr -k8 -w5 -A2 -B2 -O5,24 -E4,2 -t ${task.cpus} $reference $reads > ${subject}.sam
  """
}


/*
 * Sort SAM files and convert them to BAM
 */
process sortAndConvert {
  tag "$subject"
  publishDir "output/bam/${subject}", mode: 'copy'
  cpus 4

  input:
  tuple val(subject), path(sam)
  
  output:
  tuple val(subject), path("${subject}.bam"), path("${subject}.bam.bai")
  
  script:
  """
  samtools sort -@ ${task.cpus} $sam -o ${subject}.bam
  samtools index ${subject}.bam
  """
}


/*
 * Make consensus from BAM
 */
process makeConsensus {
  tag "$subject"
  label 'julia'
  publishDir "output/consensus/${subject}", mode: 'copy'
  cpus 8
  memory { 16.GB + (16.GB * task.attempt) }
  errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
  tuple val(subject), path(bam), path(bam_index), path(reference), path(barcodes)

  output:
  path('*.csv')

  script:
  """
  export JULIA_NUM_THREADS=${task.cpus}
  bcr_consensus.jl -o ${subject}_HC.csv $bam $barcodes $reference ${subject}_HC
  bcr_consensus.jl -o ${subject}_LC.csv $bam $barcodes $reference ${subject}_LC
  """
}


/*
 * Filter created consensus by a set of given rules
 */
process filterConsensus {
  tag "$subject_chain"
  label 'julia'
  publishDir "output/filtered_consensus/${subject}", mode: 'copy'
  cpus 1

  input:
  path consensus_file
  path region_positions
  
  output:
  path "${subject_chain}.fasta"
  
  script:
  subject_chain = consensus_file.baseName
  subject = subject_chain - ~/_[HL]C/
  """
  export JULIA_NUM_THREADS=${task.cpus}
  filter_consensus.jl \
  --min-vdj-coverage ${params.consensus_rules.vdj_coverage} \
  --min-cdr-coverage ${params.consensus_rules.cdr_coverage} \
  --min-reads ${params.consensus_rules.reads} \
  --reference-name ${consensus_file.baseName} \
  -o ${subject_chain}.fasta $consensus_file $region_positions
  """
}


/*
 * Compare sequences against external ones
 */
process searchSequences {
  tag "$subject_chain"
  label 'julia'
  publishDir "output/alignment/${subject}", mode: 'copy'
  cpus 16
  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }

  input:
  tuple path(own_consensus), path(external_consensus)
  
  output:
  path "${subject_chain}.csv"
  
  script:
  subject_chain = own_consensus.baseName
  subject = subject_chain - ~/_[HL]C/
  """
  export JULIA_NUM_THREADS=${task.cpus}
  search_cells.jl -o ${subject_chain}.csv -d ${params.consensus_max_dist} \
  $own_consensus $external_consensus 
  """
}
