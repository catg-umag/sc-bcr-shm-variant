nextflow.enable.dsl = 2

include { mappingMinimap2; sortAndConvert } from './common.nf'


workflow BuildConsensus {
  take:
    splitted_reads        // channel [subject_chain, refname, FASTQ files]
    splitted_references   // channel [subject_chain, refname, FASTA file]
    refs_by_subject_chain // channel [subject_chain, FASTA file]
    reference_regions     // channel [subject_chain, CSV file]

  main:
    splitted_reads
      | join(splitted_references, by: [0, 1])
      | map { ["${it[0]}--${it[1]}", it[2], it[3] ] }
      | multiMap {
          reads: it
          step: 'final'
        }
      | mappingMinimap2
      | sortAndConvert
      | map { [it[0].split("--")[0]] + it[1..2] }
      | groupTuple
      | mergeBams
      | join(refs_by_subject_chain)
      | join(reference_regions)
      | makeConsensus

  emit:
    makeConsensus.out
}


/*
 * Merge BAM files from same subject/chain
 */
process mergeBams {
  tag "$subject_chain"
  label 'samtools'
  publishDir "output/bam/${subject}", mode: 'copy'
  cpus 4

  input:
  tuple val(subject_chain), path(bam_files), path(bam_indexes)

  output:
  tuple val(subject_chain), path("${subject_chain}.bam"), path("${subject_chain}.bam.bai")

  script:
  subject = subject_chain - ~/-[HL]C/
  """
  samtools merge -@ ${task.cpus} ${subject_chain}.bam $bam_files
  samtools index ${subject_chain}.bam
  """
}


/*
 * Make consensus from BAM
 */
process makeConsensus {
  tag "$name"
  label 'julia'
  publishDir "output/consensus/${subject}", mode: 'copy'
  cpus 8
  memory { 32.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
  maxRetries 3

  input:
  tuple val(subject_chain), path(bam), path(bam_index), path(references), path(regions)

  output:
  tuple val(name), path("${name}.csv")

  script:
  name = subject_chain
  subject = name - ~/-[HL]C/
  """
  export JULIA_NUM_THREADS=${task.cpus}
  bcr_consensus.jl --input-bamfile $bam --references $references \
    --reference-regions $regions --output ${name}_unsrt.csv

  # sort
  export TMPDIR=.
  head -n 1 ${name}_unsrt.csv > ${name}.csv
  tail -n +2 ${name}_unsrt.csv | sort --field-separator=, --key=1,2 >> ${name}.csv
  """
}