nextflow.enable.dsl = 2

include { mappingMinimap2; sortAndConvert } from './common.nf'


workflow PreselectReferences {
  take:
    reads                   // channel [subject, FASTA files]
    refs_by_subject         // channel [subject, FASTA file]
    refs_by_subject_chain   // channel [subject-chain, FASTA file]

  main:
    reads
      | join(refs_by_subject)
      | multiMap {
          reads: it
          ncpus: 8
        }
      | mappingMinimap2
      | sortAndConvert
      | selectReference

    reads
      | combine(selectReference.out.mix(), by: 0)
      | map { ["${it[0]}-${it[2]}", it[1], it[3]] }
      | splitFastqByReference
      | flatMap {
          [it[1], it[2]]
            .transpose()
            .collect { x -> [it[0], x[0].simpleName.split("--")[1], x] }
        }
      | set { splitted_reads }

    selectReference.out
      | mix
      | map { ["${it[0]}-${it[1]}", it[2]] }
      | set { cell_references }

  emit:
    reads = splitted_reads
    cell_references = cell_references
}


/*
 * Select best reference 
 */
process selectReference {
  tag "$subject"
  label 'julia'
  publishDir "output/cell_references/", mode: 'copy'

  input:
  tuple val(subject), path(bam), path(bam_index)

  output:
  tuple val(subject), val('HC'), path("${subject}-HC.csv"), emit: hc
  tuple val(subject), val('LC'), path("${subject}-LC.csv"), emit: lc
  
  script:
  """
  select_reference.jl -i $bam -H ${subject}-HC.csv -L ${subject}-LC.csv
  """
}


/*
 * Splits FASTQ file depending the reference assigned to each cell
 */
process splitFastqByReference {
  tag "$subject_chain"
  label 'julia'

  input:
  tuple val(subject_chain), path(reads), path(reference_index)

  output:
  tuple \
    val(subject_chain), \
    path("${subject_chain}_split/*.r1.fastq"), \
    path("${subject_chain}_split/*.r2.fastq")

  script:
  (r1, r2) = reads.sort { it.name }
  """
  split_fastq_by_reference.jl --r1 $r1 --r2 $r2 \
    --reference-index $reference_index --output-dir ${subject_chain}_split
  """
}