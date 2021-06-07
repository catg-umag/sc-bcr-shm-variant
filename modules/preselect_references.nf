nextflow.enable.dsl = 2

include { mappingMinimap2; sortAndConvert } from './common.nf'

def splitSubjectChain(subject_chain) {
  return (subject_chain =~ /^(.*)-([HL]C)$/).findAll()[0][1..-1]
}


workflow PreselectReferences {
  take:
    reads                   // channel [subject, FASTA files]
    refs_by_subject_chain   // channel [subject-chain, FASTA file]
    refs_single             // channel [subject-chain]
    refs_multi              // channel [subject-chain]

  main:
    refs_by_subject_chain
      | map { splitSubjectChain(it[0]) + [it[1]] }
      | combine(reads, by: 0)
      | map { ["${it[0]}-${it[1]}", it[3], it[2]] }
      | multiMap {
          reads: it
          step: 'pre'
        }
      | mappingMinimap2
      | sortAndConvert
      | selectCellReference

    // splitted reads for single reference (dummy)
    refs_single
      | join(refs_by_subject_chain)
      | map { splitSubjectChain(it[0]) + [(it[1].text =~ /(?<=>).*/)[0]] }
      | combine(reads, by:0)
      | map { ["${it[0]}-${it[1]}", it[2], it[3]] }
      | set { splitted_reads_single }

    // splitted reads for multi reference
    selectCellReference.out
      | join(refs_multi)
      | map { splitSubjectChain(it[0]) + [it[1]] }
      | combine(reads, by: 0)
      | map { ["${it[0]}-${it[1]}", it[3], it[2]] }
      | splitFastqByReference
      | flatMap {
          [it[1], it[2]]
            .transpose()
            .collect { x -> [it[0], x[0].simpleName.split("--")[1], x] }
        }
      | set { splitted_reads_multi }

  emit:
    reads = splitted_reads_single.mix(splitted_reads_multi)
    cell_references = selectCellReference.out
}


/*
 * Select best reference 
 */
process selectCellReference {
  tag "$subject_chain"
  label 'julia'
  publishDir "output/cell_references/", mode: 'copy'

  input:
  tuple val(subject_chain), path(bam), path(bam_index)

  output:
  tuple val(subject_chain), path("${subject_chain}.csv")
  
  script:
  """
  select_reference.jl -i $bam -o ${subject_chain}.csv
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