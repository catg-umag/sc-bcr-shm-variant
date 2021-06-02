nextflow.enable.dsl = 2


workflow AnalyzeConsensus {
  take:
    consensus           // channle [subject_chain, CSV file]
    cell_references     // channel [subject_chain, CSV file]
    reference_regions   // channel [subject_chain, CSV file]

  main:
    consensus
      | filterConsensus

    consensus
      | join(cell_references)
      | join(reference_regions)
      | getShmPlaces
}


/*
 * Filter created consensus by a set of given rules
 */
process filterConsensus {
  tag "$subject_chain"
  publishDir "output/filtered_consensus/${subject}", mode: 'copy'
  label 'julia'
  cpus 2

  input:
  tuple val(name), path(consensus_file)
  
  output:
  path "${subject_chain}.fasta"
  
  script:
  subject_chain = consensus_file.baseName
  subject = subject_chain - ~/-[HL]C/
  """
  export JULIA_NUM_THREADS=${task.cpus}
  filter_consensus.jl \
    --input-consensus $consensus_file --output ${subject_chain}.fasta \
    --min-vdj-coverage ${params.consensus_rules.vdj_coverage} \
    --min-cdr-coverage ${params.consensus_rules.cdr_coverage} \
    --min-reads ${params.consensus_rules.reads}
  """
}


/*
 * Get possible SHM places
 */
process getShmPlaces {
  tag "$subject_chain"
  publishDir "output/shm/${subject}", mode: 'copy'
  label 'julia'

  input:
  tuple val(subject_chain), path('consensus_summary.csv'), \
    path('cell_consensus.csv'), path('reference_regions.csv')

  output:
  tuple val(subject_chain), path("${subject_chain}.shm.csv")

  script:
  subject = subject_chain - ~/-[HL]C/
  """
  get_shm_places.jl --input consensus_summary.csv --output ${subject_chain}.shm.csv \
    --cell-references cell_consensus.csv --reference-regions reference_regions.csv
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

  input:
  tuple path(own_consensus), path(external_consensus)
  
  output:
  path "${subject_chain}.csv"
  
  script:
  subject_chain = own_consensus.baseName - ~/_sorted/
  subject = subject_chain - ~/_[HKL]/
  """
  export JULIA_NUM_THREADS=${task.cpus}
  search_cells.jl -o ${subject_chain}.csv -d ${params.consensus_max_dist} \
    $own_consensus $external_consensus 
  """
}
