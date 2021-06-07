nextflow.enable.dsl = 2


workflow AnalyzeConsensus {
  take:
    consensus           // channle [subject_chain, CSV file]
    cell_references     // channel [subject_chain, CSV file]
    reference_regions   // channel [subject_chain, CSV file]
    aligned_v_regions   // channel [subject_chain, FASTA file]

  main:
    consensus
      | filterConsensus
      | getConsensusProductivity

    consensus
      | join(cell_references)
      | join(reference_regions)
      | join(aligned_v_regions)
      | join(getConsensusProductivity.out)
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
  tuple val(name), path("${subject_chain}.fasta")
  
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
  tuple \
    val(subject_chain), \
    path('consensus_summary.csv'), \
    path('cell_consensus.csv'), \
    path('reference_regions.csv'), \
    path('alined_v_genes.fasta'), \
    path('consensus_productivity.csv')

  output:
  tuple val(subject_chain), path("${subject_chain}.shm.csv")

  script:
  subject = subject_chain - ~/-[HL]C/
  """
  get_shm_places.jl \
    --input consensus_summary.csv \
    --output ${subject_chain}.shm.csv \
    --cell-references cell_consensus.csv \
    --reference-regions reference_regions.csv \
    --aligned-v-regions alined_v_genes.fasta \
    --consensus-productivity consensus_productivity.csv
  """
}


/*
 * Run IgBlast on each consensus and extract its productive status
 */
process getConsensusProductivity {
  tag "$subject_chain"
  label 'immcantation'
  cpus 8

  input:
  tuple val(subject_chain), path(filtered_consensus)
  
  output:
  tuple val(subject_chain), path("${name}_productivity.csv")
  
  script:
  name = subject_chain
  """
  AssignGenes.py igblast -s $filtered_consensus -b /usr/local/share/igblast \
    --outdir . --outname $name --format airr --nproc ${task.cpus}
  
  awk '
  NR==1 { for (i=1; i<=NF; i++) { f[\$i] = i } }
  { print \$(f["sequence_id"]), \$(f["productive"]) }
  ' ${name}_igblast.tsv | sed 's/ /,/g' > ${name}_productivity.csv
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
