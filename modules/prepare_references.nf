nextflow.enable.dsl = 2


workflow PrepareReferences {
  take:
    references  // channel [subject, chain, FASTA file]

  main:
    // collect references by chain and subject/chain
    references
      | collectFile() {
          ["${it[0]}-${it[1]}.fasta", it[2].text]
        }
      | map { [it.baseName, it] }
      | set { references_by_subject_chain }

    // write references in output directory
    references_by_subject_chain
      | map { it[1] }
      | collectFile(storeDir: 'output/references/fasta')

    references_by_subject_chain
      | getReferenceRegions

    references_by_subject_chain
      | separateReferences
      | branch {
          multi: it[1] instanceof List
            return it[1].collect { x -> [it[0], x.baseName, x] }
          single: true
            return [[it[0], it[1].baseName, it[1]]]
        }
      | mix
      | flatMap
      | set { splitted_references }

    references_by_subject_chain
      | map { [it[0], it[1].text.findAll(">").size()] }
      | branch {
          single: it[1] == 1
          multi:  it[1] > 1
        }
      | set { reference_counts }

  emit:
    by_subject_chain = references_by_subject_chain
    splitted = splitted_references
    regions = getReferenceRegions.out
    single = reference_counts.single.map { it[0] }
    multi = reference_counts.multi.map { it[0] }
}


/*
 * Do MSA on references
 */
process referenceAlignment {
  tag "$name"
  label 'clustalo'
  publishDir "output/references/msa", mode: 'copy'
  
  input:
  tuple val(name), path(references)

  output:
  tuple val(name), path("${name}_msa.fasta")

  script:
  """
  clustalo -i $references --wrap 10000 > ${name}_msa.fasta 
  """
}


/*
 * Gets reference VDJ and CDR positions
 */
process getReferenceRegions {
  tag "$name"
  label 'immcantation'
  publishDir "output/references/regions", mode: 'copy'
  cpus 4

  input:
  tuple val(name), path(references)

  output:
  tuple val(name), path("${name}_regions.csv")
  
  script:
  """
  AssignGenes.py igblast -s $references -b /usr/local/share/igblast \
    --outdir . --outname $name --format airr --nproc ${task.cpus}
  
  get_regions.py -i ${name}_igblast.tsv -o ${name}_regions.csv
  """
}


process separateReferences {
  tag "$name"
  label 'fasplit'

  input:
  tuple val(name), path(references)
  
  output:
  tuple val(name), path("${name}/*.fa")

  script:
  """
  mkdir ${name}/
  faSplit byname $references ${name}/ 
  """
}