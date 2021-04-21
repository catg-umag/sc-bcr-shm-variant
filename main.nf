#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { build_cutadapt_args } from './lib/nextflow/utils'


/*
 * Prepare inputs
 */
cell_index = file(params.cell_index)
cell_whitelist = file(params.cell_whitelist)

channel
  .fromPath(params.references)
  .map { [(it.baseName =~ /(.+)_([HKL])/)[0][1 .. -1], it].flatten() }
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
  // references_by_subject_chain
  //   | (referenceAlignment & getReferenceRegions)

  Preprocessing(reads, cell_index, cell_whitelist)
  PrepareReferences(Preprocessing.out, references)
  BuildConsensus(
    PrepareReferences.out.reads,
    PrepareReferences.out.references,
    PrepareReferences.out.cell_references,
    references
  )

  // filterConsensus.out
  //   .map { [it.baseName.replaceAll(/_.*/, ""), it] }
  //   .combine(external_consensus_with_name, by: 0)
  //   .map { it[1..2] }
  //   .set { cons_with_ext }

  // searchSequences(cons_with_ext)
}


workflow Preprocessing {
  take:
    reads
    cell_index
    cell_whitelist

  main:
    // correct and extract barcodes
    getBarcodeCorrections(reads, cell_whitelist)
    extractBarcodes(reads.join(getBarcodeCorrections.out))
    
    // trim adapters
    if (params.trim_adapters) {
      trimAdaptersFastp(extractBarcodes.out)
      splitFastqBySample(
        trimAdaptersFastp.out.reads,
        cell_index
      )  
    } else {
      splitFastqBySample(
        extractBarcodes.out.map { [it[0], [it[1], it[2]]] },
        cell_index
      )
    }
    
    // split FASTQ files by sample and merge lanes
    splitFastqBySample.out
      .transpose()
      .map { [it[0].simpleName.replaceAll(/_L\d+$/, '')] + it }
      .groupTuple()
      .set { grouped_lanes }

    mergeLanes(grouped_lanes)    

  emit:
    mergeLanes.out.reads
}


workflow PrepareReferences {
  take:
    reads
    references

  main:
    // save reference in output
    references
      .collectFile(storeDir: 'output/references/fasta') {
        [it[2].name, it[2].text]
      }

    // prepare references
    references
      .collectFile(newLine: true) {
        ["${it[0]}.fasta", it[2].text]}
      .map { [it.baseName, it] }
      .set { references_by_subject }

    references
      .map { ["${it[0]}_${it[1]}", it[2]] }
      .set { references_by_subject_chain }

    // align the reads against all references
    reads
      .filter { it[0] =~ /K?B_S\d+/ }
      .join(references_by_subject)
      .set { reads_with_ref }

    sortAndConvert(mapping(reads_with_ref))

    // select best reference for each cell based on previous alignment
    sortAndConvert.out
      .join(
        references
          .map { [it[0], it[2].simpleName] }
          .groupTuple())
      .set { bam_with_reference_names }

    selectReference(bam_with_reference_names)

    reads
      .combine(
        selectReference.out
          .flatMap { it[1].collect { x -> [it[0], x] } },
        by: 0)
      .map { [it[2].simpleName, it[1], it[2]] }
      .set { reads_with_ref_index }

    // split FASTQ files by their best reference
    splitFastqByReference(reads_with_ref_index)

    splitFastqByReference.out
      .flatMap {
        [it[1], it[2]].transpose()
          .collect { x -> [it[0], x[0].simpleName.split("--")[1], x] } }
      .set { splitted_reads }

    separateReferences(references_by_subject_chain)
      .flatMap { it[1].collect { x -> [it[0], x.baseName, x] } }
      .set { splitted_references }

  emit:
    reads = splitted_reads
    references = splitted_references 
    cell_references = selectReference.out.flatMap { it[1].collect { x -> [x.baseName, x] } }
}


workflow BuildConsensus {
  take:
    splitted_reads
    splitted_references
    cell_references
    references

  main:
    references
      .map { ["${it[0]}_${it[1]}", it[2]] }
      .set { references_by_subject_chain }

    getReferenceRegions(references_by_subject_chain)

    splitted_reads.join(splitted_references, by: [0, 1])
      .map { ["${it[0]}--${it[1]}", it[2], it[3] ] }
      .set { splitted_reads_with_ref }


    sortAndConvert(mapping(splitted_reads_with_ref))

    sortAndConvert.out
      .map { [it[0].split("--")[0]] + it[1..2] }
      .groupTuple()
      .set { grouped_bams }

    mergeBams(grouped_bams)

    mergeBams.out
      .join(references_by_subject_chain)
      .join(getReferenceRegions.out)
      .set { consensus_pre_info }

    makeConsensus(consensus_pre_info)
      | filterConsensus

    getShmPlaces(
      makeConsensus.out
        .join(cell_references)
        .join(getReferenceRegions.out)
    )
}



/*
 * Does MSA on references
 */
process referenceAlignment {
  tag "$name"
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


/*
 * Corrects barcodes and saves them in a JSON file.
 */
process getBarcodeCorrections {
  tag "$name"
  publishDir "output/bc_corrections/", pattern: '*.json', mode: 'copy'
  label 'julia'
  cpus 1
  memory 10.GB

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

  input:
  tuple val(name), path(reads), path(corrections)

  output:
  tuple val(name), path("${name}.r1.fastq.gz"), path("${name}.r2.fastq.gz")

  script:
  (r1, r2) = reads.sort { it.simpleName }
  """
  export JULIA_NUM_THREADS=${task.cpus}
  extract_barcodes.jl \
    -o ${name}.r1.fastq.gz -O ${name}.r2.fastq.gz \
    -i $r1 -I $r2 \
    -c $corrections
  """
}


/*
 * Clean the sequence of adapter sequences
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
    -e 2 \
    --times 3 \
    --overlap 5 \
    -m 20 \
    -j ${task.cpus} \
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
    --disable_quality_filtering \
    --trim_poly_g \
    --trim_front1 15 --trim_tail1 2 --trim_front2 2 --trim_tail2 2 \
    --cut_front --cut_tail --cut_mean_quality 30 \
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
  memory 2.GB

  input:
  tuple val(name), path(reads)
  path cell_index

  output:
  tuple path("*.r1.fastq.gz*"), path("*.r2.fastq.gz"), path("*_info.csv")

  script:
  (r1, r2) = reads.sort { it.simpleName }
  exp = name.replaceAll(/_L\d+$/, '')
  """
  export JULIA_NUM_THREADS=${task.cpus}
  split_fastq_by_sample.jl --r1 $r1 --r2 $r2 -i $cell_index -n $exp
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
  memory 4.GB

  input:
  tuple val(subject), path(reads), path(reference)

  output:
  tuple val(subject), file("${subject}.sam")

  script:
  """
  minimap2 -ax sr -t ${task.cpus} $reference $reads > ${subject}.sam
  """
}


/*
 * Sort SAM files and convert them to BAM
 */
process sortAndConvert {
  tag "$subject"
  label 'samtools'
  cpus 8

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
 * Select best reference 
 */
process selectReference {
  tag "$subject"
  label 'julia'
  publishDir "output/cell_references/", mode: 'copy'

  input:
  tuple val(subject), path(bam), path(bam_index), val(reference_names)

  output:
  tuple val(subject), path("*.csv")
  
  script:
  // H should be before either K or L
  reference_names = reference_names.sort()
  """
  select_reference.jl -i $bam -H ${reference_names[0]}.csv -L ${reference_names[1]}.csv
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
  split_fastq_by_reference.jl $r1 $r2 $reference_index -o ${subject_chain}_split
  """
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
  subject = subject_chain - ~/_[HKL]/
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
  tuple val(name), path(bam), path(bam_index), path(references), path(regions)

  output:
  tuple val(name), path("${name}.csv")

  script:
  subject = name - ~/_[HKL]$/
  """
  export JULIA_NUM_THREADS=${task.cpus}
  bcr_consensus.jl -o ${name}_unsrt.csv $bam $references $regions

  # sort
  export TMPDIR=.
  head -n 1 ${name}_unsrt.csv > ${name}.csv
  tail -n +2 ${name}_unsrt.csv | sort --field-separator=, --key=1,2 >> ${name}.csv
  """
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
  subject = subject_chain - ~/_[HKL]/
  """
  export JULIA_NUM_THREADS=${task.cpus}
  filter_consensus.jl \
    --min-vdj-coverage ${params.consensus_rules.vdj_coverage} \
    --min-cdr-coverage ${params.consensus_rules.cdr_coverage} \
    --min-reads ${params.consensus_rules.reads} \
    -o ${subject_chain}.fasta $consensus_file
  """
}


/*
 * Get possible SHM places
 */
process getShmPlaces {
  tag "$name"
  publishDir "output/shm/${subject}", mode: 'copy'
  label 'julia'

  input:
  tuple val(name), path('consensus_summary.csv'), \
    path('cell_consensus.csv'), path('reference_regions.csv')

  output:
  tuple val(name), path("${name}_shm.csv")

  script:
  subject = name - ~/_[HKL]/
  """
  get_shm_places.jl -i consensus_summary.csv -o ${name}_shm.csv \
    -c cell_consensus.csv -g reference_regions.csv
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
