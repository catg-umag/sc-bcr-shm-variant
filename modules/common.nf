nextflow.enable.dsl = 2


/*
 * Maps reads against reference
 */
process mappingMinimap2 {
  tag "$name"
  label 'minimap2'
  cpus "$ncpus"
  memory 4.GB

  input:
  tuple val(name), path(reads), path(reference)
  val(step)

  output:
  tuple val(name), file("${name}.sam"), emit: sam
  val(step), emit: step

  script:
  ncpus = step === 'pre' ? params.mapping_cpus_pre : params.mapping_cpus
  """
  minimap2 -ax sr -t ${task.cpus} $reference $reads > ${name}.sam
  """
}


/*
 * Sort SAM files and convert them to BAM
 */
process sortAndConvert {
  tag "$name"
  label 'samtools'
  cpus "$ncpus"

  input:
  tuple val(name), path(sam)
  val(step)
  
  output:
  tuple val(name), path("${name}.bam"), path("${name}.bam.bai")
  
  script:
  ncpus = step === 'pre' ? params.mapping_cpus_pre : params.mapping_cpus
  """
  samtools sort -@ ${task.cpus} $sam -o ${name}.bam
  samtools index ${name}.bam
  """
}