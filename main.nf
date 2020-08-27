cell_index = file(params.cell_index)
cell_whitelist = file(params.cell_whitelist)

channel
  .fromPath(params.references)
  .splitFasta( record: [id: true, sequence: true])
  .collectFile( storeDir: 'output/references ') {
    ["${it.id.split("_")[0]}.fasta", ">${it.id}\n${it.sequence}"]
  }
  .map { [it.baseName, it] }
  .set { references_ch }

channel
  .fromList(params.data)
  .map { [it.name, file("$HOME" + it.fastq.replace('$HOME', ''))]}
  .into { reads_ch; reads_ch2 }


// /*
//  * Corrects barcodes and saves them in a JSON file.
//  */
// process getBarcodeCorrections {
//   tag "$name"
//   publishDir "output/bc_corrections/", pattern: '*.json', mode: 'copy'

//   input:  
//   tuple val(name), path(reads) from reads_ch
//   path cell_whitelist

//   output:
//   tuple val(name), path(reads), path("corrections_${name}.json") into reads_w_corrections_ch

//   script:
//   (r1, r2) = reads.sort { it.simpleName }
//   """
//   get_barcode_corrections.jl $r1 $cell_whitelist -o corrections_${name}.json
//   """
// }


// /*
//  * Splits FASTQ files by subject
//  */
// process splitFastqBySample {
//   tag "$name"
//   publishDir "output/fastq_by_sample/$name/", mode: 'copy'
//   cpus 4

//   input:
//   tuple val(name), path(reads), path(corrections) from reads_w_corrections_ch
//   path cell_index

//   output:
//   path("*") into sample_fastqs_ch

//   script:
//   (r1, r2) = reads.sort { it.simpleName }
//   """
//   split_fastq_by_sample.jl --r1 $r1 --r2 $r2 -i $cell_index -n $name -c $corrections
//   pigz -p ${task.cpus} *.fastq
//   """
// }


// sample_fastqs_ch
//   .flatMap { it.collate(3) }
//   .map { [it[0].simpleName.replaceAll(/.*_/, ''), [it[0], it[1]], it[2]] }
//   .filter { it[0] ==~ /S\d+/ }
//   .set { subject_data_ch }


// subject_data_ch
//   .map { [it[0], it[1]] }
//   .join(references_ch)
//   .set { reads_withref_ch }


// process mapping {
//   tag "$subject"

//   input:
//   tuple val(subject), path(reads), path(reference) from reads_withref_ch

//   output:
//   tuple val(subject), file("${subject}.sam") into mapped_sequences

//   script:
//   """
//   minimap2 -a -xsr -k8 -w5 -A2 -B2 -O5,24 -E4,2 -t2 $reference $reads > ${subject}.sam
//   """
// }


// process sortAndConvert {
//   tag "$subject"
//   publishDir "output/bam/${subject}", mode: 'copy'

//   input:
//   tuple val(subject), file(sam) from mapped_sequences
  
//   output:
//   file "${subject}.bam"
//   file "${subject}.bam.bai"
  
//   script:
//   """
//   samtools sort $sam -o ${subject}.bam
//   samtools index ${subject}.bam
//   """
// }


process cutAdapt {
  tag "$name"
  publishDir "output/cutadapt/${name}", mode: 'copy'
  cpus 8

  input:
  tuple val(name), path(reads) from reads_ch2

  output:
  path("*")
  
  script:
  (r1, r2) = reads.sort { it.simpleName }
  """
  cutadapt \
    -e 0.12 \
    --times 3 \
    --overlap 5 \
    -j ${task.cpus} \
    -f fastq \
    -o "${name}_clean.R1.fastq.gz" \
    -p "${name}_clean.R2.fastq.gz" \
    -g spacer=TTTCTTATATGGG \
    -a R2_rc=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -a P7_rc=ATCTCGTATGCCGTCTTCTGCTTG \
    -a polyA=AAAAAAAAAAAAAAAAAAAA \
    -a rt_primer=AAGCAGTGGTATCAACGCAGAGTACAT \
    -A spacer_rc=CCCATATAAGAAA \
    -A R1_rc=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -A P5_rc=AGATCTCGGTGGTCGCCGTATCATT \
    -A polyA=AAAAAAAAAAAAAAAAAAAA \
    -A rt_primer=AAGCAGTGGTATCAACGCAGAGTACAT \
    ${r1} ${r2} > report.txt
  """
} 