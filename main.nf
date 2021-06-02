#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { Preprocessing } from './modules/preprocessing.nf'
include { PrepareReferences} from './modules/prepare_references.nf'
include { PreselectReferences } from './modules/preselect_references.nf'
include { BuildConsensus } from './modules/consensus.nf'
include { AnalyzeConsensus } from './modules/consensus_analysis.nf'


/*
 * Prepare inputs
 */
cell_index = file(params.cell_index)
cell_whitelist = file(params.cell_whitelist)

channel
  .fromPath(params.reference_list)
  .splitCsv(header: true)
  .map { [it.subject, it.chain, file(it.file)] }
  .set { references }

channel
  .fromPath(params.data_list)
  .splitCsv(header: true)
  .map { ["${it.experiment}_L${it.lane}", [file(it.r1), file(it.r2)]]}
  .set { reads }

channel
  .fromPath(params.external_consensus_path)
  .map { [it.baseName, it] }
  .set { external_consensus_with_name }


/*
 * Workflow
 */
workflow {
  PrepareReferences(references)

  Preprocessing(reads, cell_index, cell_whitelist)

  PreselectReferences(
    Preprocessing.out,
    PrepareReferences.out.by_subject,
    PrepareReferences.out.by_subject_chain
  )

  BuildConsensus(
    PreselectReferences.out.reads,
    PrepareReferences.out.splitted,
    PrepareReferences.out.by_subject_chain,
    PrepareReferences.out.regions
  )

  AnalyzeConsensus(
    BuildConsensus.out,
    PreselectReferences.out.cell_references,
    PrepareReferences.out.regions
  )
}