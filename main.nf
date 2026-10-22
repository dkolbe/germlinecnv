#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dkolbe/germlinecnv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/dkolbe/germlinecnv
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GERMLINECNV  } from './workflows/germlinecnv'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_germlinecnv_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_germlinecnv_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow DKOLBE_GERMLINECNV {

    take:
    batch_meta
    samplesheet // channel: samplesheet read in from --input
    fasta
    fai
    dict
    ucsc_genome
    gcnv_regions
    gcnv_model_interval_list // channel: [meta, gcnv_target_interval_list]
    gcnv_ploidy_model
    gcnv_caller_model
    graph_all
    graph_some
    txdb
    download_txdb

    main:
    //samplesheet.view()

    //
    // WORKFLOW: Run pipeline
    //
    GERMLINECNV (
        batch_meta,
        samplesheet,
        fasta,
        fai,
        dict,
        ucsc_genome,
        gcnv_regions,
        gcnv_model_interval_list,
        gcnv_ploidy_model,
        gcnv_caller_model,
        graph_all,
        graph_some,
        txdb,
        download_txdb,
        params.outdir,
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    batch_meta                  = [id:params.batch]
    ch_fasta                    = Channel.fromPath(params.fasta                   ).map { it -> [[id:it.simpleName], it] }.collect()
    ch_fai                      = Channel.fromPath(params.fai                     ).map { it -> [[id:it.simpleName], it] }.collect()
    ch_dict                     = Channel.fromPath(params.dict                    ).map { it -> [[id:it.simpleName], it] }.collect()
    ucsc_genome                 = params.ucsc_genome
    ch_gcnv_regions             = Channel.fromPath(params.gcnv_regions            ).map { it -> [[id:it.simpleName], it] }.collect()
    ch_gcnv_model_interval_list = Channel.fromPath(params.gcnv_model_interval_list).map { it -> [[id:it.simpleName], it] }.collect()
    ch_gcnv_ploidy_model        = Channel.fromPath(params.gcnv_ploidy_model       ).map { it -> [[id:it.simpleName], it] }.collect()
    ch_gcnv_caller_model        = Channel.fromPath(params.gcnv_caller_model       ).map { it -> [[id:it.simpleName], it] }.collect()
    graph_all                   = params.graph_all
    graph_some                  = params.graph_some
    if (params.txdb) {
                        ch_txdb = Channel.fromPath(params.txdb                    ).map { it -> [[id:it.simpleName], it] }.collect()
                        download_txdb = false
    }
    else {
                        ch_txdb = Channel.empty()
                        download_txdb = true
    }
    //ch_txdb                     = params.txdb ?
    //                              Channel.fromPath(params.txdb                    ).map { it -> [[id:it.simpleName], it] }.collect()
    //                              : Channel.empty()

    //
    // WORKFLOW: Run main workflow
    //
    DKOLBE_GERMLINECNV (
        batch_meta,
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_fasta,
        ch_fai,
        ch_dict,
        ucsc_genome,
        ch_gcnv_regions,
        ch_gcnv_model_interval_list,
        ch_gcnv_ploidy_model,
        ch_gcnv_caller_model,
        graph_all,
        graph_some,
        ch_txdb,
        download_txdb
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
