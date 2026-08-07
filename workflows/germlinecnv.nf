/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_germlinecnv_pipeline'

include { CSVTK_CONCAT                        } from '../modules/nf-core/csvtk/concat/main'
include { GATK4_COLLECTREADCOUNTS             } from '../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../modules/nf-core/gatk4/determinegermlinecontigploidy/main'
include { GATK4_GERMLINECNVCALLER             } from '../modules/nf-core/gatk4/germlinecnvcaller/main'
include { GATK4_POSTPROCESSGERMLINECNVCALLS   } from '../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main'
include { SAMTOOLS_INDEX                      } from '../modules/nf-core/samtools/index'

include { AGGREGATE_COUNTS                    } from '../modules/local/aggregate_counts/main'
include { CNV_VCF_TO_TABLE                    } from '../modules/local/cnv_vcf_to_table/main'
include { DOWNLOAD_DB                         } from '../modules/local/download_db/main'
include { VISUALIZE_CNV                       } from '../modules/local/visualize_cnv/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GERMLINECNV {

    take:
    batch_meta
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_fai
    ch_dict
    ucsc_genome
    ch_regions
    ch_interval_list
    ch_ploidy_model
    ch_gcnv_caller_model
    ch_txdb
    outdir

    main:

    // Check input samples for missing indices and create them
    ch_samplesheet
        .branch { meta, bam, bai ->
            alignment_with_index: bai.size() > 0
            return [meta, bam, bai]
            alignment_without_index: bai.size() == 0
            return [meta, bam]
        }
        .set { ch_for_mix }

    SAMTOOLS_INDEX(ch_for_mix.alignment_without_index)

    ch_index = SAMTOOLS_INDEX.out.index

    // Collect alignment files and their indices
    ch_for_mix.alignment_without_index
        .join(ch_index)
        .mix(ch_for_mix.alignment_with_index)
        //.combine(GATK4_PREPROCESSINTERVALS.out.interval_list.map { it -> it[1] })
        .combine(ch_interval_list.map { it -> it[1] })
        .set { ch_readcounts_in }

    GATK4_COLLECTREADCOUNTS(
        ch_readcounts_in,
        ch_fasta,
        ch_fai,
        ch_dict
    )

    GATK4_COLLECTREADCOUNTS.out.tsv
        .map { meta, tsv -> return [meta, tsv, [], []] }
        .set { ch_readcounts_out }

    // ch_readcounts_out.view()
    GATK4_DETERMINEGERMLINECONTIGPLOIDY ( ch_readcounts_out, ch_ploidy_model, [] )

    GATK4_COLLECTREADCOUNTS.out.tsv
        .join(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
        .combine(ch_gcnv_caller_model)
        .map({ meta, tsv, calls, meta2, model -> return [meta, tsv, [], calls, model ]})
        .set {ch_gcnvc_in}

    GATK4_GERMLINECNVCALLER ( ch_gcnvc_in )

    GATK4_GERMLINECNVCALLER.out.casecalls
        .map { meta, model_calls ->
            return [meta.id, meta, model_calls]
        }
        .groupTuple(by: 0)
        .map { id, metas, model_calls ->
            def meta = metas[0] // All metas should be the same for a given sample
            def models = model_calls.collect { it.toString() }
            return [meta, models]
        }
        .combine(ch_gcnv_caller_model.collect{it[1]}.toList())
        .join(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
        .set {ch_postproc_in}

    GATK4_POSTPROCESSGERMLINECNVCALLS ( ch_postproc_in )

    CNV_VCF_TO_TABLE( GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments )

    CNV_VCF_TO_TABLE.out.table
        .map { meta, table -> return [table]}
        .collect()
        .map { table -> return [[id:batch_meta.id], table] }
        .set {ch_concat_in}

    CSVTK_CONCAT( ch_concat_in, "tsv", "tsv" )

    GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised
         .collect(flat:false)
         .map { it.transpose() }
         .set {ch_aggregate_in}

    AGGREGATE_COUNTS( batch_meta, ch_aggregate_in )

    ch_txdb | ifEmpty { ['NEED_TXDB', ucsc_genome] } | DOWNLOAD_DB
    // DOWNLOAD_DB( ucsc_genome )
    // DOWNLOAD_DB.out.txdb

    CNV_VCF_TO_TABLE.out
        .map { meta, calls -> return [meta.id, meta, calls] }
        .join(GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised.map({ meta, counts -> return [meta.id, meta, counts] }))
        .map { id, meta1, calls, meta2, counts -> return [meta1, calls, counts] }
        .set {ch_visualize_in}

    VISUALIZE_CNV( ch_visualize_in,
                   AGGREGATE_COUNTS.out.table,
                   ch_regions.map { meta, file -> return file },
                   ch_txdb,
                   ucsc_genome)

    def ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'germlinecnv_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    indices        = SAMTOOLS_INDEX.out.index
    readcounts     = ch_readcounts_out
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
