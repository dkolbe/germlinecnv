/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_germlinecnv_pipeline'
include { GATK4_COLLECTREADCOUNTS } from '../modules/nf-core/gatk4/collectreadcounts/main'
include { SAMTOOLS_INDEX          } from '../modules/nf-core/samtools/index'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GERMLINECNV {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_fai
    ch_dict
    ch_interval_list
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
        .mix(GATK4_COLLECTREADCOUNTS.out.hdf5)
        .collect { _meta, file -> [file] }
        //.map { tsv -> [[id: val_pon_name], tsv] }
        .set { ch_readcounts_out }

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
