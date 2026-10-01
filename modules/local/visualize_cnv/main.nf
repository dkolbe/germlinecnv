process VISUALIZE_CNV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomeinfodbdata_bioconductor-genomicfeatures_bioconductor-genomicranges_pruned:2256b92d8b38ab36':
        'community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomeinfodbdata_bioconductor-genomicfeatures_bioconductor-genomicranges_pruned:280edc874ea7fec2' }"

    input:
    tuple val(meta), path(calls_table), path(counts)
    path(aggregated_counts)
    path(regions)
    path(txdb)
    val(ucsc_genome)

    output:
    tuple val(meta), path("*.png"), emit: png, optional: true

    script:
    """
    echo "Before Rscript..."
    touch dummy.png
    """
    template 'visualize-CNV.R'

    stub:
    """
    echo "Commented out Rscript..."
    """

}
