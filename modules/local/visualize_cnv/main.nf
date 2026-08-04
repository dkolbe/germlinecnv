process VISUALIZE_CNV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomicfeatures_bioconductor-genomicranges_bioconductor-gviz:6b706937e78bb65b':
        'community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomicfeatures_bioconductor-genomicranges_bioconductor-gviz:4ef161684da1a400' }"

    input:
    tuple val(meta), path(calls_table), path(counts)
    path(aggregated_counts)
    path(regions)
    path(txdb)
    val(ucsc_genome)

    output:
    tuple val(meta), path("*.png"), emit: png

    script:
    //template 'visualize-CNV.R'
    """
    echo "Commented out Rscript..."
    """

    stub:
    """
    echo "Commented out Rscript..."
    """

}
