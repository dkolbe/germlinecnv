process DOWNLOAD_DB {
    tag "txdb:$ucsc_genome"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomeinfodbdata_bioconductor-genomicfeatures_bioconductor-genomicranges_pruned:2256b92d8b38ab36':
        'community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomeinfodbdata_bioconductor-genomicfeatures_bioconductor-genomicranges_pruned:280edc874ea7fec2' }"

    input:
    val(ucsc_genome)

    output:
    tuple val(id:txdb_name), path(txdb_name), emit: txdb

    script:
    txdb_name = "txdb.refseq.${ucsc_genome}.sqlite"
    template 'download-DB.R'

    stub:
    txdb_name = "txdb.refseq.${ucsc_genome}.sqlite"
    """
    touch "$txdb_name"
    """
}
