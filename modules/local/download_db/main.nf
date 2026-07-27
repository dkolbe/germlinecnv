process DOWNLOAD_DB {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomicfeatures_bioconductor-genomicranges_bioconductor-gviz:6b706937e78bb65b':
        'community.wave.seqera.io/library/bioconductor-annotationdbi_bioconductor-genomicfeatures_bioconductor-genomicranges_bioconductor-gviz:4ef161684da1a400' }"

    input:
    val(ucsc_genome)

    output:
    path(txdb), emit: txdb

    script:
    template 'download-DB.R'

    stub:
    touch "txdb.refseq.$ucsc_genome.sqlite"
}
