process VISUALIZE_CNV {
    tag "$meta.id"
    label 'process_medium'

    # FIXME find and pull an appropriate container

    input:


    output:
    tuple val(meta), path("*.png"), emit: png

    script:
    template 'visualize-CNV.R'

    stub:

}
