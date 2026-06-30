process CNV_VCF_TO_TABLE {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.CNV_calls.tsv"), emit: table

    script:
    """
    echo "SAMPLE\tCONTIG\tSTART\tEND\tREF\tALT\tGenotype\tCopyNumber\tQUAL\n" > ${meta.id}.CNV_calls.tsv
    bcftools query \\
         -i 'GT="alt"' \\
         -f "[%SAMPLE]\t%CHROM\t%POS\t%INFO/END\t%REF\t%ALT\t[%GT]\t[%CN]\t[%NP]\t%QUAL\n" \\
         ${vcf} \\
         >> ${meta.id}.CNV_calls.tsv 
    """

    stub:
    """
    touch ${meta.id}.CNV_calls.tsv
    """
}
