process AGGREGATE_COUNTS {
    tag "$meta.id"
    label "process_single"

    input:
    tuple val(meta), path(counts)

    output:
    path("aggregate_counts.tsv"), emit: table

    script:
    def out_file = "aggregate_counts.tsv"
    def first_file = counts.collect()[0]
    def header = meta ? meta.id.collect().join('\t') : ""
    def paste_argument = counts ? counts.collect{"<( grep -v \"^@\" $it | cut -f 4 | tail -n +2)"}.join(' ') : ""
    """
    grep -v "^@" ${first_file} | cut -f 1-3 | tail -n +2 > row_labels.txt
    echo "CONTIG	START	END	${header}" > ${out_file}
    paste row_labels.txt ${paste_argument} >> ${out_file}
    """

    stub:
    """
    touch aggregate_counts.tsv
    """
}
