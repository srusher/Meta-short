process QUAST {
    tag "${meta.id}"
    errorStrategy 'ignore'
    conda "bioconda::quast=5.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2' :
        'biocontainers/quast:5.0.2--py37pl526hb5aa323_2' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("QUAST/*report_rawassemblies.tsv")   , emit: qc
    path "QUAST/*report_rawassemblies.tsv", emit: report
    path "versions.yml"                  , emit: versions

    script:
    """
    quast.py --threads "${task.cpus}" -l "${meta.id}" "${assembly}" -o "QUAST"
    cp QUAST/report.tsv QUAST/${meta.id}_report_rawassemblies.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """
}