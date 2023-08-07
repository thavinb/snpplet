process CREATE_SAMPLE_MAP {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path valid_samplesheet

    output:
    path "all.sample_map", emit: sample_map

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in cenmig/snpplet/bin/
    """
    tail +2 ${valid_samplesheet} | \
        cut -f1 -d ',' | \
        awk '{print \$1,\$1".vcf.gz"}' OFS="\\t" > all.sample_map

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
