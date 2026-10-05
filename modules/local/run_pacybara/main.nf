// RUN_PACYBARA: one library through pacybara_simplex.sh, Pacybara's
// single-machine mode (Nextflow, not clusterutil, does the scheduling).
//
// The generated parameter file names the fastq by its host path; the run
// uses a copy pointed at the staged fastq and at a workspace named after the
// TITLE. The original file is then copied into that workspace, which
// pacybara.sh (but not pacybara_simplex.sh) would have done itself.
// Publishes <output_dir>/<TITLE>/: the alignment, <fastq>_extract/,
// <fastq>_clustering/ (clusters_transl*.csv.gz, qc/) and <TITLE>.txt.
//
// errorStrategy 'ignore': one failed library doesn't stop the others; its
// fastq is then never archived.
process RUN_PACYBARA {
    tag "${title}"
    label 'process_high'
    errorStrategy 'ignore'
    container "${params.container_image}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(title), path(fastq), path(param_file)

    output:
    tuple val(title), path(title), emit: workspace

    script:
    """
    sed -e 's|^INFASTQ=.*|INFASTQ=${fastq}|' \\
        -e 's|^WORKSPACE=.*|WORKSPACE=${title}|' \\
        '${param_file}' > run_parameters.txt
    pacybara_simplex.sh --cpus ${task.cpus} run_parameters.txt
    cp '${param_file}' '${title}/'
    rm -rf tmp
    """
}
