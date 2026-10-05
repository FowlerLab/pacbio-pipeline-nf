// ARCHIVE_FASTQ: moves a successfully mapped library's fastq.gz from the
// data directory into completed_fastq_dir. Only runs when
// params.archive_inputs is true, and only for libraries RUN_PACYBARA
// finished. It moves the host file itself (a val, not a staged path), so it
// runs on the head node with no container (nextflow.config).
//
// Refuses to overwrite: a same-named fastq already in completed_fastq_dir
// fails this task and leaves both files where they are.
process ARCHIVE_FASTQ {
    tag "${title}"
    executor 'local'

    input:
    tuple val(title), val(fastq)

    output:
    val title, emit: archived

    script:
    def dest = "${file(params.completed_fastq_dir)}/${file(fastq).name}"
    """
    if [ -e '${dest}' ]; then
        echo "ERROR: ${dest} already exists; not moving ${fastq}" >&2
        exit 1
    fi
    mkdir -p '${file(params.completed_fastq_dir)}'
    mv '${fastq}' '${dest}'
    """
}
