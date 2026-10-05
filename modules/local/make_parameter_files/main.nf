// MAKE_PARAMETER_FILES: wraps `pacybara-make-params`. Runs once: validates
// the whole sample sheet (failing the run before any Pacybara job starts)
// and writes one <TITLE>.txt per row plus manifest.csv (title, fastq,
// param_file), which the workflow splits into one RUN_PACYBARA task per row.
//
// run_date is passed only when set. Left null, the CLI uses today's date --
// and since the script text then doesn't change from day to day, a
// -resume on a later day reuses this task's cached output, keeping every
// TITLE (and therefore every RUN_PACYBARA cache entry) stable.
process MAKE_PARAMETER_FILES {
    label 'process_single'
    container "${params.container_image}"

    input:
    path sample_sheet

    output:
    path "manifest.csv", emit: manifest
    path "*.txt", emit: param_files

    script:
    def defaults = groovy.json.JsonOutput.toJson(params.pacybara_defaults)
    def data_dir_arg = params.data_dir ? "--data-dir '${file(params.data_dir)}'" : ''
    def run_date_arg = params.run_date ? "--run-date '${params.run_date}'" : ''
    """
    pacybara-make-params \\
        --sample-sheet '${sample_sheet}' \\
        --defaults '${defaults}' \\
        ${data_dir_arg} \\
        ${run_date_arg} \\
        --output-dir .
    """
}
