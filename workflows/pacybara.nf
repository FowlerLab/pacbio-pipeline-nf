// PacybaraWorkflow:
//
//   sample_sheet ──► MAKE_PARAMETER_FILES ──► manifest.csv (one row per library)
//                                                   │
//                                                   ▼
//                                      RUN_PACYBARA (per library) ──► <output_dir>/<TITLE>/
//                                                   │
//                                                   ▼  (archive_inputs: true; successful runs only)
//                                      ARCHIVE_FASTQ ──► <completed_fastq_dir>/<fastq>
//
// Validation and parameter-file rendering live in Python
// (src/pacybara_workflow/); this file only fans the manifest out.

include { MAKE_PARAMETER_FILES } from '../modules/local/make_parameter_files/main.nf'
include { RUN_PACYBARA         } from '../modules/local/run_pacybara/main.nf'
include { ARCHIVE_FASTQ        } from '../modules/local/archive_fastq/main.nf'

workflow PacybaraWorkflow {
    def archive = params.archive_inputs.toString().toBoolean()
    if (!params.sample_sheet) {
        error "params.sample_sheet is required"
    }
    if (!params.output_dir) {
        error "params.output_dir is required"
    }
    if (archive && !params.completed_fastq_dir) {
        error "params.completed_fastq_dir is required when archive_inputs is true"
    }

    MAKE_PARAMETER_FILES(file(params.sample_sheet, checkIfExists: true))

    // checkIfExists: a missing fastq fails here, on the host, before any
    // Pacybara job is submitted.
    libraries = MAKE_PARAMETER_FILES.out.manifest
        .splitCsv(header: true)
        .map { r -> tuple(r.title, file(r.fastq, checkIfExists: true), file(r.param_file)) }

    RUN_PACYBARA(libraries)

    if (archive) {
        ARCHIVE_FASTQ(
            RUN_PACYBARA.out.workspace
                .join(libraries.map { title, fastq, _param_file -> tuple(title, fastq.toString()) })
                .map { title, _workspace, fastq -> tuple(title, fastq) }
        )
    }
}
