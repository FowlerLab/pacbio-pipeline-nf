#!/usr/bin/env nextflow
// pacybara-workflow. PacBio sample sheet -> one Pacybara parameter file per
// library -> pacybara_simplex.sh per library -> (optionally) completed
// fastqs moved out of the data directory.
//
//   nextflow run . -params-file params.yaml \
//       --sample_sheet /path/to/sample_sheet.csv \
//       --data_dir /path/to/to_process/data \
//       --output_dir /path/to/pacybara_output [-profile apptainer]
//
// See README.md.

include { PacybaraWorkflow } from './workflows/pacybara'

workflow {
    PacybaraWorkflow()
}
