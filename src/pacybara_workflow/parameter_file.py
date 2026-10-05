"""Render a Pacybara parameter file for one library.

The template is `pacybara_prep`'s `make_parameter.py` (by Dan Holmes,
github.com/dlholmes2117/pacybara_prep, BSD-3-Clause) byte for byte, except
that INFASTQ is filled in and the Pacybara arguments come from the
library's parameters instead of being hard-coded. With pacybara_prep's
values as the defaults, the output is identical to that script's apart from
INFASTQ (tests/unit/test_parameter_file.py checks this against every row of
its example sheet).
"""

from __future__ import annotations

from pacybara_workflow.sample_sheet import Library

TEMPLATE = """#########################
#Pacybara parameter file
#########################

#BEGIN ARGUMENTS
#Experiment title
TITLE={title}
# Long read sequencing input file (fastq.gz)
INFASTQ={infastq}
# Workspace directory
WORKSPACE={title}
# The barcode degeneracy sequence used in the amplicon below
BARCODE={barcode}
# The start position of the ORF in the amplicon below (1-based)
ORFSTART={orfstart}
# The end position of the ORF in the amplicon below (1-based)
ORFEND={orfend}
# The maximum number of low-quality bases allowed in a given barcode.
MAXQDROPS={maxqdrops}
# The minimum average quality score allowed in a given barcode.
MINBCQ={minbcq}
#The minimum Jaccard coefficient (relative overlap of variants) for a cluster merge
MINJACCARD={minjaccard}
#The minimum number of variants two cluster need to have in common to merge
MINMATCHES={minmatches}
#The maximum number of errors allowed between two barcode reads
MAXDIFF={maxdiff}
#The minimum Q-score for a variant basecall to be considered real based on a single read alone
MINQUAL={minqual}
#Cluster based on which barcodes? uptag, downtag, or virtual
CLUSTERMODE={clustermode}
#END ARGUMENTS

#BEGIN AMPLICON SEQUENCE
>{seq_title}
{amplicon}
#END AMPLICON SEQUENCE
"""


def render_parameter_file(library: Library, run_date: str) -> str:
    """The parameter file text for `library`, titled with `run_date`."""
    return TEMPLATE.format(
        title=library.title(run_date),
        infastq=library.fastq,
        orfstart=library.orfstart,
        orfend=library.orfend,
        seq_title=library.name,
        amplicon=library.amplicon,
        **library.parameters,
    )
