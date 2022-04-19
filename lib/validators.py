"""Validators for assorted files involved in this pipeline."""

import warnings

import numpy as np
import pandas as pd


class FamFileValidator:
    """Class for validating trio files for evaluation of phase accuracy."""

    def __init__(self, fam_file):
        """Initialize the fam-file object."""
        self.fam_file = fam_file
        self.fam_df = None

    def read_fam_file(self):
        """Read in the trio validation."""
        if self.fam_file != "":
            self.fam_df = pd.read_csv(self.fam_file, sep="\s+", engine="python")  # noqa
            self.fam_df.columns = ["sampleID", "fatherID", "motherID"]

    def check_names(self, vcf_ids):
        """Check that all of the trio names are inside the VCF files."""
        assert self.fam_df is not None
        for _, row in self.fam_df.iterrows():
            assert np.isin(row["sampleID"], vcf_ids)
            assert np.isin(row["fatherID"], vcf_ids)
            assert np.isin(row["motherID"], vcf_ids)

    def validate_fam(self, vcf_ids):
        """Full validation of a fam file."""
        if self.fam_file != "":
            self.read_fam_file()
            self.check_names(vcf_ids)
            return True
        else:
            return False
