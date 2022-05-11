"""Test suite for validators for snakemake rules."""

import pandas as pd
import pytest

from lib.validators import FamFileValidator


@pytest.fixture
def good_fam_file():
    """Good fam file dataframe example."""
    return pd.DataFrame(
        {
            "sampleID": ["A", "B"],
            "motherID": ["D", "E"],
            "fatherID": ["G", "H"],
        }
    )


def test_good_fam_file(good_fam_file):
    """Testing that a good fam file."""
    validator = FamFileValidator("")
    validator.fam_df = good_fam_file
    validator.check_names(["A", "B", "D", "E", "G", "H"])
    validator.check_names(["A", "B", "C", "D", "E", "F", "G", "H", "I"])
    # The empty string should also validate
    validator = FamFileValidator("")
    assert not validator.validate_fam(["A", "B"])


def test_good_fam_file_from_file(tmp_path, good_fam_file):
    """Testing that validation directly from a file works well."""
    p = tmp_path / "test.tsv"
    good_fam_file.to_csv(p, index=None, sep="\t")
    validator = FamFileValidator(p)
    assert validator.validate_fam(["A", "B", "D", "E", "G", "H"])


@pytest.mark.parametrize("vcfids", [["A", "E", "H"], ["X"], []])
def test_incomplete_vcfids(good_fam_file, vcfids):
    """Test that all families are contained in the VCF files."""
    validator = FamFileValidator("")
    validator.fam_df = good_fam_file
    with pytest.raises(AssertionError):
        validator.check_names(vcfids)
