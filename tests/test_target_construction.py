"""Test suite for target_construction for snakemake rules."""

import pandas as pd
import pytest

from lib import target_construction as tc


@pytest.fixture
def good_config_single():
    """Good test configuration."""
    return {
        "test1": {
            "tools": {
                "shapeit4": {"enabled": True},
                "eagle": {"enabled": False},
            },
            "genome_build": "b38",
        },
    }


@pytest.fixture
def bad_config_notools():
    """Bad test configuration with no tools."""
    return {
        "test1": {
            "genome_build": "b38",
        },
    }


@pytest.mark.parametrize(
    "chroms,exp_output",
    [
        (
            ["chr1"],
            [
                "results/shapeit4/test1.b38.chr1.vcf.gz",
                "results/shapeit4/test1.b38.chr1.vcf.gz.tbi",
            ],
        )
    ],
)
def test_good_config(good_config_single, chroms, exp_output):
    """Test that the good configuration works."""
    targets = tc.construct_phasing_targets(good_config_single, chroms, "test1")
    for t in targets:
        assert t in exp_output


def test_notools_config(bad_config_notools):
    """Test for key error with bad keys."""
    with pytest.raises(KeyError):
        tc.construct_phasing_targets(bad_config_notools, ["chroms"], "test1")
    with pytest.raises(KeyError):
        tc.construct_phasing_targets(bad_config_notools, ["chroms"], "test2")
