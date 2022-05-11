"""Test suite for pattern conversions for snakemake rules."""

import pandas as pd
import pytest

from lib import pattern_conversions as pc


@pytest.mark.parametrize("test_in, exp_out", [("b38", "hg38")])
def test_build_string(test_in, exp_out):
    """Test that the expected build string makes sense."""
    assert exp_out == pc.convert_build_to_hg(test_in)


@pytest.mark.parametrize("test_in", ["hg19", "b37", "foolofatook", "", 37])
def test_build_string_error(test_in):
    """Test that we actually get an error on non b38 builds."""
    with pytest.raises(ValueError):
        pc.convert_build_to_hg(test_in)


@pytest.mark.parametrize(
    "input_file,exp_output",
    [("test.vcf.gz", "test.vcf.gz.tbi"), ("test.bcf", "test.bcf.csi")],
)
def test_index_determination(tmp_path, input_file, exp_output):
    """Test indexes for gzipped vcf/bcf."""
    p = tmp_path / input_file
    p.touch(exist_ok=True)
    index_path = pc.determine_index(p)
    out = tmp_path / exp_output
    # We convert to string here due to pathlib conflicts.
    assert index_path == str(out)


@pytest.mark.parametrize(
    "input_file", ["test.vcf", "a.txt", "x.bcf.gz", "hello.m3vcf.gz"]
)
def test_index_bad_ext(tmp_path, input_file):
    """Test with real files, but bad extensions."""
    p = tmp_path / input_file
    p.touch(exist_ok=True)
    with pytest.raises(ValueError):
        pc.determine_index(p)


@pytest.mark.parametrize("input_file", ["test.vcf", "test.vcf.gz", "x.bcf.gz"])
def test_index_nofile(tmp_path, input_file):
    """Test with not actual files."""
    with pytest.raises(FileNotFoundError):
        pc.determine_index(input_file)


@pytest.fixture
def good_hapmap_df():
    """Correct fields for a Hapmap formatted genetic map."""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "position": [1, 2],
            "rate (cM/Mb)": [0.01, 0.001],
            "map (cM)": [0.01, 0.02],
        }
    )


@pytest.fixture
def good_hapmap_df_shapeit4():
    """["pos", "chr", "cM"]."""
    return pd.DataFrame({"pos": [1, 2], "chr": ["chr1", "chr2"], "cM": [0.01, 0.02]})


@pytest.fixture
def good_hapmap_df_eagle():
    """["chr", "position", "Rate(cM/Mb)", "Genetic_Map(cM)"]."""
    return pd.DataFrame(
        {
            "chr": ["chr1", "chr2"],
            "position": [1, 2],
            "Rate(cM/Mb)": [0.01, 0.001],
            "Genetic_Map(cM)": [0.01, 0.02],
        }
    )


@pytest.fixture
def lower_columns_hapmap_df():
    """Fewer than expected columns should raise an error."""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "position": [1, 2],
        }
    )


@pytest.mark.parametrize("algo", ["shapeit4", "eagle"])
def test_convert_hapmap_genmap_good(
    tmp_path, good_hapmap_df, good_hapmap_df_eagle, good_hapmap_df_shapeit4, algo
):
    """Test for conversion of hapmap files into genetic maps."""
    p = tmp_path / f"test.{algo}.tsv"
    good_hapmap_df.to_csv(p, index=None, sep="\t")
    res_df = pc.convert_hapmap_genmap(p, algo)
    if algo == "eagle":
        pd.testing.assert_frame_equal(res_df, good_hapmap_df_eagle)
    elif algo == "shapeit4":
        pd.testing.assert_frame_equal(res_df, good_hapmap_df_shapeit4)


@pytest.mark.parametrize("algo", ["shapeit4", "eagle"])
def test_convert_hapmap_bad_cases(tmp_path, lower_columns_hapmap_df, algo):
    """Run conversion on bad cases."""
    p = tmp_path / "test.lower_column.tsv"
    lower_columns_hapmap_df.to_csv(p, index=None, sep="\t")
    with pytest.raises(ValueError):
        pc.convert_hapmap_genmap(p, algo)


@pytest.mark.parametrize("algo", ["shapeit", "eagle2", "x", ""])
def test_bad_algorithm(tmp_path, good_hapmap_df, algo):
    """Test that bad algorithm definitions get errors."""
    p = tmp_path / f"test.{algo}.tsv"
    good_hapmap_df.to_csv(p, index=None, sep="\t")
    with pytest.raises(ValueError):
        pc.convert_hapmap_genmap(p, algo)
