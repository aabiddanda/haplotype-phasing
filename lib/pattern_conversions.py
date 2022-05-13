"""Module for some pattern conversions."""

import typing
from pathlib import Path

import pandas as pd


def convert_build_to_hg(build_str: str) -> str:
    """Convert build versions to hg-version notation.

    NOTE: we have to do this because EAGLE refers to files with hg notation.
    """
    if build_str == "b38":
        return "hg38"
    else:
        raise ValueError(
            "build -> hg conversion not supported for {}".format(build_str)
        )


def determine_index(filename: str) -> str:
    """Determine the appropriate index for a VCF/BCF file."""
    path = Path(filename)
    if path.is_file():
        # print(path.suffixes, path.suffixes[-1], path.suffixes[-2])
        if path.suffix == ".bcf":
            return path.as_posix() + ".csi"
        elif path.suffixes[-1] == ".gz" and path.suffixes[-2] == ".vcf":
            return path.as_posix() + ".tbi"
        else:
            raise ValueError(f"{filename} is not a valid VCF/BCF file!")
    else:
        raise FileNotFoundError(f"{filename} is not a file!")


def convert_hapmap_genmap(filename: str, algo: str) -> pd.DataFrame:
    """Convert HapMap formatted genetic maps to alternative formats."""
    hapmap_df = pd.read_csv(filename, sep="\t")
    if len(hapmap_df.columns) != 4:
        raise ValueError(
            "HapMap Genetic Maps should only have [chrom, position, rate (cM/Mb), map (cM)] columns!"
        )
    hapmap_df.columns = ["chr", "pos", "rate", "cM"]
    if algo == "shapeit4":
        return hapmap_df[["pos", "chr", "cM"]]
    elif algo == "eagle":
        hapmap_df.columns = ["chr", "position", "Rate(cM/Mb)", "Genetic_Map(cM)"]
        return hapmap_df
    else:
        raise ValueError("Not a valid algorithm choice!")
