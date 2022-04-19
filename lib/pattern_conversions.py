"""Module for some pattern conversions."""

import typing


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
