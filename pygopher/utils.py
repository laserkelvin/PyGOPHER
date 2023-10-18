from __future__ import annotations

from pathlib import Path
from shutil import which
from subprocess import run
from typing import Any

import pandas as pd
import yaml


def is_pgopher_in_path() -> bool:
    """
    Simple function that checks if PGopher is present in PATH.

    Returns
    -------
    bool
        True if found, False if not
    """
    is_present = which("pgo") is not None
    if not is_present:
        raise FileNotFoundError(
            "pgo binary not found in PATH;"
            " please download the PGopher linux distribution"
            " and ensure the `pgo` executable is in PATH."
        )
    return is_present


def run_pgopher(filepath: str | Path, *args):
    """
    Wrapper to run PGopher from Python.

    Optional arguments can be supplied, which will be appended
    to the list of flags to run PGopher with.

    Parameters
    ----------
    filepath : str
        Path to the PGopher file to be run

    Returns
    -------
    CompletedProcess object

    """
    assert is_pgopher_in_path()
    if isinstance(filepath, str):
        filepath = Path(filepath)
    assert filepath.exists(), f"{filepath} not found."
    keywords = ["pgo"]
    if args:
        keywords.extend(args)
    keywords.append(filepath)
    process = run(keywords, capture_output=True)
    return process


def parse_linelist(text_stream, delimiter=",") -> pd.DataFrame:
    """
    Takes the standard output of the PGopher run that calculates
    a linelist.

    Parameters
    ----------
    text_stream : str
        Text stream containing line list information

    Returns
    -------
    DataFrame
        Pandas DataFrame containing the linelist
    """
    # Make sure to convert into text, not bytes
    text_stream = str(text_stream)
    lines = text_stream.split("\\n")
    for index, line in enumerate(lines):
        if "Line list" in line:
            break
    lines = lines[index + 1 :]
    # First line is not needed, and second line contains the header
    # _ = lines.pop(0)
    labels = lines.pop(0).split(delimiter)
    for symbol in ["\\", '"']:
        labels = [name.replace(symbol, "") for name in labels]
    data = list()
    for line in lines:
        line = line.split(",")
        # The branch is not read correctly
        line[-3:-1] = ["".join(line[-3:-1])]
        data.append(line)
    df = pd.DataFrame(data[:-2], columns=labels)
    return df


def parse_partition_func(text_stream) -> pd.DataFrame:
    text_stream = str(text_stream)
    lines = text_stream.split("\\n")
    data = list()
    read = False
    for line in lines:
        if read is True:
            if "levels" in line:
                pass
            else:
                data.append(line.split())
        if "T/K" in line:
            labels = line.split()[:3]
            read = True
    df = pd.DataFrame(data, columns=labels)
    return df


def read_yaml(filepath) -> dict[str, Any]:
    with open(filepath, "r") as stream:
        try:
            data = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise yaml.YAMLError(f"Cannot read {filepath} as YAML.") from exc
    return data
