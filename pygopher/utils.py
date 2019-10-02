
from subprocess import run
import pandas as pd
import yaml


def run_pgopher(filepath: str, *args):
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
    keywords = ["pgo"]
    if args:
        keywords.extend(args)
    keywords.append(filepath)
    process = run(
        keywords,
        capture_output=True
    )
    return process


def parse_linelist(text_stream, delimiter=","):
    """
    Takes the standard output of the PGopher run that calculates
    a linelist.
    
    Parameters
    ----------
    text_stream : str
        [description]
    
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
    lines = lines[index + 1:]
    # First line is not needed, and second line contains the header
    #_ = lines.pop(0)
    labels = lines.pop(0).split(delimiter)
    labels = [name.replace('\\', "") for name in labels]
    labels = [name.replace('"', "") for name in labels]
    data = list()
    for line in lines:
        line = line.split(",")
        # The branch is not read correctly
        line[-3:-1] = ["".join(line[-3:-1])]
        data.append(line)
    df = pd.DataFrame(data[:-2], columns=labels)
    return df


def parse_partition_func(text_stream):
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
    

def read_yaml(filepath):
    with open(filepath, "r") as stream:
        try:
            data = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return data
