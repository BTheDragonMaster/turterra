from sys import argv
from pathlib import Path
import os
import numpy as np
from shutil import copyfile

def get_file_parts(input_filename):
    """
    Gets directory path, name, and extension from a filename
    Parameters
    ----------
    input_filename
    Returns
    -------
    (path, name, extension)
    """
    input_filename = Path(input_filename)
    path = str(input_filename.parent)
    extension = input_filename.suffix
    name = input_filename.stem
    return path, name, extension

def parse_log_file(input_model_dir, log_file, ignore_none=True):
    """
    Parses a log file to get the NDOPE scores of each pdb file
    If the score is None, it tries to get te score directly from the PDB file
    If this doesn't work and ignore_none is False
        then it runs MODELLER on the PDB file to get the NDOPE scores (this takes longer)
    Parameters
    ----------
    input_model_dir
    log_file
    ignore_none
    Returns
    -------
    """
    filenames, scores = [], []
    with open(log_file) as f:
        for line in f:
            name, score = line.strip().split("\t")
            if score == 'None':
                filename = Path(input_model_dir) / f"{name}.pdb"
                score = get_ndope_score(filename)
                if score is None:
                    if ignore_none:
                        continue
                    else:
                        try:
                            score = get_ndope_score_expensive(filename)
                            with open(filename) as f1:
                                lines = f1.readlines()
                            with open(filename, 'w') as f1:
                                for i, l in enumerate(lines):
                                    if i == 2:
                                        f1.write(f"REMARK   6 Normalized DOPE score: {score}\n")
                                    f1.write(l)
                        except modeller.FileFormatError:
                            continue
            filenames.append(name)
            scores.append(float(score))
    sorted_scores = np.argsort(scores)
    return filenames, scores, sorted_scores



def copy_best_models(input_model_dir, output_model_dir, num_best = 1, ignore_none=True):
    """
    Copies the top num_best models (based on their NDOPE scores) of each accession in input_model_dir to output_model_dir

    Parameters
    ----------
    input_model_dir
    output_model_dir
    num_best
    ignore_none
        if True, ignores PDB files without NDOPE information
        if False, calculates NDOPE on these files by running MODELLER
    """
    log_files = Path(input_model_dir).glob("*.log")
    for log_file in log_files:
        _, name, _ = get_file_parts(log_file)
        output_log_file = Path(output_model_dir) / f"{name}.log"
        if not output_log_file.exists():
            copyfile(log_file, output_log_file)
            filenames, scores, sorted_scores = parse_log_file(input_model_dir, log_file, ignore_none=ignore_none)
            for i in range(num_best):
                copyfile(Path(input_model_dir) / f"{filenames[sorted_scores[i]]}.pdb",
                         Path(output_model_dir) / f"{name}_{i}.pdb")




if __name__ == "__main__":
    copy_best_models(argv[1], argv[2])
