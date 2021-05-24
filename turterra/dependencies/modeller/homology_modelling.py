import sys
from pathlib import Path

import modeller
import numpy as np
from modeller.automodel import automodel, dopehr_loopmodel, refine

import typing
path_type = typing.Union[str, Path]

import subprocess

def run_muscle(in_file, out_file):
    subprocess.call(["muscle", "-in", in_file, "-out", out_file])
    
def get_sequences_from_fasta_yield(fasta_file: typing.Union[str, Path], prune_headers: bool = True) -> tuple:
    """
    Returns (accession, sequence) iterator
    Parameters
    ----------
    fasta_file
    prune_headers
        only keeps accession upto first /
    Returns
    -------
    (accession, sequence)
    """
    with open(fasta_file) as f:
        current_sequence = ""
        current_key = None
        for line in f:
            if not len(line.strip()):
                continue
            if "==" in line:
                continue
            if ">" in line:
                if current_key is None:
                    if "/" in line and prune_headers:
                        current_key = line.split(">")[1].split("/")[0].strip()
                    else:
                        current_key = line.split(">")[1].strip()
                    if "|" in current_key and prune_headers:
                        current_key = current_key.split("|")[1].strip()
                else:
                    yield (current_key, current_sequence)
                    current_sequence = ""
                    if "/" in line and prune_headers:
                        current_key = line.split(">")[1].split("/")[0].strip()
                    else:
                        current_key = line.split(">")[1].strip()
                    if "|" in current_key and prune_headers:
                        current_key = current_key.split("|")[1].strip()
            else:
                current_sequence += line.strip()
        yield (current_key, current_sequence)
        
def get_sequences_from_fasta(fasta_file: typing.Union[str, Path], prune_headers: bool = True) -> dict:
    """
    Returns dict of accession to sequence from fasta file
    Parameters
    ----------
    fasta_file
    prune_headers
        only keeps accession upto first /
    Returns
    -------
    {accession:sequence}
    """
    return {key: sequence for (key, sequence) in get_sequences_from_fasta_yield(fasta_file, prune_headers)}

def get_file_parts(input_filename: path_type) -> tuple:
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


def align_templates_to_sequence(sequence_file: path_type, templates_file: path_type, alignment_folder: path_type) -> path_type:
    """
    clustalo align sequence to template sequences

    Parameters
    ----------
    sequence_file
        file containing a single sequence
    templates_file
        subsequence file of templates
    alignment_folder
        where to store the output file
    hmm_file
        use hmm to align

    Returns
    -------
    alignment filename
    """
    _, sequence_name, _ = get_file_parts(sequence_file)
    _, templates_name, _ = get_file_parts(templates_file)
    sequences = get_sequences_from_fasta(sequence_file)
    template_sequences = get_sequences_from_fasta(templates_file)
    unaligned_file = Path(alignment_folder) / f"{sequence_name}_{templates_name}.fasta"
    with open(unaligned_file, "w") as f:
        for key in sequences:
            f.write(f">{key}\n{sequences[key]}\n")
        for key in template_sequences:
            f.write(f">{key}\n{template_sequences[key]}\n")
    aligned_file = Path(alignment_folder) / f"{sequence_name}_{templates_name}_aln.fasta"
    run_muscle(str(unaligned_file), aligned_file)
    return aligned_file


def convert_fasta_to_ali(sequence_name: str, fasta_file: path_type, ali_file: path_type=None) -> path_type:
    """
    Converts fasta of sequence and templates to ali format

    Parameters
    ----------
    sequence_name
        name of sequence, rest are assumed as template structures
    fasta_file
    ali_file
    templates_start_end
        dict with start and end residue identifiers and chains for each template

    Returns
    -------
    ali filename
    """
    if ali_file is None:
        path, name, _ = get_file_parts(fasta_file)
        ali_file = Path(path) / f"{name}.ali"
    sequences = get_sequences_from_fasta(fasta_file)
    with open(ali_file, "w") as f:
        for key in sequences:
            if key == sequence_name:
                f.write(f">P1;{key}\nsequence:{key}:::::::0.00: 0.00\n{sequences[key].upper()}*\n")
            else:
                start_chain, start_pos, end_chain, end_pos = '.', '', '.', ''
                f.write(f">P1;{key}\nstructureX:{key}:{start_pos}:{start_chain}:{end_pos}:{end_chain}:::0.00: 0.00\n{sequences[key].upper()}*\n")
    return ali_file


def make_env(pdb_directory: path_type, rand_seed: int=5):
    """
    Make a modeller environ with a random seed.
    Change seed if parallelizing.

    Parameters
    ----------
    pdb_directory
    rand_seed

    Returns
    -------
    environ

    """
    env = modeller.environ(rand_seed=rand_seed)
    env.io.atom_files_directory = [str(pdb_directory)]
    env.io.hetatm = True
    return env


def get_ndope_score(pdb_file: path_type) -> float:
    """
    Retrieves n-DOPE score from PDB file

    Parameters
    ----------
    pdb_file

    Returns
    -------
    score
    """
    with open(pdb_file) as f:
        for line in f:
            if "Normalized DOPE score" in line:
                score = line.split(":")[-1].strip()
                return float(score)


def model(name: str, templates: list, ali_file: path_type, pdb_dir: path_type, model_path: path_type, rand_seed: int, num_models: int=500, loop_refine: bool=False):
    """
    Makes models for name.
    DELETES all intermediate files
    MOVES final model pdbs and log file to model_path

    Parameters
    ----------
    name
    templates
    ali_file
    model_path
    pdb_dir
    rand_seed
    num_models
    loop_refine
    """
    env = make_env(pdb_dir, rand_seed)
    if loop_refine:
        a = dopehr_loopmodel(env,
                             alnfile=str(ali_file),
                             knowns=tuple(templates),
                             sequence=name,
                             assess_methods=(modeller.automodel.assess.DOPE,
                                             modeller.automodel.assess.normalized_dope))
        a.md_level = None
        a.loop.starting_model = 1
        a.loop.ending_model = 5
        a.loop.md_level = refine.fast
    else:
        a = automodel(env,
                      alnfile=str(ali_file),
                      knowns=tuple(templates),
                      sequence=name,
                      assess_methods=(modeller.automodel.assess.DOPE,
                                      modeller.automodel.assess.normalized_dope))
    a.starting_model = 1
    a.ending_model = num_models
    a.make()
    files = Path.cwd().glob(name + ".*")
    log_file = Path(model_path) / f"{name}.log"
    with open(log_file, "w") as f:
        for file_path in files:
            _, filename, ext = get_file_parts(file_path)
            if ext == ".pdb":
                f.write(f"{filename}\t{get_ndope_score(file_path)}\n")
                new_file = Path(model_path) / f"{filename}{ext}"
                file_path.rename(new_file)
            else:
                file_path.unlink()


def parse_log_file(input_model_dir: path_type, log_file: path_type, ignore_none=True):
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


def main():
    """
    Makes homology models for a single key
    Takes command line arguments
    key, index, template_dir, ali_dir, model_dir
    RUN USING parallel_functions.py

    """
    arguments = sys.argv
    print(arguments)
    assert len(arguments) == 7
    key, index, template_file, template_dir, ali_dir, model_dir = arguments[1:]
    template_dir = Path(template_dir)
    ali_dir = Path(ali_dir)
    model_dir = Path(model_dir)
    _, templates_name, _ = get_file_parts(template_file)
    templates = [get_file_parts(filename)[1] for filename in template_dir.glob("*.pdb")]
    ali_file = ali_dir / f"{key}_{templates_name}_aln.ali"
    if not ali_file.exists():
        print(f"{ali_file} missing")
    if ali_file.exists() and not (model_dir / f"{key}.log").exists():
        model(key, templates, ali_file, template_dir, model_dir, int(index), num_models=250, loop_refine=False)


if __name__ == '__main__':
    main()
