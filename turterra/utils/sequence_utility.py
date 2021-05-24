from collections import Counter
import typing
from pathlib import Path
import subprocess
import os


def run_muscle(guide_alignment, in_file, out_file):

    command = ['muscle', '-quiet', '-profile',  '-in1', guide_alignment, '-in2', in_file, '-out', out_file]

    subprocess.check_call(command)

def add_sequences_to_alignment(new_fasta_file, old_fasta_file, out_file):
    run_muscle(old_fasta_file, new_fasta_file, out_file)
    alignment = get_sequences_from_fasta(out_file)
    return alignment


def get_sequences_from_fasta_yield(fasta_file: typing.Union[str, Path]) -> tuple:
    """
    Returns (accession, sequence) iterator
    Parameters
    ----------
    fasta_file

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
            if ">" in line:
                if current_key is None:
                    current_key = line.split(">")[1].strip()
                else:
                    yield current_key, current_sequence
                    current_sequence = ""
                    current_key = line.split(">")[1].strip()
            else:
                current_sequence += line.strip()
        yield current_key, current_sequence


def get_sequences_from_fasta(fasta_file: typing.Union[str, Path]) -> dict:
    """
    Returns dict of accession to sequence from fasta file
    Parameters
    ----------
    fasta_file

    Returns
    -------
    {accession:sequence}
    """
    return {
        key: sequence for (key, sequence) in get_sequences_from_fasta_yield(fasta_file)
    }


def get_alignment_subselection(
    alignment: typing.Dict[str, str], keys: typing.List[str], error: bool = True
):
    """
    Obtain a subselection of full sequence alignment

    Parameters
    ----------
    alignment
        dict of key: aligned sequence
    keys
        set of keys to include in the output alignment
    error
        if True, errors if key not in alignment
    Returns
    -------
    dictionary of aligned sequences
    """
    # compile fasta of selected headers
    if error:
        for key in keys:
            assert key in alignment, f"{key} not found in alignment"
    else:
        keys = [k for k in keys if k in alignment]
    alignment_length = len(alignment[keys[0]])

    # find positions where all sequences have a break
    breaks = set()
    for i in range(alignment_length):
        if all(alignment[key][i] == "-" for key in keys):
            breaks.add(i)

    # keep only char's where not all sequences have a break
    sub_alignment = {}
    for key in keys:
        sub_alignment[key] = "".join(
            [char for idx, char in enumerate(alignment[key]) if idx not in breaks]
        )
    return sub_alignment


def alignment_to_fasta(alignment: typing.Dict[str, str]) -> str:
    """
    Convert alignment dict to fasta format

    Parameters
    ----------
    alignment
        dict of key: aligned sequence

    Returns
    -------
    fasta-formatted string
    """
    fasta = []
    for key, sequence in alignment.items():
        fasta.append(f">{key}\n{sequence}")
    return "\n".join(fasta)


def alignment_conservation(alignment: typing.Dict[str, str]) -> typing.List[float]:
    """
    Obtain residue conservations of each position in a sequence alignment

    Parameters
    ----------
    alignment
        dict of key: aligned sequence

    Returns
    -------
    list of conservation values
    """
    conservations = []
    length = len(alignment[list(alignment.keys())[0]])
    for pos in range(length):
        aligned_residues_count = Counter(
            [alignment[key][pos] for key in alignment if alignment[key][pos] != "-"]
        )
        conservation = aligned_residues_count.most_common(1)[0][1] / len(alignment)
        conservations.append(conservation)
    return conservations


def format_as_fasta(alignment):
    fasta = []
    for key in alignment:
        fasta.append(f">{key}\n{alignment[key]}\n")
    return "".join(fasta)
