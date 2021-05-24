from sys import argv
from pathlib import Path, PurePath
import os


import prody as pd

from turterra.dependencies.modeller.homology_modelling import align_templates_to_sequence, convert_fasta_to_ali


def read_fasta(fasta_dir):
    with open(fasta_dir, 'r') as fasta_file:
        fasta_dict = {}
        sequence = []
        for line in fasta_file:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    fasta_dict[ID] = ''.join(sequence)

                    ID = line[1:]
                    sequence = []
                else:
                    ID = line[1:]

            else:
                sequence.append(line)

        fasta_dict[ID] = ''.join(sequence)
        fasta_file.close()
    return fasta_dict


def separate_fastas(input_fasta, out_folder):
    fasta_dict = read_fasta(input_fasta)
    for ID in fasta_dict:
        seq_ID = ID.split()[0]
        seq_ID = seq_ID.replace('|', '_')
        out_file = os.path.join(out_folder, seq_ID + ".faa")
        out = open(out_file, 'w')
        out.write(">%s\n%s\n" % (seq_ID, fasta_dict[ID]))
        out.close()


def get_sequence_from_pdb(templates_folder, out_dir):

    sequences = {}
    for pdb in os.listdir(templates_folder):

        if pdb[-4:] == '.pdb':
            directory = templates_folder + pdb
            sequence = pd.parsePDB(directory).select("name CA and not hetatm").getSequence()
            sequences[pdb[:-4]] = sequence

    with open(out_dir, 'w') as out_file:

        for ID, sequence in sequences.items():
            out_file.write(f">{ID}\n{sequence}\n")

def make_key_file(fasta, key_file_dir):
    fasta_dict = read_fasta(fasta)
    with open(key_file_dir, 'w') as key_file:
        for id, sequence in fasta_dict.items():
            seq_id = id.split()[0]
            seq_id = seq_id.replace('|', '_')

            key_file.write(f'{seq_id}\n')



def make_modeller_files(sequence_file, templates_folder, alignment_folder, fasta_folder, data_folder):
    separate_fastas(sequence_file, fasta_folder)
    templates_name = PurePath(templates_folder).name
    templates_fasta = os.path.join(data_folder, f'{templates_name}.faa')
    get_sequence_from_pdb(templates_folder, templates_fasta)
    key_file = os.path.join(data_folder, 'key_file.txt')
    make_key_file(sequence_file, key_file)

    for sequence in os.listdir(fasta_folder):
        if sequence[-4:] == '.faa':
            sequence_name = sequence[:-4]
            sequence_file = os.path.join(fasta_folder, sequence)

            aligned_file = align_templates_to_sequence(sequence_file, templates_fasta, alignment_folder)
            convert_fasta_to_ali(sequence_name, aligned_file)

    return templates_fasta, key_file

