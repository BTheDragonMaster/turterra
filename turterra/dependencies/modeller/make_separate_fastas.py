#!/usr/bin/env python

from sys import argv

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
        out_file = out_folder + seq_ID + ".faa"
        out = open(out_file, 'w')
        out.write(">%s\n%s\n" % (seq_ID, fasta_dict[ID]))
        out.close()


if __name__ == "__main__":
    separate_fastas(argv[1], argv[2])
        
