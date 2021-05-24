from homology_modelling import *
import os
import numpy as np
from sys import argv


if __name__ == "__main__":
    in_dir = argv[1]
    out_dir = argv[2]
    for file in os.listdir(in_dir):
        if file[-4:] == '.pdb':
            file_dir = in_dir + file
            mdl = get_ndope_profile(file_dir)
            mdl = mdl.get_normalized()
            mdl.write_to_file(out_dir + file[:-4] + '_ndope.txt')

    
