import prody as pd
import os
from sys import argv

if __name__ == "__main__":
    sequences = {}
    for pdb in os.listdir(argv[1]):
        
        if pdb[-4:] == '.pdb':
            directory = argv[1] + pdb
            sequence = pd.parsePDB(directory).select("name CA and not hetatm").getSequence()
            sequences[pdb[:-4]] = sequence


    out_file = open(argv[2], 'w')
    for ID in sequences:
        out_file.write(">%s\n%s\n" % (ID, sequences[ID]))

    out_file.close()
            
    
