import sys
import os
from urllib.request import urlretrieve

def download_pdb(pdbcode, datadir=os.getcwd(), downloadurl="https://files.rcsb.org/download/"): 
    pdb = pdbcode + ".pdb"
    url = downloadurl + pdb
    outfnm = os.path.join(datadir, pdb)
    try:
        urlretrieve(url, outfnm)
        return pdb
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None

pdbsfile = open(sys.argv[2], "r") # List of PDB file to download
pdblist = []
cpt = 0

for line in pdbsfile:
    pdblist.append(line[:4])

pdblist = set(pdblist)   # Avoids duplication

for elem in pdblist:
    download_pdb(elem, sys.argv[1])
    cpt += 1
    print(str(cpt), end = "\r")