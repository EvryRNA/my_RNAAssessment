# my_RNAAssessment

## Programs

### Downloader.py

**Dependency : urllib**

Allows the user to download a list of PDB files, if he has one, in the format:
```
1BTKA
2F46A 
4TXRA
```
or at least:
```
1BTK
2F46
4TXR
```
using the command :
```
python3 Downloader.py your_dataset_repository_pathway/ PDB_list.txt
```
### angle_calc.cpp

Since this is a C++ file, it must be compiled first by entering the command:
```
g++ angle_calc.cpp -o angle_calc
```
Than, it can be run without any argument or with the option `-h` to get help.

## How to use

Basically, if you want to have the list of all the torsion angle values from a list of protein PDB files (Psi/Phi), you can simply run the command :
```
./angle_calc -d dataset_pathway/ -l list_of_PDB.txt -o output_file_name
```

To apply this program to a list of RNA PDB files (eta/theta), you can use the previous command line adding to it at the end the option `-R` :
```
./angle_calc -d dataset_pathway/ -l list_RNA_PDB.txt -o output_RNA_file -R
```
You can also access to the pseudotorsion angle values eta' and theta' using options `-Ra` or `-RA` instead of option `-R`.

There are also few more options that we leave it to the user to discover them by himself with the help option `-h`, like `-O` option to add the Omega torsion angle values for proteins, , just to name one. 
