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
:warning: `your_dataset_repository_pathway/` must be an existing repository !
### Install angle_calculation & distance_calculation

Since they are C++ files, they must be compiled. Enter the following command, when you are in `my_RNAAssessment/`, to do it:
```
make
```
Than, run them without any argument or with the option `-h` to get their help.

## How to use
### angle_calculation

Basically, if you want to have the list of all the torsion angle values from a list of protein PDB files (Psi/Phi), you can simply run the command :
```
./angle_calculation -d dataset_pathway/ -l list_of_PDB.txt -o output_file_name
```

To apply this program to a list of RNA PDB files (eta/theta), you can use the previous command line adding to it at the end the option `-R` :
```
./angle_calculation -d dataset_pathway/ -l list_RNA_PDB.txt -o output_RNA_file -R
```
You can also access to the pseudotorsion angle values **eta'** and **theta'** using options `-Ra` or `-RA` instead of option `-R`.

There are also few more options that we leave it to the user to discover them by himself with the help option `-h`, like `-O` option to add the Omega torsion angle values for proteins, just to name one of them. 
