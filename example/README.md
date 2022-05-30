# Example

Guide to learn how to use the different tools present in this repository step by step

## Step 1 - Acquisition of raw data
**For training set :**
```
./distance_calculation -d Training_set/ -l train_rna_set.txt -R -o Distance_trainset
```
**For test set :**
```
./distance_calculation -d 1Z43/ -l test_rna_set.txt -R -o Distance_1Z43
```

These 2 commands will generate 2 differents output files, `Distance_trainset.txt` and `Distance_1Z43.txt`. These files contain all the distances for each pair of residues for each PDB in the dataset repositories (respectively `Training_set/` and `1Z43/`).

## Step 2 - Data preprocessing
**For training set :**
```
python3 Preprocessing.py -D Distance_trainset.txt -o rna_pp_train -rna -set_train_test
```
**For test set :**
```
python3 Preprocessing.py -D Distance_1Z43.txt -o rna_pp_test -rna -no_mean_struct
```


These 2 commands will generate 2 files ".npz" (files of numpy array). 
