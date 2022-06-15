# my_rnaa.sh

You can execute all the programs (with their basic options), from the raw structures to scored structures (with ML_CNN.py), by using one of the following command lines.

**For proteins**

```
./my_rnaa.sh training_structures_directory/ training_list.txt test_structures_directory/ test_list.txt
```

**For RNA**

Add `rna` as the last argument.
```
./my_rnaa.sh training_structures_directory/ training_list.txt test_structures_directory/ test_list.txt rna
```

:pencil: **Note:** `training_list.txt` and `test_list.txt` are text files containing each one a list of PDB files, which are supposed to be respectively in `training_structures_directory/` and `test_structures_directory/`.
