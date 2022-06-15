# my_rnaa.sh

You can execute all the programs (with their basic options), from the raw structures to scored structures (with ML_CNN.py by default), by using one of the following command lines.

**For proteins**

```
./my_rnaa.sh training_structures_directory/ training_list.txt test_structures_directory/ test_list.txt
```

**For RNA**

Add `rna` as the $5^{th}$ argument.
```
./my_rnaa.sh training_structures_directory/ training_list.txt test_structures_directory/ test_list.txt rna
```

:pencil: **Notes:** 
- `training_list.txt` and `test_list.txt` are text files containing each one a list of PDB files, which are supposed to be respectively in `training_structures_directory/` and `test_structures_directory/`.
- You can also change the script used to score the structures by adding at the end of the command line the name of the script :
```
./my_rnaa.sh training_structures_directory/ training_list.txt test_structures_directory/ test_list.txt ML_SVM.py
```
**OR**
```
./my_rnaa.sh training_structures_directory/ training_list.txt test_structures_directory/ test_list.txt rna ML_SVM.py
```
