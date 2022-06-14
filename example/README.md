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

These 2 commands will generate 2 files ".npz" (files of numpy array), containing the data tables for the training set and the test set.

:warning: **Important :** For the example, we used datasets of RNA structures. If you want to do the same with protein structures, you just need to erase the flags `-R` and `-rna` for the 2 programs, `./distance_calculation` and `Preprocessing.py` respectively.

It is also possible to see the profile of your data, by plotting the heatmap, using the option `-hmap`. This will generate in output a PNG file containing the heatmap of your data, in addition to the file '.npz'. Tipically, you can have this using :

**Input**
```
python3 Preprocessing.py -D Distance_trainset.txt -o rna_pp_train -rna -set_train_test -hmap
```
**Output:**
![heatmap](https://github.com/FranGASTRIN/my_RNAAssessment/blob/main/img/Heatmap_rna_pp_train.png)

## Step 3 - Machine Learning scripts

### Training and prediction

Basically, the 3 programs for Machine Learning (`ML_SVM.py`, `ML_RandomForest.py` and `ML_CNN.py`) can be use the same way. To train a model and directly assign scores using this model, you can use a command like:
```
python3 ML_SVM.py -D rna_pp_train.npz -pred rna_pp_test.npz
```
This will generate in output a file like `svm_output.txt`, containing the score for each structure in the test set.

### Save and Load trained models

- For `ML_SVM.py` and `ML_RandomForest.py`, which use Scikit-learn package, assign scores by using trained models can be done using 2 command lines.

**Save :**
```
python3 ML_SVM.py -D rna_pp_train.npz -save svm_trained_model.job
```
This will save your trained model in the file `svm_trained_model.job`. By default, the model will always be saved in a file like `svm_model.job`.

**Load :**
```
python3 ML_SVM.py -load svm_trained_model.job -pred rna_pp_test.npz
```
The trained model will be load, giving to you the possibility to assign scores with this model.

- For `ML_CNN.py`, which use Tensorflow/Keras package, the way to save a model will be the same as above, but it give you 1 more option for the loading part. By default, the architecture of a model and ist weights after training are save in files `CNN_model.json` and `CNN_model.h5`.

**Loading with weights (pre-trained model) :**
```
python3 ML_CNN.py -load CNN_model.json -load_weights CNN_model.h5 -pred rna_pp_test.npz
```

Option `-load_weights` will load the weights from a H5 file for the corresponding model's architecture from a JSON file. Than, this pre-trained model can be used to assign scores.

**Loading without weights (train a model) :**

If you remove option `-load_weights`, you will load the model's architecture (JSON file) and start a new training, using this architecture, before assign scores.
```
python3 ML_CNN.py -D rna_pp_train.npz -load CNN_model.json -pred rna_pp_test.npz
```
## More informations
- To customize settings for training models with the Scikit-learn package ([SVM](https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html) and [Random Forest](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)), you can go on the web page of the 2 tools to have all the information on them.
- To have information about the [CNN](https://keras.io/api/layers/convolution_layers/convolution1d/) used in the script, or more, you can go on the website of Keras [here](https://keras.io/api/).
