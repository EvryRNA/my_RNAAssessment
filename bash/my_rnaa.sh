#!/bin/sh

if [ "$5" = "rna" ];then
	echo "Build trainset"
	./../distance_calculation -d $1 -l $2 -R -o Distance_trainset
	echo "\nBuild testset"
	./../distance_calculation -d $3 -l $4 -R -o Distance_testset

	echo "\nPreprocessing trainset"
	python3 ../src/Preprocessing.py -D Distance_trainset.txt -o rna_pp_train -rna -set_train_test -hmap
	echo "\nPreprocessing testset"
	python3 ../src/Preprocessing.py -D Distance_testset.txt -o rna_pp_test -rna -no_mean_struct

	echo "\nMachine learning training and scoring"
	if [ -z $6 ];then
		python3 ../src/ML_scripts/ML_CNN.py -D rna_pp_train.npz -pred rna_pp_test.npz
	else
		python3 ../src/ML_scripts/$6 -D rna_pp_train.npz -pred rna_pp_test.npz
	fi
else
	echo "Build trainset"
	./../distance_calculation -d $1 -l $2 -o Distance_trainset
	echo "\nBuild testset"
	./../distance_calculation -d $3 -l $4 -o Distance_testset

	echo "\nPreprocessing trainset"
	python3 ../src/Preprocessing.py -D Distance_trainset.txt -o prot_pp_train -set_train_test -hmap
	echo "\nPreprocessing testset"
	python3 ../src/Preprocessing.py -D Distance_testset.txt -o prot_pp_test -no_mean_struct

	echo "\nMachine learning training and scoring"
	if [ -z $5 ];then
		python3 ../src/ML_scripts/ML_CNN.py -D prot_pp_train.npz -pred prot_pp_test.npz
	else
		python3 ../src/ML_scripts/$5 -D prot_pp_train.npz -pred prot_pp_test.npz
	fi
fi	
