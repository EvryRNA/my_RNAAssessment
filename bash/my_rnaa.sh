#!/bin/sh

if [ $5 = "rna" ];then
	echo "Build trainset"
	./distance_calculation -d $1 -l $2 -R -o Distance_trainset
	echo "\nBuild testset"
	./distance_calculation -d $3 -l $4 -R -o Distance_testset

	echo "\nPreprocessing trainset"
	python3 Preprocessing.py -D Distance_trainset.txt -o rna_pp_train -rna -set_train_test -hmap
	echo "\nPreprocessing testset"
	python3 Preprocessing.py -D Distance_testset.txt -o rna_pp_test -rna -no_mean_struct

	echo "\nMachine learning training and scoring"
	python3 ML_CNN.py -D rna_pp_train.npz -pred rna_pp_test.npz
else
	echo "Build trainset"
	./distance_calculation -d $1 -l $2 -o Distance_trainset
	echo "\nBuild testset"
	./distance_calculation -d $3 -l $4 -o Distance_testset

	echo "\nPreprocessing trainset"
	python3 Preprocessing.py -D Distance_trainset.txt -o prot_pp_train -set_train_test -hmap
	echo "\nPreprocessing testset"
	python3 Preprocessing.py -D Distance_testset.txt -o prot_pp_test -no_mean_struct

	echo "\nMachine learning training and scoring"
	python3 ML_CNN.py -D prot_pp_train.npz -pred prot_pp_test.npz
fi	
