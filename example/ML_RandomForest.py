import argparse, os
from joblib import dump as job_dump
from joblib import load as job_load
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split


def isfile(path):
    """Check if path is an existing file"""
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program"""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-D', dest='datas', type=isfile, required=True,
                        help = "File (*.npz) containing RNA or Protein structures (output of 'Preprocessing.py')")
    parser.add_argument('-load', dest='load', type=isfile,
                        help = "File containing pretrained Random Forest model")
    parser.add_argument('-save', dest='save_model', type=str, default='rf_model.job',
                        help = "File where Random Forest model will be saved (default : rf_model.job)")
    parser.add_argument('--test', action='store_true',
                        help = "Print the accuracy after the training")
    parser.add_argument('-pred', dest='predict', type=isfile,
                        help = "File (*.npz) containing structures to process for scoring")
    parser.add_argument('-o', dest='output', type=str, default='rf_output.txt',
                        help = "Output file obtained after scoring")

    return parser.parse_args()


if __name__ == '__main__':
    # Get arguments
    args = get_arguments()

    arr_data = np.load(args.datas)

    # Load NPZ file data
    if len(arr_data.files) == 4: 
        X_train = arr_data["x_train"]
        Y_train = arr_data["y_train"]
        X_test = arr_data["x_test"]
        Y_test = arr_data["y_test"]
    else:
        Data = arr_data["Data"]
        Target = arr_data["Target"]
        X_train, X_test, Y_train, Y_test = train_test_split(Data, Target, test_size=0.3,random_state=109)

    # Reshape data to be correctly used by Random Forest model
    T1, T2, T3 = X_train.shape
    t1, t2, t3 = X_test.shape
    xtrain = X_train.reshape(T1, T2*T3)
    xtest = X_test.reshape(t1, t2*t3)

    # Load a pretrained model
    if args.load:
        clf = job_load(args.load)
        if args.test == True:
            scor = clf.score(xtest, Y_test)
            print("Training accuracy : {}".format(scor))

    # Train a new model
    else:
        """ Parameters of the Random Forest object can be changed at any moment by user to optimize learning """ 
        clf = RandomForestClassifier()  # Initialize a Random Forest classifier

        if args.test == False:    # Training with 100% of the dataset
            X_train = np.concatenate((xtrain, xtest))
            Y_train = np.concatenate((Y_train, Y_test))

            clf.fit(X_train, Y_train)
        else:                     # Training with 70%, check accuracy with 30%
            clf.fit(xtrain, Y_train)
            scor = clf.score(xtest, Y_test)
            print("Training accuracy : {}".format(scor))

        # Save the model in a string file
        job_dump(clf, args.save_model)
        print("Model save in file : "+args.save_model)

    if args.predict:
        arr_pred = np.load(args.predict)
        data_pred = arr_pred["Data"]  # ajouter une condition si pas 'Data'

        d1, d2, d3 = data_pred.shape
        my_data = data_pred.reshape(d1, d2*d3)

        y_pred = clf.predict_proba(my_data)
        with open(args.output, 'w') as out:
            for i in range(0,len(y_pred),1):
                out.write(str(y_pred[i][1])+'\n')

        print("Scores saved in file : "+args.output)





