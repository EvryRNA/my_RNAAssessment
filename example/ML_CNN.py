import argparse, os
import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow import keras
from tensorflow.keras import layers


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

    parser.add_argument('-D', dest='datas', type=isfile,
                        help = "File (*.npz) containing RNA or Protein structures (output of 'Preprocessing.py')")
    parser.add_argument('-load', dest='load', type=isfile,
                        help = "File containing architecture of a CNN model")
    parser.add_argument('-load_weights', dest='load_w', type=isfile,
                        help = "File (*.h5) containing weights of a CNN model")
    parser.add_argument('-save', dest='save_model', type=str, default='CNN_model.json',
                        help = "File where CNN model will be saved (default : CNN_model.json)")
    parser.add_argument('--no_weights', action='store_false',
                        help = "Do not save the weights of the model after training")
    parser.add_argument('--test', action='store_true',
                        help = "Print the accuracy after the training")
    parser.add_argument('-pred', dest='predict', type=isfile,
                        help = "File (*.npz) containing structures to process")
    parser.add_argument('-o', dest='output', type=str, default='cnn_output.txt',
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

    # Load a pretrained model or just the architecture of a model
    if args.load:
        with open(args.load,"r") as json:
            for line in json:
                file_json = line
        model = keras.models.model_from_json(file_json)
        model.summary()
        print("\n")
        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
        if args.load_w:
            model.load_weights(args.load_w)
            if args.test == True:
                scor = model.evaluate(X_test, Y_test)
                print("Training accuracy : {}".format(scor[1]))

    # Train a new model
    else:
        """ Layers/parameters of the CNN model can be changed at any moment by user to optimize learning """ 
        T1, T2, T3 = X_train.shape
        
        inputs = layers.Input(shape=(T2, T3))

        conv1 = layers.Conv1D(15, 3, activation = 'relu')(inputs)
        conv2 = layers.Conv1D(10, 5, activation = 'relu')(conv1)
        flatt = layers.Flatten()(conv2)
        dense1 = layers.Dense(30)(flatt)
        dense2 = layers.Dense(15)(dense1)
        dense3 = layers.Dense(7)(dense2)
        outputs = layers.Dense(1, activation='sigmoid')(dense3)

        model = keras.Model(inputs, outputs)

        print(model.summary())
        print("\n")

        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

        # Save the model architecture in a json file
        if ".json" not in args.save_model:
            name = args.save_model+".json"
        else:
            name = args.save_model
        j_son = model.to_json()
        with open(name,"w") as json:
            json.write(j_son)
        print("Model architecture save in file : "+name)

    if args.load_w == None:  # No weights loaded --> new training with an existing architecture
        if args.test == False:    # Training with 100% of the dataset
            X_train = np.concatenate((X_train, X_test))
            Y_train = np.concatenate((Y_train, Y_test))

            history = model.fit(X_train, Y_train, validation_split= 0.3, epochs= 10)
            print("\n")
        else:                     # Training with 70%, check accuracy with 30%
            history = model.fit(X_train, Y_train, validation_split= 0.3, epochs= 10)
            print("\n")
            scor = model.evaluate(X_test, Y_test)
            print("Training accuracy : {}".format(scor[1]))
        
    # Save the model weights after training in H5 file
    if Load_w == None and args.no_weights != False:
        if ".json" in args.save_model:
            nameW = args.save_model[0:-5]+".h5"
        else:
            nameW = args.save_model+".h5"
        model.save_weights(nameW)
        print("Model weights save in file : "+nameW)     

    if args.predict:
        arr_pred = np.load(args.predict)
        my_data = arr_pred["Data"]

        y_pred = model.predict(my_data)
        with open(args.output, 'w') as out:
            for i in range(0,len(y_pred),1):
                out.write(str(y_pred[i])+'\n')

        print("Scores saved in file : "+args.output)




