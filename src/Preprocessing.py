import argparse, os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

__doc__ = "Preprocessing program for data obtained at the output of the 'distance_calculation' program for Machine Learning methods use"


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
    """Retrieves the arguments of the program."""
    parser = argparse.ArgumentParser(description=__doc__)


    parser.add_argument('-D', dest='distance_data', type=isfile, required=True,
                        help = "File containing interatomic distances (Å)")
    parser.add_argument('-o', dest='output', type=str, default='out_prepro',
                        help = "File name of the output (file '.npz')")
    parser.add_argument('-rna', action = 'store_true',
                        help = "Preprocessing for values which come ​​from an RNA structure dataset")
    parser.add_argument('-max', dest='max_length', type=int, default=15,
                        help = "Maximum distance (default : 15)")
    parser.add_argument('-min', dest='min_length', type=int, default=0,
                        help = "Minimum distance (default : 0)")
    parser.add_argument('-bins', dest='bins', type=float, default=1.0,
                        help = "Interval between values (default : 1.0)")
    parser.add_argument('-set_train_test', action = 'store_true',
                        help = "Get training and test sets in output from the program. Default : only data and their targets")
    parser.add_argument('-no_mean_struct', action = 'store_true',
                        help = "Build a set without mean structures")
    parser.add_argument('-hmap', action = 'store_true',
                        help = "Heatmap of the data")

    return parser.parse_args()


def get_AtomPair(RNA = False):
    prot_atom = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    rna_atom = ["A", "C", "G", "U"]

    atom_pair = []

    if RNA == False:
        for i in range(0, len(prot_atom), 1):
            for j in range(i, len(prot_atom), 1):
                atom_pair.append(prot_atom[i]+prot_atom[j])
    else:
        for i in range(0, len(rna_atom), 1):
            for j in range(i, len(rna_atom), 1):
                atom_pair.append(rna_atom[i]+rna_atom[j])

    return atom_pair


def get_bins(Min = 0, Max = 15, interv = 1):
    atom_pair = []
    j = Min
    while j <= Max-interv:
        atom_pair.append(str(j)+"-"+str(j+interv))
        j += interv
    return atom_pair


def reinit_dict(ele, k):
    if type(ele) is list:
        return [reinit_dict(val, k) for val in ele]
    elif type(ele) is dict:
        return {key : reinit_dict(ele[key], k) for key in ele}
    return k


def check_pair(pair,dico_pair):
	if pair not in dico_pair:
		pair = pair[-1::-1]
	return pair


def get_mean_struct(dico_pair, interv):
    mn_struct = []
    len_dico = len(dico_pair.keys())
    for i in range(0, interv, 1):
        pos = 0
        for key in dico_pair:
            pos += dico_pair[key][i]
        mn_struct.append(round(pos/len_dico))
    mean_struct = [mn_struct]*len_dico
    return mean_struct


if __name__ == '__main__':
    # Get arguments
    args = get_arguments()

    if args.min_length < 0 or args.max_length < 0:
        sys.stderr.write("Error (argument): Values 'bins', 'max' and 'min' must be positive\n")
        sys.exit()
    if args.bins < 0:
        sys.stderr.write("Error (argument): Values 'bins', 'max' and 'min' must be positive\n")
        sys.exit()
    if args.min_length > args.max_length:
        sys.stderr.write("Error (argument): 'min' value must be less than 'max' value !\n")
        sys.exit()

    tot_len = args.max_length - args.min_length
    if tot_len % args.bins != 0:
        sys.stderr.write("Error (bins): 'bins' value must be a divisor of the difference into the minimum distance (min) and the maximum distance (max)\n")
        sys.exit()

    print("Reading file "+args.distance_data+"... ", end="\r")

    if args.rna == True:
        atompair = get_AtomPair(True)
    else:
        atompair = get_AtomPair()

    BINS = int((tot_len) / args.bins)

    dic_pair = {}
    for pair in atompair:
        if pair not in dic_pair:
            dic_pair[pair] = [0]*BINS
    
    # Get file in a list
    file2vec = []
    with open(args.distance_data, "r") as fl:
        for line in fl:
            ln = line.strip("\n")
            file2vec.append(ln.split())

    print("Reading file "+args.distance_data+"... Done")
    print("Distances preprocessing in progress...")

    count_list = []
    y_list = []

    crt_pdb = file2vec[0][-1]

    for line in file2vec:
        if line[-1] == crt_pdb:
            val = float(line[0]) 
            if val > args.max_length or val < args.min_length:
                continue
            if val == args.max_length: 
                value = int((val- args.min_length) / args.bins) - 1
            else:
                value = int((val- args.min_length) / args.bins)
            crt_pair = check_pair(line[-2], dic_pair)
            if crt_pair in dic_pair:
                dic_pair[crt_pair][value] +=1
            else:
                sys.stderr.write("\tWarning : '"+crt_pair+"' is not in the list of residue pairs ("+line[-1]+")\n")
        else:
            altlist = [dic_pair[key] for key in dic_pair]
            count_list.append(altlist)
            y_list.append(0)
            if args.no_mean_struct == False:
                count_list.append(get_mean_struct(dic_pair, BINS))
                y_list.append(1)
            crt_pdb = line[-1]
            dic_pair = reinit_dict(dic_pair, 0)

            val = float(line[0]) 
            if val > args.max_length or val < args.min_length:
                continue
            if val == args.max_length: 
                value = int((val- args.min_length) / args.bins) - 1
            else:
                value = int((val- args.min_length) / args.bins)
            crt_pair = check_pair(line[-2], dic_pair)
            if crt_pair in dic_pair:
                dic_pair[crt_pair][value] +=1
            else:
                sys.stderr.write("Warning : '"+crt_pair+"' is not in the list of residue pairs ("+line[-1]+")\n")

    altlist = [dic_pair[key] for key in dic_pair]
    count_list.append(altlist)
    y_list.append(0)
    if args.no_mean_struct == False:
        count_list.append(get_mean_struct(dic_pair, BINS))
        y_list.append(1)
    
    
    print("Preprocessing done !")
    X = np.asarray(count_list)
    Y = np.asarray(y_list)

    print("Dimension of the data table : "+str(X.shape))

    if args.hmap == True:
        my_bins = get_bins(args.min_length, args.max_length, args.bins)
        
        mytable = []
        T1, T2, T3 = X.shape

        for i in range(0, T2, 1):
          aamean = []
          for j in range(0, T3, 1):
            aaI = []
            for k in range(0, T1, 1):
              if Y[k] != 1:
                aaI.append(X[k][i][j])
            aamean.append(sum(aaI))
          mytable.append(aamean)

        AA_table = pd.DataFrame(mytable, columns=my_bins, index=atompair)
        
        if args.rna != True:
            sns.set(rc={"figure.figsize":(20, 20)})
        else:
            sns.set(rc={"figure.figsize":(15, 5)})
        sns.heatmap(AA_table, center=AA_table.max().max()/2)
        plt.xticks(rotation=45)
        plt.xlabel('Distance (Å)')
        plt.ylabel('Atom pairs')
        plt.title('Heatmap '+args.output, fontsize = 20)
        plt.tight_layout()   # To not have trimmed plot
        plt.savefig('Heatmap_'+args.output+'.png', dpi=300)
        
        print("Heatmap figure saved in file : Heatmap_"+args.output+".png")

    if args.set_train_test == True:
        X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3,random_state=109) # 70% training and 30% test
        np.savez(args.output, x_train = X_train, x_test = X_test, y_train = Y_train, y_test = Y_test)
    else:
        np.savez(args.output, Data = X, Target = Y)
    
    print("Numpy array saved in file : "+args.output+".npz")




