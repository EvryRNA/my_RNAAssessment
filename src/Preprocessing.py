import argparse, os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

__doc__ = "Preprocessing program for data obtained at the output of the 'distance_calculation' program or/and the 'angle_calculation' program for Machine Learning methods use"


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


    parser.add_argument('-D', dest='distance_data', type=isfile,
                        help = "File containing interatomic distances (Å)")
    parser.add_argument('-A', dest='angle_data', type=isfile,
                        help = "File containing dihedral angles ([0°;360°])")
    parser.add_argument('-o', dest='output', type=str, default='out_prepro',
                        help = "File name of the output (file '.npz')")
    parser.add_argument('-rna', action = 'store_true',
                        help = "Preprocessing for interatomic distances values which come ​​from an RNA structure dataset")
    parser.add_argument('-repr', choices = ['single', 'CA+CB', 'allatom'], default='single', 
                        help = "Choice of model representation: 'single' and 'allatom' (RNA/protein), 'CA+CB' (only protein). Default: 'single'")  ## Add name for new representation
    parser.add_argument('-max', dest='max_length', type=int, default=15,
                        help = "Maximum distance (default: 15)")
    parser.add_argument('-min', dest='min_length', type=int, default=0,
                        help = "Minimum distance (default: 0)")
    parser.add_argument('-bins', dest='bins', type=float, default=1.0,
                        help = "Interval between distance values (default: 1.0)")
    parser.add_argument('-set_train_test', action = 'store_true',
                        help = "Get training and test sets in output from the program. Default: only data and their targets")
    parser.add_argument('-no_mean_struct', action = 'store_true',
                        help = "Build a set without mean structures")
    parser.add_argument('-hmap', action = 'store_true',
                        help = "Heatmap of the data (for distances only)")
    parser.add_argument('-v', '--verbose', action = 'store_true',
                        help = "Verbosity of the program")

    return parser.parse_args()


def get_AtomPair(RNA = False, representation = 'single'):
    all_prot_atom = ["AC", "ACA", "ACB", "AN", "AO", "CC", "CCA", "CCB", "CN", "CO", "CSG", "DC", "DCA", "DCB", "DCG", "DN", "DO", "DOD1", "DOD2", "EC", "ECA", "ECB", "ECD", "ECG", "EN", "EO", "EOE1", "EOE2", "FC", "FCA", "FCB", "FCD1", "FCD2", "FCE1", "FCE2", "FCG", "FCZ", "FN", "FO", "GC", "GCA", "GN", "GO", "HC", "HCA", "HCB", "HCD2", "HCE1", "HCG", "HN", "HND1", "HNE2", "HO", "IC", "ICA", "ICB", "ICD1", "ICG1", "ICG2", "IN", "IO", "KC", "KCA", "KCB", "KCD", "KCE", "KCG", "KN", "KNZ", "KO", "LC", "LCA", "LCB", "LCD1", "LCD2", "LCG", "LN", "LO", "MC", "MCA", "MCB", "MCE", "MCG", "MN", "MO", "MSD", "NC", "NCA", "NCB", "NCG", "NN", "NND2", "NO", "NOD1", "PC", "PCA", "PCB", "PCD", "PCG", "PN", "PO", "QC", "QCA", "QCB", "QCD", "QCG", "QN", "QNE2", "QO", "QOE1", "RC", "RCA", "RCB", "RCD", "RCG", "RCZ", "RN", "RNE", "RNH1", "RNH2", "RO", "SC", "SCA", "SCB", "SN", "SO", "SOG", "TC", "TCA", "TCB", "TCG2", "TN", "TO", "TOG1", "VC", "VCA", "VCB", "VCG1", "VCG2", "VN", "VO", "WC", "WCA", "WCB", "WCD1", "WCD2", "WCE2", "WCE3", "WCG", "WCH2", "WCZ2", "WCZ3", "WN", "WNE1", "WO", "YC", "YCA", "YCB", "YCD1", "YCD2", "YCE1", "YCE2", "YCG", "YCZ", "YN", "YO", "YOH"]

    all_rna_atom = ["AC1'", "AC2", "AC2'", "AC3'", "AC4", "AC4'", "AC5", "AC5'", "AC6", "AC8", "AN1", "AN3", "AN6", "AN7", "AN9", "AO2'", "AO3'", "AO4'", "AO5'", "AOP1", "AOP2", "AP", "CC1'", "CC2", "CC2'", "CC3'", "CC4", "CC4'", "CC5", "CC5'", "CC6", "CN1", "CN3", "CN4", "CO2", "CO2'", "CO3'", "CO4'", "CO5'", "COP1", "COP2", "CP", "GC1'", "GC2", "GC2'", "GC3'", "GC4", "GC4'", "GC5", "GC5'", "GC6", "GC8", "GN1", "GN2", "GN3", "GN7", "GN9", "GO2'", "GO3'", "GO4'", "GO5'", "GO6", "GOP1", "GOP2", "GP", "UC1'", "UC2", "UC2'", "UC3'", "UC4", "UC4'", "UC5", "UC5'", "UC6", "UC8", "UN1", "UN3", "UO2", "UO2'", "UO3'", "UO4", "UO4'", "UO5'", "UOP1", "UOP2", "UP"]

    #''' First - Add new representation(s) - Create a list with all the reference atoms, if there is more than 1 atom for reference '''#
    ## list_representation = []   # 1_letter_residue + atom_name (e.g. "LCB" for Leucine + Carbon Beta)
    
    AB_prot_atom = ["ACA", "ACB", "CCA", "CCB", "DCA", "DCB", "ECA", "ECB", "FCA", "FCB", "GCA", "HCA", "HCB", "ICA", "ICB", "KCA", "KCB", "LCA","LCB", "MCA", "MCB", "NCA", "NCB", "PCA", "PCB", "QCA", "QCB", "RCA", "RCB", "SCA", "SCB", "TCA","TCB", "VCA", "VCB", "WCA", "WCB", "YCA", "YCB"]
    
    prot_atom = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    rna_atom = ["A", "C", "G", "U"]

    def pair_loop(atom_repr, over2 = False):
        at_pair = []
        for i in range(0, len(atom_repr), 1):
            for j in range(i, len(atom_repr), 1):
                if over2 == True:
                    at_pair.append(atom_repr[i]+"/"+atom_repr[j])
                else:
                    at_pair.append(atom_repr[i]+atom_repr[j])
        return at_pair

    ## Choice of the model representation ##
    if RNA == False:  # Protein
        if representation == 'allatom':
            atom_pair = pair_loop(all_prot_atom, True)
        elif representation == 'CA+CB':
            atom_pair = pair_loop(AB_prot_atom, True)
        #''' Second - Add your representation as an option '''#
        ## elif representation == 'representation_name':  # Name of the option in the command line
        ##     atom_pair = pair_loop(list_representation, True)  # /!\ Always 'True' if representation != 'single'
        else:
            atom_pair = pair_loop(prot_atom)  # representation == 'single'
    else:  # RNA
        if representation == 'allatom':
            atom_pair = pair_loop(all_rna_atom, True)
        ## elif representation == 'representation_name':  # Same process if the representation is for RNA structures
        ##     atom_pair = pair_loop(list_representation, True) 
        else:
            atom_pair = pair_loop(rna_atom)
        
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
    if len(pair) <= 2:    ## e.g. CA => AC
        if pair not in dico_pair:
            pair = pair[-1::-1]
        return pair
    else:                 ## e.g. CCB/ACA => ACA/CCB
        if pair not in dico_pair:
            sep = pair.index("/")
            pair = pair[sep+1:]+"/"+pair[:sep]
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
    
    
def get_Ag_mean_struct(dico_pair, interv, angle_sup):
    if angle_sup == "OMEGA":
        intv1 = int(interv/3)
        intv2 = int(interv/3*2)
        intv3 = interv
    elif angle_sup == "THETA'":
        intv1 = int(interv/4)
        intv2 = int(interv/4*2)
        intv3 = int(interv/4*3)
        intv4 = interv
    else: 
        intv1 = int(interv/2)
        intv2 = interv
        
    mean1 = []
    len_dico = len(dico_pair.keys())
    for i in range(0, intv1, 1):
        pos = 0
        for key in dico_pair:
            pos += dico_pair[key][i]
        mean1.append(round(pos/len_dico))
    
    mean2 = []
    for j in range(intv1, intv2, 1):
        pos = 0
        for key in dico_pair:
            pos += dico_pair[key][j]
        mean2.append(round(pos/len_dico))
    
    if angle_sup == "OMEGA" or angle_sup == "THETA'":
        mean3 = []
        for k in range(intv2, intv3, 1):
            pos = 0
            for key in dico_pair:
                pos += dico_pair[key][k]
            mean3.append(round(pos/len_dico))
            
    if angle_sup == "THETA'":
        mean4 = []
        for l in range(intv3, intv4, 1):
            pos = 0
            for key in dico_pair:
                pos += dico_pair[key][l]
            mean4.append(round(pos/len_dico))
    
    if angle_sup == "OMEGA":
        mean_struct = [mean1 + mean2 + mean3]*len_dico
    elif angle_sup == "THETA'":
        mean_struct = [mean1 + mean2 + mean3 + mean4]*len_dico
    else:
        mean_struct = [mean1 + mean2]*len_dico
    return mean_struct

if __name__ == '__main__':
    # Get arguments
    args = get_arguments()
    
    if not args.distance_data and not args.angle_data:
        sys.stderr.write("Error: You need to have at least 1 file (distances [-D] or dihedral angles [-A]) in input to process\n")
        sys.exit()

    if args.distance_data:

        # Check bin size
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
        
        if args.rna == True and args.repr == 'CA+CB':
            sys.stderr.write("Error (option): No CA/CB representation for RNA structure\n")
            sys.exit()
    
        # General variables
        print("Reading file "+args.distance_data+"... ", end="\r")

        if args.rna == True:
            atompair = get_AtomPair(True, args.repr)
        else:
            atompair = get_AtomPair(representation = args.repr)

        BINS = int((tot_len) / args.bins)

        dic_pair = {}
        for pair in atompair:
            if pair not in dic_pair:
                dic_pair[pair] = [0]*BINS
    

        count_list = []
        y_list = []

        fl = open(args.distance_data, "r")
        lin1 =  next(fl).strip("\n")

        crt_pdb = lin1.split()[-1]  # PDB file
        del fl
    
        print("Reading file "+args.distance_data+"... Done")

    
        if args.rna == True:
            MODE = "RNA"
        else:
            MODE = "Protein"
            
        print("Distances preprocessing in progress ("+MODE+" mode)...")

        with open(args.distance_data, "r") as fl:
            for ln in fl:
                ln1 = ln.strip("\n")
                line = ln1.split()
                if line[-1] == crt_pdb:
                    # Counting phase
                    val = float(line[0]) 
                    if val > args.max_length or val < args.min_length:
                        continue
                    if val == args.max_length: 
                        value = int((val- args.min_length) / args.bins) - 1  #
                    else:                                                    # Adapted value for bins != 1
                        value = int((val- args.min_length) / args.bins)      #
                    crt_pair = check_pair(line[-2], dic_pair)
                    if crt_pair in dic_pair:
                        dic_pair[crt_pair][value] +=1
                    else:
                        if args.verbose == True:
                            sys.stderr.write("\tWarning: '"+crt_pair+"' is not in the list of residue pairs ("+line[-1]+")\n")
                else:
                    # Save data of native and generate non-native with it
                    altlist = [dic_pair[key] for key in dic_pair]
                    count_list.append(altlist)
                    y_list.append(0)
                    if args.no_mean_struct == False:
                        count_list.append(get_mean_struct(dic_pair, BINS))
                        y_list.append(1)
                    crt_pdb = line[-1]
                    dic_pair = reinit_dict(dic_pair, 0)
    
                    # Start next native's counting phase
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
                        if args.verbose == True:
                            sys.stderr.write("Warning: '"+crt_pair+"' is not in the list of residue pairs ("+line[-1]+")\n")

        # For the last native and non-native 
        altlist = [dic_pair[key] for key in dic_pair]
        count_list.append(altlist)
        y_list.append(0)
        if args.no_mean_struct == False:
            count_list.append(get_mean_struct(dic_pair, BINS))
            y_list.append(1)
    
    
        print("Distances preprocessing done !")
        X = np.asarray(count_list)
        Y = np.asarray(y_list)

        print("Dimension of the data table: "+str(X.shape))

        # Heatmap to see the profile of the input data
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
        
            print("Heatmap figure saved in file: Heatmap_"+args.output+".png")

        if args.set_train_test == True:
            X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3,random_state=109) # 70% training and 30% test
            np.savez(args.output, x_train = X_train, x_test = X_test, y_train = Y_train, y_test = Y_test)
        else:
            np.savez(args.output, Data = X, Target = Y)
    
        print("Numpy array saved in file: "+args.output+".npz", end="\n\n")


    ### ANGLES [0°;360°]###
    if args.angle_data:
        
        # General variables           
        print("Reading file "+args.angle_data+"... ", end="\r")
        

        count_list = []
        y_list = []
    
        fl = open(args.angle_data, "r")
        head1 = next(fl).strip("\n")
        header = head1.split()
        
        apair = header.index("PAIR")
        apdb = header.index("PDB_FILE")
        Asuppl = header[2]  # Additional angles
        
        lin1 = next(fl).strip("\n")
        crt_pdb = lin1.split()[apdb]  # PDB file
        del fl
        
        if Asuppl == "OMEGA":
            MODE = "Protein"
            atompair = get_AtomPair()
            dic_pair = {}
            interv = 1080
            for pair in atompair:
                if pair not in dic_pair:
                    dic_pair[pair] = [0]*interv  # 360x3
        elif Asuppl == "THETA'":
            MODE = "RNA"
            atompair = get_AtomPair(True)
            dic_pair = {}
            interv = 1440
            for pair in atompair:
                if pair not in dic_pair:
                    dic_pair[pair] = [0]*interv  # 360x4
        else:
            if header[0] == "THETA" or header[0] == "THETA'":
                MODE = "RNA"
                atompair = get_AtomPair(True)
            else:
                MODE = "Protein"
                atompair = get_AtomPair()
                
            dic_pair = {}
            interv = 720
            for pair in atompair:
                if pair not in dic_pair:
                    dic_pair[pair] = [0]*interv  # 360x2
        
        print("Reading file "+args.angle_data+"... Done")
    
    
        print("Dihedral angles preprocessing in progress ("+MODE+" mode)...")
        
        myfile = open(args.angle_data, "r")
        next(myfile)
       
        for ln in myfile:
            ln1 = ln.strip("\n")
            line = ln1.split()
            if line[apdb] == crt_pdb:  
                # Counting phase
                if line[0] != "NA":
                    val1 = float(line[0]) # angle 1: position 0 to 359
                if line[1] != "NA":
                    val2 = float(line[1]) + 360 # angle 2: position 360 to 719 # 360+360
                if Asuppl == "OMEGA" and line[2] != "NA":
                    val3 = float(line[2]) + 720 # angle 3: position 720 to 1079 # 360+720
                if Asuppl == "THETA'":
                    if line[2] != "NA":
                        val3 = float(line[2]) + 720 # angle 3: position 720 to 1079 # 360+720
                    if line[3] != "NA":
                        val4 = float(line[3]) + 1080 # angle 4: position 1080 to 1440 # 360+1080
                
                if val1 == 360: 
                    val1 = int(val1 - 1)
                else:
                    val1 = int(val1) 
                if val2 >= 720: 
                    val2 = int(val2 - 1)
                else:
                    val2 = int(val2)
                if Asuppl == "OMEGA" or Asuppl == "THETA'":
                    if val3 >= 1080: 
                        val3 = int(val3 - 1)
                    else:
                        val3 = int(val3)
                if Asuppl == "THETA'":
                    if val4 >= 1440: 
                        val4 = int(val4 - 1)
                    else:
                        val4 = int(val4)
                crt_pair = check_pair(line[apair], dic_pair)  
                if crt_pair in dic_pair:
                    dic_pair[crt_pair][val1] +=1
                    dic_pair[crt_pair][val2] +=1
                    if Asuppl == "OMEGA" or Asuppl == "THETA'":
                        dic_pair[crt_pair][val3] +=1
                    if Asuppl == "THETA'":
                        dic_pair[crt_pair][val4] +=1
                else:
                    if args.verbose == True:
                        sys.stderr.write("\tWarning: '"+crt_pair+"' is not in the list of residue pairs ("+line[apdb]+")\n")  
            else:
                # Save data of native and generate non-native with it
                altlist = [dic_pair[key] for key in dic_pair]
                count_list.append(altlist)
                y_list.append(0)
                if args.no_mean_struct == False:
                    count_list.append(get_Ag_mean_struct(dic_pair, interv, Asuppl)) 
                    y_list.append(1)
                crt_pdb = line[-1]
                dic_pair = reinit_dict(dic_pair, 0)
    
                # Start next native's counting phase
                if line[0] != "NA":
                    val1 = float(line[0])
                if line[1] != "NA":
                    val2 = float(line[1]) + 360
                if Asuppl == "OMEGA" and line[2] != "NA":
                    val3 = float(line[2]) + 720
                if Asuppl == "THETA'":
                    if line[2] != "NA":
                        val3 = float(line[2]) + 720
                    if line[3] != "NA":
                        val4 = float(line[3]) + 1080
                
                if val1 == 360: 
                    val1 = int(val1 - 1)  
                else:
                    val1 = int(val1)  
                if val2 >= 720: 
                    val2 = int(val2 - 1)
                else:
                    val2 = int(val2)
                if Asuppl == "OMEGA" or Asuppl == "THETA'":
                    if val3 >= 1080: 
                        val3 = int(val3 - 1)
                    else:
                        val3 = int(val3)
                if Asuppl == "THETA'":
                    if val4 >= 1440: 
                        val4 = int(val4 - 1)
                    else:
                        val4 = int(val4)
                crt_pair = check_pair(line[apair], dic_pair)  
                if crt_pair in dic_pair:
                    dic_pair[crt_pair][val1] +=1
                    dic_pair[crt_pair][val2] +=1
                else:
                    if args.verbose == True:
                        sys.stderr.write("Warning: '"+crt_pair+"' is not in the list of residue pairs ("+line[apdb]+")\n")

        # For the last native and non-native 
        altlist = [dic_pair[key] for key in dic_pair]
        count_list.append(altlist)
        y_list.append(0)
        if args.no_mean_struct == False:
            count_list.append(get_Ag_mean_struct(dic_pair, interv, Asuppl)) 
            y_list.append(1)

        
        print("Dihedral angles preprocessing done !")
        Xa = np.asarray(count_list)
        Ya = np.asarray(y_list)

        print("Dimension of the data table: "+str(Xa.shape))

        if args.set_train_test == True:
            X_train, X_test, Y_train, Y_test = train_test_split(Xa, Ya, test_size=0.3,random_state=109) # 70% training and 30% test
            np.savez(args.output+"_angle", x_train = X_train, x_test = X_test, y_train = Y_train, y_test = Y_test)
        else:
            np.savez(args.output+"_angle", Data = Xa, Target = Ya)
    
        print("Numpy array saved in file: "+args.output+"_angle.npz", end="\n\n")
        
        if args.hmap == True and not args.distance_data:
            sys.stderr.write("Warning: Heatmap only available for interatomic distances\n")



