#include <iostream>
#include <fstream>
#include <string>
#include <map>      // for dictionary
#include <cstring>  // stof()
#include <vector>
#include <cmath>
#include <getopt.h>  // add options

using namespace std;

void add_coords(vector<vector<string>> &xyz, string &lines, string &residus, string &atoms, float &occup, int &iter, bool nx_step = true){
	string X = lines.substr(30,8);  //
	string Y = lines.substr(38,8);  // Atomic coordinates x, y, z
	string Z = lines.substr(46,8);  //
	string pos = lines.substr(22,4);  // Residue position
	string Chain = lines.substr(21,1);
	try {occup = stof(lines.substr(54,6));}  // For residue alternate location (RAL)
	catch (...){ occup = 1;}                 //
	if (nx_step){xyz.push_back(vector<string>(7));
	iter += 1;}
	xyz[iter][0] = X;        // 
	xyz[iter][1] = Y;        //
	xyz[iter][2] = Z;        // Stock the coordinates, residue name,
	xyz[iter][3] = residus;  // atom and other informations in an
	xyz[iter][4] = atoms;    // intermediate vector
	xyz[iter][5] = pos;      //
	xyz[iter][6] = Chain;    //
}


vector<vector<vector<string> >> sch_coord_pdb(string pdbfile, string chain, bool rna = false){
	vector<vector<vector<string> >> tableau;         // Atomic coordinates ('main' dimension)
	vector<vector<string> > interm_tab;              // Atomic coordinates for current chain
	string line;
	string Atom;
	string previous_letter = " ";  // For alternates locations that don't start with A
	int plc1; int vrf1; int plc2; int vrf2; int plc3; // Atom and residue parameters
	float occupancy;
	string refer;  // Atoms that interest us
	int i = -1;

	if (!rna){
		refer = "CA";  // CA and C4' for RNA
		plc1 = 2; vrf1 = 17; plc2 = 3; vrf2 = 16; plc3 = 4; // Atom and residue parameters for protein
	} else {
		refer = "C4'";
		plc1 = 3; vrf1 = 19; plc2 = 1; vrf2 = 18; plc3 = 2; // Atom and residue parameters for RNA
	}
	
	ifstream fl(pdbfile);
	while(getline(fl, line)){
		if (line.substr(0,4) == "ATOM"){
			if (chain == line.substr(21,1))  // Chain selected
			{
				Atom = line.substr(13,plc1);
				if (refer == Atom)
				{
					string residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
					if ((line.substr(vrf2,plc3) == " "+residu)  || (line.substr(vrf2,plc3) == "A"+residu))
					{
						add_coords(interm_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(vrf2,1);
					}else if (previous_letter == " ")
					{
						add_coords(interm_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(vrf2,1);
					} else { try { if (stof(line.substr(54,6)) > occupancy)  // Compare the occupancy if there is RAL
					{
						add_coords(interm_tab, line, residu, Atom, occupancy, i, false);
					}} catch (...) { add_coords(interm_tab, line, residu, Atom, occupancy, i);}}
				}
			}
		}
	}
	if (!interm_tab.empty())             // Stock atomic coordinates for the last chain
	{                                    // in the main vector
		tableau.push_back(interm_tab);   //
	}
	return tableau;
	}


/**** Store the coordinates of the atoms we are interested in for 1 PDB file ****/
vector<vector<vector<string> >> coord_pdb(string pdbfile, bool rna = false){
	vector<vector<vector<string> >> tableau;         // Atomic coordinates for all chains
	vector<vector<string> > interm_tab;              // Atomic coordinates for current chain
	vector<string> tab_chain = {"A"};               // To know in which chain we are currently (start with "A")
	string line;
	string Atom;
	string previous_letter = " ";  // For alternates locations that don't start with A
	int plc1; int vrf1; int plc2; int vrf2; int plc3; // Atom and residue parameters
	float occupancy;
	string refer;  // Atoms that interest us
	int i = -1;

	if (!rna){
		refer = "CA";  // P and C4' for RNA
		plc1 = 2; vrf1 = 17; plc2 = 3; vrf2 = 16; plc3 = 4; // Atom and residue parameters for protein
	} else {
		refer = "C4'";
		plc1 = 3; vrf1 = 19; plc2 = 1; vrf2 = 18; plc3 = 2; // Atom and residue parameters for RNA
	}
	

	ifstream fl(pdbfile);
	while(getline(fl, line)){
		if (line.substr(0,4) == "ATOM"){
			Atom = line.substr(13,plc1);
			if (tab_chain.back() == line.substr(21,1))  // If still in the same chain
			{
				if (refer == Atom)
				{
					string residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
					if ((line.substr(vrf2,plc3) == " "+residu)  || (line.substr(vrf2,plc3) == "A"+residu))
					{
						add_coords(interm_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(vrf2,1);
					} else if (previous_letter == " ")
					{
						add_coords(interm_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(vrf2,1);
					} else { try { if (stof(line.substr(54,6)) > occupancy)  // Compare the occupancy if there is RAL
					{
						add_coords(interm_tab, line, residu, Atom, occupancy, i, false);
					}} catch (...) { add_coords(interm_tab, line, residu, Atom, occupancy, i);}}
				}
			} else {
				tab_chain.push_back(line.substr(21,1));       // Chain changing
				if (!interm_tab.empty())             // Avoid to stock an empty vector if the
				{                                    // first chain is not "A"
					tableau.push_back(interm_tab);   // Otherwise stock atomic coordinates of
					interm_tab.clear();              // previous chain in the main vector
				}
				i = -1;
				if (refer == Atom)
				{
					string residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
					add_coords(interm_tab, line, residu, Atom, occupancy, i);
					previous_letter = line.substr(vrf2,1);
				}
			}		
		}
	}
	if (!interm_tab.empty())             // Stock atomic coordinates for the last chain
	{                                    // in the main vector
		tableau.push_back(interm_tab);   //
	}
	return tableau;
}


// Distance between 2 atoms (A)
float distance(vector<string> &atom1, vector<string> &atom2){
	return sqrt(pow(stof(atom1[0])-stof(atom2[0]), 2)+pow(stof(atom1[1])-stof(atom2[1]), 2)+pow(stof(atom1[2])-stof(atom2[2]), 2));
}


string ftsround(float num, int deci){  // For round correctly distance values
	if (deci == 0){
		int NUM = round(num);
		return to_string(NUM);
	}
	int precision = pow(10, deci);
	int NUM = round(num*precision);
	string sNUM = to_string(NUM);
	string entier = sNUM.substr(0, sNUM.size()-deci);
	if (entier.empty()){ entier = "0";}
	if (entier == "-"){ entier = "-0";}
	string rnum = entier+"."+sNUM.substr(sNUM.size()-deci, deci);
	return rnum;
}


int main(int argc, char** argv)
{
    string optlist =
		"   Usage:\n"
		"   ./distance_calculation [-d PATHWAY_DATASET] [-l INPUT_LIST] [-o OUTPUT_FILE_NAME] [-m MAX_DIST] [-R] \n\n"
		"   Options:\n"
		"   -d   string   Pathway of the repository where PDB files you interested of are\n"
		"   -l   string   List of all PDB files you want to be processed\n"
		"   -o   string   Name of your file in output (ex.: YourOuputName_pdbcode.txt). The extension '.txt'\n"
		"                 is automatically add to the output file name that you chose\n"
		"   -R            Calculation of the distances for RNA structures (Using atoms C4')\n"
        	"   -m   float      Maximum interatomic distance (Å) (default=15.0)\n"          
		"   -h            Help\n\n";

    string in_dir;    // Pathway of the repository
    string listpdb;   // File containing PDB files list
    string output = "Distances";    // Output file name without extension (".txt",".out",etc...)
    float maxdist = 15.0;
    bool Rna = false;

    int opt;
    while ((opt = getopt(argc,argv, "hRd:l:o:m:")) != EOF){
        switch(opt){
            case 'd': in_dir = optarg; break;
            case 'l': listpdb = optarg; break;
            case 'o': output = optarg; break;
            case 'm': maxdist = stof(optarg); break;
            case 'R': Rna = true; break;
            case 'h': fprintf(stderr, "%s", optlist.c_str()); return 0;
        }
    }

    if (argc == 1){ fprintf(stderr, "%s", optlist.c_str()); return 1; }
    if (listpdb.empty() and output.empty()){ cerr << "Error (argument) : Missing -l and -o arguments\n" << endl; return 1;}
    if (listpdb.empty()){ cerr << "Error (argument) : Missing -l argument\n" << endl; return 1;}
    if (output.empty()){ cerr << "Error (argument) : Missing -o argument\n" << endl; return 1;}

    if (!Rna)
    {
        map<string,string> code3to1;
        code3to1["ALA"]="A"; code3to1["CYS"]="C"; code3to1["ASP"]="D"; code3to1["GLU"]="E"; code3to1["PHE"]="F"; code3to1["GLY"]="G"; code3to1["HIS"]="H"; code3to1["ILE"]="I"; 
        code3to1["LYS"]="K"; code3to1["LEU"]="L"; code3to1["MET"]="M"; code3to1["ASN"]="N"; code3to1["PRO"]="P"; code3to1["GLN"]="Q"; code3to1["ARG"]="R"; code3to1["SER"]="S"; 
        code3to1["THR"]="T"; code3to1["VAL"]="V"; code3to1["TRP"]="W"; code3to1["TYR"]="Y";

    /**** Processing for each file in PDB files list ****/

        ifstream my_pdbs(listpdb);
        string line;

        string ffile = output+".txt";  // Output file name
        ofstream file_out;
        file_out.open(ffile);           // Open a new file
        if (file_out.is_open())
        {
            int cpt = 0; int cptot = 0; bool cutoff = false;
	    while(getline(my_pdbs, line)) // 1 line = 1 PDB file
	    {
            string fl = line;
            vector<vector<vector<string> >> Coords;

            if ((fl.size() == 9) || (fl.size() == 5))
            {
                Coords = sch_coord_pdb(in_dir+fl.substr(0,4)+".pdb", fl.substr(4,1));
            } else if (fl.size() == 4){
                Coords = coord_pdb(in_dir+fl+".pdb");
            } else {
                Coords = coord_pdb(in_dir+fl);}  // 3D vector {Chain[Atom[informations]]}

            if (Coords.empty())
            {
                cerr << "\n" << fl << "\nError: There is no protein sequence in this PDB file" << endl;
                cptot += 1;
            } else {

            for (int k = 0; k < Coords.size(); k++)
            {
                if (Coords[k].size() >= 5)
                {
                    for (int i = 0; i < Coords[k].size()-4; i+=1)
                    {
                        for (int j = i+4; j < Coords[k].size(); j+=1)
                        {
                            float Dist = distance(Coords[k][i], Coords[k][j]);
                            if (Dist < maxdist)
                            {
                                file_out << Dist << "    " << code3to1[Coords[k][i][3]]+code3to1[Coords[k][j][3]] 
                                << "    " << fl << endl;
                            }
                            
                        }
                    }
                } else {
				cutoff = true;}  // Check the length of the sequence, return an error if it is under the cutoff
            }

            if (cutoff) {
                cerr << "\n" << fl << "\nWarning (chain length too short): Protein residues in insufficient number for at least 1 chain" <<endl;
                cpt +=1; cptot += 1;
            } else {
                cpt +=1; cptot += 1;
                cout.flush();
                cout << "\r" << "Processed PDB files :" << cpt;  // To see the evolution of the processing
            }
        }
    }
    cout << "\nTotal processed files : " << cpt << " on " << cptot << " given" << endl; // To see how many file was processed at the end
    file_out.close();}
    return 0;}


    else {   /* If you want to work with RNA sequences */

	/**** Processing for each file in PDB files list for RNA structures ****/

        ifstream my_pdbs(listpdb);
        string line;

        string ffile = output+".txt";  // Output file name
        ofstream file_out;
        file_out.open(ffile);           // Open a new file
        if (file_out.is_open())
        {
            int cpt = 0; int cptot = 0; bool cutoff = false;
	    while(getline(my_pdbs, line)) // 1 line = 1 PDB file
	    {
            string fl = line;
            vector<vector<vector<string> >> Coords;

            if ((fl.size() == 9) || (fl.size() == 5))
            {
                Coords = sch_coord_pdb(in_dir+fl.substr(0,4)+".pdb", fl.substr(4,1), true);
            } else if (fl.size() == 4){
                Coords = coord_pdb(in_dir+fl+".pdb", true);
            } else {
                Coords = coord_pdb(in_dir+fl, true);}  // 3D vector {Chain[Atom[informations]]}

            if (Coords.empty())
            {
                cerr << "\n" << fl << "\nError: There is no RNA sequence in this PDB file" << endl;
                cptot += 1;
            } else {

            for (int k = 0; k < Coords.size(); k++)
            {
                if (Coords[k].size() >= 5)
                {
                    for (int i = 0; i < Coords[k].size()-4; i+=1)
                    {
                        for (int j = i+4; j < Coords[k].size(); j+=1)
                        {
                            float Dist = distance(Coords[k][i], Coords[k][j]);
                            if (Dist < maxdist)
                            {
                                file_out << Dist << "    " << Coords[k][i][3]+Coords[k][j][3] << "    " << fl << endl;
                            }
                        }
                    }
                } else {
                    cutoff = true;}  // Check the length of the sequence, return an error if it is under the cutoff
            }

            if (cutoff) {
                cerr << "\n" << fl << "\nWarning (chain length too short): Protein residues in insufficient number for at least 1 chain" <<endl;
                cptot += 1;
            } else {
                cpt +=1; cptot += 1;
                cout.flush();
                cout << "\r" << "Processed PDB files :" << cpt;  // To see the evolution of the processing
            }
        }
    }
    cout << "\nTotal processed files : " << cpt << " on " << cptot << " given" << endl; // To see how many file was processed at the end
    file_out.close();}
    return 0;}
}
