#include <iostream>
#include <fstream>
#include <string>
#include <map>      // for dictionary
#include <cstring>  // stof()
#include <vector>
#include <cmath>
#include <algorithm> // find()
#include <getopt.h>  // add options

using namespace std;

// Protein all atoms

vector<string> set_prot_allatom(){
    vector<string> atypes = {"C  ", "CA ", "CB ", "N  ", "O  ", "SG ", "CG ", "OD1", "OD2", "CD ", "OE1", "OE2", "CD1", "CD2", "CE1", "CE2", "CZ ", "ND1", "NE2", "CG1", "CG2", "CE ", "NZ ", "SD ", "ND2", "NE ", "NH1", "NH2", "OG ", "OG1", "CE3", "CH2", "CZ2", "CZ3", "NE1", "OH "};
    return atypes;
}

// RNA all atoms 

vector<string> set_rna_allatom(){
    vector<string> atypes = {"C1'", "C2'", "C3'", "C4'", "C5'", "C2 ", "C4 ", "C5 ", "C6 ", "C8 ", "N1 ", "N2 ", "N3 ", "N4 ", "N6 ", "N7 ", "N9 ", "O2'", "O3'", "O4'", "O5'", "O2 ", "O4 ", "O6 ", "P  ", "OP1", "OP2"};
    return atypes;
}

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


vector<vector<vector<vector<string> >>> sch_coord_pdb(string pdbfile, string chain, string Atref, bool rna = false){
	vector<vector<vector<vector<string> >>> tableau;         // Atomic coordinates ('main' dimension)
	vector<vector<vector<string> >> interm_tab;              // Atomic coordinates for current chain
	vector<vector<string> > residue_tab;                      // Atomic coordinates for current residue
	string line;
	string Atom;
	string previous_letter = " ";  // For alternates locations that don't start with A
	string space = "";
	int plc1; int vrf1; int plc2;  // Atom and residue parameters
	float occupancy;
	vector<string> Reference;
	vector<string> refer;  // Atoms that interest us
	string residu;
	int i = -1;

	if (!rna){
		if (Atref == "allatom"){ Reference = set_prot_allatom();}   // Choice of the model
		else if (Atref == "CA+CB"){ Reference = {"CA ", "CB "};}    // representation
		else { Reference.push_back(Atref+" ");}  // e.g. CA and C4' (protein/RNA)
		plc1 = 3; vrf1 = 17; plc2 = 3; // Atom and residue parameters for protein
	} else {
		if (Atref == "allatom"){ Reference = set_rna_allatom();}
		else { Reference.push_back(Atref.substr(0,2)+"'");}
		plc1 = 3; vrf1 = 19; plc2 = 1; space = "  "; // Atom and residue parameters for RNA
	}
	
	ifstream fl(pdbfile);
	while(getline(fl, line)){
        refer = Reference;
		if (line.substr(0,4) == "ATOM"){
			if (chain == line.substr(21,1))  // Chain selected
			{
				Atom = line.substr(13,plc1);
				residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
				if ((Atref != "allatom") and (residu == "GLY")){ refer = {"CA "};}
				if (*find(refer.begin(), refer.end(), Atom) == Atom)
				{
					string ver_res = space+residu;  // for RNA alter. residue (ex: 'B  G')
					if ((line.substr(16,4) == " "+ver_res)  || (line.substr(16,4) == "A"+ver_res))
					{
						add_coords(residue_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(16,1);
					}else if (previous_letter == " ")
					{
						add_coords(residue_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(16,1);
					} else { try { if (stof(line.substr(54,6)) > occupancy)  // Compare the occupancy if there is RAL
					{
						add_coords(residue_tab, line, residu, Atom, occupancy, i, false);
					}} catch (...) { add_coords(residue_tab, line, residu, Atom, occupancy, i);}}					
				}
			}
		}
	}

	
	if (!residue_tab.empty()){ 
	vector<vector<string> > join;
	join.push_back(residue_tab[0]);
	for (int i = 1; i < residue_tab.size(); i++)
	{
		if (residue_tab[i][5] == residue_tab[i-1][5]){ join.push_back(residue_tab[i]);}    // 4D vector
		else { interm_tab.push_back(join); join.clear(); join.push_back(residue_tab[i]);}  // {Chain[Residue[Atom[informations]]]}

	}interm_tab.push_back(join);}
	
	if (!interm_tab.empty())             // Stock atomic coordinates for the last chain
	{                                    // in the main vector
		tableau.push_back(interm_tab);   //
	}
	return tableau;
	}


/**** Store the coordinates of the atoms we are interested in for 1 PDB file ****/
vector<vector<vector<vector<string> >>> coord_pdb(string pdbfile, string Atref, bool rna = false){
	vector<vector<vector<vector<string> >>> tableau;         // Atomic coordinates ('main' dimension)
	vector<vector<vector<string> >> interm_tab;              // Atomic coordinates for current chain
	vector<vector<string> > residue_tab;                      // Atomic coordinates for current residue
	vector<string> tab_chain = {"A"};               // To know in which chain we are currently (start with "A")
	string line;
	string Atom;
	string previous_letter = " ";  // For alternates locations that don't start with A
	string space = "";
	int plc1; int vrf1; int plc2;  // Atom and residue parameters
	float occupancy;
	vector<string> Reference;
	vector<string> refer;  // Atoms that interest us
	string residu;
	int i = -1;

	if (!rna){
		if (Atref == "allatom"){ Reference = set_prot_allatom();}   // Choice of the model
		else if (Atref == "CA+CB"){ Reference = {"CA ", "CB "};}    // representation
		else { Reference.push_back(Atref+" ");}  // e.g. CA and C4' (protein/RNA)
		plc1 = 3; vrf1 = 17; plc2 = 3; // Atom and residue parameters for protein
	} else {
		if (Atref == "allatom"){ Reference = set_rna_allatom();}
		else { Reference.push_back(Atref.substr(0,2)+"'");}
		plc1 = 3; vrf1 = 19; plc2 = 1; space = "  "; // Atom and residue parameters for RNA
	}
	

	ifstream fl(pdbfile);
	while(getline(fl, line)){
		refer = Reference;
		if (line.substr(0,4) == "ATOM"){
			Atom = line.substr(13,plc1);
			residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
			if ((Atref != "allatom") and (residu == "GLY")){ refer = {"CA "};}
			if (tab_chain.back() == line.substr(21,1))  // If still in the same chain
			{
				if (*find(refer.begin(), refer.end(), Atom) == Atom)
				{
					string ver_res = space+residu;  // for RNA alter. residue (ex: 'B  G')
					if ((line.substr(16,4) == " "+ver_res)  || (line.substr(16,4) == "A"+ver_res))
					{
						add_coords(residue_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(16,1);
					} else if (previous_letter == " ")
					{
						add_coords(residue_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(16,1);
					} else { try { if (stof(line.substr(54,6)) > occupancy)  // Compare the occupancy if there is RAL
					{
						add_coords(residue_tab, line, residu, Atom, occupancy, i, false);
					}} catch (...) { add_coords(residue_tab, line, residu, Atom, occupancy, i);}}
				}
			} else {
				tab_chain.push_back(line.substr(21,1));       // Chain changing
				
				if (!residue_tab.empty()){ 
				vector<vector<string> > join;
				join.push_back(residue_tab[0]);
				for (int i = 1; i < residue_tab.size(); i++)
				{
					if (residue_tab[i][5] == residue_tab[i-1][5]){ join.push_back(residue_tab[i]);}    // 4D vector
					else { interm_tab.push_back(join); join.clear(); join.push_back(residue_tab[i]);}  // {Chain[Residue[Atom[informations]]]}
					
				}interm_tab.push_back(join); residue_tab.clear();}
				
				if (!interm_tab.empty())             // Avoid to stock an empty vector if the
				{                                    // first chain is not "A"
					tableau.push_back(interm_tab);   // Otherwise stock atomic coordinates of
					interm_tab.clear();              // previous chain in the main vector
				}
				i = -1;
				if (*find(refer.begin(), refer.end(), Atom) == Atom)
				{
					string residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
					add_coords(residue_tab, line, residu, Atom, occupancy, i);
					previous_letter = line.substr(16,1);
				}
			}		
		}
	}
	
	if (!residue_tab.empty()){ 
	vector<vector<string> > join;
	join.push_back(residue_tab[0]);
	for (int i = 1; i < residue_tab.size(); i++)
	{
		if (residue_tab[i][5] == residue_tab[i-1][5]){ join.push_back(residue_tab[i]);}
		else { interm_tab.push_back(join); join.clear(); join.push_back(residue_tab[i]);}
		
	}interm_tab.push_back(join);}
	
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


string noblank(string atom){
	if (atom.substr(1,2) == "  "){ return atom.substr(0,1);}
	else if (atom.substr(2,1) == " "){ return atom.substr(0,2);}
	else { return atom;}
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
        "   ./distance_calculation [-d PATHWAY_DATASET] [-l INPUT_LIST] [-o OUTPUT_FILE_NAME] [-R] [-m MIN_DIST]\n"
        "                          [-M MAX_DIST] [-i MIN_SEQ_SEPARATION] [-j MAX_SEQ_SEPARATION] [-c REPRESENTATION]\n\n"
        "   Options:\n"
        "   -d   string   Pathway of the repository where PDB files you interested of are\n"
        "   -l   string   List of all PDB files you want to be processed\n"
        "   -o   string   Name of your file in output (ex.: YourOuputName_pdbcode.txt). The extension '.txt'\n"
        "                 is automatically add to the output file name that you chose\n"
        "   -R            Calculation of the distances for RNA structures (Using atoms C4')\n\n"
        "   Additional options:\n"
        "   -i   int      Minimum number of positions separating the residue pair (default=4)\n"
        "   -j   int      Maximum number of positions separating the residue pair\n"
        "   -m   float    Minimum interatomic distance (Å) (default=0.0)\n"
        "   -M   float    Maximum interatomic distance (Å) (default=15.0)\n"
        "   -c   string   Carbon of Reference : CA (Calpha), CB (Cbeta), CA+CB, allatom [Protein, default=CA]\n"
        "                                       C4p (C4'), C1p (C1'), allatom [RNA, default=C4p]\n"
        "   -v            Verbosity of the program\n"        
        "   -h            Help\n\n";

    string in_dir;    // Pathway of the repository
    string listpdb;   // File containing PDB files list
    string output = "Distances";    // Output file name without extension (".txt",".out",etc...)
    string carbref = "CA";   // Reference atom
    int Min = 4;   // Minimum positions between 2 residues
    int Max = 6000;   //  Maximum positions between 2 residues (Max >> any chain size)
    float mindist = 0.0;
    float maxdist = 15.0;
    bool Rna = false;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc,argv, "hRvc:d:l:o:i:j:m:M:")) != EOF){
        switch(opt){
            case 'd': in_dir = optarg; break;
            case 'l': listpdb = optarg; break;
            case 'o': output = optarg; break;
            case 'c': carbref = optarg; break;
            case 'i': Min = stoi(optarg); break;
            case 'j': Max = stoi(optarg); break;
            case 'm': mindist = stof(optarg); break;
            case 'M': maxdist = stof(optarg); break;
            case 'R': Rna = true; break;
            case 'v': verbose = true; break;
            case 'h': fprintf(stderr, "%s", optlist.c_str()); return 0;
        }
    }

    vector<string> protref = {"CA","CB","CA+CB","allatom"};   // Command line representation names
    vector<string> rnaref = {"C4p","C1p","allatom"};          // 

    if (argc == 1){ fprintf(stderr, "%s", optlist.c_str()); return 1; }
    if (listpdb.empty() and output.empty()){ cerr << "\nError (argument) : Missing -l and -o arguments\n" << endl; return 1;}
    if (listpdb.empty()){ cerr << "\nError (argument) : Missing -l argument\n" << endl; return 1;}
    if (output.empty()){ cerr << "\nError (argument) : Missing -o argument\n" << endl; return 1;}
    if (Rna)
    {
    	if (carbref == protref[0]) { carbref = "C4p";}
    	else if (carbref == protref[1])
    	{ cerr << "\nError (user) : -c argument (CB) not available for RNA mode \n" << endl; return 1;}
        else if (!(*find(rnaref.begin(), rnaref.end(), carbref) == carbref)){ cerr << "\nError (argument) : Invalid -c argument for RNA mode\n" << endl; return 1;}
    }
    if (!Rna)
    {
    	if (carbref == rnaref[0]) { carbref = "CA"; cerr << "\nError : -c argument (C4p) not available for Protein mode. Default use : CA\n" << endl;}
    	else if (carbref == rnaref[1])
    	{ cerr << "\nError (user) : -c argument (C1p) not available for Protein mode \n" << endl; return 1;}
        else if (!(*find(protref.begin(), protref.end(), carbref) == carbref)){ cerr << "\nError (argument) : Invalid -c argument Protein mode\n" << endl; return 1;}
    }
    if (mindist > maxdist){ cerr << "\nError (user) : -M argument must be greater than -m argument\n" << endl; return 1; }
    if (mindist < 0){ cerr << "\nError (user) : -m argument cannot be negative\n" << endl; return 1; }
    if (maxdist < 0){ cerr << "\nError (user) : -M argument cannot be negative\n" << endl; return 1; }
    if (Min > Max){ cerr << "\nError (user) : -j argument must be greater than -i argument\n" << endl; return 1; }
    if (Min < 0){ cerr << "\nError (user) : -i argument cannot be negative\n" << endl; return 1; }
    if (Max < 0){ cerr << "\nError (user) : -J argument cannot be negative\n" << endl; return 1; }

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
            vector<vector<vector<vector<string> >>> Coords;

            if ((fl.size() == 9) || (fl.size() == 5))
            {
                Coords = sch_coord_pdb(in_dir+fl.substr(0,4)+".pdb", fl.substr(4,1), carbref);
            } else if (fl.size() == 4){
                Coords = coord_pdb(in_dir+fl+".pdb", carbref);
            } else {
                Coords = coord_pdb(in_dir+fl, carbref);}  // 4D vector {Chain[Residue[Atom[informations]]]}


            if (Coords.empty())
            {
                cptot += 1;
                cerr << "\n" << cptot  << " : " << "Error: There is no protein sequence in this PDB file ("+fl+")" << endl;
            } else {

            for (int k = 0; k < Coords.size(); k++)   // Chains
            {
                if (Coords[k].size() >= Min+1)
                {
                    int lenseq;
                    if (Coords[k].size() < Max){
                    	lenseq = Coords[k].size();
                    } else { lenseq = Max;}
                    for (int i = 0; i < lenseq-Min; i+=1)   // Residues
                    {
                        for (int j = i+Min; j < lenseq; j+=1)
                        {
                            for (int l1 = 0; l1 < Coords[k][i].size(); l1++)   // Atoms
                            {
                                for (int l2 = 0; l2 < Coords[k][j].size(); l2++)
                                {
                                    float Dist = distance(Coords[k][i][l1], Coords[k][j][l2]);
                                    if ((Dist > mindist) && (Dist < maxdist))
                                    {
                                        string at1; string at2;
                                        if ((carbref != "CA") && (carbref != "CB")){  // Avoid blank space
                                        at1 = noblank(Coords[k][i][l1][4])+"/"; at2 = noblank(Coords[k][j][l2][4]);}
                                        file_out << Dist << "    " << code3to1[Coords[k][i][l1][3]]+at1+code3to1[Coords[k][j][l2][3]]+at2 << "    " << fl << endl;
                                    }
                                }
                            }
                        }
                    }
                } else {
				cutoff = true;}  // Check the length of the sequence, return an error if it is under the cutoff
            }


            if ((verbose) and (cutoff)) {
                cpt +=1; cptot += 1;
                cerr << "\n" << cpt << " : "  << "Warning (chain length too short): Protein residues in insufficient number for at least 1 chain ("+fl+")" << endl;
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
            vector<vector<vector<vector<string> >>> Coords;

            if ((fl.size() == 9) || (fl.size() == 5))
            {
                Coords = sch_coord_pdb(in_dir+fl.substr(0,4)+".pdb", fl.substr(4,1), carbref, true);
            } else if (fl.size() == 4){
                Coords = coord_pdb(in_dir+fl+".pdb", carbref, true);
            } else {
                Coords = coord_pdb(in_dir+fl, carbref, true);}  // 4D vector {Chain[Residue[Atom[informations]]]}


            if (Coords.empty())
            {
                cptot += 1;
                cerr << "\n" << cptot  << " : " << "Error: There is no RNA sequence in this PDB file ("+fl+")" << endl;
            } else {

            for (int k = 0; k < Coords.size(); k++)   // Chains
            {
                if (Coords[k].size() >= Min+1)
                {
                    int lenseq;
                    if (Coords[k].size() < Max){
                    	lenseq = Coords[k].size();
                    } else { lenseq = Max;}
                    for (int i = 0; i < lenseq-Min; i++)   // Residues
                    {
                        for (int j = i+Min; j < lenseq; j++)
                        {
                            for (int l1 = 0; l1 < Coords[k][i].size(); l1++)   // Atoms
                            {
                                for (int l2 = 0; l2 < Coords[k][j].size(); l2++)
                                {
                                    float Dist = distance(Coords[k][i][l1], Coords[k][j][l2]);
                                    if ((Dist > mindist) && (Dist < maxdist))
                                    {
                                        string at1; string at2;
                                        if (carbref == "allatom"){  // Avoid blank space
                                        at1 = noblank(Coords[k][i][l1][4])+"/"; at2 = noblank(Coords[k][j][l2][4]);}
                                        file_out << Dist << "    " << Coords[k][i][l1][3]+at1+Coords[k][j][l2][3]+at2 << "    " << fl << endl;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    cutoff = true;}  // Check the length of the sequence, return an error if it is under the cutoff
            }


            if ((verbose) and (cutoff)) {
                cpt +=1; cptot += 1;
                cerr << "\n" << cpt << " : "  << "Warning (chain length too short): RNA residues in insufficient number for at least 1 chain ("+fl+")" << endl;
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
