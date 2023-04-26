#include <iostream>
#include <fstream>
#include <string>
#include <map>      // for dictionary
#include <cstring>  // stof()
#include <vector>
#include <cmath>
#include <algorithm> // find()
#include <filesystem>
#include <getopt.h>  // add options

namespace fs = std::filesystem;
using namespace std;

vector<string> get_list_dir(string dir_path) {
    vector<string> vec;
    if (!fs::is_directory(dir_path)){
        if (fs::is_regular_file(dir_path)) {
            vec.push_back(dir_path) ;
            return vec ;
        }
    }
    for (const auto& entry : fs::directory_iterator(dir_path))
    {
        if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb")
        {
            string f_name = entry.path().filename();
            string all_path = dir_path + '/' + f_name;
            vec.push_back(all_path);
        }
    }
    return vec;
}
string removeWhitespace(string str) {
    string result = "";
    for (int i = 0; i < str.length(); i++) {
        if (str[i] != ' ') {
            result += str[i];
        }
    }
    return result;
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


vector<vector<vector<string> >> sch_coord_pdb(string pdbfile, string chain, bool rna = false){
	vector<vector<vector<string> >> tableau;         // Atomic coordinates ('main' dimension)
	vector<vector<string> > interm_tab;              // Atomic coordinates for current chain
	vector<string> tab_atom;
	string line;
	string Atom;
	string previous_letter = " ";  // For alternates locations that don't start with A
	string space = "";
	int plc1; int vrf1; int plc2;  // Atom and residue parameters
	float occupancy;
	vector<string> refer;  // Atoms that interest us
	int i = -1;

	if (!rna){
		refer = {"N ","CA","C "};  // P and C4' for RNA
		plc1 = 2; vrf1 = 17; plc2 = 3;  // Atom and residue parameters for protein
	} else {
		refer = {"P  ","C4'","C1'"};
		plc1 = 3; vrf1 = 19; plc2 = 1; space = "  ";  // Atom and residue parameters for RNA 18 2
	}
	
	ifstream fl(pdbfile);
	while(getline(fl, line)){
		if (line.substr(0,4) == "ATOM"){
			if (chain == line.substr(21,1))  // Chain selected
			{
				Atom = line.substr(13,plc1);
				tab_atom.push_back(Atom);               //
				if (tab_atom.front() != refer.front())  // Avoids possible mismatches if the N-ter part
				{                                       // is missing from the PDB file 
					tab_atom.pop_back();                //
				} else {
					if (*find(refer.begin(), refer.end(), Atom) == Atom)
						{
							string residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
							string ver_res = space+residu;  // for RNA alter. residue (ex: 'B  G')
							if ((line.substr(16,4) == " "+ver_res)  || (line.substr(16,4) == "A"+ver_res))
							{
								add_coords(interm_tab, line, residu, Atom, occupancy, i);
								previous_letter = line.substr(16,1);
							}else if (previous_letter == " ")
							{
								add_coords(interm_tab, line, residu, Atom, occupancy, i);
								previous_letter = line.substr(16,1);
							} else { try { if (stof(line.substr(54,6)) > occupancy)  // Compare the occupancy if there is RAL
							{
								add_coords(interm_tab, line, residu, Atom, occupancy, i, false);
							}} catch (...) { add_coords(interm_tab, line, residu, Atom, occupancy, i);}}
						}
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
	vector<string> tab_atom;
	vector<string> tab_chain = {"A"};               // To know in which chain we are currently (start with "A")
	string line;
	string Atom;
	string previous_letter = " ";  // For alternates locations that don't start with A
	string space = "";
	int plc1; int vrf1; int plc2;  // Atom and residue parameters
	float occupancy;
	vector<string> refer;  // Atoms that interest us
	int i = -1;

	if (!rna){
		refer = {"N ","CA","C "};  // P and C4' for RNA
		plc1 = 2; vrf1 = 17; plc2 = 3;  // Atom and residue parameters for protein
	} else {
		refer = {"P  ","C4'","C1'"};
		plc1 = 3; vrf1 = 19; plc2 = 1; space = "  ";  // Atom and residue parameters for RNA
	}
	

	ifstream fl(pdbfile);
	while(getline(fl, line)){
		if (line.substr(0,4) == "ATOM"){
			Atom = line.substr(13,plc1);
			tab_atom.push_back(Atom);               //
			if (tab_atom.front() != refer.front())  // Avoids possible mismatches if the N-ter part
			{                                       // is missing from the PDB file 
				tab_atom.pop_back();            //
			} else {
				if (tab_chain.back() == line.substr(21,1))  // If still in the same chain
				{
					if (*find(refer.begin(), refer.end(), Atom) == Atom)
					{
						string residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
						string ver_res = space+residu;  // for RNA alter. residue (ex: 'B  G')
						if ((line.substr(16,4) == " "+ver_res)  || (line.substr(16,4) == "A"+ver_res))
						{
							add_coords(interm_tab, line, residu, Atom, occupancy, i);
							previous_letter = line.substr(16,1);
						} else if (previous_letter == " ")
						{
							add_coords(interm_tab, line, residu, Atom, occupancy, i);
							previous_letter = line.substr(16,1);
						} else { try { if (stof(line.substr(54,6)) > occupancy)  // Compare the occupancy if there is RAL
						{
							add_coords(interm_tab, line, residu, Atom, occupancy, i, false);
						}} catch (...) { add_coords(interm_tab, line, residu, Atom, occupancy, i);}}
					}
				}else {
					tab_chain.push_back(line.substr(21,1));       // Chain changing
					if (!interm_tab.empty())             // Avoid to stock an empty vector if the
					{                                    // first chain is not "A"
						tableau.push_back(interm_tab);   // Otherwise stock atomic coordinates of
						interm_tab.clear();              // previous chain in the main vector
					}
					i = -1;
					if (refer.front() == Atom)
					{
						string residu = line.substr(vrf1,plc2);  // Residue name (3 letters -> protein ; 1 letter -> RNA)
						add_coords(interm_tab, line, residu, Atom, occupancy, i);
						previous_letter = line.substr(16,1);
					} else {
						tab_atom.clear();
					}
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
float distance(vector<float> atom1, vector<float> atom2){
	return sqrt(pow(atom1[0]-atom2[0], 2)+pow(atom1[1]-atom2[1], 2)+pow(atom1[2]-atom2[2], 2));
}


// Norme of a vector
float norme(vector<float> vecteur){
	vector<float> Vec0 = {0.0, 0.0, 0.0};
	float D = distance(vecteur, Vec0);
	return D;
}


// Return the vector product of 2 vectors
vector<float> vector_prod(vector<float> vect1, vector<float> vect2){
	vector<float> product;
	product.push_back(vect1[1]*vect2[2]-vect1[2]*vect2[1]);
	product.push_back(vect1[2]*vect2[0]-vect1[0]*vect2[2]);
	product.push_back(vect1[0]*vect2[1]-vect1[1]*vect2[0]);
	return product;
}


// Return the scalar product of 2 vectors
float scalar_prod(vector<float> vect1, vector<float> vect2){
	float product = vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2];
	return product;
}


// Return the torsion angle between atom2 and atom3
float torsion_angle(vector<string> atom1, vector<string> atom2, vector<string> atom3, vector<string> atom4, bool to360 = false){
	vector<float> vecteur12 = {stof(atom1[0])-stof(atom2[0]),stof(atom1[1])-stof(atom2[1]),stof(atom1[2])-stof(atom2[2])};  //
    vector<float> vecteur23 = {stof(atom2[0])-stof(atom3[0]),stof(atom2[1])-stof(atom3[1]),stof(atom2[2])-stof(atom3[2])};  //  stof() : convert string to float
    vector<float> vecteur34 = {stof(atom3[0])-stof(atom4[0]),stof(atom3[1])-stof(atom4[1]),stof(atom3[2])-stof(atom4[2])};  //

    vector<float> vecteur_normal1 = vector_prod(vecteur12, vecteur23);
    vector<float> vecteur_normal2 = vector_prod(vecteur23, vecteur34);

    if (scalar_prod(vecteur23,vector_prod(vecteur_normal1,vecteur_normal2))<0)  // Know if sign of the value will be '+' or '-'
    {
	    float angle = acos(scalar_prod(vecteur_normal1,vecteur_normal2)/
	    	          (norme(vecteur_normal1)*norme(vecteur_normal2)))*180/M_PI; // Formula to calculate a torsion angle ]0;180°]
	    return angle;
    } else {
    	float angle = -acos(scalar_prod(vecteur_normal1,vecteur_normal2)/
	    	          (norme(vecteur_normal1)*norme(vecteur_normal2)))*180/M_PI; // [-180°;0[
    	if (to360){ return angle+360;} else { return angle;}  // For ]180°;360°]	
    }
}


string ftsround(float num, int deci){  // For round correctly angle values
	try {
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
	if (rnum.size() > 5+deci){                   // Control
		string sNUM = "180.";
		if (deci == 0){ sNUM = "180";} else {
		for (int i = 0; i < deci; ++i)
		{
			sNUM += "0";
		}}
		return sNUM;}
	else {return rnum;}} catch (...){             // If num is a NaN (because of the acos() in torsion_angle()) 
		string sNUM = "180.";
		if (deci == 0){ sNUM = "180";} else {
		for (int i = 0; i < deci; ++i)
		{
			sNUM += "0";
		}}
		return sNUM;
	}
}


void get_theta_eta(vector<vector<vector<string> >> &pdbcoord, string &order, string &theta, string &eta, 
                   int &k, int &i, int &j, bool &pdbmistake,bool &to360, int &deci, bool adjust = true){
	if (order == "P  C4'P  C4'"){
		theta = ftsround(torsion_angle(pdbcoord[k][i], pdbcoord[k][i+1], pdbcoord[k][i+3], pdbcoord[k][i+4], to360), deci);    // ATOMS : P-C4'-P-C4'
		if (pdbcoord[k][i+6][4] == "P  "){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+2], pdbcoord[k][j+3], pdbcoord[k][j+5], to360), deci);}    // ATOMS : C4'-P-C4'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+2], pdbcoord[k][j+3], pdbcoord[k][j+7], to360), deci);}
		else { eta = "NA";}}
	else if (order == "P  C4'C1'C4'"){
		theta = ftsround(torsion_angle(pdbcoord[k][i], pdbcoord[k][i+1], pdbcoord[k][i+5], pdbcoord[k][i+4], to360), deci);  // ATOMS : P-C4'-C4'-P --> P-C4'-P-C4'
		if (pdbcoord[k][i+6][4] == "P  "){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+4], pdbcoord[k][j+3], pdbcoord[k][j+5], to360), deci);}    // ATOMS : C4'-C4'-P-P --> C4'-P-C4'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+4], pdbcoord[k][j+3], pdbcoord[k][j+7], to360), deci);}
		else { eta = "NA";}}
	else if (order == "C1'C4'P  C4'"){
		theta = ftsround(torsion_angle(pdbcoord[k][i+2], pdbcoord[k][i+1], pdbcoord[k][i+3], pdbcoord[k][i+4], to360), deci);  // ATOMS : C4'-P-P-C4' --> P-C4'-P-C4'
		if (pdbcoord[k][i+6][4] == "P  "){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+2], pdbcoord[k][j+3], pdbcoord[k][j+5], to360), deci);}  // ATOMS : C4'-(P-)P-C4'-P--> (P-)C4'-P-C4'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+2], pdbcoord[k][j+3], pdbcoord[k][j+7], to360), deci);}
		else { eta = "NA";}}
	else if (order == "C1'C4'C1'C4'"){
		theta = ftsround(torsion_angle(pdbcoord[k][i+2], pdbcoord[k][i+1], pdbcoord[k][i+5], pdbcoord[k][i+4], to360), deci);  // ATOMS : C4'-P-C4'-P --> P-C4'-P-C4'
		if (pdbcoord[k][i+6][4] == "P  "){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+4], pdbcoord[k][j+3], pdbcoord[k][j+5], to360), deci);} // ATOMS : C4'-(P-)C4'-P-P --> (P-)C4'-P-C4'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			eta = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+4], pdbcoord[k][j+3], pdbcoord[k][j+7], to360), deci);}
		else { eta = "NA";}}
	else {
		theta = "NA";  // Returns NA if the atoms in the backbone are 
		eta = "NA";    // not well referenced for 1 pair of residue
		pdbmistake = true;
		if (adjust){
		if (pdbcoord[k][i+1][4] == "P  ")        //
		{                                        //
			i -= 2;                              //
		} else if (pdbcoord[k][i+2][4] == "P  ") // Try to find another pair of residue
		{                                        // with all their backbone atoms for the
			i -= 1;                              // next step
		} else {                                 //
		while(pdbcoord[k][i+3][4] != "P  "){     //
			i += 1;}                             //
		}
		}
	}
}


void get_thetaP_etaP(vector<vector<vector<string> >> &pdbcoord, string &order, string &thetaP, string &etaP, 
                   int &k, int &i, int &j, bool &pdbmistake,bool &to360, int &deci){
	if (order == "P  C1'P  C1'"){
		thetaP = ftsround(torsion_angle(pdbcoord[k][i], pdbcoord[k][i+2], pdbcoord[k][i+3], pdbcoord[k][i+5], to360), deci);    // ATOMS : P-C1'-P-C1'
		if (pdbcoord[k][i+6][4] == "P  "){
			etaP = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+1], pdbcoord[k][j+3], pdbcoord[k][j+4], to360), deci);}    // ATOMS : C1'-P-C1'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			etaP = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+1], pdbcoord[k][j+3], pdbcoord[k][j+6], to360), deci);}
		else { etaP = "NA";}}
	else if (order == "P  C1'C1'P  "){
		thetaP = ftsround(torsion_angle(pdbcoord[k][i], pdbcoord[k][i+2], pdbcoord[k][i+5], pdbcoord[k][i+3], to360), deci);  // ATOMS : P-C1'-C1'-P --> P-C1'-P-C1'
		if (pdbcoord[k][i+6][4] == "P  "){
			etaP = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+3], pdbcoord[k][j+1], pdbcoord[k][j+4], to360), deci);}    // ATOMS : C1'-C1'-P-P --> C1'-P-C1'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			etaP = ftsround(torsion_angle(pdbcoord[k][j], pdbcoord[k][j+3], pdbcoord[k][j+1], pdbcoord[k][j+6], to360), deci);}
		else { etaP = "NA";}}
	else if (order == "C1'P  P  C1'"){
		thetaP = ftsround(torsion_angle(pdbcoord[k][i+2], pdbcoord[k][i], pdbcoord[k][i+3], pdbcoord[k][i+5], to360), deci);  // ATOMS : C1'-P-P-C1' --> P-C1'-P-C1'
		if (pdbcoord[k][i+6][4] == "P  "){
			etaP = ftsround(torsion_angle(pdbcoord[k][j-2], pdbcoord[k][j+1], pdbcoord[k][j+3], pdbcoord[k][j+4], to360), deci);}  // ATOMS : C1'-(P-)P-C1'-P--> (P-)C1'-P-C1'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			etaP = ftsround(torsion_angle(pdbcoord[k][j-2], pdbcoord[k][j+1], pdbcoord[k][j+3], pdbcoord[k][j+6], to360), deci);}
		else { etaP = "NA";}}
	else if (order == "C1'P  C1'P  "){
		thetaP = ftsround(torsion_angle(pdbcoord[k][i+2], pdbcoord[k][i], pdbcoord[k][i+5], pdbcoord[k][i+3], to360), deci);  // ATOMS : C1'-P-C1'-P --> P-C1'-P-C1'
		if (pdbcoord[k][i+6][4] == "P  "){
			etaP = ftsround(torsion_angle(pdbcoord[k][j-2], pdbcoord[k][j+3], pdbcoord[k][j+1], pdbcoord[k][j+4], to360), deci);} // ATOMS : C1'-(P-)C1'-P-P --> (P-)C1'-P-C1'-P
		else if ((pdbcoord[k][i+6][4] == "C1'") && (pdbcoord[k][i+8][4] == "P  ")){
			etaP = ftsround(torsion_angle(pdbcoord[k][j-2], pdbcoord[k][j+3], pdbcoord[k][j+1], pdbcoord[k][j+4], to360), deci);}
		else { etaP = "NA";}}
	else {
		thetaP = "NA";  // Returns NA if the atoms in the backbone are 
		etaP = "NA";    // not well referenced for 1 pair of residue
		pdbmistake = true;
		if (pdbcoord[k][i+1][4] == "P  ")        //
		{                                        //
			i -= 2;                              //
		} else if (pdbcoord[k][i+2][4] == "P  ") // Try to find another pair of residue
		{                                        // with all their backbone atoms for the
			i -= 1;                              // next step
		} else {                                 //
		while(pdbcoord[k][i+3][4] != "P  "){     //
			i += 1;}                             //  
		}
	}
}

string getHeaderPdbOutput(string mol, bool Omega, bool PosChain, bool ShowFile, bool alterC1, bool C4andC1, string separator){
    vector<string> names;
    string pos = "POSITION";
    string chain = "CHAIN";
    string res = "RESIDUE";
    string pdb_file = "PDB_FILE";
    string output;
    if (mol == "prot") {
        output = "PHI" + separator + "PSI";
    } else if (mol == "rna") {
        string head_alterC1 = "ETA'"+ separator + "THETA'" ;
        string head_C1 = "ETA" + separator + "THETA";
        if (alterC1){
            output = head_alterC1;
        }
        if (C4andC1){
            output = head_C1 + separator + head_alterC1;
        }
        else {
            output = head_C1;
        }
    }
    if (Omega){
        output += separator + "OMEGA";
    }
    if (PosChain){
        output += separator + res + separator + chain + separator + pos;
    }
    if (ShowFile){
        output += separator + pdb_file;
    }
    output += "\n";
    return output;
}

string addValuesToOutput(string angles, string angle_omega,
                    string coords, string positions, string residus, string pdb_file,
                    bool Omega, bool PosChain, bool ShowFile, string separator){
        string output = angles;
        if (Omega){
            output += separator + angle_omega;
        }
        if (PosChain){
            output += separator + coords + separator + residus + separator + positions;
        }
        if (ShowFile){
            output+= separator + pdb_file;
        }
        output += "\n" ;
        return output;
}

int main(int argc, char** argv)
{
	string optlist =
		"   Usage:\n"
		"   ./angle_calculation [-d PATHWAY_DATASET] [-o OUTPUT_FILE_NAME] [-O] [-R] [-a|-A]\n"
		"                       [-t] [-i DECIMAL_PLACE] [-p] [-f] [-v]\n\n"
		"   Options:\n"
		"   -d   string   Pathway of the repository where PDB files you interested of are\n"
		"   -o   string   Name of your file in output (ex.: YourOuputName_pdbcode.txt).\n"
		"   -O            Add Omega angles values to those of Psi and Phi angles\n\n"
		"   Access to the RNA processing mode and its options:\n"
		"   -R            Calculation of Theta and Eta pseudotorsion angles (Using atoms P and C4')\n"
		"   -Ra           Alternative calculation (Theta'/Eta') using atoms P and C1'\n"
		"   -RA           Calculation of pseudotorsion angles using both methods\n\n"
		"   Additional options:\n"
		"   -p            Allow the user to see the residues sequence number, the type of residues \n"
		"                   and their chain identifiers on the output file\n"
		"   -f            Allow the user to see from which PDB the observed angle values come from\n"
		"   -i   int      Choosing the number of decimal places (default : 3)\n"
		"   -t            For users who want angle values between 0° and 360° (default : [-180°;180°])\n"          
		"   -v            Output in the terminal the values.\n\n"
		"   -h            Help\n\n";

	string in_dir;    // Pathway of the repository
	string output;    // Output file name without extension (".txt",".out",etc...)
	bool PosChain = false;  // To add residue sequence number and chain identifier
	bool ShowFile = false;  // To add PDB file code
	bool Omega = false;   // Omega angle values are in option (default : false)
	bool Rna = false;     // Turn on RNA mode for pseudotorsion angles (default : protein psi-phi)
	bool alterC4 = true;  // Only pseudotorsion angles with C4' (RNA)
	bool alterC1 = false; bool check_a = false;  // Only pseudotorsion angles with C1' (RNA)
	bool C4andC1 = false; bool check_A = false;  // Pseudotorsion angles with C4' and C1' (RNA)
	int decimal = 3;  // Decimal places
	bool to360 = false;  // [-180°;180°] or [0°;360°]
	bool verbose= false;
	bool toTerminal = false;

	int opt;
	while ((opt = getopt(argc,argv, "hpfORaAtivd:o:")) != EOF){
		switch(opt){
			case 'd': in_dir = optarg; break;
			case 'o': output = optarg; break;
			case 'p': PosChain = true; break;
			case 'f': ShowFile = true; break;
			case 'i': decimal = stoi(optarg); break;
			case 't': to360 = true; break;
			case 'O': Omega = true; break;
			case 'R': Rna = true; break;
			case 'v': verbose=true; break;
			case 'a': alterC4 = false; alterC1 = true; check_a = true; break;
			case 'A': alterC4 = false; C4andC1 = true; check_A = true; break;
			case 'h': fprintf(stderr, "%s", optlist.c_str()); return 0;
		}
	}
	if (argc == 1){ fprintf(stderr, "%s", optlist.c_str()); return 1; }
	if (Rna and Omega){ cerr << "Error (option) : There is no omega pseudotorsion angle [-O] in RNA structure\n" << endl; return 1;}
	if (!Rna and alterC1){ cerr << "Error (option) : Turn into the RNA mode [-R] to use option [-a]\n" << endl; return 1;}
	if (!Rna and C4andC1){ cerr << "Error (option) : Turn into the RNA mode [-R] to use option [-A]\n" << endl; return 1;}
	if (check_a and check_A){ cerr << "Error (option) : -a and -A are incompatible options\n" << endl; return 1;}
	if (decimal < 0) {cerr << "Error (user) : The number of decimal places [-i] must be greater than or equal to 0\n" << endl; return 1;}
	if (output.empty()){ output = "output.csv"; toTerminal = true; }

	string position1; string position2; string position3; string res_chain; string pdb_file;
	string head3; string head4; string head_res = "PAIR"; string head_pos; string head_ch; string head_f;

	if (!Rna)
	{
	map<string,string> code3to1;
	code3to1["ALA"]="A"; code3to1["CYS"]="C"; code3to1["ASP"]="D"; code3to1["GLU"]="E"; code3to1["PHE"]="F"; code3to1["GLY"]="G"; code3to1["HIS"]="H"; code3to1["ILE"]="I"; 
	code3to1["LYS"]="K"; code3to1["LEU"]="L"; code3to1["MET"]="M"; code3to1["ASN"]="N"; code3to1["PRO"]="P"; code3to1["GLN"]="Q"; code3to1["ARG"]="R"; code3to1["SER"]="S"; 
	code3to1["THR"]="T"; code3to1["VAL"]="V"; code3to1["TRP"]="W"; code3to1["TYR"]="Y";

	/**** Processing for each file in PDB files list ****/

	string line;

	string ffile = output;  // Output file name
	ofstream file_out;
	file_out.open(ffile);           // Open a new file for angle values
    string separator = ",";
	if (file_out.is_open())
	{
        string head_file = getHeaderPdbOutput("prot", Omega, PosChain, ShowFile, false, false, separator);
        if (!toTerminal){
            file_out << head_file;
        } else {
            cout << head_file;
        }

	int cpt = 0; int cptot = 0;
    vector<string> pdb_names = get_list_dir(in_dir);
    for(int n_file = 0; n_file < pdb_names.size(); n_file++)
	{
		bool pdbmistake = false; bool cutoff = false;
		string fl = pdb_names[n_file];
		vector<vector<vector<string> >> Coords;
		if ((fl.size() == 9) || (fl.size() == 5))
		{
			Coords = sch_coord_pdb(fl.substr(0,4)+".pdb", fl.substr(4,1));
			if (ShowFile){ pdb_file = "          "+fl.substr(0,4);}
		} else if (fl.size() == 4){
			Coords = coord_pdb(fl+".pdb");
			if (ShowFile){ pdb_file = "          "+fl;}
		} else {
			Coords = coord_pdb(fl);  // 3D vector {Chain[Atom[informations]]}
			if (ShowFile){ pdb_file = "          "+fl.substr(0,fl.size()-4);}}

		if (Coords.empty())
		{
			cptot += 1;
			cerr << "\n" << cptot  << " : " << "Error: There is no protein sequence in this PDB file ("+fl+")" << endl;
		} else {

			for (int k = 0; k < Coords.size(); k++)
			{
				if (Coords[k].size() >= 6) {
					for (int i = 0; i < Coords[k].size()-5; i+=3)
					{
                        string coords;
                        string positions;
                        string angle_omega ;
                        string angle_psi;
                        string angle_phi;

                        if (PosChain){
                              position1 = removeWhitespace(Coords[k][i][5]);
                              position2 = removeWhitespace(Coords[k][i+3][5]);
                              res_chain = removeWhitespace(Coords[k][i][6]);  // Add position, chain and PDB code
                              coords = code3to1[Coords[k][i][3]] + code3to1[Coords[k][i+3][3]];
                              positions = position1 + "/" + position2;
                        }
                        string order = Coords[k][i][4]+Coords[k][i+1][4]+Coords[k][i+2][4]+Coords[k][i+3][4];
                        if (order == "N CAC N "){
                            int j = i+2; angle_psi = ftsround(torsion_angle(Coords[k][i], Coords[k][i+1], Coords[k][i+2], Coords[k][i+3], to360), decimal);    // ATOMS : N-CA-C-N
                            angle_phi = ftsround(torsion_angle(Coords[k][j], Coords[k][j+1], Coords[k][j+2], Coords[k][j+3], to360), decimal);    // ATOMS : C-N-CA-C
                            if (Omega){
                                int a = i+1;
                                angle_omega = ftsround(torsion_angle(Coords[k][a], Coords[k][a+1], Coords[k][a+2], Coords[k][a+3], to360), decimal);   // ATOM : CA-C-N-CA
                            }
                        } else if (order == "N C CAN "){

                              int j = i+1;
                              angle_psi = ftsround(torsion_angle(Coords[k][i], Coords[k][i+2], Coords[k][i+1], Coords[k][i+3], to360), decimal);    // ATOMS : N-C-CA-N --> N-CA-C-N
                              angle_phi = ftsround(torsion_angle(Coords[k][j], Coords[k][j+2], Coords[k][j+1], Coords[k][j+3], to360), decimal);    // ATOMS : C-CA-N-C --> C-N-CA-C

                              if (Omega)
                              {
                                  int a = i+2;
                                  angle_omega = ftsround(torsion_angle(Coords[k][a], Coords[k][a+2], Coords[k][a+1], Coords[k][a+3], to360), decimal);   // ATOM : CA-N-C-CA --> CA-C-N-CA
                              }
                        } else {
                              angle_psi = "NA";    // Returns NA if the atoms in the backbone are
                              angle_phi = "NA";    // not well referenced for 1 pair of residue
                              pdbmistake = true;
                        }
                            string angles = angle_phi + separator + angle_psi ;
                            string output_to_file = addValuesToOutput(angles, angle_omega,
                                            coords, positions,res_chain, removeWhitespace(pdb_file),
                                            Omega, PosChain, ShowFile, separator);
                            if (!toTerminal){
                                file_out << output_to_file;
                            } else {
                                cout << output_to_file;
                            }


							if (Coords[k][i+1][4] == "N ")         //
							{                                      //
								i -= 2;                            //
							} else if (Coords[k][i+2][4] == "N ")  // Try to find another pair of residue
							{                                      // with all their backbone atoms for the
								i -= 1;                            // next step
							} else {                               //
							while(Coords[k][i+3][4] != "N "){      //
								i += 1;}                           //
							}
						}

				} else {
				cutoff = true;}  // Check the length of the sequence, return an error if it is under the cutoff
			}
		
	if (pdbmistake) {
		cpt += 1; cptot += 1;
		cerr << "\n" << cptot  << " : " << "Warning: Potential badly written text in the PDB file ("+fl+")" << endl;  // Insert an error message if there is at least 1 written mistake in the PDB file
	} else if (cutoff) {
		cpt += 1; cptot += 1;
		cerr << "\n" << cptot  << " : " << "Warning (chain length too short): Presence of protein residues, but in insufficient number for at least 1 chain ("+fl+")" << endl;	
	} else {
	cpt +=1; cptot += 1;
	if (verbose){
        cout.flush();
        cout << "\r" << "Processed PDB files :" << cptot;
	}
    }  // To see the evolution of the processing
		}
	}
	if (verbose){
        cout << "\nTotal processed files : " << cpt << " on " << cptot << " given" << endl; // To see how many file was processed at the end
	}
	file_out.close();}
	return 0;}


	else {   /* If you want to work with RNA sequences */

	/**** Processing for each file in PDB files list for RNA structures ****/
	string line;

	string ffile = output;  // Output file name
	ofstream file_out;
	file_out.open(ffile);           // Open a new file for angle values
	string separator = ",";
	if (file_out.is_open())
	{
        string head_file = getHeaderPdbOutput("rna", false, PosChain, ShowFile, alterC1, C4andC1, separator);
        if (!toTerminal){
            file_out << head_file;
        } else {
            cout << head_file;
        }


	int cpt = 0; int cptot = 0;
    vector<string> pdb_names = get_list_dir(in_dir);
     for(int n_file = 0; n_file < pdb_names.size(); n_file++)
	{
		bool pdbmistake = false; bool cutoff = false;
		string fl = pdb_names[n_file];
		vector<vector<vector<string> >> Coords;
		if ((fl.size() == 9) || (fl.size() == 5))
		{
			Coords = sch_coord_pdb(fl.substr(0,4)+".pdb", fl.substr(4,1), true);
			if (ShowFile){ pdb_file = "          "+fl.substr(0,4);}
		} else if (fl.size() == 4){
			Coords = coord_pdb(fl+".pdb", true);
			if (ShowFile){ pdb_file = "          "+fl;}
		} else {
			Coords = coord_pdb(fl, true);  // 3D vector {Chain[Atom[informations]]}
			if (ShowFile){ pdb_file = "          "+fl.substr(0,fl.size()-4);}}

		if (Coords.empty())
		{
			cptot += 1;
			cerr << "\n" << cptot  << " : " << "Error: There is no RNA sequence in this PDB file ("+fl+")" << endl;
		} else {

			for (int k = 0; k < Coords.size(); k++)
			{
				if (Coords[k].size() >= 7)
				{
					for (int i = 0; i < Coords[k].size()-6; i+=3)
					{
                        string angles = "";
                        string coords ;
                        string positions;
						if (PosChain)
						{
							position1 = removeWhitespace(Coords[k][i][5]);
							position2 = removeWhitespace(Coords[k][i+3][5]);
							position3 = removeWhitespace(Coords[k][i+6][5]);
							res_chain = removeWhitespace(Coords[k][i][6]);
                            coords = Coords[k][i][3] + Coords[k][i+3][3] + Coords[k][i+6][3];
                            positions = position1 + "/" + position2 + "/" + position3;
						}

						string order1 = Coords[k][i][4]+Coords[k][i+1][4]+Coords[k][i+3][4]+Coords[k][i+4][4];
						string order2 = Coords[k][i][4]+Coords[k][i+2][4]+Coords[k][i+3][4]+Coords[k][i+5][4];
						if (alterC4)
						{
							int j = i+1; string angle_theta; string angle_eta;
							get_theta_eta(Coords, order1, angle_theta, angle_eta, k, i, j, pdbmistake, to360, decimal);
							angles+= angle_eta + separator + angle_theta ;
						} else if (alterC1)
						{
							int j = i+2; string angle_thetaP; string angle_etaP;
							get_thetaP_etaP(Coords, order2, angle_thetaP, angle_etaP, k, i, j, pdbmistake, to360, decimal);
							angles+= angle_etaP + separator + angle_thetaP ;
						} else if (C4andC1)
						{
							int j1 = i+1; int j2 = i+2; string angle_theta; string angle_eta; string angle_thetaP; string angle_etaP;
							get_theta_eta(Coords, order1, angle_theta, angle_eta, k, i, j1, pdbmistake, to360, decimal, false);  // Reagjustement of i during get_thetaP_etaP()
							get_thetaP_etaP(Coords, order2, angle_thetaP, angle_etaP, k, i, j2, pdbmistake, to360, decimal);
							angles += angle_eta + separator + angle_theta + separator + angle_etaP + separator + angle_thetaP;
						}
                    string output = addValuesToOutput(angles, "", coords, positions, res_chain, removeWhitespace(pdb_file), false, PosChain, ShowFile, separator);
                    if (!toTerminal){
                        file_out << output;
                    } else {
                        cout << output;
                    }
					}
				} else {
				cutoff = true;}  // Check the length of the sequence, return an error if it is under the cutoff
			}

	if (pdbmistake) {
		cpt += 1; cptot += 1;
		cerr << "\n" << cptot  << " : " << "Warning: Potential badly written text in the PDB file ("+fl+")" << endl;  // Insert an error message if there is at least 1 written mistake in the PDB file
	} else if (cutoff) {
		cpt += 1; cptot += 1;
		cerr << "\n" << cptot  << " : " << "Warning (chain length too short): Presence of RNA residues, but in insufficient number for at least 1 chain ("+fl+")" << endl;
	} else {
	cpt +=1; cptot += 1;
	if (verbose){
        cout.flush();
        cout << "\r" << "Processed PDB files :" << cptot;
	}
    }
		}
	}
	if (verbose){
        cout << "\nTotal processed files : " << cpt << " on " << cptot << " given" << endl;
	}
	file_out.close();}
	return 0;}
}
