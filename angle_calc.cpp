#include <iostream>
#include <fstream>
#include <string>
#include <map>      // for dictionary
#include <cstring>  // stof()
#include <vector>
#include <cmath>
#include <algorithm> // find()
#include <getopt.h>

using namespace std;

/**** Store the coordinates of the atoms we are interested in for 1 PDB file ****/
vector<vector<vector<string> >> coord_pdb(string pdbfile){
	vector<vector<vector<string> >> tableau;         // Atomic coordinates for all chains
	vector<vector<string> > interm_tab;              // Atomic coordinates for current chain
	vector<string> tab_atom;
	vector<string> tab_chain = {"A"};               // To know in which chain we are currently (start with "A")
	string line;
	string Atom;
	float occupancy;
	const vector<string> refer = {"N ","CA","C "};  // Atoms that interest us (P and C4' for RNA)
	int i = -1;

	ifstream fl(pdbfile);
	while(getline(fl, line)){
		if (line.substr(0,4) == "ATOM"){
			Atom = line.substr(13,2);
			tab_atom.push_back(Atom);               //
			if (tab_atom.front() != refer.front())  // Avoids possible mismatches if the N-ter part
			{                                       // is missing from the PDB file 
				tab_atom.pop_back();                //
			} else {
				if (tab_chain.back() == line.substr(21,1))  // If still in the same chain
				{
					if (*find(refer.begin(), refer.end(), Atom) == Atom)
					{
						string residu = line.substr(17,3);  // Residu name (3 letters)
						if ((line.substr(16,4) == " "+residu)  || (line.substr(16,4) == "A"+residu))
						{
							string X = line.substr(30,8);  //
							string Y = line.substr(38,8);  // Atomic coordinates x, y, z
							string Z = line.substr(46,8);  //
							occupancy = stof(line.substr(54,6));  // For residu alternate location (RAL)
							interm_tab.push_back(vector<string>(5));
							i += 1;
							interm_tab[i][0] = X;       // 
							interm_tab[i][1] = Y;       //
							interm_tab[i][2] = Z;       // Stock the coordinates, residu name
							interm_tab[i][3] = residu;  // and atom in an intermediate vector
							interm_tab[i][4] = Atom;    //
						} else if (stof(line.substr(54,6)) > occupancy) // Compare the occupancy if there is RAL
						{
							string X = line.substr(30,8);  //
							string Y = line.substr(38,8);  // Atomic coordinates x, y, z
							string Z = line.substr(46,8);  //
							occupancy = stof(line.substr(54,6));  // For residu alternate location
							interm_tab[i][0] = X;       // 
							interm_tab[i][1] = Y;       //
							interm_tab[i][2] = Z;       // Stock the coordinates, residu name
							interm_tab[i][3] = residu;  // and atom in an intermediate vector
							interm_tab[i][4] = Atom;    //
						}
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
						string X = line.substr(30,8);  //
						string Y = line.substr(38,8);  // Atomic coordinates x, y, z
						string Z = line.substr(46,8);  //
						string residu = line.substr(17,3);  // Residu name (3 letters)
						occupancy = stof(line.substr(54,6));  // For residu alternate location
						interm_tab.push_back(vector<string>(5));
						i += 1;
						interm_tab[i][0] = X;       // 
						interm_tab[i][1] = Y;       //
						interm_tab[i][2] = Z;       // Stock the coordinates, residu name
						interm_tab[i][3] = residu;  // and atom in an intermediate vector
						interm_tab[i][4] = Atom;    //
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
float torsion_angle(vector<string> atom1, vector<string> atom2, vector<string> atom3, vector<string> atom4){
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
    	return angle;
    }   
}


int main(int argc, char** argv)
{
	string in_dir;    // Pathway of the repository
	string listpdb;   // File containing PDB files list
	string output;    // Output file name without extension (".txt",".out",etc...)
	bool Omega = false;   // Omega angle values are in option (default : false)

	int opt;
	while ((opt = getopt(argc,argv, "Od:l:o:")) != EOF){
		switch(opt){
			case 'd': in_dir = optarg; break;
			case 'l': listpdb = optarg; break;
			case 'o': output = optarg; break;
			case 'O': Omega = true; break;
		}
	}

	//string in_dir = argv[1];      // Pathway of the repository

	map<string,string> code3to1;
	code3to1["ALA"]="A"; code3to1["CYS"]="C"; code3to1["ASP"]="D"; code3to1["GLU"]="E"; code3to1["PHE"]="F"; code3to1["GLY"]="G"; code3to1["HIS"]="H"; code3to1["ILE"]="I"; 
	code3to1["LYS"]="K"; code3to1["LEU"]="L"; code3to1["MET"]="M"; code3to1["ASN"]="N"; code3to1["PRO"]="P"; code3to1["GLN"]="Q"; code3to1["ARG"]="R"; code3to1["SER"]="S"; 
	code3to1["THR"]="T"; code3to1["VAL"]="V"; code3to1["TRP"]="W"; code3to1["TYR"]="Y";

	/**** Processing for each file in PDB files list ****/

	ifstream my_pdbs(listpdb);    // File containing PDB files list //argv[2]
	string line;
	while(getline(my_pdbs, line)) // 1 line = 1 PDB file
	{
		bool pdbmistake = false;
		string fl = line;
		cout << fl << endl;
		vector<vector<vector<string> >> Coords;
		//string namef = argv[3];                   // Output file name without extension (".txt",".out",etc...) 
		string ffile = output+"_"+fl.substr(0,fl.size()-4)+".txt";  // Output file name + processed PDB code
		
		ofstream file_out;

		Coords = coord_pdb(in_dir+fl);  // 3D vector {Chain[Atom[informations]]}

		file_out.open(ffile);           // Open a new file for angle values
		if (file_out.is_open())
		{
			for (int k = 0; k < Coords.size(); k++)
			{
		        file_out << "Chain " << k+1 << endl;
				cout << "Chain " << k+1 << endl;
				for (int i = 0; i < Coords[k].size()-5; i+=3)
				{
					string order = Coords[k][i][4]+Coords[k][i+1][4]+Coords[k][i+2][4];
					if (order == "N CAC ")
					{
						int j = i+2;
						float angle_psi = torsion_angle(Coords[k][i], Coords[k][i+1], Coords[k][i+2], Coords[k][i+3]);    // ATOMS : N-CA-C-N 
						float angle_phi = torsion_angle(Coords[k][j], Coords[k][j+1], Coords[k][j+2], Coords[k][j+3]);    // ATOMS : C-N-CA-C
						if (Omega)     // If user wants also the omega angle values
						{
							int a = i+1;
							float angle_omega = torsion_angle(Coords[k][a], Coords[k][a+1], Coords[k][a+2], Coords[k][a+3]);   // ATOM : CA-C-N-CA
							file_out << angle_psi << "      \t" << angle_phi << "      \t" << angle_omega <<"      \t" << code3to1[Coords[k][i][3]] << code3to1[Coords[k][i+3][3]] << endl;
						} else {
						file_out << angle_psi << "      \t" << angle_phi << "      \t" << code3to1[Coords[k][i][3]] << code3to1[Coords[k][i+3][3]] << endl;}   // Residue concerned for each angle
					} else if (order == "N C CA")
					{
						int j = i+1;
						float angle_psi = torsion_angle(Coords[k][i], Coords[k][i+2], Coords[k][i+1], Coords[k][i+3]);    // ATOMS : N-C-CA-N --> N-CA-C-N
						float angle_phi = torsion_angle(Coords[k][j], Coords[k][j+2], Coords[k][j+1], Coords[k][j+3]);    // ATOMS : C-CA-N-C --> C-N-CA-C
						if (Omega)
						{
							int a = i+2;
							float angle_omega = torsion_angle(Coords[k][a], Coords[k][a+2], Coords[k][a+1], Coords[k][a+3]);   // ATOM : CA-N-C-CA --> CA-C-N-CA
							file_out << angle_psi << "      \t" << angle_phi << "      \t" << angle_omega <<"      \t" << code3to1[Coords[k][i][3]] << code3to1[Coords[k][i+3][3]] << endl;
						} else {
						file_out << angle_psi << "      \t" << angle_phi << "      \t" << code3to1[Coords[k][i][3]] << code3to1[Coords[k][i+3][3]] << endl;}
					} else {

						string angle_psi = "  NA  ";    // Returns NA if the atoms in the backbone 
						string angle_phi = "  NA  ";    // are not referenced in a common way
						pdbmistake = true;
						file_out << angle_psi << "      \t" << angle_phi << "      \t" << code3to1[Coords[k][i][3]] << code3to1[Coords[k][i+3][3]] << endl;
						while(Coords[k][i][4] != "N "){
							i += 1;}   // Try to go to another pair of residu with all their backbone atoms
					}
				}
			}
		}
	if (pdbmistake) {
		cerr << "\nError: Potential badly written text in the PDB file\n" << endl;  // Insert an error message if there is at least 1 written mistake in the PDB file
	} else {
	cout << "Done !" << endl;}
	file_out.close();
	}
	return 0;
}
