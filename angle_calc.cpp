#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// Store the coordinates of the atoms we are interested in for 1 PDB file
vector<vector<float> > coord_pdb(string pdbfile){
	vector<vector<float> > tableau;
	vector<string> tab_atom;
	ifstream fl(pdbfile);
	string line;
	string Atom;
	const vector<string> refer = {"N ","CA","C "};
	int i = -1;

	while(getline(fl, line)){
		if (line.substr(0,4) == "ATOM"){
			Atom = line.substr(13,2);
			tab_atom.push_back(Atom);
			if (tab_atom.front() != refer.front())
			{
				tab_atom.pop_back();
			} else {
				if (*find(refer.begin(), refer.end(), Atom) == Atom)
				{
					float X = stof(line.substr(30,8));  //
					float Y = stof(line.substr(38,8));  // Convert string to float
					float Z = stof(line.substr(46,8));  //
					tableau.push_back(vector<float>(3,0.0));
					i += 1;
					tableau[i][0] = X;  //
					tableau[i][1] = Y;  // Atomic coordinates x, y, z
					tableau[i][2] = Z;  //
				}
			}		
		}
	}
	return tableau;
}


// Distance between 2 atoms (A)
float distance(vector<float> atom1, vector<float> atom2){
	return sqrt(pow(atom1[0]-atom2[0], 2)+pow(atom1[1]-atom2[1], 2)+pow(atom1[2]-atom2[2], 2));
}


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
float torsion_angle(vector<float> atom1, vector<float> atom2, vector<float> atom3, vector<float> atom4){
	vector<float> vecteur12 = {atom1[0]-atom2[0],atom1[1]-atom2[1],atom1[2]-atom2[2]};
    vector<float> vecteur23 = {atom2[0]-atom3[0],atom2[1]-atom3[1],atom2[2]-atom3[2]};
    vector<float> vecteur34 = {atom3[0]-atom4[0],atom3[1]-atom4[1],atom3[2]-atom4[2]};

    vector<float> vecteur_normal1 = vector_prod(vecteur12, vecteur23);
    vector<float> vecteur_normal2 = vector_prod(vecteur23, vecteur34);

    if (scalar_prod(vecteur23,vector_prod(vecteur_normal1,vecteur_normal2))<0)
    {
	    float angle = acos(scalar_prod(vecteur_normal1,vecteur_normal2)/
	    	          (norme(vecteur_normal1)*norme(vecteur_normal2)))*180/M_PI;
	    return angle;
    } else {
    	float angle = -acos(scalar_prod(vecteur_normal1,vecteur_normal2)/
	    	          (norme(vecteur_normal1)*norme(vecteur_normal2)))*180/M_PI;
    	return angle;
    }   
}


int main(int argc, char** argv)
{
	string in_dir = argv[1];      // Pathway of the repository

	/* Processing for each file in PDB files list */

	ifstream my_pdbs(argv[2]);    // File containing PDB files list
	string line;
	while(getline(my_pdbs, line)) // 1 line = 1 PDB file
	{
		string fl = line;
		cout << fl << endl;
		vector<vector<float> > Coords;
		string namef = argv[3];         // Output file name without extension (".txt",".out",etc...)
		string ffile = namef+"_"+fl.substr(0,fl.size()-4)+".txt";  // Output file name + processed PDB code
		
		ofstream file_out;

		Coords = coord_pdb(in_dir+fl);
		file_out.open(ffile);           // Open a new file for angle values
		if (file_out.is_open())
		{
			for (int i = 0; i < Coords.size()-6; i += 4)
			{
				int j = i+2;
				float angle_psi = torsion_angle(Coords[i], Coords[i+1], Coords[i+2], Coords[i+3]);
				float angle_phi = torsion_angle(Coords[j], Coords[j+1], Coords[j+2], Coords[j+3]);
				file_out << angle_psi << "   " << angle_phi << endl;
			}
			cout << "Done !" << endl;
			file_out.close();
		}
	}
	return 0;
}
