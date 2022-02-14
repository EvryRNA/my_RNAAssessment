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
	ifstream fl(pdbfile);
	string line;
	const vector<string> refer = {"N ","CA","C "};
	int i = -1;

	while(getline(fl, line)){
		if (line.substr(0,4) == "ATOM"){
			if (*find(refer.begin(), refer.end(), line.substr(13,2)) == line.substr(13,2))
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

    float angle = sin(scalar_prod(vecteur_normal1,vecteur_normal2)/(norme(vecteur_normal1)
    			*norme(vecteur_normal2)))*180;
    return angle;
}


int main(int argc, char** argv)
{
	string in_dir = argv[1];			// Pathway of the repository

// Stokage of the file names of the directory given in input
	DIR *dir; 				// Pointeur du répertoire, fonctionne avec readdirr (package <dirent.h>)
	struct dirent *diread;	// Accès aux noms des fichiers dans le répertoire par la structure dirent (package <dirent.h>)
    vector<char *> files;								//
	if ((dir = opendir(argv[1])) != nullptr) {			//
		while ((diread = readdir(dir)) != nullptr) {	// Listes des fichiers dans le répertoire
	               files.push_back(diread->d_name);		//
		}
		closedir (dir);
	} else {
		perror ("opendir");
		        return EXIT_FAILURE;
	}

// File processing if it is a PDB file
	ofstream flo;
	flo.open("liste_pdb.txt"); // Open a new file to stock all the PDB file names in the repository
	if (flo.is_open())
	{
		for (auto file : files)
		{
			int lenf = strlen(file);
			string fl = file;
			if ((lenf > 5) && (fl.substr(lenf-4,4) == ".pdb"))
			{
				flo << file << endl;
			}
		}
	}
	flo.close();

	ifstream my_pdb("liste_pdb.txt");
	string line;
	while(getline(my_pdb, line)) // 1 line = 1 PDB file
	{
		string fl = line;
		cout << fl << endl;
		vector<vector<float> > Coords;
		string namef = argv[2];
		string ffile = namef+"_"+fl.substr(0,fl.size()-4)+".txt";  // Output file name + processed PDB code
		
		ofstream file_out;

		Coords = coord_pdb(in_dir+fl);
		file_out.open(ffile); // Open a new file
		if (file_out.is_open())
		{
			for (int i = 1; i < Coords.size()-4; i += 4)
			{
				float angle = torsion_angle(Coords[i], Coords[i+1], Coords[i+2], Coords[i+3]);
				file_out << angle << endl;
			}
			cout << "Done !" << endl;
			file_out.close();
		}
	}
	return 0;
}
