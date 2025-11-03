// Reading dump file from LAMMPS-21May08
// Author: Jordi Andreu
// Date: 17-October-2010

//
// L'objecte utilitzat per a llegir cadascuna de les configuracions
// d'un arxiu dump del LAMMPS és de la classe CDump
// Les funcions definides dins d'aquesta classe ens permeten
// manipular el Dump (veure la definicio d'aquestes dins la classe)
//


#include<cstdio>
#include<sstream>
#include<iostream>
#include <cstring>
#include<cstdlib>
#include"dump.h"


using namespace std;

void CDump::Open(){

//	cout << "Old Infile pointer: " << inFile << endl;
	inFile = new ifstream;
//	cout << "New Infile pointer: " << inFile << endl;
	
	do {

		cout << "Enter dump filename: \n";
		cin >> filename;

		inFile->open(filename);

		if (!inFile->is_open()){cerr << "Unable to open file " << filename << endl;}
		
	} while (!inFile->is_open());


//----------- Select input file pattern -------

	option="0";
		
	while (option!="1" && option!="2" && option!="3"){
		cout << "Select fields pattern:\n";
		cout << "1) | id | tag | x | y | z |\n";
		cout << "2) | id | tag | x | y | z | vx | vy | vz |\n";
		cout << "3) | id | tag | x | y | z | vx | vy | vz | q |\n";
		cin >> option;
	}
//------------ Choose rescaled/simulation coordinates ------
	normalized="0";
	while (normalized!="y" && normalized!="n"){
		cout << "Coordinates rescaled to unity by box length?(y/n)";
		cin >> normalized;
	}
	
	nflag = false;
	
	if (normalized == "y"){nflag=true;}

	particle = NULL;

}

void CDump::Close(){

//		cout << "Particle pointer after delete: " << particle << endl;
		delete [] particle;
//		cout << "Particle pointer after delete: " << particle << endl;

		inFile->close();
//		cout << "Infile pointer before delete: " << inFile << endl;		
		delete inFile;
//		cout << "Infile pointer after delete: " << inFile << endl;		
}

//- loadDump ----------------------------------
// Load a complete configuration from a dump
// file. This function returns a boolean, which
// returns 0 when eof is reached.
//---------------------------------------------

bool CDump::Load(){

	string sLine;
	//double xbox=1.0;
	//double ybox=1.0;
	//double zbox=1.0;
	xbox=1.0;
	ybox=1.0;
	zbox=1.0;


//	cout << "inFile memory address: " << inFile << endl;

	getline(*inFile, sLine);
	if (inFile->eof()) {return inFile->eof();}
//	*inFile >> header.timestep; inFile->ignore();

	getline(*inFile, sLine); header.timestep = atof(sLine.c_str());
	
	cout << "new timestep:\t" << header.timestep << endl; 

	getline(*inFile, sLine);//ITEM: NUMBER OF ATOMS
	*inFile >> header.natoms; inFile->ignore();

	getline(*inFile, sLine);//ITEM: BOX BOUNDS
	*inFile >> header.xlo >> header.xhi;
	*inFile >> header.ylo >> header.yhi;
	*inFile >> header.zlo >> header.zhi;
	inFile->ignore();

// Calculate box length for normalized coordinates
	if(normalized=="y"){
		xbox=(header.xhi-header.xlo);
		ybox=(header.yhi-header.ylo);
		zbox=(header.zhi-header.zlo);
	}
	
	if(particle!=NULL){delete [] particle;}
	
	//Memory allocation for each particle. Deallocated in previous line if already allocated previously!
//	cout << "Old particle pointer: " << particle << endl;
	particle =new t_line[header.natoms];
//	cout << "New particle pointer: " << particle << endl;
	
	getline(*inFile, sLine);//ITEM: ATOMS
	for (int i=0; i<header.natoms; i++){

	string token,endline;
	
				getline(*inFile, token,' ');
				particle[i].id = atoi(token.c_str());
				
				getline(*inFile, token,' ');  
				particle[i].type = atoi(token.c_str());
				
				getline(*inFile, token,' ');
				particle[i].r[0] = atof(token.c_str());
								
				getline(*inFile, token,' ');  
				particle[i].r[1] = atof(token.c_str());

				if (option == "1"){
				
					getline(*inFile, token);
					particle[i].r[2] = atof(token.c_str());
					
				}
				else {
					
					getline(*inFile, token, ' ');
					particle[i].r[2] = atof(token.c_str());
				
					getline(*inFile, token,' ');
					particle[i].v[0] = atof(token.c_str());
								
					getline(*inFile, token,' ');  
					particle[i].v[1] = atof(token.c_str());

					}
								
				if (option == "2"){
				
					getline(*inFile, token);
					particle[i].v[2] = atof(token.c_str());
					
				}
				else if (option > "2") {
					
					getline(*inFile, token, ' ');
					particle[i].v[2] = atof(token.c_str());

					getline(*inFile, token);
					particle[i].q = atof(token.c_str());
				}

// Rescaling coordinates from box length normalization to simulation units
		
		particle[i].r[0] = particle[i].r[0]*xbox;
		particle[i].r[1] = particle[i].r[1]*ybox;
		particle[i].r[2] = particle[i].r[2]*zbox;

	}

	return inFile->eof();
}


//- PrintHeader -------------------------------
// Write the header corresponding to the loaded
// configuration in screen.
//---------------------------------------------

void CDump::PrintHeader(){

	cout << "timestep: " << header.timestep << endl;
	cout << "natoms: " << header.natoms << endl;
	cout << "x range: " << header.xlo << "," << header.xhi << endl;
	cout << "y range: " << header.ylo << "," << header.yhi << endl;
	cout << "z range: " << header.zlo << "," << header.zhi << endl;
	cout << "------------------------------------------" << endl;

}

//- PrintConfig -------------------------------
// Print the whole configuration corresponding
// to the loaded configuration in screen.
//---------------------------------------------
void CDump::PrintConfig(){

	for (int i=0; i<header.natoms; i++){
	
		cout << particle[i].id << "\t" << particle[i].type << "\t" << particle[i].r[0] << "\t" << particle[i].r[1] << "\t" << particle[i].r[2] ;
		
		if (option=="2"){
			cout << "\t" << particle[i].v[0] << "\t" << particle[i].v[1] << "\t" << particle[i].v[2] ;
		}
		
		if (option=="3"){
			cout << "\t" << particle[i].q;	
		}
		
		cout << endl;	
	}

}

//- Print -------------------------------------
// Prints the whole dump in screen
//---------------------------------------------
void CDump::Print(){

	PrintHeader();
	PrintConfig();

}

//-// DEPRECATED!! Write -------------------------------------
// Write the loaded dump into a file
//---------------------------------------------

void CDump::Shift(ofstream& out,int shift){

	cout << "Warning: Deprecated Shift function in use" << endl;

	out << "ITEM: TIMESTEP" << endl;
	out << header.timestep + shift << endl;
	out << "ITEM: NUMBER OF ATOMS" << endl;
	out << header.natoms << endl;
	out << "ITEM: BOX BOUNDS" << endl;
	out << header.xlo << " " << header.xhi << endl;
	out << header.ylo << " " << header.yhi << endl;
	out << header.zlo << " " << header.zhi << endl;
	out << "ITEM: ATOMS" << endl;

	for (int i=0; i<header.natoms; i++){
	
		out << particle[i].id << " " << particle[i].type << " " << particle[i].r[0] << " " << particle[i].r[1] << " " << particle[i].r[2] ;
		
		if (option=="2"){
			out << " " << particle[i].v[0] << " " << particle[i].v[1] << " " << particle[i].v[2] ;
		}
		
		if (option=="3"){
			out << " " << particle[i].q;	
		}
		
		out << endl;	
	}
}

//- Write -------------------------------------
// Write the loaded dump into a file
//---------------------------------------------

void CDump::Write(ofstream& out, bool norm, int shift){

	out << "ITEM: TIMESTEP" << endl;
	out << header.timestep + shift << endl;
	out << "ITEM: NUMBER OF ATOMS" << endl;
	out << header.natoms << endl;
	out << "ITEM: BOX BOUNDS" << endl;
	out << header.xlo << " " << header.xhi << endl;
	out << header.ylo << " " << header.yhi << endl;
	out << header.zlo << " " << header.zhi << endl;
	out << "ITEM: ATOMS" << endl;

	if ((nflag && norm) || (!nflag && !norm)){ xbox = 1.0; ybox = 1.0; zbox = 1.0; }
	if (nflag && !norm){ xbox = 1.0/xbox; ybox = 1.0/ybox; zbox = 1.0/zbox;	}
//	if (!nflag && norm){}

	for (int i=0; i<header.natoms; i++){
	
		out << particle[i].id << " " << particle[i].type << " " << particle[i].r[0]/xbox << " " << particle[i].r[1]/ybox << " " << particle[i].r[2]/zbox ;		
		if (option=="2"){ out << " " << particle[i].v[0] << " " << particle[i].v[1] << " " << particle[i].v[2] ;}
		if (option=="3"){ out << " " << particle[i].q; }
		
		out << endl;	
	}
}



// DEPRECATED!! Write the loaded dump WITH COORDINATES NORMALIZED into a file
//---------------------------------------------

void CDump::WriteNormalized(ofstream& out){

	double xbox,ybox,zbox;

	cout << "Warning: Deprecated WriteNormalized function in use" << endl;

	out << "ITEM: TIMESTEP" << endl;
	out << header.timestep << endl;
	out << "ITEM: NUMBER OF ATOMS" << endl;
	out << header.natoms << endl;
	out << "ITEM: BOX BOUNDS" << endl;
	out << header.xlo << " " << header.xhi << endl;
	out << header.ylo << " " << header.yhi << endl;
	out << header.zlo << " " << header.zhi << endl;
	out << "ITEM: ATOMS" << endl;

	xbox=(header.xhi-header.xlo);
	ybox=(header.yhi-header.ylo);
	zbox=(header.zhi-header.zlo);


	for (int i=0; i<header.natoms; i++){
	
		out << particle[i].id << " " << particle[i].type << " " << particle[i].r[0]/xbox << " " << particle[i].r[1]/ybox << " " << particle[i].r[2]/zbox ;
		
		if (option=="2"){
			out << " " << particle[i].v[0] << " " << particle[i].v[1] << " " << particle[i].v[2] ;
		}
		
		if (option=="3"){
			out << " " << particle[i].q;	
		}
		
		out << endl;	
	}
}

