// Reading dump file from LAMMPS-21May08
// Author: Jordi Andreu
// Date: 17-Octobber-2010

/*! \file dump.h
L'objecte utilitzat per a llegir cadascuna de les configuracions
d'un arxiu dump del LAMMPS és de la classe CDump
Les funcions definides dins d'aquesta classe ens permeten
manipular el Dump (veure la definicio d'aquestes dins la classe)
*/

#ifndef LIBDUMP_H
#define LIBDUMP_H

#include<cstdlib>
#include<cstdio>
#include<fstream>
#include<iostream>
//#include <cstring>

using namespace std;

//bool FileExists(const string& filename);

// ------ Structs --------------------------------------------------
struct t_header
{
	int natoms;		/** < number of atoms. */
	double timestep;	/** < configuration timestep. */
	double xlo;		/** < lowest boundary value in x. */
	double xhi;		/** < highest boundary value in x. */
	double ylo;		/** < lowest boundary value in y. */
	double yhi;		/** < highest boundary value in y. */
	double zlo;		/** < lowest boundary value in z. */
	double zhi;		
	
};

struct t_line
{
	int id;			/** < particle identifier. */
	int type;		/** < particle type. */
	double r[3];	/** < position vector. */
	double v[3];	/** < velocity vector. */
	double q;		/** < electrical charge. */
};

// ------ Classes --------------------------------------------------
class CDump{

public:

	char filename[20];	/** < dump filename. */
	string option;		/** < pattern mode. */
	string normalized;	/** < normalization. */
	ifstream* inFile;	/** < stream object mapping input dump file. */
	t_header header;	/** < header struct. */
	t_line* particle;	/** < particle record. */

	double xbox;	/** < x-axis simulation box length. */
	double ybox;	/** < y-axis simulation box length. */
	double zbox;	/** < z-axis simulation box length. */
	bool nflag;		/** < true indicates normalized coordinates (Default false). */
	
	/// Prints the loaded configuration
	/**	Prints the header and the body of a configuration dump
	*/
	void Print();

	/// Prints the header
	/**	Prints the header af a given configuration
	*/
	void PrintHeader();

	/// Prints the body of the configuration
	/**	Prints the body af a given configuration
	*/
	void PrintConfig();
	
	/// Open a file associated to a dump
	/**	Open a file associated to a dump
	*/	
	void Open();

	/// Close a file associated to a dump
	/**	Close a file associated to a dump
	*/	
	void Close();
	
	/// Shift the timestep of a given loaded dump
	/**	Shift the timestep of a given loaded dump file by a certain amount given by
	* \param outfile ofstream stream object associated to the output file
	* \param shift number of configurations to be shifted in timesteps	
	*/	
	void Shift(ofstream& outfile,int shift);
	
	/// Write dump to file (Deprecated)
	/**	Writes the given dump to a file (Deprecated)
	* \param outfile ofstream stream object associated to the output file
	*/		
	void Write(ofstream& outfile, bool rescale=true, int shift=0);
	
	/// Write dump to file with normalized coordinates (Deprecated)
	/**	Writes the given dump to a file in normalized coordinates (Deprecated)
	* \param outfile ofstream stream object associated to the output file
	* \param rescale (optional) true sets rescaled coordinates as output. Default: true
	* \param rescale (optional) shifts the timestep counter of each configuration this amount. Default: 0
	*/		
	void WriteNormalized(ofstream&);
	
	/// Load a given configuration
	/**	Loads a given configuration by reading from file
	*/		
	bool Load();

// -------------------------------------
	CDump(){}
	
	~CDump(){}
	
};

#endif	//LIBDUMP_H
