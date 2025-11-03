/*! \file aggregate.cpp 
\brief Evaluate the aggregates distribution along z axis and the kinetic constants (Author: Jordi Andreu).
Dump file from LAMMPS-May-2008 required.
Last Modification: 28-Març-2010.
Developing State: Optimization.
*/ 

#include<iostream>
//#include<cstdlib>
#include<cstdio>
#include<cmath>
#include <vector>
#include <ctime>
#include"bead.h"
#include"aggregate.h"

#define _MAX_AGG_SIZE 160 ///< Maximum number of allowed beads in an aggregate.
#define _ts_offset 0 ///< Set the timestep offset for computing averages.
#define _DEBUG false ///< Switch on/off several debugging functions in the code.

using namespace std;


CAggregate::CAggregate(){};

	/// Constructor.
	/** Initializes an aggregate with a given id and with size=0.
	* \param i Id number assigned to the aggregate.
	*/
CAggregate::CAggregate(int i){ id=i; size=0; }
	
	/// Adds a bead to this aggregate.
	/**	This memeber function adds a new bead object to the aggregate and updates its size.
	* \param a A CBead element to be added to this aggregate.
	*/
void CAggregate::AddBead(CBead& a){

		bead.push_back(a);
		size=bead.size();
		
	}

	/// Prints aggregate information on screen.
void CAggregate::Print(){

		cout << "aggregate " << id << " with size " << bead.size() << " and composed by bead(s)";
		for (int i=0; i<bead.size(); i++){ cout << " " << bead[i].id; }
		cout << endl;

	}
	///  Search a given bead in this aggregate.
	/** Checks if bead a belongs to the aggregate by comparing the corresponding id.
	* \param a Bead to be searched in the aggregate.
	* \return boolean expression.
	*/
bool CAggregate::CheckBead(CBead a){
	
		bool found=false;
		
		for (int i=0; i<bead.size(); i++){
		
			if (bead[i].id == a.id){found=true;}
		
		}
		
		return found;	
		
	}

