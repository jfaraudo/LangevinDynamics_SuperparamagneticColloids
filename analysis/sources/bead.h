/*! \file bead.h 
\brief Evaluate the aggregates distribution along z axis and the kinetic constants (Author: Jordi Andreu).
Dump file from LAMMPS-May-2008 required.
Last Modification: 28-Març-2010.
Developing State: Optimization.
*/ 

#ifndef BEAD_H
#define BEAD_H


#include<iostream>
//#include<cstdlib>
#include<cstdio>
#include<cmath>
#include <vector>
#include <ctime>
#include"dump.h"

#define _MAX_AGG_SIZE 160 ///< Maximum number of allowed beads in an aggregate.
#define _ts_offset 0 ///< Set the timestep offset for computing averages.
#define _DEBUG false ///< Switch on/off several debugging functions in the code.

using namespace std;


/// This class is used to map each particle (bead) in the dump file.
/* Each bead can be linked to other beads as neighbors.
*/
class CBead{
	
public:

        int id;         		/**< Bead id corresponding to the id from LAMMPS dump file. */
	int aggid;			/**< Parent aggregate id. This id is the same for all beads being neighbor or belonging to the same aggregate. */
	vector <CBead*> neighbor;	/**< Vector (pointer) list of neighbor beads.*/

	/// Add a neighbor to this bead.
	/**	This memeber function adds \a a as neighbor by copying its address as a new pointer in the neighbor vector.
	* \param a A CBead element to be added as neighbor.
	*/
	void AddNeighbor(CBead* a);

	/// Print bead id on screen.
	/** Display bead id on screen.
	*/
	void Print ();

};


#endif	//BEAD_H
