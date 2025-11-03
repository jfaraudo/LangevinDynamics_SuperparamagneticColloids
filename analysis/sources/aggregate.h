/*! \file aggregate.h 
\brief Evaluate the aggregates distribution along z axis and the kinetic constants (Author: Jordi Andreu).
Dump file from LAMMPS-May-2008 required.
Last Modification: 28-Març-2010.
Developing State: Optimization.
*/ 

#ifndef AGGREGATE_H
#define AGGREGATE_H


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


/// This class is used to represent an aggregate made of unlimmited beads.
/* Each aggregate has a unique id and a list of beads.
*/
class CAggregate{
	
public:
	
	int id;			/**< Aggregate id. */
	int size;		/**< Total aggregate size. */
	vector <CBead> bead;	/**< Vector list of beads belonging to this aggregate. */

	CAggregate();

	/// Constructor.
	/** Initializes an aggregate with a given id and with size=0.
	* \param i Id number assigned to the aggregate.
	*/
	CAggregate(int i);
	
	/// Adds a bead to this aggregate.
	/**	This memeber function adds a new bead object to the aggregate and updates its size.
	* \param a A CBead element to be added to this aggregate.
	*/
	void AddBead(CBead& a);

	/// Prints aggregate information on screen.
	void Print();

	///  Search a given bead in this aggregate.
	/** Checks if bead a belongs to the aggregate by comparing the corresponding id.
	* \param a Bead to be searched in the aggregate.
	* \return boolean expression.
	*/
	bool CheckBead(CBead a);

};

#endif	//AGGREGATE_H
