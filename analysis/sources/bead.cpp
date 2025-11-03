/*! \file bead.cpp 
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

#define _MAX_AGG_SIZE 160 ///< Maximum number of allowed beads in an aggregate.
#define _ts_offset 0 ///< Set the timestep offset for computing averages.
#define _DEBUG false ///< Switch on/off several debugging functions in the code.

using namespace std;


void CBead::AddNeighbor(CBead* a){ neighbor.push_back(a); }

void CBead::Print(){ cout << "bead id: " << id; }
