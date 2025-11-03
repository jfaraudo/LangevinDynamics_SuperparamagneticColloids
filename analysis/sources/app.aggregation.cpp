/*! \file app.aggregation.cpp 
\brief Evaluate the aggregates distribution along z axis and the kinetic constants from a Dump file from LAMMPS-May-2008.(Author: Jordi Andreu).
The format of the dump file must follow one of the options defined in CDump class.
It makes use of classes CBead, CAggregate and CDump.
Last Modification: 28-Març-2010.
Developing State: Optimization.
*/ 

#include<iostream>
//#include<cstdlib>
#include<cstdio>
#include<cmath>
#include <vector>
#include <ctime>
#include"dump.h"
#include"bead.h"
#include"aggregate.h"

#define _MAX_AGG_SIZE 160 ///< Maximum number of allowed beads in an aggregate.
#define _ts_offset 0 ///< Set the timestep offset for computing averages.
#define _DEBUG false ///< Switch on/off several debugging functions in the code.

using namespace std;


///// Global VARIABLES Declarations /////
CBead* bead;	//< Array of beads corresponding to each particle in dump file.
vector <CAggregate> aggregate;	//< List of aggregates at a given time t.
vector <CAggregate> prev_aggregate; //< List of aggregates at at a given time t-dt.
CDump dump; //< System configuration info at a given time.
	
int total_beads;	/*< Total number of beads at time t. */
double total_agg;	/*< Total number of aggregates at time t. */
double contact;		/*< Constact distance definition. */
double dt;			/*< Time interval (in seconds) between configurations. */
float first_ts;		/*< First configuration that will be processed. */
float last_ts;		/*< Last configuration that will be processed. */
ofstream out_average, out_histogram, out_norm_histogram,out_mean_cluster_size;
ofstream out_events,out_kinetics;

int hist[_MAX_AGG_SIZE];	/*< aggregate's distribution histogram: hist[n]==#aggregates of n+1 particles. */

time_t start,end;
struct tm * timeinfo;
double elapsed_time;

///// FUNCTION Prototypes/////

void Initialize();
void GetConnections();
void CountAggregates();
void GetSize(CBead* bead,int i, int &size, int aggid);
void Kinetics();
void ShowOutput();

int Event_11_2(int i);
int Event_2_11(int i);
int Event_3_111(int i);
int Event_111_3(int i);
int Event_3_12(int i);
int Event_12_3(int i);

bool SearchAggregateInAggregateList(CAggregate target,vector <CAggregate> aggregate);
bool SearchBeadInBeadList(int id,vector <CBead> beadlist);

double AvoidZeroDivision(int n, int d);

void DbgInitialAggId(bool debug);
void DbgNeighborList(bool debug);
void DbgAggregatesList(bool debug);
void DbgPrintAggregatesList(vector <CAggregate> aggregate);
//void DbgPrintSelectedBead(vector <CAggregate> aggregate,int id);

/////////////// MAIN FUNCTION //////////
int main(int argc, char * const argv[]){

time (&start);
		
	Initialize();

	dump.Open();

	while(!dump.Load()){

		if((double(dump.header.timestep) >= first_ts && double(dump.header.timestep) <= last_ts) || first_ts==-1 ){				
				
			// inicialitzem l'array de beads
//			dump.Print();
			//
			
			bead = new CBead[dump.header.natoms];
					
			for (int i=0; i<dump.header.natoms; i++){
	
				bead[i].id=dump.particle[i].id;
				bead[i].aggid=0;
				
				
//				cout << "id: " << bead[i].id << endl;
			}

//			ofstream logout;
//			logout.open("log");
//			dump.Write(logout);


//			DbgInitialAggId(_DEBUG);
			GetConnections();
//			DbgInitialAggId(_DEBUG);

//			DbgNeighborList(_DEBUG);
			CountAggregates();
//			DbgNeighborList(_DEBUG);		

			Kinetics();
		
			ShowOutput();

//			DbgAggregatesList(_DEBUG);
			
			// Reiniciem les variables a calcular a per a cada configuracio
			for(int h=0; h<_MAX_AGG_SIZE;h++){ hist[h]=0; }
			
			total_beads=0;
			total_agg=0;
		
			delete [] bead;

		}
	
	}

	dump.Close();
	
	out_average.close();
	out_histogram.close();
	out_norm_histogram.close();
	out_mean_cluster_size.close();

	out_events.close();
	out_kinetics.close();
	
	time (&end);
	elapsed_time = difftime (end,start);
	cout << "Elapsed time " << elapsed_time << " second(s)" << endl;
	
	timeinfo= localtime ( &start );		
	cout << "process started on " << asctime (timeinfo);
	timeinfo = localtime ( &end );
	cout << "process finished on " << asctime (timeinfo);
}
/////////////// END MAIN //////////

////////////////////////////////////////////////////////////////////////////////
/** \fn void Initialize() 
* \brief Reads STDIN and stores initial values for global variables.
*
* All global variables, input/output streams are initialized here.
*/
////////////////////////////////////////////////////////////////////////////////
void Initialize(){

	total_beads=0;
	total_agg=0;

	cout << "Insert contact distance: " << endl;
	cin >> contact;
	cout << "Insert time interval (in seconds) between configurations: " << endl;
	cin >> dt;

	cout << "First timestep (type -1 to process all canfigurations): " << endl;
	cin >> first_ts;

	if (first_ts!=-1){
		cout << "Last timestep: " << endl;
		cin >> last_ts;
	}


	out_average.open("stdmean.txt");
	out_histogram.open("hist.txt");
	out_norm_histogram.open("norm.txt");
	out_mean_cluster_size.open("s-mean.txt");

	out_events.open("events.txt");	
	out_kinetics.open("kinetics.txt");
	
	out_average << "timestep\t " << "average_size" << endl;
	out_histogram << "timestep\t";
	out_norm_histogram << "timestep\t";
	out_mean_cluster_size << "timestep\t" << "S(t)" << endl;
	
	for(int i=0; i<_MAX_AGG_SIZE;i++){

		hist[i]=0;
		out_histogram << "\t" << "aggregate-" << i+1 ;
		out_norm_histogram << "\t" << "aggregate-" << i+1 ;
	}

	out_histogram << endl;
	out_norm_histogram << endl;
	
	// Set output format

	out_kinetics.setf(ios::scientific,ios::floatfield);
	out_kinetics.precision(12);

	out_events.setf(ios::scientific,ios::floatfield);
	out_events.precision(12);

	out_histogram.setf(ios::scientific,ios::floatfield);
	out_histogram.precision(12);

	out_norm_histogram.setf(ios::scientific,ios::floatfield);
	out_norm_histogram.precision(12);
	
	out_mean_cluster_size.setf(ios::scientific,ios::floatfield);
	out_mean_cluster_size.precision(12);
	
	out_average.setf(ios::scientific,ios::floatfield);
	out_average.precision(12);

}
////////////////////////////////////////////////////////////////////////////////
/** \fn void GetConnections() 
* \brief Creates a complete list of neighbors.
*
* This function searches all the neighbors of each bead in the dump file from LAMMPS.
* Two beads are neighbors if their distance is less than the contact distance defined.
* Each bead can have unlimited neighbors.
*/
////////////////////////////////////////////////////////////////////////////////
void GetConnections(){

	double xbox,ybox,zbox;
	double r=0.0;	// Distance between beads
	double dist[3];		// Vector distance

	xbox=(dump.header.xhi-dump.header.xlo);
	ybox=(dump.header.yhi-dump.header.ylo);
	zbox=(dump.header.zhi-dump.header.zlo);

		for(int i=0; i<dump.header.natoms-1; i++){

//			cout << "particle "<< i << " has:" << endl;

			for(int j=i+1; j<dump.header.natoms; j++){

				dist[0]=dump.particle[j].r[0]-dump.particle[i].r[0];
				dist[1]=dump.particle[j].r[1]-dump.particle[i].r[1];
				dist[2]=dump.particle[j].r[2]-dump.particle[i].r[2];				
			
				//MINIMUM IMAGE CONVENTION (Computer Simulation of Liquids. pag30)

				double ix,iy,iz;

				modf(dist[0]/(xbox/2.0),&ix);
				modf(dist[1]/(ybox/2.0),&iy);
				modf(dist[2]/(zbox/2.0),&iz);

				dist[0]=dist[0]-xbox*ix;
				dist[1]=dist[1]-ybox*iy;
				dist[2]=dist[2]-zbox*iz;

				r=sqrt(pow(dist[0],2.0)+pow(dist[1],2.0)+pow(dist[2],2.0));

				if (r< contact){
											
							bead[i].AddNeighbor(&bead[j]);
							bead[j].AddNeighbor(&bead[i]);

				}

			} // close 'j' iteration
		} // close 'i' iteration
}

////////////////////////////////////////////////////////////////////////////////
/** \fn void CountAggregates() 
* \brief Creates a list of aggregates once the bead list has been created.
*
* Starting by the first bead, an aggregate is created by checking all the neighbors
* of this bead, and all neighbors are searched recursively. Once one bead belongs to
* an aggregate, the bead is check as a used neighbor and cannot be included in any other aggregate.
*/
////////////////////////////////////////////////////////////////////////////////
void CountAggregates(){

		int aggid=1;
		
		// per a tots els beads...
		for (int i=0; i<dump.header.natoms; i++){

		// reservem memoria pels aggregats
			aggregate.reserve(dump.header.natoms);
			prev_aggregate.reserve(dump.header.natoms);

		//...que encara no formin part d'un agregat
			if (bead[i].aggid==0){
			int size=1;
			bead[i].aggid=aggid;

		// Si trobem un bead sense aggregar (aggid=0) creem un aggreagat nou i l'incloem en la llista de bead de l'aggregat
			aggregate.push_back(aggid);
			aggregate[aggid-1].AddBead(bead[i]);
			
			GetSize(&bead[i],i,size,aggid);		
		
//if(_DEBUG){	cout << "Particle " << i << " with bead number " << bead[i].id << " has length: " << size << endl; }

//cout << "Size: " << size << endl;

//			cout << "size: " << size << endl;
//			cout << "hist" << hist[size-1] << endl;
			hist[size-1]++;
			aggid++;

			}
		}
		
// aggid-1 es el nombre 'aggregats que hem contabilitzat
// cout << "Nombre aggregats:" << aggid-1 << endl;		

//DbgPrintAggregatesList(aggregate);

}

////////////////////////////////////////////////////////////////////////////////
/** \fn void GetSize(CBead *bead,int id, int &size, int aggid) 
* \brief Searches for all common neighbors belonging to the same aggregate.
*
* \param bead Pointer to the next bead to check.
* \param id Bead id being checked.
* \param size Actual size of the aggregate being checked.
* \param aggid Aggregate id that will be assigned to any bead belonging to the aggregate.
* This recursive function is used to get all the beads belonging to the same neighbor
* and to set is total size. Once finished, a unique id is set to the aggregate made.
*/
////////////////////////////////////////////////////////////////////////////////

void GetSize(CBead *bead, int id, int &size, int aggid){

//cout << "Number of neighbors: " << bead->neighbor.size() << endl;
	for(int k=0; k< bead->neighbor.size(); k++){

//cout << "Aggregate id: " << bead->neighbor[k]->aggid << endl;

		if(bead->neighbor[k]->aggid==0 ){
		
			bead->neighbor[k]->aggid=aggid;
			size++;
			
			aggregate[aggid-1].AddBead(*bead->neighbor[k]);
			
			GetSize(bead->neighbor[k],k,size,aggid);
				
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
/** \fn void Kinetics() 
* \brief Calculate the kinetic constants
*
* Any kinetic constant to be calculated is placed here. The results are stored in
* two different files named "events.txt" and "kinetics.txt"
* ATTENTION: we do need the histogram in order to compute the kinetic constants and
* histogram is computed in CountAggregates function.
*/
////////////////////////////////////////////////////////////////////////////////
void Kinetics(){

int events_11_2=0;
int events_2_11=0;
int events_12_3=0;
int events_3_12=0;
int events_111_3=0;
int events_3_111=0;

		if(!prev_aggregate.empty()){

//		cout << "pev_agg not empty" << endl;

			// Creacio d'aggregats
			for (int i=0;i<aggregate.size();i++){

				events_11_2 += Event_11_2(i);
				events_111_3 += Event_111_3(i);
				events_12_3 += Event_12_3(i);
			}
			
			// Destruccio d'aggregats
			for (int i=0;i<prev_aggregate.size();i++){
			
				events_2_11 += Event_2_11(i);
				events_3_111 += Event_3_111(i);
				events_3_12 += Event_3_12(i);
			}

			cout << "event:1+1->2 " << events_11_2  << " times" << endl;
			cout << "event:2->1+1 " << events_2_11  << " times" << endl;
			cout << "event:1+1+1->3 " << events_111_3  << " times" << endl;
			cout << "event:3->1+1+1 " << events_3_111  << " times" << endl;
			cout << "event:1+2->3 " << events_12_3  << " times" << endl;
			cout << "event:3->1+2 " << events_3_12  << " times" << endl;
			
			// Print events to file
			out_events << dump.header.timestep-_ts_offset << "\t";
			out_events << events_11_2 << "\t";
			out_events << events_2_11 << "\t";
			out_events << events_111_3 << "\t";
			out_events << events_3_111 << "\t";
			out_events << events_12_3 << "\t";
			out_events << events_3_12 << endl;

			double k_11_2=AvoidZeroDivision(events_11_2,(hist[0]*hist[0]));
			double k_2_11=AvoidZeroDivision(events_2_11,hist[1]);
			double k_111_3=AvoidZeroDivision(events_111_3,hist[0]*hist[0]*hist[0]);
			double k_3_111=AvoidZeroDivision(events_3_111,hist[2]);
			double k_12_3=AvoidZeroDivision(events_12_3,(hist[0]*hist[1]));
			double k_3_12=AvoidZeroDivision(events_3_12,hist[2]);

			// Print kinetic constants to file
			out_kinetics << dump.header.timestep-_ts_offset << "\t";
			out_kinetics << k_11_2/dt << "\t";
			out_kinetics << k_2_11/dt << "\t";
			out_kinetics << k_111_3/dt << "\t";
			out_kinetics << k_3_111/dt << "\t";
			out_kinetics << k_12_3/dt << "\t";
			out_kinetics << k_3_12/dt << endl;

		} else {
		
			out_events << "timestep\t" << "e(1+1->2)\t" << "e(2->1+1)\t" << "e(1+1+1->3)\t" << "e(3->1+1+1)\t" << "e(1+2->3)\t" << "e(3->1+2)\t" << endl;
			out_kinetics << "timestep\t" << "k(1+1->2)\t" << "k(2->1+1)\t" << "k(1+1+1->3)\t" << "k(3->1+1+1)\t" << "k(1+2->3)\t" << "k(3->1+2)\t" << endl;
						
		}

// Store list of aggregates as previous list
		prev_aggregate=aggregate;
		
// Erase aggregate list
		aggregate.erase (aggregate.begin(),aggregate.begin()+aggregate.size()-1);

}


////////////////////////////////////////////////////////////////////////////////
/** \fn int Event_11_2(int i) 
* \brief Calculate number of events for the process 1+1 -> 2.
*
* \param i Index of the aggregate that could be involved in the process (as a product).
* \return Number of events found.
*/
////////////////////////////////////////////////////////////////////////////////
int Event_11_2(int i){

int e_11_2=0;

//for (int i=0;i<aggregate.size();i++){

	if (aggregate[i].size==2){

// We split the aggregate in one set 1+1+1

		CAggregate aggsize1[2];

		// we DO initialze the auxiliar aggregate id to 0
		
		aggsize1[0].id=0;
		aggsize1[1].id=0;

		aggsize1[0].AddBead(aggregate[i].bead[0]);
		aggsize1[1].AddBead(aggregate[i].bead[1]);

//		cout << "auxiliar aggregates: " << endl;
//		aggsize1[0].Print();
//		aggsize1[1].Print();

		if (SearchAggregateInAggregateList(aggsize1[0],prev_aggregate) && SearchAggregateInAggregateList(aggsize1[1],prev_aggregate)){
		
			e_11_2++;

		}

	}
//}
	return e_11_2;
}


////////////////////////////////////////////////////////////////////////////////
/** \fn int Event_2_11(int i) 
* \brief Calculate number of events for the process 2 -> 1+1.
*
* \param i Index of the aggregate that could be involved in the process (as a reactant).
* \return Number of events found.
*/
////////////////////////////////////////////////////////////////////////////////
int Event_2_11(int i){

int e_2_11=0;

//for (int i=0;i<prev_aggregate.size();i++){

	if (prev_aggregate[i].size==2){

// We split the aggregate in one set 1+1+1

		CAggregate aggsize1[2];

		// we DO initialze the auxiliar aggregate id to 0
		
		aggsize1[0].id=0;
		aggsize1[1].id=0;

		aggsize1[0].AddBead(prev_aggregate[i].bead[0]);
		aggsize1[1].AddBead(prev_aggregate[i].bead[1]);
				
		if (SearchAggregateInAggregateList(aggsize1[0],aggregate) && SearchAggregateInAggregateList(aggsize1[1],aggregate)){
		
			e_2_11++;

		}

	}
//}
	return e_2_11;
}

////////////////////////////////////////////////////////////////////////////////
/** \fn int Event_111_3(int i) 
* \brief Calculate number of events for the process 1+1+1 -> 3.
*
* \param i Index of the aggregate that could be involved in the process (as a product).
* \return Number of events found.
*/
////////////////////////////////////////////////////////////////////////////////
int Event_111_3(int i){

int e_111_3=0;

//for (int i=0;i<aggregate.size();i++){

	if (aggregate[i].size==3){

// We split the aggregate in one set 1+1+1

		CAggregate aggsize1[3];

		// we DO initialze the auxiliar aggregate id to 0
		
		aggsize1[0].id=0;
		aggsize1[1].id=0;
		aggsize1[2].id=0;
		
		aggsize1[0].AddBead(aggregate[i].bead[0]);
		aggsize1[1].AddBead(aggregate[i].bead[1]);
		aggsize1[2].AddBead(aggregate[i].bead[2]);
						
		if (SearchAggregateInAggregateList(aggsize1[0],prev_aggregate) && SearchAggregateInAggregateList(aggsize1[1],prev_aggregate) && SearchAggregateInAggregateList(aggsize1[2],prev_aggregate)){
		
			e_111_3++;

		}

	}
//}

	return e_111_3;
}


////////////////////////////////////////////////////////////////////////////////
/** \fn int Event_3_111(int i) 
* \brief Calculate number of events for the process 3 -> 1+1+1.
*
* \param i Index of the aggregate that could be involved in the process (as a reactant).
* \return Number of events found.
*/
////////////////////////////////////////////////////////////////////////////////
int Event_3_111(int i){

int e_3_111=0;

//for (int i=0;i<prev_aggregate.size();i++){

	if (prev_aggregate[i].size==3){

// We split the aggregate in one set 1+1+1

		CAggregate aggsize1[3];

		// we DO initialze the auxiliar aggregate id to 0
		
		aggsize1[0].id=0;
		aggsize1[1].id=0;
		aggsize1[2].id=0;

		aggsize1[0].AddBead(prev_aggregate[i].bead[0]);
		aggsize1[1].AddBead(prev_aggregate[i].bead[1]);
		aggsize1[2].AddBead(prev_aggregate[i].bead[2]);
				
		if (SearchAggregateInAggregateList(aggsize1[0],aggregate) && SearchAggregateInAggregateList(aggsize1[1],aggregate) && SearchAggregateInAggregateList(aggsize1[2],aggregate)){
		
			e_3_111++;
			
		}

	}
//}

	return e_3_111;
}


////////////////////////////////////////////////////////////////////////////////
/** \fn int Event_12_3(int i) 
* \brief Calculate number of events for the process 1+2 -> 3.
*
* \param i Index of the aggregate that could be involved in the process (as a product).
* \return Number of events found.
*/
////////////////////////////////////////////////////////////////////////////////
int Event_12_3(int i){
// Calculate the rate for the processes 3 -> 1+2

int e_12_3=0;

//for (int i=0;i<aggregate.size();i++){

	if (aggregate[i].size==3){

// We split the aggregate in three sets of 1+2 aggregates

		CAggregate aggsize1[3];
		CAggregate aggsize2[3];

		// we DO initialze the auxiliar aggregate id to 0
		
		aggsize1[0].id=0;
		aggsize1[1].id=0;
		aggsize1[2].id=0;

		aggsize2[0].id=0;
		aggsize2[1].id=0;
		aggsize2[2].id=0;
				
		aggsize1[0].AddBead(aggregate[i].bead[0]);
		aggsize1[1].AddBead(aggregate[i].bead[1]);
		aggsize1[2].AddBead(aggregate[i].bead[2]);


		aggsize2[0].AddBead(aggregate[i].bead[1]);
		aggsize2[0].AddBead(aggregate[i].bead[2]);
		aggsize2[1].AddBead(aggregate[i].bead[0]);
		aggsize2[1].AddBead(aggregate[i].bead[2]);
		aggsize2[2].AddBead(aggregate[i].bead[0]);
		aggsize2[2].AddBead(aggregate[i].bead[1]);

		for(int k=0; k<3;k++){
		
			if (SearchAggregateInAggregateList(aggsize1[k],prev_aggregate) && SearchAggregateInAggregateList(aggsize2[k],prev_aggregate)){

				e_12_3++;
		
			}
		}
		
	}
//}

	return e_12_3;
}


////////////////////////////////////////////////////////////////////////////////
/** \fn int Event_3_12(int i) 
* \brief Calculate number of events for the process 3 -> 1+2.
*
* \param i Index of the aggregate that could be involved in the process (as a reactant).
* \return Number of events found.
*/
////////////////////////////////////////////////////////////////////////////////
int Event_3_12(int i){

int e_3_12=0;

//for (int i=0;i<prev_aggregate.size();i++){

	if (prev_aggregate[i].size==3){

// We split the aggregate in three sets of 1+2 aggregates

		CAggregate aggsize1[3];
		CAggregate aggsize2[3];

		// we DO initialze the auxiliar aggregate id to 0
		
		aggsize1[0].id=0;
		aggsize1[1].id=0;
		aggsize1[2].id=0;

		aggsize2[0].id=0;
		aggsize2[1].id=0;
		aggsize2[2].id=0;
				
		aggsize1[0].AddBead(prev_aggregate[i].bead[0]);
		aggsize1[1].AddBead(prev_aggregate[i].bead[1]);
		aggsize1[2].AddBead(prev_aggregate[i].bead[2]);


		aggsize2[0].AddBead(prev_aggregate[i].bead[1]);
		aggsize2[0].AddBead(prev_aggregate[i].bead[2]);
		aggsize2[1].AddBead(prev_aggregate[i].bead[0]);
		aggsize2[1].AddBead(prev_aggregate[i].bead[2]);
		aggsize2[2].AddBead(prev_aggregate[i].bead[0]);
		aggsize2[2].AddBead(prev_aggregate[i].bead[1]);

		for(int k=0; k<3;k++){
		
			if (SearchAggregateInAggregateList(aggsize1[k],aggregate) && SearchAggregateInAggregateList(aggsize2[k],aggregate)){

				e_3_12++;
		
			}
		}
		
	}
//}

	return e_3_12;
}


////////////////////////////////////////////////////////////////////////////////
/** \fn bool SearchBeadInBeadList(int id,vector <CBead> beadlist)
* \brief Check if bead with \a id is in the beads list \a beadlist
*
* \param id Target bead.
* \param beadlist List of possible beads.
* \return boolean value.
*/
////////////////////////////////////////////////////////////////////////////////
bool SearchBeadInBeadList(int id,vector <CBead> beadlist){

	bool a=false;

	if (!beadlist.empty()){

	for (int i=0; i<beadlist.size();i++){

//		cout << "SBIBL " << "beadlist size: " << beadlist.size() << endl;
//		cout << "SBIBL " << "beadlist is: " << beadlist[i]->id << endl;
					
		if (beadlist[i].id==id){a=true;}

	}
	}
	return a;

}

////////////////////////////////////////////////////////////////////////////////
/** \fn bool SearchAggregateInAggregateList(CAggregate target,vector <CAggregate> aggregate)
* \brief Check if aggregate \a taget is in the aggregates list \a aggregate
*
* \param target Target aggregate.
* \param aggregate List of possible aggregates.
* \return boolean value.
*/
////////////////////////////////////////////////////////////////////////////////
bool SearchAggregateInAggregateList(CAggregate target,vector <CAggregate> aggregate){

	bool a=false;
	int counter=0;
	
//	cout << "target: " ;
//	target.Print();
	
	for (int i=0; i<aggregate.size();i++){
//	cout << "aggregate in list ";

////////////////////////////aqui peta!!	
//	aggregate[i].Print();
////////////////////////////
//		if (aggregate[i].size == target.size && a==false){
		if (aggregate[i].bead.size() == target.size && a==false){
		
			for (int j=0; j<target.bead.size();j++){

//				cout << "SAIAL: " << "target bead id: " <<  target.bead[j]->id << endl
//				cout << "SAIAL: " << "bead size in aggregate: " <<  aggregate[i].bead.size() << endl;
//				cout << "SAIAL: " << "bead id in aggregate: " <<  aggregate[i].bead[j]->id << endl;
				
				if(SearchBeadInBeadList(target.bead[j].id, aggregate[i].bead)){counter++;}			
			
			}
			
			if(counter==target.bead.size()){
			
				a=true;
			
			}
		
		}
			
	}
		
	return a;

}



////////////////////////////////////////////////////////////////////////////////
/** \fn void ShowOutput()
* \brief Calculates and prints averaged and total results to output files.
*
*/
////////////////////////////////////////////////////////////////////////////////
void ShowOutput(){

// LES PARTICULES INDIVIDUALS SI/NO ES CONSIDEREN AGGREGATS!!!

		cout << "Timestep: " << dump.header.timestep << " Atoms: " << dump.header.natoms << endl;
//		out_average << "Dump Info -> Timestep: " << dump.header.timestep << " Atoms: " << dump.header.natoms << endl;

		out_histogram << dump.header.timestep << " ";
		out_norm_histogram << dump.header.timestep << " ";

// 1- Contem el nombre total de beads i creem l'histograma			
		for(int h=0; h<_MAX_AGG_SIZE ;h++){

			total_beads+=hist[h]*(h+1);
//			total_agg+=hist[h];
//			cout << "aggregate length " << h+1 << " : " << hist[h] << endl;			
//			out_average << h+1 << " " << hist[h] << endl;
			out_histogram << "\t"<< hist[h];
		}

// 2- Contem el nombre total d'agregats (minim 1 beads)
		for(int h=0; h<_MAX_AGG_SIZE ;h++){

			total_agg+=hist[h];

		}


// 3- Normalitzem l'histograma
		//Normalitzed aggregates population
		for(int h=0; h<_MAX_AGG_SIZE ;h++){
		
			out_norm_histogram << "\t"<< double(hist[h])/double(total_agg);
		
		}

		out_histogram << endl;
		out_norm_histogram << endl;

// 4- Calculem el promig standard de la mida dels agregats (minim 1 beads)
		// Average aggregate size
		double ave_agg_size=0;
		
		for(int h=0; h<_MAX_AGG_SIZE ;h++){
		
			ave_agg_size+=hist[h]*(h+1);
		
		}
		
		ave_agg_size=ave_agg_size/total_agg;
		
		cout << "Counted atoms: " << total_beads << endl << "Total aggregates: " << total_agg << endl << "Average aggregate size: " << ave_agg_size << endl;
		cout << " " << endl;
		out_average << dump.header.timestep-_ts_offset << " " << " " << ave_agg_size << endl;

//		out_average << " " << endl;

// 5- Calculem el s-promig
		// Mean Cluster Size at each timestep (M.Carmen Miguel, R.Pastor-Satorras, PRE 59, 826-834, 1999)
		
		double s2n, sn;
		
		s2n=0;
		sn=0;
		
		for(int s=0; s<_MAX_AGG_SIZE ;s++){
		
			s2n+=(s+1)*(s+1)*hist[s];
			sn+=(s+1)*hist[s];
		}
		
		out_mean_cluster_size << dump.header.timestep-_ts_offset << " " << s2n/sn << endl;
			
}

////////////////////////////////////////////////////////////////////////////////
/** \fn double AvoidZeroDivision(int n, int d){
* \brief Calculated the division avoiding non finite results.
*
* \param n numerator.
* \param d denominator.
*/
////////////////////////////////////////////////////////////////////////////////
double AvoidZeroDivision(int n, int d){

double result=0;

	if (d == 0){result=n;} else {result=double(n)/double(d);}

return result;
}

////////// DEBUGGING FUNCTIONS //////////

// Check that aggid of each bead is equal to 0 once has been initialized
void DbgInitialAggId(bool debug){
	
	if(debug){

			cout << "Id check-1..." << endl;
			for (int i=0; i<dump.header.natoms; i++){
			
				if (bead[i].aggid!=0){cout << bead[i].id << " is not 0." << endl;exit(0);}
			
			}
	}
}

// Check that neighbors list has been created successfully
void DbgNeighborList(bool debug){

if(debug){

	for (int i=0; i<dump.header.natoms; i++){
	
		cout << "Particle " << i << " with bead number " << bead[i].id << " and aggregate number " << bead[i].aggid << " is connected to beads: ";
	
		for (int k=0; k< bead[i].neighbor.size(); k++){
		
			cout << bead[i].neighbor[k]->id << " ";
		
		}
		
		cout << endl;
	}
}


}
// Write list of beads in each aggregate
void DbgAggregatesList(bool debug){

// Show aggregates
if(debug){
	
	int i=1;
	
	while (i<= total_agg){
	
		cout << "Aggregate " << i << " composed by: ";
	
		for (int k=0; k<dump.header.natoms; k++){
		
			if (bead[k].aggid==i){cout << bead[k].id << " ";}
					
		}
		
		cout << endl;
		i++;
		
	}
}

}

void DbgPrintAggregatesList(vector <CAggregate> aggregate){
// DEBUG Imprimim la llista d'agregats creats, juntament amb els beads que els formen:

		for (int i=0; i<aggregate.size();i++){
		
			aggregate[i].Print();
		
		}
}


////////// UNUSED FUNCTIONS //////////

/*
// DEBUG Imprimeix l'agregat que conte el bead amb id "id"
void DbgPrintSelectedBead(vector <CAggregate> aggregate,int id){
		
	aggregate[SearchBeadInAggregateList2(id, aggregate)-1].Print();
		
}
*/

/*
// returns the aggregate id if the bead (id) from aggregate vector
int SearchBeadInAggregateList2(int id,vector <CAggregate> aggregate){

	int a=0;

	for (int i=0; i<aggregate.size();i++){
			
		for (int j=0;j<aggregate[i].bead.size();j++){
		
			if (aggregate[i].bead[j]->id==id){a=aggregate[i].id;}
		
		}	
	}
		
	return a;

}
*/
