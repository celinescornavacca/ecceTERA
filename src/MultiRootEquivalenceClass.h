#ifndef MULTIROOTADJ_CLASS_H_
#define MULTIROOTADJ_CLASS_H_


/*

This file contains a class for adjacency classes

Created the: 17-04-2016
by: Wandrille Duchemin

Last modified the: 17-04-2016
by: Wandrille Duchemin

*/

#include "EquivalenceClass.h"
#include "MultiRootAdjMatrix.h"


class MultiRootEquivalenceClass : public EquivalenceClass
{
//protected:
//	MultiRootAdjMatrix * Amat;

public:
	MultiRootEquivalenceClass(int fam1, int fam2) : EquivalenceClass(fam1, fam2)
	{//cout << "MultiRootEquivalenceClass creation" << endl;
	};

	~MultiRootEquivalenceClass()
	{ //cout << "plop MREC"<<endl;
	};

	vector<EquivalenceClass *> refineEqClass(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool forceLTrefining, bool verbose);


	vector<EquivalenceClass *> refineEqClassWhole(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool verbose);


	//void createAdjMatrix(double Gcost, double Bcost, ReconciledTree * rtree1, ReconciledTree * rtree2, bool VERBOSE, bool boltzmann , double temp, double absencePenalty );
};

#endif