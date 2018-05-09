#ifndef CLADE_RECONCILIATION_H_
#define CLADE_RECONCILIATION_H_



/*

This file contains a class for reconciliated clades

Created the: 28-10-2015
by: Wandrille Duchemin

Last modified the: 10-02-2016
by: Wandrille Duchemin

*/
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <Bpp/Exceptions.h>
#include <fstream>


#include "CladesAndTripartitions.h"
#include "DTLGraph.h"
#include "MySpeciesTree.h"
#include "ReconciliationEvent.h"

#include "XMLUtils.h"

using namespace std;
using namespace bpp;

class CladeReconciliation
{
public:

	int idU;//clade id
	int idUl;//clade id of the left child of idU
	int idUr;//clade id of the right child of idU
	ReconciliationEvent previousEvent;//split event that lead to this clade. Left empty if this is the root clade
	vector<ReconciliationEvent> Events;
	string name;


	bool isRoot();
	bool isLeaf();

	CladeReconciliation(){Events.clear();}

	//constructor for ecceTERA classes
	CladeReconciliation(CladesAndTripartitions * CandT, DTLGraph *graph, MySpeciesTree * speciesTree, int recIndex, vector< vector<DTLGraph::MyGraph::Vertex> > *reconciliation,ReconciliationEvent previous, bool VERBOSE);


	CladeReconciliation(string recline, map<string, int> speciesLeafNameToId, map<string, int> geneLeafNameToId);


	~CladeReconciliation(){}

	void printMe();

};



#endif