#ifndef DECO_OUTPUT_MANAGER_H_
#define DECO_OUTPUT_MANAGER_H_


/*

This file contains a class that manages DeCo's outputs

Created the: 06-01-2016
by: Wandrille Duchemin

Last modified the: 21-04-2016
by: Wandrille Duchemin

*/

#include <fstream>

#include "ReconciledTree.h"
#include "AdjTree.h"
#include "CoEvent.h"

using namespace std;

class DeCoOutputManager
{
	protected:
		void beginLine(ofstream& OUT, int indent_level);
		void WritePhyloXMLRecEvent(ofstream& OUT, ReconciledTree * Rtree, int nodeid, int indent_level, bool hasLoss);
		void WritePhyloXMLRecTreeAux(ofstream& OUT, ReconciledTree * Rtree, int nodeid, int indent_level);

		void WritePhyloXMLAdjEvent(ofstream& OUT, AdjTree * Atree, int nodeid, int indent_level);
		void WritePhyloXMLAdjTreeAux(ofstream& OUT, AdjTree * Atree, int nodeid, int indent_level);

		void WritePhyloXMLSpeTreeAux(ofstream& OUT, MySpeciesTree * Stree, int nodeid, int indent_level);

		void WriteAdjTreeAdjacencies(ofstream& OUT, AdjTree * ATree);

		void WriteCoEvent(ofstream& OUT, CoEvent coevent);

	public:
		DeCoOutputManager(){}
		~DeCoOutputManager(){}

		void WritePhyloXMLRecTree(ofstream& OUT, ReconciledTree * Rtree);
		void WritePhyloXMLAdjTree(ofstream& OUT, AdjTree * Atree);
		void WritePhyloXMLSpeTree(ofstream& OUT, MySpeciesTree * Stree);
		void WritePhyloXMLAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest);

		void WriteNewickRecTree(ofstream& OUT, ReconciledTree * Rtree, bool hideLosses = false);
		void WriteNewickAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, bool hideLosses = false);

		void WriteRecTree(ofstream& OUT, ReconciledTree * Rtree, bool newick, bool hideLosses = false);
		void WriteAdjForest(ofstream& OUT, vector<AdjTree *> * Aforest, bool newick, bool hideLosses = false);

		void WriteAdjForestAdjacencies(ofstream& OUT, vector<AdjTree *> * Aforest);

		void WriteCoEventSet(ofstream& OUT, vector <CoEvent> * CoEventSet);

};

#endif