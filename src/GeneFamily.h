#ifndef GENE_FAMILY_H_
#define GENE_FAMILY_H_


/*

This file contains a class for Gene Family as a Clade and Tripartition instance associated to a reconciled tree

Created the: 18-11-2015
by: Wandrille Duchemin

Last modified the: 15-06-2016
by: Wandrille Duchemin

*/

#include "MyGeneTree.h"
#include "MySpeciesTree.h"
#include "MyCladesAndTripartitions.h"
#include "ReconciliationEvent.h"
#include "CladeReconciliation.h"
#include "ReconciledTree.h"

#include "CladesAndTripartitions.h"
#include "DTLMatrix.h"
#include "DTLGraph.h"

#include <Bpp/Exceptions.h>

class GeneFamily
{
protected:

	MyCladesAndTripartitions * CCPDistrib;
	MyGeneTree UnrootedTree;
	DTLMatrix * ReconciliationMatrix;
	ReconciledTree RecTree;

	//utility
	char CharSep;
	vector <string> leafNames;

	//scores
	double TreeLikelihood; // distinct from the associated score which is the -log of this normalized by maximum ccp score
	double MaxLkh;
	double RecScore; 

	//booleans about the different states of the Gene Family
	bool CanChangeTree;//implies CanChangeReconciliation. If false, this means that CCPDistrib does not point to a valid object
	bool CanChangeReconciliation;
	bool TreeSet;
	bool ReconciliationSet;//implies TreeSet
	bool DTLMatSet;//implies TreeSet


	//booleans about the reconciliation options
//	bool ultrametricSpeciesTree;
//	bool AllowTransfer;//implies ultrametricSpeciesTree?
//	bool AllowTransferFromDead;//implies AllowTransfer

	//verbosity
	bool VERBOSE;
	bool SUPERVERBOSE;//implies VERBOSE

	//costs? (cost of events, relative cost of ccp and rec components)
	//maps?



	/*
	Takes:
	 - geneTree (MyGeneTree *) : tree to check
	 - errStr (string &): string to get the nature of the error

	Returns:
		(bool): true if the tree is valid for CladesAndTripartions: more than leaves and no identical leaf names
	*/
	bool CheckGeneTree(MyGeneTree * geneTree, string &errStr);

	void drawReconciliation(MySpeciesTree * speciesTree, bool random, int dated);

	void resetReconciliation();
	void resetUnrootedTree();//Also reset theDTLMatrix and ReconciledTree

public:


	//constructors and destructors
	GeneFamily(bool verbose = false, bool superverbose=false)
	{
		VERBOSE = verbose;
		SUPERVERBOSE = superverbose;
	}

	~GeneFamily()
	{
		//cout << "~GeneFamily Utree's nb of nodes: " << UnrootedTree.getNumberOfNodes() << endl;
		delete CCPDistrib;
		if(DTLMatSet)
			delete ReconciliationMatrix;
	}


	GeneFamily(vector<MyGeneTree*> &geneTrees, char charsep,bool verbose, bool superverbose);//with a bunch of trees
	
	GeneFamily(string aleFileName, char charsep,bool verbose, bool superverbose);//ale file reader

	GeneFamily(MyGeneTree UnrootedgeneTree, char charsep, bool verbose, bool superverbose);//only one tree -> won't use a CladesAndTripartitions instance

	GeneFamily(ReconciledTree RTree, bool verbose, bool superverbose);//directly the reconciled tree -> won't use a CladesAndTripartitions instance nor a DTLMatrix instance


	//also set the tree likelihood
	void makeUnrootedTree(bool random); // draw an unrooted tree from the CCPDistrib; at random if random is true; otherwise uses best split
	void setUnrootedTree(MyGeneTree UnrootedgeneTree); // version where the tree is given directly. CAUTION: no particular test are done to ensure that a tree of the correct family was given


	void makeReconciliation(MySpeciesTree * speciesTree, bool computeT, bool computeTL, double DupCost, double HgtCost, double LossCost, int maxTS, double SplitWeight,bool SubdividedSpTree, bool tryAllAmalgamation);
	void setReconciliation(ReconciledTree RTree);

	void changeReconciliation(MySpeciesTree * speciesTree, bool SubdividedSpTree, bool BTS); // just re-draw a tree in the already computed DTLMatrix. Won't work if no DTLMatrix have been computed

	//getter and setters
	int getReconciledTreeTimeSliceStatus(); // 0 if no time slices (NTS) , 1 if precise time slices (TS), 2 if bounded time slices (BTS)
	void setReconciledTreeTimeSlicesToTS(MySpeciesTree * spTree); // set the tree to a complete time sliced one (TS)
	void setReconciledTreeTimeSlicesToBTS(MySpeciesTree * spTree); // set the tree to a bounded time slice one (BTS)
	void setReconciledTreeTimeSlicesToNoTS();

	ReconciledTree * getRecTree();
	ReconciledTree getRecTreeCopy();

	double getTreeLikelihood();
	double getNormalizedTreeLikelihood();
	void setTreeLikelihood(double lkh);
	double getRecScore();
	void setRecScore(double rscore);
	void setRecScore(double DupCost, double HgtCost, double LossCost); // Counts events in the RecTree and uses the costs to compute the score

	void setMaxLikelihood() 
	{ MaxLkh = CCPDistrib->getMaxLikelihood(); };

	void printRecTree();

	bool hasLeaf( string lName );

	void dumpStuff(); // free CCPDistrib and ReconciliationMatrix from the memory



};



#endif