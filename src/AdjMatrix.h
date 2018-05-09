#ifndef ADJ_MATRIX_H_
#define ADJ_MATRIX_H_


/*

This file contains a class for a matrix for adjacency history

Created the: 30-11-2015
by: Wandrille Duchemin

Last modified the: 12-06-2016
by: Yoann Anselmetti

*/

#define MYINFINIY numeric_limits<double>::max();
#define MYMINFINIY -numeric_limits<double>::max();


#include <stdlib.h>

#include "MyMatrix.h"
#include "ReconciledTree.h"
#include "AdjTree.h"
#include "MyMatrixAdjSolution.h"

const double BOLTZMANN_K = 1; // not true but for now this will do

using namespace std;
using namespace bpp;

class AdjMatrix
{
protected:

	///////////////WMODIF
	MyMatrix RootMatrix;
	vector < pair<int,int> > rootList;
	bool rootComputed;

	void backtrackAux( vector< AdjTree *> * AdjacencyTrees, bool stochastic, bool alwaysGainAtTop , double c1proba , int Root1, int Root2);


	bool isImpossibleCase(vector<AdjSolution> aSolC1, vector<AdjSolution> aSolC0);

	void initRootMatrix();
	void scanSolutionForRooting(int id1, int id2);

	void scanForRoots();
	vector < pair<int,int> > getMatrixRoots();
	///////////////WMODIF


	double GainCost; // cost of a single gain
	double BreakCost; // cost of a single break

	double WeightedDupCost ;//Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score. This is defaulting at 0 and set only if requested when computing the matrix. Used to subtract a reconciliation event to the adjacency scores of co-events.
	double WeightedLossCost;//Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score. This is defaulting at 0 and set only if requested when computing the matrix. Used to subtract a reconciliation event to the adjacency scores of co-events.
	double WeightedHgtCost ;//Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score. This is defaulting at 0 and set only if requested when computing the matrix. Used to subtract a reconciliation event to the adjacency scores of co-events.	

	double worstScore; // worst possible score (infinity for the DeCo case; 0 for the DeCoBoltzmann case)
	double bestScore; // best possible score (0 for the DeCo case; 1 for the DeCoBoltzmann case)
	double defaultScore; // default score for cases that haven't been computed yet ( -1 in both case)
	double worstAbsenceScore; // worst possible score for an absent adjacency presence( == worstScore unless some option is set at the AdjMatrix creation)

	map<int,vector<float> > speciesC0C1;
	map<int, map<string,int> > speGeneAdjNb;

	ReconciledTree Rtree1;
	ReconciledTree Rtree2;


	MyMatrix MatrixC0;
	MyMatrix MatrixC1;

	MyMatrixAdjSolution SolutionMatrixC0;
	MyMatrixAdjSolution SolutionMatrixC1;

	//matrix ids go from 0 to Dim without interuption; while tree ids have no particular obligations baout that -> we need map linking one to the other

	map <int,int> TreeToMatrixId1;
	map <int,int> TreeToMatrixId2;

	map <int,int> MatrixToTreeId1;
	map <int,int> MatrixToTreeId2;

	bool matrixComputed;

	bool verbose;

	bool useBoltzmann;
	double Temperature;

	bool decoLTalgo;

//// methods that will be used by the score algebra
	double addition(double const& a, double const& b);
	double multiplication(double const& a, double const& b);

	double getminimum(vector <double> const& v);
	double getsum(vector <double> const& v);

/////pointers to these methods
	typedef double (AdjMatrix::*ScoreAggregator)(double const& a,double const& b);
	typedef double (AdjMatrix::*ScoreComparator)(vector <double> const& v);  	

	ScoreAggregator scoreAggregatorfunc;
	ScoreComparator scoreComparatorfunc;


///methods
	void AdjMatrixAux(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb,
						 double Gcost, double Bcost, 
						 ReconciledTree * rtree1, ReconciledTree * rtree2, 
						 vector< pair <int,int> > adjacencies, 
						 bool VERBOSE, bool boltzmann =false, double temp = 1, double absencePenalty = -1);


	void BuildIdMaps();

	void addAdjacency(pair <int,int> adjacency);


	void initMatrix();
	void initAdjAbsence(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb);
	void initMatrix(vector< pair <int,int> > adjacencies, map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb);


	//void setC0(int id1, int id2, double value);
	//void setC1(int id1, int id2, double value);

	void resetCase(int id1, int id2);

///// Boltzmann specific function /////
	void setComputationToBoltzmann();

//// Backtrack Functions ////
	AdjSolution chooseRandomSolution(vector <AdjSolution> Vsolution);
	AdjSolution chooseBestSolution(vector <AdjSolution> Vsolution);
	AdjSolution chooseRandomSolutionWeighted(vector <AdjSolution> Vsolution);

	Node * backtrackAuxC1(int id1, int id2, vector< AdjTree *> * AdjacencyTrees, bool stochastic);
	void backtrackAuxC0(int id1, int id2, vector< AdjTree *> * AdjacencyTrees, bool stochastic);
	pair <int,int> getAdjEventAndSpecies(int NodeId1, int NodeId2);

///////////////////////////////////
/////// Cost functions //////////// in another file
///////////////////////////////////


/////////////////////////////////// OLD FUNCTIONS ////////////////////////////////////////////////////////

	double AggregateScore(vector <double> scores, int nbGain = 0, int nbBreak = 0);
/*
	pair <double, double> computeScore(int NodeId1, int NodeId2);

// simple DeCo cases
	double C1ExtantWithExtant(int NodeId1, int NodeId2);
	double C0ExtantWithExtant(int NodeId1, int NodeId2);

	double C1LossWithOther(int NodeId1, int NodeId2);//Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
	double C0LossWithOther(int NodeId1, int NodeId2);//Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss

	double C1LossWithLoss(int NodeId1, int NodeId2); //add something later to account for possible co-events
	double C0LossWithLoss(int NodeId1, int NodeId2); //add something later to account for possible co-events


	double D1(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert ); //auxiliary of C1DupWithOther 
	double D0(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert ); //auxiliary of C0DupWithOther 
	double D12(int NodeId1, int NodeId2);

	double C1DupWithOther(int NodeId1, int NodeId2, bool firstDup );
	double C0DupWithOther(int NodeId1, int NodeId2, bool firstDup );

	double C1SpecWithSpec(int NodeId1, int NodeId2);
	double C0SpecWithSpec(int NodeId1, int NodeId2);

	double C1DupWithDup(int NodeId1, int NodeId2); //add something later to account for possible co-events
	double C0DupWithDup(int NodeId1, int NodeId2); //add something later to account for possible co-events


//DeCoLT specific cases

	double C1NullWithNull(int NodeId1, int NodeId2);
	double C0NullWithNull(int NodeId1, int NodeId2);

	pair <int,int> CostSoutWithExtantOrSpecOrNullAux(int NodeId1, int NodeId2, bool firstSout);

	double C1SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout);
	double C0SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout);

	double C1SoutWithSoutSynchronous(int NodeId1, int NodeId2);
	double C1SoutWithSoutaSynchronous(int NodeId1, int NodeId2);

	double C1SoutWithSout(int NodeId1, int NodeId2); //add something later to account for possible co-events
	double C0SoutWithSout(int NodeId1, int NodeId2); //add something later to account for possible co-events

	double C1RecWithRec(int NodeId1, int NodeId2); //add something later to account for possible co-events
	double C0RecWithRec(int NodeId1, int NodeId2); //add something later to account for possible co-events

	double B1(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert); // auxiliary of C1BoutWithBout and C1BoutWithOther
	double B0(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert); // auxiliary of C0BoutWithBout and C0BoutWithOther
	double B12(int NodeId1, int NodeId2);// auxiliary of C1BoutWithBout and C1BoutWithOther

	double C1BoutWithBout(int NodeId1, int NodeId2); //add something later to account for possible co-events
	double C0BoutWithBout(int NodeId1, int NodeId2); //add something later to account for possible co-events

	double C1BoutWithRec(int NodeId1, int NodeId2, bool firstBout);
	double C0BoutWithRec(int NodeId1, int NodeId2, bool firstBout);
/////////////////////////////////// END OLD FUNCTIONS ////////////////////////////////////////////////////////
*/

///score functions with AdjSolutions ///
	void computeSolution(int NodeId1, int NodeId2, vector <AdjSolution> &VsolutionC1, vector <AdjSolution> &VsolutionC0 );
	double compareScore( vector <AdjSolution> Vsolution);

	vector<AdjSolution> SolutionC1DefaultImpossibleCase();
	vector <AdjSolution> SolutionC1ExtantWithExtant(int NodeId1, int NodeId2);
	vector <AdjSolution> SolutionC0ExtantWithExtant(int NodeId1, int NodeId2);
	vector <AdjSolution> SolutionC1LossWithOther(int NodeId1, int NodeId2);
	vector <AdjSolution> SolutionC0LossWithOther(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC1LossWithLoss(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0LossWithLoss(int NodeId1, int NodeId2);

	vector<AdjSolution> SolutionC1synchronousOneChild(int NodeId1, int NodeId2, bool inTheDead, bool iscoevent);
	vector<AdjSolution> SolutionC1synchronousTwoChildren(int NodeId1, int NodeId2, bool inTheDead, bool iscoevent, bool isSout);
	vector<AdjSolution> SolutionC1synchronousTwoAndOneChildren(int NodeId1, int NodeId2, bool inTheDead);
	vector<AdjSolution> SolutionC1asynchronousOneChild(int NodeId1, int NodeId2, bool NodeId2First, bool inTheDead, bool isRec);
	vector<AdjSolution> SolutionC1asynchronousTwoChildren(int NodeId1, int NodeId2, bool NodeId2First, bool inTheDead);

	vector<AdjSolution> SolutionC0DefaultImpossibleCase();
	vector<AdjSolution> SolutionC0DefaultImpossibleCaseNew(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0DefaultImpossibleCaseAux(int NodeId1, int NodeId2, bool invert);

	vector<AdjSolution> SolutionC0synchronousOneChild(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0asynchronousOneChild(int NodeId1, int NodeId2, bool NodeId2First);
	vector<AdjSolution> SolutionC0synchronousTwoAndOneChildren(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0asynchronousTwoChildren(int NodeId1, int NodeId2, bool NodeId2First);
	vector<AdjSolution> SolutionC0synchronousTwoChildren(int NodeId1, int NodeId2);


	//vector<AdjSolution> SolutionC1General(int NodeId1, int NodeId2);
	//vector<AdjSolution> SolutionC0General(int NodeId1, int NodeId2);

/*
	vector<AdjSolution> SolutionC1ExtantWithExtant(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0ExtantWithExtant(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC1LossWithOther(int NodeId1, int NodeId2);//Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
	vector<AdjSolution> SolutionC0LossWithOther(int NodeId1, int NodeId2);//Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
	vector<AdjSolution> SolutionC1LossWithLoss(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionC0LossWithLoss(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionD1(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert );
	vector<AdjSolution> SolutionD0(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert );
	vector<AdjSolution> SolutionD12(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC1DupWithOther(int NodeId1, int NodeId2, bool firstDup );
	vector<AdjSolution> SolutionC0DupWithOther(int NodeId1, int NodeId2, bool firstDup );
	vector<AdjSolution> SolutionC1SpecWithSpec(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0SpecWithSpec(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC1DupWithDup(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0DupWithDup(int NodeId1, int NodeId2);
	/// decolt ///
	vector<AdjSolution> SolutionC1NullWithNull(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0NullWithNull(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC1SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout);
	vector<AdjSolution> SolutionC0SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout);
	vector<AdjSolution> SolutionC1SoutWithSoutSynchronous(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC1SoutWithSoutaSynchronous(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC1SoutWithSout(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionC0SoutWithSout(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionC1RecWithRec(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionC0RecWithRec(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionC1RecWithOther(int NodeId1, int NodeId2, bool firstRec); 
	vector<AdjSolution> SolutionC0RecWithOther(int NodeId1, int NodeId2, bool firstRec); 
	vector<AdjSolution> SolutionB1(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert); // auxiliary of C1BoutWithBout and C1BoutWithOther
	vector<AdjSolution> SolutionB0(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert); // auxiliary of C0BoutWithBout and C0BoutWithOther
	vector<AdjSolution> SolutionB12(int NodeId1, int NodeId2);// auxiliary of C1BoutWithBout and C1BoutWithOther
	vector<AdjSolution> SolutionC1BoutWithBout(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionC0BoutWithBout(int NodeId1, int NodeId2); //add something later to account for possible co-events
	vector<AdjSolution> SolutionC1BoutWithRec(int NodeId1, int NodeId2, bool firstBout);
	vector<AdjSolution> SolutionC0BoutWithRec(int NodeId1, int NodeId2, bool firstBout);
	vector<AdjSolution> SolutionC1DefaultImpossibleCase();
	vector<AdjSolution> SolutionC0DefaultImpossibleCase();


	vector<AdjSolution> SolutionC1RecWithSout(int NodeId1, int NodeId2, bool firstSout);
	vector<AdjSolution> SolutionC0RecWithSout(int NodeId1, int NodeId2, bool firstSout);
	vector<AdjSolution> SolutionC1DupWithSout(int NodeId1, int NodeId2, bool firstDup );
	vector<AdjSolution> SolutionC0DupWithSout(int NodeId1, int NodeId2, bool firstDup );
	vector<AdjSolution> SolutionC1DupWithRec(int NodeId1, int NodeId2, bool firstDup );
	vector<AdjSolution> SolutionC0DupWithRec(int NodeId1, int NodeId2, bool firstDup );


	vector<AdjSolution> SolutionC0DefaultImpossibleCaseNew(int NodeId1, int NodeId2);
	vector<AdjSolution> SolutionC0DefaultImpossibleCaseAux(int NodeId1, int NodeId2, bool invert);
*/

public:

//	AdjMatrix(double Gcost, double Bcost, bool VERBOSE);
	AdjMatrix(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb,
				double Gcost, double Bcost,
				ReconciledTree * rtree1, ReconciledTree * rtree2, 
				vector< pair <string,string> > adjacencies, 
				bool VERBOSE, bool boltzmann = false, double temp = 1 , double absencePenalty = -1);
	AdjMatrix(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb,
				double Gcost, double Bcost, 
				ReconciledTree * rtree1, ReconciledTree * rtree2, 
				vector< pair <int,int> > adjacencies, 
				bool VERBOSE, bool boltzmann = false, double temp = 1  , double absencePenalty = -1);


	AdjMatrix(){}

	~AdjMatrix()
	{
		//cout << "ploup" << endl;
		/*
		if(TreeToMatrixId1.size() >0)
			TreeToMatrixId1.clear();
		if(TreeToMatrixId2.size() >0)
			TreeToMatrixId2.clear();

		if(MatrixToTreeId1.size() > 0)
			MatrixToTreeId1.clear();
		if(MatrixToTreeId2.size() > 0)
			MatrixToTreeId2.clear();
		*/
		/*
		if(Rtree1 !=NULL)
			delete Rtree1;

		if(Rtree2 != NULL)
			delete Rtree2;
		*/
		//the reconciled tree are referenced by the gene familiy too, so don't delete them <- false because we use cloneSubtree which create a new instance
	}

	void resetMatrix();
	void partiallyResetMatrix(int id1, int id2);

	double getC1(int id1, int id2, bool invert = false);
	double getC0(int id1, int id2, bool invert = false);

	vector <AdjSolution> getSolutionC1(int id1, int id2);
	vector <AdjSolution> getSolutionC0(int id1, int id2);
	void setSolutionC1(int id1, int id2, vector <AdjSolution> &Vsolution);
	void setSolutionC0(int id1, int id2, vector <AdjSolution> &Vsolution);


	bool issetC1(int id1, int id2);
	bool issetC0(int id1, int id2);

	void setdecoLTalgo(bool decolt);

	bool isComputed()
	{return matrixComputed;}

	void setComputed()
	{matrixComputed = true;}

	double getGainCost();
	double getBreakCost();
	double getBoltzmannGainBreakCost(int nbGain, int nbBreak);


	void setRtree1(ReconciledTree * rtree);
	void setRtree2(ReconciledTree * rtree);

	void addAdjacencies(vector < pair <int,int> > adjacencies);

	void computeMatrix();
	void computeMatrix(double WDupCost, double WLossCost, double WHgtCost);


	void printC0();
	void printC1();
	void printMe();

	void backtrack( vector< AdjTree *> * AdjacencyTrees, bool stochastic, bool alwaysGainAtTop = true, double c1proba = 0.5 );

	///////////////WMODIF
	bool getRootMatrix(int id1, int id2);
	void setRootMatrixPossible(int id1, int id2);
	void setRootMatrixNotPossible(int id1, int id2);

	AdjMatrix* getClone();

	void setC0(int id1, int id2, double value);
	void setC1(int id1, int id2, double value);


	///////////////WMODIF



};


#endif