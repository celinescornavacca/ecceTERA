#ifndef MULTIROOT_ADJ_MATRIX_H_
#define MULTIROOT_ADJ_MATRIX_H_


/*

This file contains a class for a matrix for adjacency history

Created the: 17-04-2016
by: Wandrille Duchemin

Last modified the: 17-04-2016
by: Wandrille Duchemin

*/

#include "AdjMatrix.h"


class MultiRootAdjMatrix: public AdjMatrix
{
protected:
	MyMatrix RootMatrix;
	vector < pair<int,int> > rootList;
	bool rootComputed;

	void backtrackAux( vector< AdjTree *> * AdjacencyTrees, bool stochastic, bool alwaysGainAtTop , double c1proba , int Root1, int Root2);


	bool isImpossibleCase(vector<AdjSolution> aSolC1, vector<AdjSolution> aSolC0);

	void initRootMatrix();
	void scanSolutionForRooting(int id1, int id2);

	void scanForRoots();
	vector < pair<int,int> > getMatrixRoots();

public:
	MultiRootAdjMatrix(double Gcost, double Bcost, 
								ReconciledTree * rtree1, ReconciledTree * rtree2, 
								vector< pair <string,string> > adjacencies, 
								bool VERBOSE, bool boltzmann = false, double temp = 1 , double absencePenalty = -1);

	MultiRootAdjMatrix(double Gcost, double Bcost, 
							ReconciledTree * rtree1, ReconciledTree * rtree2, 
							vector< pair <int,int> > adjacencies, 
							bool VERBOSE, bool boltzmann = false, double temp = 1  , double absencePenalty = -1);

	~MultiRootAdjMatrix()
	{
		//cout << "ploup" << endl;
	};

	bool getRootMatrix(int id1, int id2);
	void setRootMatrixPossible(int id1, int id2);
	void setRootMatrixNotPossible(int id1, int id2);

	void backtrack( vector< AdjTree *> * AdjacencyTrees, bool stochastic, bool alwaysGainAtTop = true, double c1proba = 0.5 );

};

#endif