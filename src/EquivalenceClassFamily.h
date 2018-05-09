#ifndef ADJ_CLASS_FAMILY_H_
#define ADJ_CLASS_FAMILY_H_


/*

This file contains a class contaning all equivalence classes between two families.

Created the: 26-02-2016
by: Wandrille Duchemin

Last modified the: 13-06-2016
by: Wandrille Duchemin

*/

#include <string>
#include <vector>

#include <Bpp/Exceptions.h>

#include "ReconciledTree.h"
#include "AdjMatrix.h"
#include "AdjTree.h"
#include "EquivalenceClass.h"

#include "MultiRootEquivalenceClass.h"

class EquivalenceClassFamily
{
protected:
	int Gfam1;
	int Gfam2;

	//vector <MultiRootEquivalenceClass * > EclassList; 
	vector <MultiRootEquivalenceClass > EclassList; 

public:
	EquivalenceClassFamily()
	{};

	EquivalenceClassFamily(int f1, int f2);

	~EquivalenceClassFamily()
	{
		//cout <<"destroy ECF"<<endl;
		EclassList.clear();
	}

	
	int getGfamily1() const;
	int getGfamily2() const;

	void setGfamily1( int fam1);
	void setGfamily2( int fam2);

	int getNbEqClasses() const;
	int getNbAdj(); // accross all eq classes
	int getNbAdjTrees();// accross all eq classes ; check if the adj trees have been computed

	int getNbAdjGain();// accross all eq classes ; check if the adj trees have been computed
	int getNbAdjBreak();// accross all eq classes ; check if the adj trees have been computed

	bool areSetAdjMatrix();
	bool areSetAdjForest();

	//hasers and finders
	bool hasLname1(string name);
	bool hasLname2(string name);
	pair <int, int> findLname1(string name);// returns -1 if the name is absent
	pair <int, int> findLname2(string name);// returns -1 if the name is absent

	//adders with checks on gfamily names
	bool CheckAddAdj(string name1, string name2, int fam1, int fam2); // adds a pair of Leaf names if fam1 and fam2 correspond to valid Gfamily names 
	bool CheckAddAdj(pair <string, string> names, pair <int, int> fams); // adds a pair of Leaf names if fams correspond to valid Gfamily names 
	bool CheckAddAdjList(vector <string> nameList1, vector<string> nameList2, pair <int, int> fams);
	bool CheckAddAdjList(vector <pair <string,string> > nameList, pair <int, int> fams); // adds several pair of leaf names


	const EquivalenceClass * getEClass(int i) const;
	vector <AdjTree * > * getAdjForest(int i);

	void setAdjForest(int i, vector <AdjTree * > * newAForest);

	void reset();

	void refine(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool useWholeClass, bool forceLTrefining, bool verbose);


	///ADJMatrixes
	void createAdjMatrix(map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, double Gcost, double Bcost, ReconciledTree * rtree1, ReconciledTree * rtree2, bool VERBOSE, bool boltzmann , double temp = 1 , double absencePenalty = -1);
	void computeAdjMatrix();
	void computeAdjMatrix(double WDupCost, double WLossCost, double WHgtCost);

	bool iscomputedAdjMatrix();
	void backtrackAdjMatrix(ReconciledTree * rtree1, ReconciledTree * rtree2, vector < vector< AdjTree *> * > * AdjacencyTreesVector, bool stochastic, bool alwaysGainAtTop = true, double c1proba = 0.5); // creates the adj forest

	void backtrackAdjMatrixForSelf(ReconciledTree * rtree1, ReconciledTree * rtree2, bool stochastic, bool alwaysGainAtTop = true, double c1proba = 0.5); // creates the adj forest
		

	void printMe(bool verbose);


	void clone(const EquivalenceClassFamily * ECF);

	EquivalenceClassFamily& operator=(const EquivalenceClassFamily& rhs)
	{
		//cout << "equal operator - EquivalenceClassFamily"<<endl;

		EclassList.clear();		

		Gfam1 = rhs.getGfamily1();
		Gfam2 = rhs.getGfamily2();
	
		for(unsigned i = 0; i  < rhs.getNbEqClasses();i++)
		{
			EclassList.push_back(MultiRootEquivalenceClass(Gfam1,Gfam2));
			EclassList.back().clone( rhs.getEClass(i) );
		}

	

		return *this;
	}

	void dumpStuff()
	{
		for(unsigned i = 0; i  < getNbEqClasses();i++)
		{
			EclassList[0].dumpStuff();
		}

	};

};


#endif