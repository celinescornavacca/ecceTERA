#ifndef CLADEFOREST_H_
#define CLADEFOREST_H_


/*

CladeForest represnent a bunch of partial trees. 
Nodes are identified by their cladeId, hence the name.


Created the: 09-06-2016
by: Wandrille Duchemin

Last modified the: 09-06-2016
by: Wandrille Duchemin

*/

#include <Bpp/Exceptions.h>


#include <vector>
#include <map>


#include "ReconciledTree.h"
#include "AdjTree.h"


using namespace bpp;
using namespace std;


struct Clade 
{ // parent is init at -1 and children as an empty vector
	int parent;
	vector <int> children;
};

class CladeForest
{
protected:

	map <int , Clade> cladeMap; //keys are clades id; values are clades

public:
	CladeForest(){};

	~CladeForest(){};

	bool hasClade(int cladeId);
	Clade getClade(int cladeId, bool check = false);

	int getCladeParent(int cladeId, bool check = false);

	int getNbCladeChildren(int cladeId, bool check = false);
	vector <int> getCladeChildren(int cladeId, bool check = false);

	int getNbClades();
	vector <int> getRoots();
	int getNbRoots();

	bool addClade(int nodeId, ReconciledTree * Rtree);
	bool addCladesInAdjTree(AdjTree * Atree, ReconciledTree * Rtree);

	vector <int> getHighestsRoots(ReconciledTree * Rtree);
};


#endif