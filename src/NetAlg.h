#ifndef NETALG_H_
#define NETALG_H_

/**

@file
@author Celine Scornavacca
@author Edwin Jacox

@section LICENCE

Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/

#include <string>
#include <vector>

#include "MyNetwork.h"
#include "MyGeneTree.h"

/**
 */
class NetAlg 
{
private:
    MyNetwork *mSpeciesNetwork;
    MyGeneTree *mGeneTree;
    char mCharSep;
    double mDupliCost;
    double mLossCost;
    double mTransferCost;


    struct NodeVal {
        bool dupe;
        int losses;
        double val;
        int c1nodeId;
        int c2nodeId;
        MySpeciesNode *node;
    };

    vector< vector< NodeVal > > mCLists;

    int countPathLosses( vector<MySpeciesNode*> &path, MySpeciesNode *lca );
    int countPathTransfers( vector<MySpeciesNode*> &path, MySpeciesNode *lca );
    MySpeciesNode *commonAncestor( MySpeciesNode *node1, MySpeciesNode *node2,
             vector<MySpeciesNode*> &spPath1, vector<MySpeciesNode*> &spPath2 );
    double LCAaux( MyGeneNode *geneNode, int &numLosses, 
                   int &numDupli, int &numTransfers, vector<MySpeciesNode*> &recon );
    double LCA( MyNetwork *speciesTree, MyGeneTree *geneTree, int &numLosses,
                 int &numDupli, int &numTransfers, vector<MySpeciesNode*> &recon );
    MySpeciesNode *makeBicompTree( MySpeciesNode *node, MySpeciesNode *parent,
                vector<bool> &seenIt, vector<int> &biCompMap );
    MyGeneNode *makeGnAux( MyGeneNode *node, vector<MySpeciesNode*> &recon, 
              vector<int> &mu, int &cnt, MySpeciesNode *parentBiCompNode=NULL );
    MyGeneTree *makeGn( vector<int> &mu, vector<int> &biCompMap );
    double biCompCost( vector<MyNetwork*> &switchings,
                       vector<MyGeneTree*> &biCompGeneTrees, bool isRootBiComp,
                       int &bestLosses, int &bestDupes, int &bestTransfers, int &bestSwitching );
    void merge( vector<NodeVal> &cList, NodeVal newPair );
    void checkComputeMin( MySpeciesNode *y, MySpeciesNode *z, 
                          vector<MySpeciesNode*> &otherParents,
                          vector<MySpeciesNode*> &toCheck );
    int computeMinAux( MySpeciesNode *node, vector<int> &yDists, 
                       vector<int> &zDists,
                       vector<int> &labels, vector<MySpeciesNode*> &min );
    void processInternalNode( NodeVal c1, NodeVal c2,
                vector<MySpeciesNode*> &otherParents, vector<NodeVal> &cList );
    vector<MySpeciesNode*> computeMin( MySpeciesNode *y, MySpeciesNode *z,
            vector<MySpeciesNode*> &otherParents, vector<int> &yDists, 
            vector<int> &zDists );
    void runMinReconAux( MyGeneNode* geneNode,
            vector<MySpeciesNode*> &otherParents, 
            vector<MySpeciesNode*> &leafMapping );

public:
    /** DTL Graph Constructor */
    NetAlg(
        MyNetwork *speciesNetwork, ///< species tree
        MyGeneTree *geneTree,   ///< gene tree
        char charSep,       ///< character separating gene names from taxa
        double dupliCost,   ///< duplication cost
        double lossCost,   ///< loss cost
        double transferCost )   ///< transfer cost
        : mSpeciesNetwork(speciesNetwork), mGeneTree(geneTree), 
          mCharSep(charSep), mDupliCost(dupliCost), mLossCost(lossCost), mTransferCost(transferCost)
    {}

    double runBrute( int &numLosses, int &numDupli,int &numTransfers);
    double runMinSwitch( int &numLosses, int &numDupli,int &numTransfers, vector<std::pair <int,int> > & edgesBestSwitchings );
    double runMinRecon();
    void printRecon( MyGeneNode* geneNode=NULL, int cListIdx=0);

};
#endif
