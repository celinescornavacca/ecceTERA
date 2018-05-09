#ifndef MYNETWORK_H_
#define MYNETWORK_H_

/*

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

#include "MySpeciesTree.h"

using namespace bpp;
using namespace std;

/**
 * Information associated with each tree node.
 */
/*
class NetworkNodeInfos  {
	typedef NodeTemplate<NetworkNodeInfos> MyNetworkNode;
	    
	protected:
    bool seenIt;
    bool secondaryEdge;
    MyNetworkNode * secondaryFather; //TO DO : Use this instead of otherParent...
	MyNetworkNode * primaryFather;

    
    public:
    
    void setSeenIt(bool value){
		seenIt = value;
	}
	
	bool getSeenIt(){
		return seenIt ;
	}
	
	void setSecondaryEdge(bool value){
		secondaryEdge = value;
	}
	
	bool getSecondaryEdge(){
		return secondaryEdge ;
	}
	
	void setPrimaryFather(MyNetworkNode * value){
		primaryFather = value;
	}
	
	MyNetworkNode * getPrimaryFather(){
		return primaryFather ;
	}
	
	void setSecondaryFather(MyNetworkNode * value){
		secondaryFather = value;
	}
	
	MyNetworkNode * getSecondaryFather(){
		return secondaryFather ;
	}
	
	
		

		
};

typedef NodeTemplate<NetworkNodeInfos> MyNetworkNode;
*/

// Implementation of BPP tree
class MyNetwork: public MySpeciesTree {
private:
    bool mEnewick;  ///< true if tree is a network

    MySpeciesNode *getPostorderTreeAux( MySpeciesNode *curNode );


    map<string,MySpeciesNode*> tags; // eNewick tags
    virtual MySpeciesNode *handleNewickTag( bool isLeaf, MySpeciesNode *node, 
                                 string tag, string &errString );
                                 
    virtual void handleSecondaryArcs( MySpeciesNode *node);
                                 

                              
    virtual void copyPrimSecParentsInfo(MySpeciesNode *node, MySpeciesNode *nodeToCopy );                  
                                                   
    void deleteNetwork();

public:

    /**
     * Constructor that takes root node.
     */
 	MyNetwork(
            MySpeciesNode &root ) ///< root node
        //: MyTreeTemplate<MyNetworkNode>(& root), mEnewick(true)
        : MySpeciesTree(root), mEnewick(true)
    {  }

    MyNetwork( string description, string &errString );

    /** 
     * Deep copy constructor
     */
    MyNetwork( const MyNetwork &tree ) : MySpeciesTree( tree ) 
    {}

    MyNetwork( MyNetwork &origNetwork );


    /**
     * Deconstrutor that correctly traverses network.
     */
    ~MyNetwork() {
        deleteNetwork();
    }

	void compute_RealPostOrder();

    // tree input io
    static MyNetwork* readMyNetwork( const char *treePathChar,
                                      string &errString );
	void makeBinarySwitching( MySpeciesNode *node=NULL,  bool onlyDeleteUnlabeledLeaves=false );
	int assignNetworkPostOrderIdsDFS();

    void assignNetworkPostOrderIds();
    vector<MySpeciesNode*> getNodes();
    bool checkNetwork( string &errStr );

    void printNetworkAux( MySpeciesNode*node, vector<bool> &seenIt, 
                        string indent="" );
    void printNetwork();
    vector<MySpeciesNode*> getSortedLeaves();
    vector<MyNetwork*> allSwitchings( MySpeciesNode *netRoot, 
            const int biCompId, const vector<int> &biCompMap );
    vector<int> getBiconnectedCompMap( vector<MySpeciesNode*> &biRoots );
    void assignBiconnectedComp( MySpeciesNode*node, const int biNum,
                                vector<int> &map );
    int getBiconnectedCompMapAux( MySpeciesNode*node, int &counter, 
           vector<int> &map, vector<int> &pre, 
           vector<MySpeciesNode*> &otherParents,
           MySpeciesNode *parent=NULL );
           
	void setCorrespondance(vector <MySpeciesNode *> allNodes);

};	

#endif /*MYNETWORK_H_*/
