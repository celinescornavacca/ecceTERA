#ifndef MYSPECIESTREE_H_
#define MYSPECIESTREE_H_

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

#include "MyTreeTemplate.h"

using namespace std;

/**
 * Information associated with each tree node.
 */
class NodeInfos  {
typedef bpp::NodeTemplate<NodeInfos> MySpeciesNode;
public:
    // ids
    int breadthFirstOrder;  ///< breadth first ordering number (for strale)
    unsigned long realPostOrder;  
        ///< post-order number assigned to real nodes (Hali's code)

    // species node information 
    //bool WGD;       ///< whole genome duplication occured
    bool isAlpha;   ///< node is alpha (for transfer to dead)
    int timeSlice;  ///< time slice for species trees
    vector<int> fakeIds;

    // species node specific costs
    double duplicationCost; ///< node specific duplication costs
    double hgtCost;     ///< node specific transfer cost
    double lossCost;    ///< node specific loss cost

    // for network
    bool seenIt;
	MySpeciesNode *primaryFather;
    MySpeciesNode *secondaryFather; //TO DO: Use this instead of otherParent...

    NodeInfos() : //WGD(false), 
        isAlpha(false), primaryFather(NULL), secondaryFather(NULL) {}
};
typedef bpp::NodeTemplate<NodeInfos> MySpeciesNode;


/** 
 * Species tree.
 */
class MySpeciesTree: public MyTreeTemplate<MySpeciesNode> {
protected:
	vector< vector<int> > mCorrespondenceTS; 
        ///< vectors of nodes for each time slice
private:
    bool mHasAlpha;	///< true if addAlphaForDeadTransfer called
    bool mSubdivision; ///< no subdivision
    MySpeciesNode *mAlpha; ///< original alpha node
    vector<bool> mILSnodes;    ///< tracks ils fake nodes
    vector<int> mTimeSlices; ///< time slices mapping by postOrder id

    vector< vector< pair<int,int> > > mSplits; 
        ///< all splits from species tree, by clade
    vector< vector< pair<int,int> > > mIlsSplits; 
        ///< all ILS splits (not in species tree), by clade
        
    bool checkBootstrapValues( MySpeciesNode *node = NULL, 
                               double parentBS = -1 );
    void nameInternalNodes();
    bool assignTimeSlices( vector< pair<int, int> > dateMap,
                           vector<int> &changedTimeSlices, string &errStr);
    MySpeciesNode *getPostorderTreeAux( MySpeciesNode *curNode );
    bool trimTreeAux( boost::unordered_map<string,int> &taxaNames,
            MySpeciesNode *node, int &removedCount,
            MySpeciesNode *&newNode, bool &partial );
    bool computeDistances( bool ultrametric, MySpeciesNode * node, 
                           string &errStr);
    double getLongestDistance( MySpeciesNode *node );


    // ILS Functions
    MySpeciesNode* LCA( vector<int> &set );
    bool isComplete( MySpeciesNode *lca, vector<int> &set );
    int processClades( vector<int> &cladeList );
    vector<int> computeSpeciesCladesAndSplitsAux( MySpeciesNode *node,
                    double ilsCutoff, int maxClusterSize, vector<int> &stats );
public:

    /** empty tree for Strale */
    MySpeciesTree() 
       : mHasAlpha(false), mSubdivision(false)
    {} 


    /**
     * Constructor that takes root node.
     */
 	MySpeciesTree(
            MySpeciesNode &root ) ///< root node
        : MyTreeTemplate<MySpeciesNode>(& root),
        mHasAlpha(false), mSubdivision(false)
    {}  


    MySpeciesTree( string description, string &errString, 
                   bool bootstrap=false );

    /** 
     * Deep copy constructor
     */
    MySpeciesTree( const MySpeciesTree &tree ) 
        : MyTreeTemplate<MySpeciesNode>( tree ) { }


    void assignPostOrderIds();
    static MySpeciesTree* readMySpeciesTree( const char *treePathChar,
                string &errString, bool bootstrap=false );
    MySpeciesTree *getPostorderTree( bool keepOutgroup=false ); 

    // prune leaves not in given list
    bool restrictTreeToLCA( boost::unordered_map<string,int> &taxaNames,
                bool keepRoot, bool verbose=false );
    bool trimTree( boost::unordered_map<string,int> &taxaNames,
                    bool verbose );


    int getArtificialFatherId( int id );
    int getNumberOfIds();
    vector< pair<int,int> > getSplits( int id );
    void printIds();

    /** get realPostOrder field */
    int getRPO( 
            int id ) ///< id
    {
        MySpeciesNode *x = getNodeById( id );
        return x->getInfos().realPostOrder;
    }

    /** return true if id is not assigned to a fake (subdivided) node */
    bool isReal( 
            int id ) ///< id
    {
        MySpeciesNode *node = getNodeById( id ); 
        if( id != node->getId() ) 
            return false;
        return true;
    }

    // time slices
    int getTimeSlice( int id );
	int setVectorTimeSlices();
 	vector<int> getVectorWithTS(int ts);
    int getAlphaIdForTS( int ts );
    void assignNoSubdivisionTimeSlices();

    // ILS functions 
    void computeSpeciesCladesAndSplits( double ilsCutoff, int maxClusterSize,
                                        bool verbose );
    bool hasILS() {
        if( mIlsSplits.size() == 0 ) 
            return false;
        return true;
    }
    vector< pair<int,int> > getIlsSplits( int id );
    bool isILS( size_t id ) {
        if( id > mILSnodes.size() )
            throw bpp::Exception( "MySpeciesTree::isILS: invalid id" );
        return mILSnodes[id];
    }

    // give all leaves the same depth
    bool computeSubdivision( vector< pair<int,int> > dateMap, 
            bool bootstrapOrdering, bool ultrametric,
            vector<int> &changedTimeStamps, string &errStr );

    // alpha for transfer from the dead
    /** Has alpha */
    bool hasAlpha() { return mHasAlpha; }
    void addAlphaForDeadTransfer( bool bootstrapOrdering,
                    double hgtCost, double lossCost );
    /** True if node with id is alpha */
    bool isAlpha( int id ) ///< node id
    {
        MySpeciesNode *node = getNodeById( id );
        if( node == NULL )
            return false;
        return node->getInfos().isAlpha;
    }

    // print tree functions
    void printSif( MySpeciesNode *node=NULL );
    void printTreeInfo( MySpeciesNode *node=NULL, int level = 0 );

    // strale associated functions
    void breadthFirstreNumber();
    void lexicographicBreadthFirstreNumber();
    vector<int> findChangedTimeSlices( MySpeciesTree *otherTree, 
                                       bool printChanges=false );
    string assignCosts( vector<double> &duplicationCosts, 
                        vector<double> &hgtCosts, vector<double> &lossCosts ); 
    void assignFixedCosts( double duplicationCost,
                           double hgtCost, double lossCost );
    bool assignCosts( string costsFileName );


    // Hali's functions
    void compute_RealPostOrder();

    //WDF
    bool isSubdivided();
    //WDFend
    
    string toRecPhyloXML();

    string toRecPhyloXML(MySpeciesNode *node);

};	

#endif /*MYSPECIESTREE_H_*/
