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


@section DESCRIPTION
MySpeciesTree implements species trees using the BPP tree templates.

The species tree specific functions are computing subdivions (assigning
internal nodes to time slices and adding artifical nodes) and adding
an alpha node for transfering from dead (extinct) species nodes.

*/


#include "MySpeciesTree.h"

#include <iostream>
#include <queue>
#include <boost/foreach.hpp>


/** Get bootstrap value from a node.
 * @return bootstrap value */
inline double getBootstrap( const MySpeciesNode *node ) {
    return dynamic_cast<const bpp::Number<double> *> 
               (node->getBranchProperty(bpp::TreeTools::BOOTSTRAP))->getValue();
}

/** 
 * Constructor from a string in Newick format.
 */
MySpeciesTree::MySpeciesTree( 
        string description, ///< Newick string 
        string &errString, ///< error description
        bool bootstrap ) ///< names are bootstrap values
    : MyTreeTemplate<MySpeciesNode>(), mHasAlpha(false), mSubdivision(false)
{
    readTree( description, errString, false, bootstrap );

    if( bootstrap && !checkBootstrapValues() ) 
        errString = "Bootstrap values incorrectly ordered.";
}






/**
 * Read newick tree from a file.
 *
 * @return a tree
 */
MySpeciesTree* MySpeciesTree::readMySpeciesTree( 
    const char *treePathChar, ///< path to file with newick trees
    string &errString, ///< error description
    bool bootstrap ) ///< names are bootstrap values
{
    vector<string> treeStrings = readNewickStrings( treePathChar );
    if( treeStrings.size() == 0 )
        return NULL;

    return new MySpeciesTree( treeStrings[0], errString, bootstrap );
}


/**
 * Assign ids using breadth-first postorder.
 *
 * This assigns ids outside of the BPP code, giving explicit
 * control of the type of ordering used.
 *
 * Also, set the mCorrespondence vector mapping ids to nodes.
 */
void MySpeciesTree::assignPostOrderIds()
{
    // initialize to -1 to indicate not assigned
    vector<MySpeciesNode*> allNodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, allNodes ) 
        node->setId( -1 );

    int pOrd = 0;

    // bottom-up breadth first method
    vector<MySpeciesNode*> nodes = getSortedLeaves(); // nodes for one level

    // Assign ids bottom up.
    while( nodes.size() > 0 ) {
        vector<MySpeciesNode*> parentNodes; // nodes for next level up
        BOOST_FOREACH( MySpeciesNode *node, nodes ) {
            if( node->getId() == -1 ) {
                mCorrespondence.push_back( node );
                node->setId( pOrd++ );
                if( node->hasFather() ) 
                    parentNodes.push_back( node->getFather() );
            }
        }
        nodes = parentNodes;
    }

    if( mSubdivision ) {
        // create fake ids and assign values to mCorrespondence and mTimeSlices

        int maxTS = getRootNode()->getInfos().timeSlice;
        vector<MySpeciesNode*> nodes = getSortedLeaves(); 
        int leafCount = nodes.size();
        mTimeSlices.resize( mCorrespondence.size() );
        // assign fake ids bottom up
        for( int ts=1; ts<=maxTS; ts++ ) {
            bool foundIt=false;
            // loop through all of the nodes in this time slice
            for( size_t idx=0; idx<nodes.size(); idx++ ) {
                MySpeciesNode *node = nodes[idx];
                if( node == NULL || !node->hasFather() )
                    continue;
                MySpeciesNode *father = node->getFather();
                int sonTimeSlice = node->getInfos().timeSlice;
                int fatherTimeSlice = father->getInfos().timeSlice;
                if( ts > sonTimeSlice && ts < fatherTimeSlice ) {
                    mCorrespondence.push_back( node );
                    mTimeSlices.push_back( ts );
                    int id = pOrd++;
                    node->getInfos().fakeIds.push_back( id );
                }

                // once the node for this time slice is reached, 
                // set one son to it and the other to NULL to 
                // avoid doing it twice
                if( fatherTimeSlice == ts ) {
                    if( foundIt ) 
                        nodes[idx] = NULL;
                    else {
                        // first son
                        foundIt = true;
                        nodes[idx] = father;
                    }
                }
            }
        }

        // redo real node ids by time slice
        BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
            int id;
            if( !node->isLeaf() ) {
                id = leafCount + node->getInfos().timeSlice-1;
                mCorrespondence[id] = node;
                node->setId( id );
            } else {
                id = node->getId();
            }
            mTimeSlices[id] = node->getInfos().timeSlice;
        }
    }
}


/**
 * Check if bootstrap values are ordered correctly, i.e.,
 * parents greater than children.
 *
 * Internal nodes must be greater than 0, which are checked
 * when reading the tree.
 *
 * @return True if ordered correctly.
 */
bool MySpeciesTree::checkBootstrapValues(
        MySpeciesNode *node, ///< node in traversal (default = NULL)
        double parentBS ) ///< parent bootstrap value (default == -1)
{
    if( node == NULL ) 
        node = getRootNode(); 
        // root defaults to bootstrap -1, meaning not there

    if( node->isLeaf() ) 
        return true;

    double bs = parentBS; // propagate parent value if missing
    if( node->hasBootstrapValue() ) {
        bs = getBootstrap( node );
        if( bs <= 0 )
            throw bpp::Exception( "MySpeciesTree::checkBootstrapValues:"
                " negative or zero bootstrap, should be checked before this." );
        if( parentBS != -1 && bs >= parentBS )
            return false;
    }

    size_t childCount = node->getNumberOfSons();
    for( size_t i=0; i<childCount; i++ ) 
        if( !checkBootstrapValues( node->getSon(i), bs ) )
            return false;

    return true;

    
}


/** 
 *  Remove nodes above the least common ancestor of all of the
 *  nodes whose taxa is in the given map.
 *
 *  @return False if tree has none of the taxa
 */
bool MySpeciesTree::restrictTreeToLCA (
         boost::unordered_map<string,int> &taxaNames,
            ///< map of leaf names
         bool keepRoot, ///< Make LCA a child of root
         bool verbose ) ///< print debugging info
{
    // count the number of descendants of each node
    int originalCount = getNumberOfNodes();
    vector<int> descendantCounts( originalCount, 0 );
    boost::unordered_map<string,int>::iterator iter; 
    vector<MySpeciesNode *> leaves = getLeaves();
    int foundLeafCount = 0;
    BOOST_FOREACH( MySpeciesNode *leaf, leaves ) {
        string name = leaf->getName();
        iter = taxaNames.find( name );
        if ( iter != taxaNames.end() ) {
            foundLeafCount++;
            MySpeciesNode *node = leaf;
            descendantCounts[node->getId()]++;
            while( node->hasFather() ) {
                node = node->getFather();
                descendantCounts[node->getId()]++;
            }
		}
	}
    if( foundLeafCount == 0 )
        return false; // nothing in common

    // find the LCA
    MySpeciesNode *LCA = getRootNode();
    MySpeciesNode *son = LCA;
    bool hasNewRoot = false;
    while( descendantCounts[son->getId()] == foundLeafCount ) {
        LCA = son;
        int sonCount = LCA->getNumberOfSons();
        if( sonCount == 0 ) {
            // leaf
            if( foundLeafCount != 1 )
                throw bpp::Exception( "MySpeciesTree::restrictTreeToLCA: "
                        " foundLeafCount not 1 but LCA is a leaf" );
            break;
        }
        for( int i=0; i<sonCount; i++ ) {
            son = LCA->getSon(i);
            if( descendantCounts[son->getId()] == foundLeafCount ) {
                hasNewRoot = true;
                break;
            }
        }
    }

    // trim tree
    if( hasNewRoot ) {
        int deleteCount = 0;

        // remove unneeded subtrees
        MySpeciesNode *ancestor = LCA->getFather();
        int keepSonId = LCA->getId();
        vector<MySpeciesNode*> ancestors;
        MySpeciesNode *LCApathSon = NULL; // for keepRoot
        while( true ) {
            int sonCount = ancestor->getNumberOfSons();
      
            if( ancestor->hasFather() || !keepRoot ) {
                for( int i=0; i<sonCount; i++ ) {
                    MySpeciesNode *sibling = ancestor->getSon(i);
                    if( sibling->getId() != keepSonId )
                        deleteCount += deleteSubtree( sibling );
                    else
                        LCApathSon = ancestor;
                }

                ancestors.push_back( ancestor ); // to delete list
            }
            if( !ancestor->hasFather() ) 
                break;

            // set up next iteration
            keepSonId = ancestor->getId(); 
            ancestor = ancestor->getFather();
        }

        if( ancestors.size() > 0 ) {
            LCA->removeFather(); // do this before deleting father

            if( keepRoot ) {
                // set LCA as root son, replacing its ancestor
                MySpeciesNode *root = getRootNode();

                vector<MySpeciesNode*> sons;
                for( size_t i=0; i<root->getNumberOfSons(); i++ )
                    sons.push_back( root->getSon(i) );
                root->removeSons();
                BOOST_FOREACH( MySpeciesNode *son, sons ) {
                    if( son == LCApathSon )
                        root->addSon( LCA );
                    else
                        root->addSon( son );
                }
            } else {
                // set LCA as root
                setRootNode( LCA );
            }

            // delete ancestors of LCA
            BOOST_FOREACH( MySpeciesNode *ancestor, ancestors ) {
                delete ancestor;
                deleteCount++;
            }

            resetNodesId();
        }

        if( !keepRoot && verbose ) {
            //cout << "Species tree trim found " << foundLeafCount
            //   << " / " << getNumberOfLeaves() << " leaves in common" << endl;
            cout << "Species tree trim removed " << deleteCount 
                << " / " << originalCount << " nodes" << endl;
        }
     }

    return true;
}


/**
 * Traverse tree, removing all but first unused (not in taxaNmes) sons.
 *
 * @return true if all sons are unused
 */
bool MySpeciesTree::trimTreeAux(
        boost::unordered_map<string,int> &taxaNames,
            ///< map of leaf names
        MySpeciesNode *node, ///< current nodes
        int &removedCount,  ///< number of nodes removed
        MySpeciesNode *&newNode, ///< returns this node if it is trimmed
        bool &partial ) ///< true if one used and one unused son
{
    newNode = node;
    partial = false;

    if( node->isLeaf() ) {
        string name = node->getName();
        boost::unordered_map<string,int>::iterator iter 
                = taxaNames.find( name );
        if ( iter == taxaNames.end() ) 
            return true; // not found in taxaNames
        else
            return false;
    }

    // two bad, delete this node, keep first son
    
    int sonCount = node->getNumberOfSons();
    if( sonCount > 2 ) // doesn't handle polytomy
        throw bpp::Exception( "MySpeciesTree::trimTreeAux: not binary" );

    vector<MySpeciesNode*> sons;
    for( int i=0; i<sonCount; i++ ) 
        sons.push_back( node->getSon(i) );

    bool allUnused = true;
    bool hasUnused = false;
    // returns son or whatever descendant replaces it
    MySpeciesNode *newSon = NULL; 
    MySpeciesNode *firstUnusedSon = NULL;
    bool sonPartial = false;
    BOOST_FOREACH( MySpeciesNode *son, sons ) {
        if( !trimTreeAux( taxaNames, son, removedCount, newSon, sonPartial ) ) 
        {
            // in taxaNames
            allUnused = false;
        } else if( hasUnused ) {
            // both sons unused, delete second (current) son
            node->removeSon( son );
            removedCount += deleteSubtree( son );

            // push node distance down
            if( node->hasDistanceToFather() ) {
                double dist = node->getDistanceToFather();
                node->deleteDistanceToFather();
                MySpeciesNode *nextNode = firstUnusedSon;
                while( nextNode->getNumberOfSons() == 1 ) {
                    nextNode = nextNode->getSon( 0 );
                }
                if( !nextNode->hasDistanceToFather() ) {
                    throw bpp::Exception( "MySpeciesTree::trimTreeAux:"
                            " next node has no distance to father" );
                }
                dist += nextNode->getDistanceToFather();
                nextNode->setDistanceToFather( dist );
            } else if( node->hasFather() ) {
                throw bpp::Exception( "MySpeciesTree::trimTreeAux: non-root has"
                        " no distance to father" );
            }
        } else {
            hasUnused = true;
            firstUnusedSon = newSon;
        }
    }

    return allUnused;
}



/** 
 *  Merge sibling nodes whose taxa are in the given map.
 *
 *  @return False if tree has none of the taxa
 */
bool MySpeciesTree::trimTree(
         boost::unordered_map<string,int> &taxaNames,
            ///< map of leaf names
         bool verbose ) ///< print debugging info
{
    int originalCount = getNumberOfNodes();
    int deletedCount = 0;
    MySpeciesNode *newRoot = NULL;
    bool partial = false;
    if( trimTreeAux( taxaNames, getRootNode(), deletedCount, newRoot, partial) )
        return false; // nothing in common

    resetNodesId();

    restrictTreeToLCA ( taxaNames, true, verbose );

    if( verbose ) {
        int removedCount = originalCount - getNumberOfNodes();
        cout << "Species tree trim removed " << removedCount
            << " / " << originalCount << " nodes" << endl;
    }
   
    if( verbose ) {
        int removedCount = originalCount - getNumberOfNodes();
        cout << "Species tree trimSL removed " << removedCount
            << " / " << originalCount << " nodes" << endl;
    }


    return true;
}


/**
 * Return the id of the father node in the subdivided tree.
 *
 * @return father id
 */
int MySpeciesTree::getArtificialFatherId( int id ) {
    MySpeciesNode *xNode = getNodeById( id );
    if( xNode->hasFather() )
        if( mSubdivision ) {
            int fakeSize = xNode->getInfos().fakeIds.size();
            if( fakeSize == 0 
                || id == xNode->getInfos().fakeIds[fakeSize-1] )
            {
                return xNode->getFather()->getId();
            }

            int nodeId = xNode->getId();
            int fatherIdx=0;
            int idx = 0;
            while( nodeId != id && idx < fakeSize ) 
                nodeId = xNode->getInfos().fakeIds[idx++];
            if( nodeId != id )
                throw bpp::Exception("MySpeciesTree::getArtficialFatherId:"
                                " Could not find node" );
            return xNode->getInfos().fakeIds[fatherIdx];
            
        } else
            return xNode->getFather()->getId();
    else
        return -1;
}


/**
 * Return time slice for this id.
 */
int MySpeciesTree::getTimeSlice(
        int id ) ///< tree id
{
    //if( mTimeSlices.size() > 0 ) 
    if( mSubdivision )
        return mTimeSlices[id];

    MySpeciesNode *node = getNodeById( id );
    return node->getInfos().timeSlice;
    
}

vector< pair<int,int> > MySpeciesTree::getSplits( 
    int id ) ///< a post order id
{
    // if mSplits is set, then species tree is possibly polytomic
    // and id could be a fake id
    if( mSplits.size() != 0 ) { 
        if( id >= (int) mSplits.size() ) 
            throw bpp::Exception( "MySpeciesTree:getSplits:"
                       " request for invalid id" );
        return mSplits[id]; 
    }

    MySpeciesNode *node = getNodeById( id );
    int sonCnt = node->getNumberOfSons();
    vector< pair<int,int> > idPairs;
    pair<int,int> idPair;
    if( mSubdivision ) {
        vector<int> sonIds;
        if( node->getId() == id ) {
            // return two sons (possibly artficial) for real nodes
            for( int i=0; i<sonCnt; i++ ) {
                MySpeciesNode *son = node->getSon(i);
                int fakeSize = son->getInfos().fakeIds.size();
                if( fakeSize == 0 )
                    sonIds.push_back( son->getId() );
                else 
                    sonIds.push_back( son->getInfos().fakeIds[fakeSize-1] );
            }
        } else {
            // return one son for artficial nodes
            int artIdx = -1;
            int fakeSize = node->getInfos().fakeIds.size();
            for( artIdx=0; artIdx<fakeSize; artIdx++ )
                if( id == node->getInfos().fakeIds[artIdx] )
                    break;
            // add child id
            if( artIdx == 0 )
                sonIds.push_back( node->getId() );
            else if( artIdx < 0 || artIdx >= fakeSize )
               throw bpp::Exception( "MySpeciesTree:getSplit:"
                       " could not find node" );
            else
                sonIds.push_back( node->getInfos().fakeIds[artIdx-1] );
        }

        if( sonIds.size() > 0 ) {
            pair<int,int> idPair;
            idPair.first = sonIds[0];
            if( sonIds.size() == 2 )
                idPair.second = sonIds[1];
            else
                idPair.second = -1;
            idPairs.push_back( idPair );
        }

    } else if( sonCnt > 0 ) {
        idPair.first = node->getSon(0)->getId();
        if( sonCnt == 2 )
            idPair.second = node->getSon(1)->getId();
        else
            idPair.second = -1;
        idPairs.push_back( idPair );
    }

    return idPairs;
}


/**
 * Return the post order ids of the sons of the node with
 * the given id. -1 is returned for ids of sons that don't
 * exist. This correponds to a single tripartition of id.
 */
/*
vector<int> MySpeciesTree::getSplit( 
    int id ) ///< a post order id
{
    vector<int> sonIds;

    MySpeciesNode *node = getNodeById( id );
    int sonCnt = node->getNumberOfSons();

    if( mSubdivision ) {
        if( node->getId() == id ) {
            // return two sons (possibly artficial) for real nodes
            for( int i=0; i<sonCnt; i++ ) {
                MySpeciesNode *son = node->getSon(i);
                int fakeSize = son->getInfos().fakeIds.size();
                if( fakeSize == 0 )
                    sonIds.push_back( son->getId() );
                else 
                    sonIds.push_back( son->getInfos().fakeIds[fakeSize-1] );
            }
        } else {
            // return one son for artficial nodes
            int artIdx = -1;
            int fakeSize = node->getInfos().fakeIds.size();
            for( artIdx=0; artIdx<fakeSize; artIdx++ )
                if( id == node->getInfos().fakeIds[artIdx] )
                    break;
            // add child id
            if( artIdx == 0 )
                sonIds.push_back( node->getId() );
            else if( artIdx < 0 || artIdx >= fakeSize )
               throw bpp::Exception( "MySpeciesTree:getSplit:"
                       " could not find node" );
            else
                sonIds.push_back( node->getInfos().fakeIds[artIdx-1] );
        }
    } else {
        for( int i=0; i<sonCnt; i++ ) 
            sonIds.push_back( node->getSon(i)->getId() );
    }

    return sonIds;
}
*/




/**
 * Calculate node time slices (tree levels).
 *
 * @return maximum time slice
 */
int MySpeciesTree::setVectorTimeSlices()
{
    // the root has the maximum time slice
    int maxTS = getRootNode()->getInfos().timeSlice;

    // Create lists (mCorresponanceTS) of all nodes
    // in each time slice.
    mCorrespondenceTS.clear();
    mCorrespondenceTS.resize(maxTS+1);

    //if( mTimeSlices.size() > 0 ) {
    if( mSubdivision ) {
        vector<int> alphaIds; // need to be last
        for( size_t id=0; id<mTimeSlices.size(); id++ ) {
            if( mHasAlpha && isAlpha( id ) ) 
                alphaIds.push_back( id );
            else {
                int ts = mTimeSlices[id];
                mCorrespondenceTS[ts].push_back( id );
            }
        }
        // now do alpha ids
        BOOST_FOREACH( int id, alphaIds ) {
            int ts = mTimeSlices[id];
            mCorrespondenceTS[ts].push_back( id );
        }
    } else {
        // mTimeSlices not set if there is no subdivision
        vector <MySpeciesNode *> nodes = getNodes();
        BOOST_FOREACH( MySpeciesNode *node, nodes ) {
            int ts = node->getInfos().timeSlice;
            mCorrespondenceTS[ts].push_back( node->getId() );
        }
    }

    return maxTS;
}


/**
 * Return all of the nodes in the given time slice by id.
 *
 * @return a vector of postorder ids in the given time slice
 */
vector<int> MySpeciesTree::getVectorWithTS( int ts ) ///< A time slice.
{
    if( ts < 0 || (size_t) ts >= mCorrespondenceTS.size() )
        throw bpp::Exception( "MySpeciesTree::getVectorWithTS-"
                " given invalid time slice" );
    return mCorrespondenceTS[ts]; 
}


/**
 * Return alpha for the given time slice.
 *
 * @return alpha (extinct/dead)  node for this time slice or -1
 *     if the time slice is for the root
 */
int MySpeciesTree::getAlphaIdForTS( int ts ) ///< A time slice.
{
    if( !mHasAlpha ) 
        throw bpp::Exception( "MySpeciesTree::getAlphaIdForTS called when tree"
                " has no alpha" );

    if( ts < 0 || (size_t) ts >= mCorrespondenceTS.size() )
        throw bpp::Exception( "MySpeciesTree::getAlphaIdForTS-"
                " given invalid time slice" );

    // root has no alpha
    if( (size_t) ts == mCorrespondenceTS.size()-1 )
        return -1;

    int alphaId;
    if( !mSubdivision ) 
        alphaId = mAlpha->getId(); // use original alpha
    else {
        alphaId = mCorrespondenceTS[ts].back(); 

        // Should be last in time slice.
        // check it, if not throw exception
        if( !isAlpha( alphaId ) )
            throw bpp::Exception( "MySpeciesTree::getAlphaIdForTS:"
                    " alpha is not the last vector element" );
    }

    return alphaId;
}











/**
 * Traverse tree, assigning distances (BOOTSTRAP branch property).
 *
 * Each node's distance is the sum of the distanceToFather and
 * bootstrap properties.
 *
 * An error is returned if a node does not have a distance
 * to father or if the tree is not ultrametric (son lengths different).
 *
 * @return True if there were no errors.
 */
bool MySpeciesTree::computeDistances(
        bool ultrametric, ///< report error if not ultrametric
        MySpeciesNode * node,  ///< current node in tbe recursion
        string &errStr) ///< description of error, if any
{

    double maxDist = 0; 
    for( size_t i=0; i<node->getNumberOfSons(); i++) {

        if( !computeDistances( ultrametric, node->getSon(i), errStr ) ) 
            return false;

        //ULTRAMETRIC
        if( node->getSon(i)->getInfos().isAlpha ) {
            // this is alpha
            continue;
        } else if( !node->getSon(i)->hasDistanceToFather() ) {
            stringstream ss;
            ss << "No distance to father for node " 
                << node->getId();
            errStr = ss.str();
            return false;
        }

        double dist  = node->getSon(i)->getDistanceToFather() 
                        + getBootstrap(node->getSon(i));

        double frac = abs(dist-maxDist)/dist;
        if( i!=0 && frac > 0.00001 ) {
            stringstream ss;
            ss << "Tree is not ultrametric at node " 
                << node->getId();
            errStr = ss.str();
            return false;
        }

        // The nodes should have the same length, however
        // due to rounding errors and branch lengths of zero, 
        // this might not be true. Take the maximum to ensure
        // correct ordering of nodes
        if( i==0 || dist > maxDist ) 
            maxDist = dist;
    }

    node->setBranchProperty( bpp::TreeTools::BOOTSTRAP, 
                             bpp::Number<double>(maxDist));

    return true;
}



/** 
* Sort vector from least distance to greatest.
* For ties, use greater breadth first ordering
*/
struct distance_sort {
/** distance sorting operator */
inline bool operator() (const MySpeciesNode * a, const MySpeciesNode * b ) 
{
    double bsA = getBootstrap(a);
    double bsB = getBootstrap(b);
    if( bsA == bsB ) 
        // otherwise post order so that children are before parents
        return b->getId() > a->getId(); 

    return bsB > bsA;
}
};

/**
 * Assign time slices. Leaves are 0. Internal nodes are ordered by
 * bootstrap distance.
 *
 * Time slice pairs in the dateMap vector are swapped. These values are
 * given in breadth first order. The corresponding changed time slices 
 * are returned in the changedTimeSlices variable.
 *
 * @return True if reordering doesn't change topology.
 */
bool MySpeciesTree::assignTimeSlices( 
        vector<pair<int, int> > dateMap, ///< pairs of time slices to swap
        vector<int> &changedTimeSlices, ///< assigned time slices swapped
        string &errStr) ///< description of error, if any
{

    vector <MySpeciesNode *> nodes = getNodes();
    vector <MySpeciesNode *> internalNodes;
    BOOST_FOREACH( MySpeciesNode *node, nodes) {
        if( node->isLeaf() )
            node->getInfos().timeSlice = 0;
        else
            internalNodes.push_back( node );
    }

    sort( internalNodes.begin(), internalNodes.end(), distance_sort() );
    int ts = 1;
    BOOST_FOREACH( MySpeciesNode *node, internalNodes ) {
        node->getInfos().timeSlice = ts;
        ts++;
    }


    // switch time slice from date map
    if( dateMap.size() != 0 ) {

        vector<MySpeciesNode*> nodes = getNodes();		
        vector<MySpeciesNode*> indexedNodes( nodes.size() );
        BOOST_FOREACH( MySpeciesNode* node, nodes ) 
            indexedNodes[node->getInfos().breadthFirstOrder] = node;

        pair<int,int> p;
        BOOST_FOREACH( p, dateMap ) {
            int ts1 = indexedNodes[p.first]->getInfos().timeSlice;
            int ts2 = indexedNodes[p.second]->getInfos().timeSlice;

//cout << "**** SWAPPING " << p.first << " (ts=" << ts1 << ",id=" 
//    << indexedNodes[p.first]->getId() << ") and " << p.second << " (ts=" 
//    <<  ts2 << ",id=" << indexedNodes[p.second]->getId() << ")" << endl;

            indexedNodes[p.first]->getInfos().timeSlice = ts2;
            indexedNodes[p.second]->getInfos().timeSlice = ts1;
                
            changedTimeSlices.push_back( ts1 );
            changedTimeSlices.push_back( ts2 );
        }

        BOOST_FOREACH( MySpeciesNode *node, internalNodes ) {
            if( node->hasFather() ) {
                MySpeciesNode *father = node->getFather();
                int sonTimeSlice = node->getInfos().timeSlice;
                int fatherTimeSlice = father->getInfos().timeSlice;
                if( sonTimeSlice >= fatherTimeSlice ) {
                    errStr = "Reordering changed topology";
                    return false;
                }
            }
        }
    }

    return true;
}



/** 
 * Compute the subdivision, assigning time slices and adding nodes 
 * to give all branches the same height.
 *
 * All leaves are assigned to time slice zero. Internal nodes are
 * sorted by distance and consecutively assigned time slices, with
 * the root at the maximal time slice. If bootstrapOrdering is true,
 * the bootstrap distances are used directly. If not, the branch length
 * from getDistanceToFather function is used to calculate the distance
 * to the leaves.
 *
 * If a dateMap is given, the internal nodes indicated by 
 * breadthFirstOrdering have their time slices swapped.
 *
 * @return false if the date map has errors
 */
bool MySpeciesTree::computeSubdivision( 
        vector<pair<int, int> > dateMap, ///< pairs of time slices to swap
        bool bootstrapOrdering, ///< Ordering given from bootstrap values.
        bool ultrametric, ///< report error if not ultrametric
        vector<int> &changedTimeSlices, ///< assigned time slices swapped
        string &errStr) ///< description of error, if any
{
    mSubdivision = true;

    if( bootstrapOrdering ) { 
        // if bootstrap, skip this step, but set all leaves to 0
        vector <MySpeciesNode *> leaves = getLeaves();		
        BOOST_FOREACH (MySpeciesNode *node, leaves) 
            node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, 
                                    bpp::Number<double>(0));
        // check for bootstrap values
        vector <MySpeciesNode *> nodes = getNodes();		
        BOOST_FOREACH (MySpeciesNode *node, nodes) {
            if( !node->hasBootstrapValue() ) {
                stringstream ss;
                ss << "No bootstrap value for node " 
                    << node->getId();
                errStr = ss.str();
                return false;
            }
        }
    } else {
        // set bootstrap branch properites with distances
	    if( !computeDistances( ultrametric, getRootNode(), errStr ) ) 
            return false;
    }


    if( !assignTimeSlices( dateMap, changedTimeSlices, errStr ) )
        return false;

//    printTreeInfo(); // debugging

    // sort changed time stamps and remove duplicates
    if( changedTimeSlices.size() > 0 ) {
        sort( changedTimeSlices.begin(), changedTimeSlices.end() );
        int uniqueIdx=0;
        for( size_t i=0; i<changedTimeSlices.size(); i++ ) {
            if( i==0 
                || changedTimeSlices[i] != changedTimeSlices[uniqueIdx-1] ) 
            {
                changedTimeSlices[uniqueIdx++] = changedTimeSlices[i];
            }
        }
        changedTimeSlices.resize( uniqueIdx );
    }

    return true;
} 


/**
 * Return longest distance from node to root.
 */
double MySpeciesTree::getLongestDistance(
        MySpeciesNode *node ) ///< a tree node
{
    double dist = 0;
    int nbrSons = node->getNumberOfSons();
    for( int i=0; i<nbrSons; i++ )  {
        double sonDist = getLongestDistance( node->getSon(i) );
        if( sonDist > dist )
            dist = sonDist;
    }
    if( node->hasDistanceToFather() )
        dist += node->getDistanceToFather();

    return dist;
}

/**
 * Add an outGroup to transfer from the dead by creating
 * a sibling of the root and then creating a new root.
 *
 * The assigned costs are used with variable costs and are
 * not necessary for fixed costs. There is no duplication cost
 * for the dead.
 */
void MySpeciesTree::addAlphaForDeadTransfer(
        bool bootstrapOrdering, ///< set bootstrap of new root if true
        double hgtCost, ///< transfer cost for alpha
        double lossCost ) ///< loss cost for alpha
{
    mHasAlpha = true;

    // create new root with old tree and alphas as sons
    MySpeciesNode *oldRootNode = getRootNode();
// What should this be?
    oldRootNode->setDistanceToFather(1);	

    MySpeciesNode *newRootS = new MySpeciesNode(); 
//    newRootS->getInfos().WGD = false;
    newRootS->getInfos().isAlpha = false;

    mAlpha = new MySpeciesNode();
    mAlpha->setName("OUTGROUP");
//    mAlpha->getInfos().WGD = false;
    mAlpha->getInfos().isAlpha = true;
    mAlpha->getInfos().duplicationCost = 0;
    mAlpha->getInfos().hgtCost = hgtCost;
    mAlpha->getInfos().lossCost = lossCost;

    double dist = getLongestDistance( oldRootNode );
    mAlpha->setDistanceToFather( dist );

    newRootS->addSon(oldRootNode);
    newRootS->addSon(mAlpha);

    setRootNode(newRootS); 
    resetNodesId();
    
    if( bootstrapOrdering ) {
        if( !oldRootNode->hasBootstrapValue() ) 
            throw bpp::Exception("bootstrap.ordering true, but root"
                    " has no bootstrap value" );

        // increment old bootstrap value by one for new root
        int oldBS = getBootstrap( oldRootNode );
        newRootS->setBranchProperty(bpp::TreeTools::BOOTSTRAP, 
                                        bpp::Number<double>(oldBS+1));
    }
}





/**
 * Assign internal nodes a name.
 *
 * Name is left most leaf of the subtree + @[level]
 * where [level] is the height of the subtree.
 */
void MySpeciesTree::nameInternalNodes() 
{
    vector <MySpeciesNode *> nodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, nodes ) 
        if( !node->isLeaf() ) 
            node->setName("");

    vector <MySpeciesNode *> leaves = getSortedLeaves();

    // assign un-named internal nodes bottom up
    BOOST_FOREACH( MySpeciesNode *leaf, leaves ) {
        int cnt = 1;
        MySpeciesNode *node = leaf;
        while( node->hasFather() ) {
            node = node->getFather();
            if( !node->hasName() || node->getName() == "" ) {
                stringstream name;
                name << leaf->getName() << "@" <<  cnt;
                node->setName( name.str() );
            } else {
                break;
            }
            cnt++;
        }
    }
}



/**
 * Return the number of assign ids.
 *
 * This includes artificial nodes (non-instantiated) if subdivided
 *
 * @return the number of assigned ids
 */
int MySpeciesTree::getNumberOfIds() {
    if( mCorrespondence.size() == 0 )
       throw bpp::Exception("MySpeciesTree::getNumberOfIds: postOrder not set");
    return mCorrespondence.size();
}



/**
 * pre-order print to cout
 */
void MySpeciesTree::printTreeInfo( 
        MySpeciesNode *node, ///< current node in the recursion
        int level ) ///< depth of recursion
{
    if( mCorrespondence.size() == 0 )
        throw bpp::Exception("MySpeciesTree::printTreeInfo: postOrder not set");

    if( node == NULL ) {
        nameInternalNodes();
        node = getRootNode();
    }
//return;
    //indent
    for( int i=0; i<level; i++ ) 
        cout << " ";

	int nSons=node->getNumberOfSons();
	int po = node->getId();
    int ts = node->getInfos().timeSlice;
    int bf = node->getInfos().breadthFirstOrder;
    cout << po << ":" << node->getName() 
         << " ts=" << ts << " bf=" << bf
         << endl;

	for(int i=0; i<nSons; i++)
		printTreeInfo( node->getSon(i), level+1 );
}



/**
 * Strale breadth-first ordering, which is used to identify
 * internal node names.
 */
void MySpeciesTree::breadthFirstreNumber() 
{
    if( 1 ) {
        lexicographicBreadthFirstreNumber();
        return;
    }

	int index = 0;
	std::map<MySpeciesNode*, int> color ;
	vector<MySpeciesNode*> nodes = getNodes();
	//All nodes white
	for (unsigned int i = 0; i< nodes.size() ; i++) {
		color.insert(std::pair<MySpeciesNode*,int>(nodes[i],0));
	}
	std::queue <MySpeciesNode *> toDo;
	toDo.push(getRootNode());
	color[getRootNode()] = 1;
	
	getRootNode()->getInfos().breadthFirstOrder = index;
	
	vector <int> v;
	index++;
	while(!toDo.empty()) {
		MySpeciesNode *u = toDo.front();
		toDo.pop();
		std::vector<MySpeciesNode*> sons;
		for (unsigned int j = 0 ; j< u->getNumberOfSons() ; j++) {
			sons.push_back(u->getSon(j));
		}
		for (unsigned int j = 0; j< sons.size() ; j++) {
			if (color[sons[j]]==0) {
				color[sons[j]]=1;
	            sons[j]->getInfos().breadthFirstOrder = index;
				index++;
				toDo.push(sons[j]);
			}
		}
		color[u]=2;
	}

}

/**************************************************************************
 * This function re-numbers nodes with a breadth-first traversal.
 * Contrary to the above function, deeper nodes have HIGHER indices than
 * later nodes.
 ***************************************************************************/
void MySpeciesTree::lexicographicBreadthFirstreNumber()
{
    vector<MySpeciesNode*> nodes = getNodes();

    std::map<string,MySpeciesNode*> name_node;
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
        if( node->isLeaf() ) {
            name_node[node->getName()] = node;
        } else {
            vector<string> leafnames = 
                    bpp::TreeTemplateTools::getLeavesNames(*node);
            sort(leafnames.begin(),leafnames.end());
            stringstream name;
            BOOST_FOREACH( string st, leafnames ) 
                name << st << ".";

            name_node[name.str()] = node;
        }
    }

    // register species
    int last_branch=0;
    vector<bool> saw( getNumberOfNodes(), false );
    std::map<MySpeciesNode*,int> node_ids;  //Map between node and its id.
    for( map <string,MySpeciesNode*>::iterator it=name_node.begin();
            it!=name_node.end(); it++ )
    {
        if( (*it).second->isLeaf() ) {
            MySpeciesNode *node = (*it).second;
            node_ids[node] = last_branch;
            last_branch++;
            node->getInfos().breadthFirstOrder = last_branch;
            saw[node->getId()] = true;
        }
    }

    //ad-hoc postorder
    vector<MySpeciesNode*> next_generation;
    for( map <string,MySpeciesNode*>::iterator it=name_node.begin();
            it!=name_node.end(); it++ )
    {
        if( (*it).second->isLeaf() ) {
            MySpeciesNode *node = (*it).second;
            next_generation.push_back(node);
        }
    }

    while( next_generation.size() ) {
        vector<MySpeciesNode*> new_generation;
        BOOST_FOREACH( MySpeciesNode *node, next_generation ) {
            if( node->hasFather() ) {
                MySpeciesNode *father = node->getFather();
                MySpeciesNode *sister = father->getSon( 0 );
                if( sister == node ) 
                    sister = father->getSon( 1 ); 

                if( not node_ids.count(father) and saw[sister->getId()] ) 
                {
                    node_ids[father] = last_branch;
                    stringstream name;
                    name << last_branch;
                    // father->setBranchProperty("ID",BppString(name.str()));
                    //father->setId(last_branch); // MY ADDITION
	                father->getInfos().breadthFirstOrder = last_branch;
                    last_branch++;

                    saw[father->getId()] = true;
                    new_generation.push_back( father );
                }
            }
        }

        next_generation.clear();
        BOOST_FOREACH( MySpeciesNode *node, new_generation ) 
            next_generation.push_back( node );
    }
}


   
/**
 * Copy the id's set by resetNodeIds to the realPostOrder property
 * (for Hali's code).
 * 
 * When the subdivision is calculated, artifical nodes will be added
 * and the nodes will be be assigned new ids. This preserve the
 * original id of the real nodes. The articial nodes will be
 * assigned the realPostOrder id of the first real descendant.
 */
void MySpeciesTree::compute_RealPostOrder() 
{
    const vector<MySpeciesNode*>& nodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
        if( node->isLeaf() || node->getNumberOfSons() == 2 )
            node->getInfos().realPostOrder = node->getId();  
    }
}


/** 
 * Helper for getPostorderTree
 *
 * @return newly created node
 */
MySpeciesNode *MySpeciesTree::getPostorderTreeAux( 
        MySpeciesNode *curNode ) ///< Current node in the recursion. 
{
    MySpeciesNode *newNode = NULL;
    if( curNode->isLeaf() ) {
        newNode = new MySpeciesNode();
        newNode->setName( curNode->getName() );
        if( curNode->hasDistanceToFather() )
            newNode->setDistanceToFather( curNode->getDistanceToFather() );
    } else { // two sons
        if( curNode->getNumberOfSons() != 2 ) 
            throw bpp::Exception("MySpeciesTree::getPostorderTree:"
                    " tree not binary");
        newNode = new MySpeciesNode();
        newNode->addSon( getPostorderTreeAux( curNode->getSon(0) ) );
        newNode->addSon( getPostorderTreeAux( curNode->getSon(1) ) );
        newNode->setBranchProperty( bpp::TreeTools::BOOTSTRAP, 
            bpp::Number<double>(curNode->getInfos().realPostOrder) );
        if( curNode->hasDistanceToFather() )
            newNode->setDistanceToFather( curNode->getDistanceToFather() );
    } 

    return newNode;
}


/**
 * Create a tree without artificial nodes with the post order id
 * in the bootstrap value.
 *
 * @return create tree
 */
MySpeciesTree* MySpeciesTree::getPostorderTree(
        bool keepOutgroup ) ///< keep outgroup if it is there
{
    MySpeciesNode *root = getRootNode();
    if( !keepOutgroup && mHasAlpha ) 
        root = root->getSon(0);
    MySpeciesNode *node = getPostorderTreeAux( root );
    if( node->hasDistanceToFather() ) // root has no father
        node->deleteDistanceToFather();
    return new MySpeciesTree( *node ); 
}



/**
 * Assign a time slice to each node, bottom up.
 *
 * Each node is assigned a higher time slice than either
 * of its children.
 */
void MySpeciesTree::assignNoSubdivisionTimeSlices() 
{
    // use -1 to indicate unassigned nodes
    vector<MySpeciesNode *> nodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, nodes )
        node->getInfos().timeSlice = -1;

    // visit nodes bottoms up to assign time slice
    int timeSlice = 0;
    if( mHasAlpha )
        timeSlice = 1;
    size_t nodeIdx = 0;
    nodes = getLeaves();
    while( nodeIdx < nodes.size() ) {
        MySpeciesNode *node = nodes[nodeIdx];
        if( node->getInfos().timeSlice != -1 ) {
            nodeIdx++;
            continue;
        }
        int sonCount = node->getNumberOfSons();
        int seenSons = 0;
        for( int i=0; i<sonCount; i++ ) 
            if( node->getSon(i)->getInfos().timeSlice != -1 )
                seenSons++;
        if( sonCount == seenSons ) {
            if( node->getInfos().isAlpha )
                node->getInfos().timeSlice = 0;
            else {
                node->getInfos().timeSlice = timeSlice++;
                // add parent to nodes to visit
                if( node->hasFather() ) 
                    nodes.push_back( node->getFather() );
            }
        }
        nodeIdx++;
    }
}


/**
 * Find changed time slices.
 *
 * Return the time slices in which the nodes don't correspond.
 *
 * @return changed time slices
 */
vector<int> MySpeciesTree::findChangedTimeSlices( 
        MySpeciesTree *otherTree, ///< tree to compare
        bool printChanges ) ///< output changes
{
    if( printChanges ) 
        cout << "=====findChangedTimeSlices====" << endl;

    if( mCorrespondence.size() == 0 )
        throw bpp::Exception("MySpeciesTree::findChangedTimeSlices:"
                        " postOrder not set");

    // check each time slice, starting with leaves
    vector<int> changedTimeSlices;
    for( int ts = 0; ts<(int)mCorrespondenceTS.size(); ts++ ) {
        bool changed = false;
        BOOST_FOREACH( int id, mCorrespondenceTS[ts] ) {
            // get node with the same post order
            MySpeciesNode *node = getNodeById( id );
            int bfOrd = node->getInfos().breadthFirstOrder;            
            MySpeciesNode *otherNode = otherTree->getNodeById( id );

            if( node->isLeaf() ) {
                // check if both are leaves and they have the same name
                if( !otherNode->isLeaf() 
                    || node->getName() != otherNode->getName() )
                {
                    if( printChanges ) 
                        cout << ts << " " << id << " (" << bfOrd << ") LEAF "
                            << otherNode->getName() << "->" << node->getName()
                            << endl;
                    changed = true;
                }
            } else {
                // check if children have same post order numbers
                size_t sonCount = node->getNumberOfSons();
                if( sonCount != otherNode->getNumberOfSons() ) {
                    if( printChanges ) 
                        cout << ts << " " << id << " (" << bfOrd 
                            << ") SON COUNT "
                            << otherNode->getNumberOfSons() << "->" 
                            << node->getNumberOfSons() << endl;
                    changed = true;
                } else if( sonCount == 2 ) {
                    for( int i=0; i<2; i++ ) {
                        MySpeciesNode *son = node->getSon(i);
                        if( son->getId()
                                != otherNode->getSon(0)->getId()
                            && son->getId()
                                != otherNode->getSon(1)->getId() )
                        {
                            if( printChanges ) 
                                cout << ts << " " << id << " (" << bfOrd 
                                   << ") 2 SON CHANGE "
                                   << otherNode->getSon(0)->getId()
                                   << ","
                                   << otherNode->getSon(1)->getId()
                                   << " -> " 
                                   << node->getSon(0)->getId()
                                   << ","
                                   << node->getSon(1)->getId()
                                   << endl;
                            changed = true;
                            break;
                        }
                    }
                } else 
                    throw bpp::Exception("MySpeciesTree::findChangedTimeSlices:"
                        " trees not binary" );
            }

            if( !changed && (node->getInfos().duplicationCost 
                    != otherNode->getInfos().duplicationCost 
                || node->getInfos().hgtCost 
                    != otherNode->getInfos().hgtCost 
                || node->getInfos().lossCost 
                    != otherNode->getInfos().lossCost ) )
            {
                if( printChanges ) {
                    if( node->getNumberOfSons() == 1 )
                        cout << "  ";
                    cout << ts << " " << id << " (" << bfOrd 
                        << ") 2 COSTS CHANGE "
                        << otherNode->getInfos().duplicationCost << ", "
                        << otherNode->getInfos().hgtCost << ", "
                        << otherNode->getInfos().lossCost << " -> " 
                        << node->getInfos().duplicationCost << ", "
                        << node->getInfos().hgtCost << ", "
                        << node->getInfos().lossCost << endl;
                    if( isAlpha( node->getId() ))
                        cout << "   =====ALPHA" << endl;
                }

                changed = true;
            }
        }
        if( changed )
            changedTimeSlices.push_back( ts );
    }
if( printChanges ) 
cout << "===================" << endl;
    return changedTimeSlices;
}



/**
 * Assign costs from the given vectors.
 *
 * This must be done before adding alpha.
 *
 * @return error string if an error, else ""
 */
string MySpeciesTree::assignCosts(
    vector<double> &duplicationCosts, ///< duplication costs
    vector<double> &hgtCosts,         ///< transfer costs
    vector<double> &lossCosts )       ///< loss costs
{
    size_t size = duplicationCosts.size();
    if( hgtCosts.size() != size )
        return "MySpeciesTree::assignCosts: given vectors of different size";
    if( lossCosts.size() != size )
        return "MySpeciesTree::assignCosts: given vectors of different size";

    vector <MySpeciesNode *> nodes = getNodes();
    for (size_t i=0; i<nodes.size(); i++) {
        if( nodes[i]->getInfos().isAlpha )
            return( "MySpeciesTree::assignCosts: " 
                    "costs need to be assigned before adding alpha" );

        size_t idx = nodes[i]->getInfos().breadthFirstOrder;
        if( idx > size )
            return( "MySpeciesTree::assignCosts: " 
                    "node index is greater than cost vector size" );
        nodes[i]->getInfos().duplicationCost = duplicationCosts[idx];
        nodes[i]->getInfos().hgtCost = hgtCosts[idx];
        nodes[i]->getInfos().lossCost = lossCosts[idx];
    }

    return "";
}


/**
 * Assign fixed costs to tree.
 */
void MySpeciesTree::assignFixedCosts(
        double duplicationCost, ///< duplication cost
        double hgtCost,         ///< transfer cost
        double lossCost )       ///< loss cost
{
    vector <MySpeciesNode *> nodes = getNodes();
    for (size_t i=0; i<nodes.size(); i++) {
        nodes[i]->getInfos().duplicationCost = duplicationCost;
        nodes[i]->getInfos().hgtCost = hgtCost;
        nodes[i]->getInfos().lossCost = lossCost;
    }
}



/**
 * Assign costs from the given file.
 *
 * The file format is [idx duplCost hgtCost lossCost].
 * Space separated. Idx is the breathFirstOrder index.
 *
 * @return true if there were no format errors
 */
bool MySpeciesTree::assignCosts( 
        string costsFileName )  ///< costs file name
{
    // costs file uses breadth first numbering
    breadthFirstreNumber();

    ifstream costsFile( costsFileName.c_str(), ios::in);  
    if( !costsFile)  
        throw bpp::IOException( "MySpeciesTree::assignCosts:"
                " failed to read from stream");
  
    vector<MySpeciesNode*> nodes = getNodes();

    string description;
    vector<double> dupCosts( nodes.size() );
    vector<double> hgtCosts( nodes.size() );
    vector<double> lossCosts( nodes.size() );
	while( !costsFile.eof() ) {
        string temp;
        getline( costsFile, temp, '\n' ); 
        if( temp == "" ) 
            continue;
        boost::char_separator<char> sep(" ");
        Tokenizer tok( temp, sep );
        Tokenizer::iterator iter=tok.begin();
        if( *iter == "" ) // blank line
            continue;
        try { 
            int idx = bpp::TextTools::toInt( *iter );
            if( idx >= (int) nodes.size() ) {
                cout << costsFileName << " index out of range: " << idx 
                    << " for line: " << temp << endl;
                return false;
            }

            iter++;
            if( iter==tok.end() ) {
                cout << costsFileName << " missing duplication cost: "
                    << temp << endl;
                return false;
            }
            dupCosts[idx] = bpp::TextTools::toDouble( *iter );
            
            iter++;
            if( iter==tok.end() ) {
                cout << costsFileName << " missing transfer cost: "
                    << temp << endl;
                return false;
            }
            hgtCosts[idx] = bpp::TextTools::toDouble( *iter );

            iter++;
            if( iter==tok.end() ) {
                cout << costsFileName << " missing loss cost: "
                    << temp << endl;
                return false;
            }
            lossCosts[idx] = bpp::TextTools::toDouble( *iter );

        } catch( bpp::Exception e ) {
              cout << costsFileName 
                   << " has a bad number format: " << temp << endl;
              return false;
        }
    }

    // assign the costs
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
        if( node->getInfos().isAlpha )
            return( "MySpeciesTree::assignCosts: " 
                    "costs need to be assigned before adding alpha" );

        int idx = node->getInfos().breadthFirstOrder;
        node->getInfos().duplicationCost = dupCosts[idx];
        node->getInfos().hgtCost = hgtCosts[idx];
        node->getInfos().lossCost = lossCosts[idx];
    }

    return true;
}



//////////////////////// ILS Functions ///////////////////////////////

/**
 * Get ILS splits for this node.a
 *
 * @return ILS splits
 */
vector< pair<int,int> > MySpeciesTree::getIlsSplits( int id ) {
    if( mIlsSplits.size() == 0 ) {
        vector< pair<int,int> > emptyVector;
        return emptyVector;
    }
    if( id >= (int) mIlsSplits.size() )
        throw bpp::Exception( "MySpeciesTree:getIlsSplits:"
                   " request for invalid id" );
    return mIlsSplits[id];
}

/**
 * Print all ids and their info.
 */
void MySpeciesTree::printIds() {

    for( size_t id=0; id<mCorrespondence.size(); id++ ) {
        cout << id;
        MySpeciesNode *node = mCorrespondence[id];
        if( node == NULL ) 
// get real ancestor id
            cout << " ils";
        else { 
            if( node->getId() != (int) id )
                cout << "  fake (" << node->getId() << ")";

            if( isAlpha( id ) ) 
                cout << " ALPHA";

            if( mTimeSlices.size() > 0 )
                cout << " ts=" << mTimeSlices[id];
            else
                cout << " ts=" << node->getInfos().timeSlice;

            vector< pair<int,int> > splits = getSplits( id );
            pair<int,int> split;
            if( splits.size() > 0 ) 
                cout << " sons:";
            BOOST_FOREACH( split, splits ) 
                cout << " " << split.first << "," << split.second;
        }

        // ils splits
        pair<int,int> ilsSplit; 
        if( mIlsSplits.size() > 0 ) {
            cout << " ils splits:";
            BOOST_FOREACH( ilsSplit, mIlsSplits[id] ) {
                cout << " " << ilsSplit.first << "/" << ilsSplit.second;
            }
        }

        cout << endl;
    }

    for( size_t ts=0; ts<mCorrespondenceTS.size(); ts++ ) {
        cout << "ts=" << ts << ":";
        BOOST_FOREACH( int id, mCorrespondenceTS[ts] ) {
            cout << " " << id;
            if( isAlpha( id ) )
                cout << "*";
        }
        cout << endl;
    }
}



/**
 * Find the common ancestor of set of species nodes (as ids).
 *
 * Also check if the set is complete, i.e. there are no missing nodes
 * in the tree rooted at the LCA to the nodes in the set. 
 *
 * @return True if complete.
 */
MySpeciesNode *MySpeciesTree::LCA(
        vector<int> &set ) ///< set of ids 
{
    if( set.size() == 1 ) 
        return mCorrespondence[set[0]];

    // get the ancestors for each node 
    vector<MySpeciesNode*> commonAncestors;
    int commonIdx = -1;
    BOOST_FOREACH( int idX, set ) {
        MySpeciesNode *node = mCorrespondence[idX];
        vector<MySpeciesNode*> ancestors;
        ancestors.push_back( node );
        while( node->hasFather() ) {
            node = node->getFather();
            ancestors.push_back( node ); 
        }

        reverse( ancestors.begin(), ancestors.end() );

        if( commonIdx == -1 ) {
            commonAncestors = ancestors;
            commonIdx = ancestors.size() - 1; 
        } else {
            int i = 0;
            while( i < (int) commonAncestors.size() && i <= commonIdx ) {
                if( commonAncestors[i] != ancestors[i] ) 
                    commonIdx = i-1;
                i++;
            }
        }
    }
    if( commonIdx == -1 ) 
        throw bpp::Exception( "MySpeciesTree::LCA: no common ancestor" );
   
    return commonAncestors[commonIdx];
}

/**
 * Check if the set is complete, i.e. there are no missing nodes
 * in the tree rooted at the LCA to the nodes in the set. 
 *
 * @return True if complete.
 */
bool MySpeciesTree::isComplete(
        MySpeciesNode *lca, ///< lca of set
        vector<int> &set ) ///< set of ids 
{
    map<int,int> sonCounts;

    vector<MySpeciesNode*> toCheck;
    BOOST_FOREACH( int id, set ) 
        toCheck.push_back( mCorrespondence[id] );

    // check if both sons are present for each father node
    while( toCheck.size() > 0 ) {
        vector<MySpeciesNode*> nextLevel;
        BOOST_FOREACH( MySpeciesNode* node, toCheck ) {
            if( !node->hasFather() )
                throw bpp::Exception( "MySpeciesTree::isComplete:"
                        "lca not right" );
            MySpeciesNode *fatherNode = node->getFather();
            if( fatherNode != lca ) {
                int fatherId = fatherNode->getId();
                sonCounts[fatherId]++;
                if( sonCounts[fatherId] == 1 )
                    nextLevel.push_back( fatherNode );
            }
        }
        toCheck = nextLevel;
    }

    pair<int,int> p;
    BOOST_FOREACH( p, sonCounts )
        if( p.second == 1 )
            return false; // node missing a child

    return true;
}

/**
 * Create ILS clades
 *
 * @return node id for the root (head) of this ils
 */
int MySpeciesTree::processClades(
        vector<int> &cladeList )
{
    if( cladeList.size() < 2 )
        throw bpp::Exception( "MySpeciesTree::processClades: only one set" );

    // find time slice of LCA
    MySpeciesNode *headNode = LCA( cladeList );
    int headTS = headNode->getInfos().timeSlice;
    if( cladeList.size() == 2 )
        headTS--; // only need headTS for all if there are ILS splits


    // for each subset
    int numSubsets = pow( (float) 2, (int) cladeList.size() );
    int idX = -1;
    vector<int> cladeTimeSlices(numSubsets);
    vector< vector<int> > cladeMap( numSubsets );
    for( int setNum=1; setNum<numSubsets; setNum++ ) {
        // create subset
        vector<int> set;
        vector<int> idxSet;
        int bit = 1;
        for( size_t i=0; i<cladeList.size(); i++ ) {
            if( bit & setNum ) {
                set.push_back( cladeList[i] ); 
                idxSet.push_back( bit );
            }
            bit *= 2;
        }

        MySpeciesNode *node;
        int timeSlice;
        bool real = true;
        if( set.size() == 1 ) {
            idX = set[0]; // existing clade
            node = mCorrespondence[idX];
            timeSlice = node->getInfos().timeSlice;
        } else {
            // find LCA to get time for this set
            node = LCA( set );
            real = isComplete( node, set );
            timeSlice = node->getInfos().timeSlice;

            if( real ) {
                // existing node
                idX = node->getId();
                mTimeSlices[idX] = timeSlice;
                mILSnodes[idX] = false;
                mCorrespondence[idX] = node;
            } else {
                // create ils node
                idX = mCorrespondence.size(); // new clade, assign id
                mTimeSlices.push_back( timeSlice );
                mILSnodes.push_back( true );
                mCorrespondence.push_back( node );
                mSplits.push_back( vector< pair<int, int> >() );
                mIlsSplits.push_back( vector< pair<int, int> >() );
            }
            // ordering of mCorrespondenceTS is important 
            //   - children before parents
            mCorrespondenceTS[timeSlice].push_back( idX );
        }
        cladeTimeSlices[setNum] = timeSlice;

        //- ADD UP TO TS V-1 if 2 clades given, else to V
        //- Save clades in structure indexed by setNum and V-timeSlice
        int childIdX = idX;
        int size = headTS-timeSlice+1;
        if( size > 0 ) {
            cladeMap[setNum].resize( size );
            cladeMap[setNum][headTS-timeSlice] = idX;
        } // else not needed

        for( int ts=timeSlice+1; ts<=headTS; ts++ ) {
            int artIdX = mCorrespondence.size(); // new clade, assign id
            mCorrespondence.push_back( node );
            mILSnodes.push_back( !real );
            mTimeSlices.push_back( ts );
            mCorrespondenceTS[ts].push_back( artIdX );
            vector< pair<int,int> > splits;
            splits.push_back( make_pair( childIdX, -1 ) );
            mSplits.push_back( splits );

            mIlsSplits.push_back( vector< pair<int, int> >() );
            cladeMap[setNum][headTS-ts] = artIdX;

            childIdX = artIdX;
        }

        // foreach possible split (no splits for set.size()==1)
        int numSplits = pow( (float) 2, (int) set.size() ) - 1;
        for( int splitNum=1; splitNum<numSplits; splitNum++ ) {
            // reconstruct set numbers for clade and anti-clade
            int bit = 1;
            int cladeSetNum = 0;
            int antiCladeSetNum = 0;
            for( size_t i=0; i<set.size(); i++ ) {
                if( bit & splitNum ) 
                    cladeSetNum += idxSet[i];
                else 
                    antiCladeSetNum += idxSet[i];
                bit *= 2;
            }
            if( cladeTimeSlices[cladeSetNum] < timeSlice
                && cladeTimeSlices[antiCladeSetNum] < timeSlice )
            {
                pair<int,int> split;
                //ADJUST first and second to be versions in previsous time slice
                split.first = cladeMap[cladeSetNum][headTS-(timeSlice-1)];
                split.second = cladeMap[antiCladeSetNum][headTS-(timeSlice-1)];
                if( split.first < split.second )  // don't add twice
                    mSplits[idX].push_back( split );
            } else {
                //COPY SPLIT UP TO ILS-HEAD TIME SLICE.
                for( int ts=timeSlice; ts<=headTS; ts++ ) {
                    pair<int,int> split;
                    split.first = cladeMap[cladeSetNum][headTS-ts];
                    split.second = cladeMap[antiCladeSetNum][headTS-ts];
                    //ADJUST first and second to be versions in SAME time slice.
                    if( split.first < split.second ) { // don't add twice
                        int tsIdX = cladeMap[setNum][headTS-ts];
                        mIlsSplits[tsIdX].push_back( split );
                    }
                }
            }
        }
    }

    return idX;
} 


/**
 * Create ILS nodes.
 *
 * @return child clades
 */
vector<int> MySpeciesTree::computeSpeciesCladesAndSplitsAux(
        MySpeciesNode *node, ///< current node in recursion
        double ilsCutoff,///< cutoff for ils branch lengths
        int maxClusterSize, ///< maximum cluster size, else abort
        vector<int> &stats ) ///< cluster sizes
{
    if( node->getNumberOfSons() == 1 ) // skip artificial
        return computeSpeciesCladesAndSplitsAux( node->getSon(0), ilsCutoff,
                        maxClusterSize, stats );

    vector<int> childClades;
    if( node->isLeaf() ) {
        int idX = node->getId(); 
        mCorrespondence[idX] = node;
        mTimeSlices[idX] = 0;
        mCorrespondenceTS[0].push_back( idX );
        mILSnodes[idX] = false;
        childClades.push_back( idX );
    } else {
        childClades =
                computeSpeciesCladesAndSplitsAux( node->getSon(0), ilsCutoff,
                                    maxClusterSize, stats );
        vector<int> childClades1 =
                computeSpeciesCladesAndSplitsAux( node->getSon(1), ilsCutoff,
                                    maxClusterSize, stats );
        BOOST_FOREACH( int id, childClades1 )
            childClades.push_back( id );

        if( !node->hasFather() || node->getDistanceToFather() > ilsCutoff ) {
            if( (int) childClades.size() > maxClusterSize ) {
                cout << "ABORT: maximum cluster size ("
                     << maxClusterSize << ") exceeded: " 
                     << childClades.size() << endl;
                exit(1);
            }
            if( childClades.size() > 2 ) 
                stats[childClades.size()]++;
            // not marked process clades
            int idX = processClades( childClades );

            // return just this clade
            childClades.clear();
            childClades.push_back( idX );
        } // else marked, return child clades
    }

    return childClades;
}



/**
 * Create ILS nodes.
 */
void MySpeciesTree::computeSpeciesCladesAndSplits(
        double ilsCutoff, ///< cutoff for ils branch lengths
        int maxClusterSize, ///< maximum cluster size, else abort
        bool verbose ) ///< print cluster sizes if true
{
    // initialize structures for leaves
    vector<MySpeciesNode*> allNodes = getNodes();
	int nodeCount = allNodes.size();
    mCorrespondence.clear();
    mCorrespondence.resize( nodeCount );
    mSplits.resize( nodeCount );
    mIlsSplits.resize( nodeCount );
    mILSnodes.resize( nodeCount );

    int maxTS = getRootNode()->getInfos().timeSlice;
    mCorrespondenceTS.clear();
    mCorrespondenceTS.resize( maxTS+1 ); 
    mTimeSlices.clear();
    mTimeSlices.resize( nodeCount );

    vector<int> stats( 100 );
    computeSpeciesCladesAndSplitsAux( getRootNode(), ilsCutoff, 
                                      maxClusterSize, stats );
    if( verbose ) {
        for( size_t i=0; i<stats.size(); i++ ) 
            if( stats[i] != 0 )
                cout << stats[i] << " ils cluster of size " << i << endl;
    }

}


///////////////////////////////////// MISC
/*
void printTree( MySpeciesNode *node,
        boost::unordered_map<string,int> &taxaNames,
        string spaces = "" )
{
    if( node->isLeaf() ) {
        string name = node->getName();
        boost::unordered_map<string,int>::iterator iter 
                = taxaNames.find( name );
        if ( iter == taxaNames.end() ) {
cout << spaces << node->getName() << " not found" << endl;
        } else {
cout << spaces << node->getName() << " good" << endl;
        }
    } else {

        cout << spaces << node->getId() << endl;
        for( size_t i=0; i<node->getNumberOfSons(); i++ )
            printTree( node->getSon(i), taxaNames, spaces+"  " );
    }
}
*/

/**
 * Print the tree in sif format to cout.
 */
void MySpeciesTree::printSif( 
        MySpeciesNode *node ) ///< current node in the recursion
{
    if( node == NULL ) 
        node = getRootNode();

    int nbrSons = node->getNumberOfSons();
    for( int i=0; i<nbrSons; i++ )  {
        cout << node->getId() << " -> " << node->getSon(i)->getId() << endl;
        printSif( node->getSon(i) );
    }
}

/**
 * Return the tree in recPhyloXML format .
 */
string MySpeciesTree::toRecPhyloXML(){
	string speciesTree_recPhyloXML_format="<recPhylo>\n<spTree>\n<phylogeny>";
	speciesTree_recPhyloXML_format=speciesTree_recPhyloXML_format + MySpeciesTree::toRecPhyloXML(this->getRootNode());
	speciesTree_recPhyloXML_format=speciesTree_recPhyloXML_format + " </phylogeny>\n</spTree>";
	return speciesTree_recPhyloXML_format;
};

string MySpeciesTree::toRecPhyloXML(MySpeciesNode *node){
	string speciesTree_recPhyloXML_format="<clade>\n<name>";
	
    if(node->getNumberOfSons()==0){
          speciesTree_recPhyloXML_format=speciesTree_recPhyloXML_format + node->getName();
    }
    else{
        speciesTree_recPhyloXML_format=speciesTree_recPhyloXML_format + bpp::TextTools::toString(node->getId());
        //cout << bpp::TextTools::toString(node->getId())<< endl;
    }    
    speciesTree_recPhyloXML_format=speciesTree_recPhyloXML_format + "</name>\n";

	for (int i=0; i< node->getNumberOfSons();i++)
		speciesTree_recPhyloXML_format=speciesTree_recPhyloXML_format + MySpeciesTree::toRecPhyloXML(node->getSon(i));
	speciesTree_recPhyloXML_format=speciesTree_recPhyloXML_format + "</clade>\n";

	return speciesTree_recPhyloXML_format;
};

//WDF
bool MySpeciesTree::isSubdivided()
{
    return mSubdivision;
}
//WDFend
