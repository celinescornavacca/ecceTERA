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
MyGeneTree implements gene trees using the BPP tree templates.

The gene trees specific function is restrictTreeToASetofTaxa, which removes
leaves that are not found in a given list (the leaves in the species tree).


*/


#include "MyGeneTree.h"


/**
 * Depth-first postorder numbering of nodes.
 */
void MyGeneTree::postorder( 
        int &pOrd, ///< current number for node
        MyGeneNode *node ) ///< current node in the recursion
{
    int nSons=node->getNumberOfSons();
    for(int i=0; i<nSons; i++)
        postorder( pOrd, node->getSon(i) );

    mCorrespondence.push_back( node );
    node->setId( pOrd++ );
}



/**
 * Assign ids using a depth-first postorder.
 *
 * This assigns ids outside of the BPP code, giving explicit
 * control of the type of ordering used.
 *
 * Also, set the mCorrespondence vector mapping ids to nodes.
 */
void MyGeneTree::assignPostOrderIds()
{
    int pOrd = 0;
    postorder( pOrd, getRootNode() );
}




/**
 * Collapse nodes to create a polytomic tree. 
 *
 * Nodes are collapsed if distances or bootstrap values are below
 * the given threshold.
 * mode =0 ==> distances;
 * mode =1 ==> BOOTSTRAP;
 *
 * @return maximum out degree of collapsed tree
 */
int MyGeneTree::collapseOnTree( 
        float threshold, ///< collapse threshold
        int mode) ///< mode (0=distances, 1=bootstrap)
{
    vector<MyGeneNode*> nodes = getNodes();
    for( size_t i=0; i<nodes.size(); i++ ) {
        if( nodes[i]->getNumberOfSons()==0 ) {
            // skip
        } else if( mode==0 ) {
            if( nodes[i]->hasDistanceToFather() 
                && nodes[i]->getDistanceToFather() < threshold )
            {
                collapseEdge( nodes[i], true, true );
            }
        } else if( mode==1 ) {
            if( nodes[i]->hasFather() 
                && nodes[i]->hasBranchProperty(bpp::TreeTools::BOOTSTRAP) ) 
            {
                double bootstrap = dynamic_cast<const bpp::Number<double> *> 
                        (nodes[i]->getBranchProperty(bpp::TreeTools::BOOTSTRAP))
                        ->getValue();
                if( bootstrap < threshold ) 
                    collapseEdge( nodes[i], true, true );
            }
        }
    }

    nodes = getNodes();
    int maxDegree = 0;
    for( size_t i=0; i<nodes.size(); i++ ) 
        if( nodes[i]->getNumberOfSons() > (size_t) maxDegree )
            maxDegree = nodes[i]->getNumberOfSons();

    return maxDegree;
}


/**
 * Check if the tree is unrooted.
 *
 * If the tree is unrooted, (root has three sons), 
 * add a root son and make it the parent of the second two original
 * root sons to create a rooted tree.
 *
 * If a node is given, use that as a root son.
 *
 * Otherwise, return false.
 *
 * @return True if tree is unrooted.
 */
bool MyGeneTree::rootTree(
        MyGeneNode *rootSon )
{
    int numSons = getRootNode()->getNumberOfSons(); 
    if( numSons == 3 ) {
        MyGeneNode *root = getRootNode();
        MyGeneNode *newSon = new MyGeneNode();

        vector<MyGeneNode*> sons;
        bool foundSon = false;
        for( int i=0; i<numSons; i++ ) {
            MyGeneNode *son = root->getSon( i );
            if( son == rootSon )
                foundSon = true;
            sons.push_back( son );
        }
        if( rootSon != NULL ) {
            if( !foundSon )
                throw bpp::Exception("MyGeneTree::rootTree: son not found");
        } else {
            rootSon = sons[0]; // randomly choose first son to keep
        }

        for( int i=0; i<numSons; i++ ) {
            if( sons[i] != rootSon ) {
                root->removeSon( sons[i] );
                newSon->addSon( sons[i] );
            }
        }
        root->addSon( newSon );

        // give the new son the same distances as the first
        if( rootSon->hasDistanceToFather() ) 
            newSon->setDistanceToFather( rootSon->getDistanceToFather() );

        // copy branch properties
        vector<string> names = rootSon->getBranchPropertyNames();
        for( size_t j=0; j < names.size(); j++ ) 
            newSon->setBranchProperty( names[j], 
                        *rootSon->getBranchProperty(names[j]));
       
        resetNodesId();
        return true;
    } else {
        return false;
    }
}

/**
 * Get node value for collapsing threshold.
 *
 * return value or -1 if no there
 */
double getNodeValue( MyGeneNode *node, int mode ) 
{
    double value = -1;
    if( mode==0 ) {
        if( node->hasDistanceToFather() )
            value = node->getDistanceToFather();
    } else if( mode==1 ) {
        if( node->hasFather() 
            && node->hasBranchProperty(bpp::TreeTools::BOOTSTRAP) ) 
        {
            value = dynamic_cast<const bpp::Number<double> *> 
               (node->getBranchProperty(bpp::TreeTools::BOOTSTRAP))->getValue();
        }
    } else 
        throw bpp::Exception( "MyGeneTree::getNodeValue: invalid mode" );
    return value;
}

/**
 * Find a node above or equal to threshold if mode != -1 and root there
 * mode =0 ==> distances;
 * mode =1 ==> BOOTSTRAP;
 *
 * If no such node, add a node for rooting.
 *
 * Note: the tree is invalid if it is already rooted, but the two
 * root sons have different values (distanceds or bootstrap).
 *
 * @return true, unless tree is invalid
 *
 */
bool MyGeneTree::rootWithThreshold(
        float threshold, ///< collapse threshold
        int mode ) ///< mode (0=distances, 1=bootstrap)
{
    if( getRootNode()->getNumberOfSons() == 2 ) {
        double value0 = getNodeValue( getRootNode()->getSon( 0 ), mode );
        double value1 = getNodeValue( getRootNode()->getSon( 1 ), mode );
        if( value0 != value1 ) {
            //return false;
            // use higher value and set other son to it
            int changeSon = 1;
            int otherSon = 0;
            if( value0 < value1 ) {
                changeSon = 0;
                otherSon = 1;
                value0 = value1;
            }
            MyGeneNode *node = getRootNode()->getSon( changeSon );
            MyGeneNode *otherNode = getRootNode()->getSon( otherSon );
            if( mode == 0 ) 
                node->setDistanceToFather( value0 );
            else {
                node->setBranchProperty( bpp::TreeTools::BOOTSTRAP,
                     *otherNode->getBranchProperty(bpp::TreeTools::BOOTSTRAP) );
            }
        }
        if( value0 >= threshold ) 
            return true; // already good
    }
   

    // find a node above threshold
    vector<MyGeneNode*> nodes = getNodes();
    MyGeneNode *newRoot = NULL;
    BOOST_FOREACH( MyGeneNode *node, nodes ) {
        double value = getNodeValue( node, mode );
        if( value != -1 && value >= threshold ) {
            newRoot = node;
            break;
        }
    }

    if( newRoot != NULL ) {
        // found a new root
        MyGeneNode *father = newRoot->getFather();
        rootAt( newRoot ); // tree unrooted at newRoot
        rootTree( father ); // root with old father as a root child 

    } else {
        // no node above threshold, root and add an outgroup
        rootTree(); // roots tree if unrooted

        MyGeneNode *newRoot = new MyGeneNode(); 
        MyGeneNode *outGroup = new MyGeneNode(); 
        outGroup->setName( "outgroup" );
        
        newRoot->addSon( getRootNode() );
        newRoot->addSon( outGroup );

        setRootNode( newRoot ); 
    }


    resetNodesId();

    return true;
}



/** 
 *  Remove the leaves (and their branches) that are not in
 *  the given map.
 *
 *  @return False if deleted all the nodes
 */
bool MyGeneTree::restrictTreeToASetOfTaxa(
         boost::unordered_map<string,int> &taxaNames, 
            ///< map of leaf names
         char charSep, ///< character seperator for leaf names
         bool verbose ) ///< print debugging info
{
    bool hasLeaves = true;

    vector <MyGeneNode *> leaves = getLeaves();
    BOOST_FOREACH( MyGeneNode *leaf, leaves ) {

        string name = leaf->getName();
        if( charSep != 'x' ) {
            // extract the gene name (part before seperator)         
            size_t pos = name.find(charSep); 
            name = name.substr(0,pos);
        }

        boost::unordered_map<string,int>::iterator iter 
            = taxaNames.find( name );
        if ( iter == taxaNames.end() ) {
			hasLeaves = collapseEdge( leaf );
            //if( verbose ) 
			    cout << "Gene not taken into account (not in S) <" 
                     << name << ">" << endl;
            if( hasLeaves ) 
                resetNodesId();
		}
	}
    
    return hasLeaves;
}




/**
 * Reroot the tree, using given set of leaves as an outgroup.
 *
 * @return False if outgoup not monophyletic.
 */
bool MyGeneTree::reroot( 
        vector<string> &outGroup, ///< set of leaves forming new out group
        int proportion ) ///< ratio of distances of root sons
{
    vector<string> leavesTree = getLeavesNames();
    vector<string> remainingTaxa;
    bpp::VectorTools::diff( leavesTree, outGroup, remainingTaxa );
    if( remainingTaxa.size() == 0 )  // already good
        return true;
   
    // move up the tree until the subtree contains all of the outgroup
    newOutGroup( getNode( remainingTaxa[0] ) );
    MyGeneNode *newRoot = getNode( outGroup[0] );
    vector<string> tempLeaves = bpp::TreeTemplateTools::getLeavesNames( 
                                        *newRoot );
    while( newRoot->hasFather() && 
           !(bpp::VectorTools::containsAll(tempLeaves, outGroup)))
    {   
        newRoot = newRoot->getFather();
        tempLeaves = bpp::TreeTemplateTools::getLeavesNames( *newRoot );
    }
      
    if( tempLeaves.size() == outGroup.size() ) {
        // outgroup is monophyletic 
        newOutGroup( newRoot );
    } else {
        // check if sons are monophyletic
        for( size_t f=0; f<newRoot->getNumberOfSons(); f++ ) {
            tempLeaves = bpp::TreeTemplateTools::getLeavesNames(
                                        *newRoot->getSon(f));
            vector<string> diff;
            bpp::VectorTools::diff(outGroup, tempLeaves, diff);
            if(diff.size() != 0 && diff.size() != tempLeaves.size() ) {
                //The proposed outgroup is not monophyletic. 
                // The analysis for this tree is interrupted.
                // No more outgroup are analysed.
                return false;
            }
        }

        tempLeaves = bpp::TreeTemplateTools::getLeavesNames(* newRoot);
        if( tempLeaves.size() != leavesTree.size() ) {
            MyGeneNode *cloneRoot 
                    = bpp::TreeTemplateTools::cloneSubtree<MyGeneNode>( 
                                *newRoot );
            MyGeneTree *low = new MyGeneTree( *cloneRoot );
            newOutGroup( newRoot );
            vector<string> tempLeaves2 = bpp::TreeTemplateTools::getLeavesNames(
                                            *(getRootNode())->getSon(0));
            std::sort( tempLeaves2.begin(), tempLeaves2.end() );
            vector<string> intersection = bpp::VectorTools::vectorIntersection(
                                        tempLeaves2, outGroup );
            MyGeneNode *sonUpper;
            if( intersection.size() != 0 ) 
                sonUpper = getRootNode()->getSon(1);
            else
                sonUpper = getRootNode()->getSon(0);

            int ident = bpp::TreeTools::getMaxId( *low, low->getRootId() );
            vector<MyGeneNode*> nodesTemp
                = bpp::TreeTemplateTools::getNodes( *sonUpper );
            for(unsigned int F = 0; F < nodesTemp.size(); F++)
                (*nodesTemp[F]).setId(ident + F + 1);
            low->getRootNode()->addSon( sonUpper );

            //tree = low; // replaced by newroot and deleting old tree
            MyGeneNode *oldRoot = getRootNode();
            setRootNode( low->getRootNode() );
            deleteSubtree( oldRoot );
        }
    }

    resetNodesId();

    // Change the ratio of the distances of the two branches 
    // coming from the root
    MyGeneNode* root = getRootNode();
    MyGeneNode* son1 = root->getSon(0); 
    MyGeneNode* son2 = root->getSon(1);
    double dist1 =son1->getDistanceToFather();
    double dist2 =son2->getDistanceToFather();
    son1->setDistanceToFather((dist1+dist2)/(proportion+1));
    son2->setDistanceToFather((dist1+dist2)/(proportion+1)*proportion);

    return true;
}



/**
 * Read newick trees from a file.
 *
 * @return a vector of trees
 */
vector<MyGeneTree*> MyGeneTree::readMyGeneTrees( 
    const char *treePathChar, ///< path to file with newick trees
    string &errString, ///< error description
    bool readBootstrap ) ///< read bootstrap values (default=false)
{
    vector<string> treeStrings = readNewickStrings( treePathChar );

    vector<MyGeneTree*> trees;
    int lineNumber=0;
    BOOST_FOREACH( string description, treeStrings ) {
        lineNumber ++;
        MyGeneTree *tree = new MyGeneTree( description, errString, 
                                           readBootstrap );
        if( errString != "" ) {
            ostringstream str;
            size_t size = trees.size()+1;
            str << " (error in tree " << size << ")";
            errString += str.str();
            return trees;
        }
        if(tree->getNumberOfLeaves()>=2)
            trees.push_back( tree );
        else    
            cout << "Gene tree number " << lineNumber << " has only one leaf,  skipped\n";
    }

            

    return trees;
}


/**
 * Construct a subtree for the given id, i.e.,
 * only nodes where biCompMap[nodeId] == biCompId are included.
 *
 * @return a switching tree for the given biCompId
 */
MyGeneTree* MyGeneTree::switching(
        MyGeneNode *root, ///< root of network (speciation node)
        const int biCompId, ///< biconnected set to use 
        const vector<int> &biCompMap ) ///< mapping of node ids to set
{
    // find all hybrid nodes and record which nodes (of all) are in network
    vector<int> hybridIndices( getNumberOfNodes() ); 

    // create swtiching
    // maps id->new
    vector<MyGeneNode*> newNodes( getNumberOfNodes() ); 
    vector<bool> deletedNodes( getNumberOfNodes(), false ); 
    allSwitchingsAux( root, 0, biCompId, biCompMap,
                    hybridIndices, newNodes, deletedNodes );
    MyGeneNode *newRoot = newNodes[root->getId()];

    return new MyGeneTree( *newRoot );
}
