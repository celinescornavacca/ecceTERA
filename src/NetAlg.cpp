/*
@file netAlg.cpp
@author Celine Scornavacca
@author Edwin Jacox
@version 1.0 
@date 15/10/2015

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

Network algorithms described in "Efficient algorithms for reconciling
gene trees and species networks via duplication and loss events."
Thu-Hien To and Celine Scornavacca. 
2015

**/

#include <boost/foreach.hpp>

#include "NetAlg.h"

const double COST_DIFF = 0.00001;
bool gDebug = false;


/**
 * Count path losses. A loss for each parent with two children.
 *
 * @return number of losses
 */
int NetAlg::countPathLosses(
        vector<MySpeciesNode*> &path, ///< a list of nodes
        MySpeciesNode *lca )          ///< lca
{
    int losses = 0;
    for( size_t i=1; i<path.size(); i++ ) 
        if( path[i]->getNumberOfSons() > 1 )
            losses++;
    
    return losses;
}

/**
 * Count path transfers.
 *
 * @return number of transfers
 */
int NetAlg::countPathTransfers(
        vector<MySpeciesNode*> &path, ///< a list of nodes
        MySpeciesNode *lca )          ///< lca
{
    int transfers = 0;
    for( size_t i=0; i<path.size()-1; i++ ) {
        if( (path[i]->getInfos().secondaryFather !=NULL)&& (path[i]->getInfos().secondaryFather->getId()== path[i+1]->getId())){
            transfers++;      
        }    
    } 
    return transfers;
}

/**
 * Find the least common ancestor of the two given nodes and return
 * the path to it from each.
 *
 * @return The least common ancestor.
 */
MySpeciesNode *NetAlg::commonAncestor(
    MySpeciesNode *node1,      ///< a node
    MySpeciesNode *node2,      ///< another node
    vector<MySpeciesNode*> &spPath1,  ///< path of node1 to lca
    vector<MySpeciesNode*> &spPath2 )  ///< path of node2 to lca
{
    spPath1.push_back( node1 );
    spPath2.push_back( node2 );
    MySpeciesNode *lca = NULL;
    bool foundIt = false;
    while( !foundIt && (spPath1.back()->hasFather() 
            || spPath2.back()->hasFather()) ) 
    {
        // get next ancestor of node1 and compare to those of node2
        if( spPath1.back()->hasFather() ) {
            lca = spPath1.back()->getFather();
            spPath1.push_back( lca );
            for( size_t i=0; i<spPath2.size(); i++ ) 
                if( lca == spPath2[i] ) {
                    spPath2.resize(i+1);
                    foundIt = true;
                    break;
                }
        }
        // get next ancestor of node2 and compare to those of node1
        if( !foundIt && spPath2.back()->hasFather() ) {
            lca = spPath2.back()->getFather();
            spPath2.push_back( lca );
            for( size_t i=0; i<spPath1.size(); i++ ) 
                if( lca == spPath1[i] ) {
                    spPath1.resize(i+1);
                    foundIt = true;
                    break;
                }
        }
    }
    if( !foundIt ) 
        throw bpp::Exception("NetAlg::LCAaux: no LCA found" );

    return lca;
}

/**
 * Recursively find the least common ancestor of internal nodes
 * for a reconciliation between a gene tree and a species tree.
 *
 * The leaves should already be reconciled before calling.
 *
 * @return cost of reconciliation for the subtree rooted at geneNode
 */
double NetAlg::LCAaux(
    MyGeneNode *geneNode,   ///< current node in the recursion
    int &numLosses,     ///< number of losses in the reconciliation
    int &numDupli,      ///< number of duplications in the reconciliation
    int &numTransfers,    ///< number of transfers in the reconciliation 
    vector<MySpeciesNode*> &recon ) ///< reconciliation to return
{
    if( geneNode->isLeaf() ) 
        return 0;

    // recursion
    MyGeneNode *geneSon1 = geneNode->getSon(0);
    MyGeneNode *geneSon2 = geneNode->getSon(1);
    double cost1 = LCAaux( geneSon1, numLosses, numDupli, numTransfers, recon );
    double cost2 = LCAaux( geneSon2, numLosses, numDupli, numTransfers, recon );

    // duplication with no losses
    if( recon[geneSon1->getId()] == recon[geneSon2->getId()] ) {
        recon[geneNode->getId()] = recon[geneSon1->getId()];
        numDupli++;
        return mDupliCost + cost1 + cost2;
    }

    // find common ancestor of species nodes
    vector<MySpeciesNode*> spPath1;
    vector<MySpeciesNode*> spPath2;
    MySpeciesNode *lca = commonAncestor( recon[geneSon1->getId()],
                                  recon[geneSon2->getId()], spPath1, spPath2 );
    recon[geneNode->getId()] = lca;

    // count duplications and losses for subtree
    int dupli = 0;
    int losses = 0;
    int transfers = 0;

    if( lca == recon[geneSon1->getId()] || lca == recon[geneSon2->getId()] ) 
    {
        dupli = 1;
        if( lca != recon[geneSon1->getId()] ){ 
            losses += countPathLosses( spPath1, lca );
            
        }
        else if( lca != recon[geneSon2->getId()] ){ 
            losses += countPathLosses( spPath2, lca );
        }
    } else {
        losses += countPathLosses( spPath1, lca ) - 1; // minus speciation 
        losses += countPathLosses( spPath2, lca ) - 1;
    }
		
	transfers	+= countPathTransfers( spPath1, lca );
	transfers	+= countPathTransfers( spPath2, lca );
		
    // add subtree values to total
    numLosses += losses;
    numDupli += dupli;
    numTransfers += transfers;
    
    //cout << numTransfers << endl;

    return dupli*mDupliCost + losses*mLossCost + transfers*mTransferCost + cost1 + cost2;
}


/**
 * Reconcile the species tree and the gene tree using the least
 * common ancestor method.
 *
 * @return cost of the reconciliation
 */
double NetAlg::LCA(
    MyNetwork *speciesNetwork,///< species tree to reconcile
    MyGeneTree *geneTree,   ///< gene tree to reconcile
    int &numLosses,     ///< number of losses in the reconciliation
    int &numDupli,      ///< number of duplications in the reconciliation
    int &numTransfers,    ///< number of transfers in the reconciliation 
    vector<MySpeciesNode*> &recon ) ///< reconciliation to return
{
    numLosses = 0;
    numDupli = 0;


    if( geneTree->getNumberOfNodes() < 2 ) 
       throw bpp::Exception( "NetAlg::LCA gene tree has less than two nodes" ); 

    // match leaves between the species tree and gene tree
    size_t speciesIdx = 0;
    vector<MySpeciesNode*> speciesLeaves = speciesNetwork->getSortedLeaves();
    vector<MyGeneNode*> geneLeaves = geneTree->getSortedLeaves();
    MySpeciesNode *prevSpecies = NULL;
    bool allSame = true;
    BOOST_FOREACH( MyGeneNode *geneLeaf, geneLeaves ) {
        string leafName = geneLeaf->getName();
        size_t pos = leafName.find( mCharSep );
        string speciesName = leafName.substr( 0, pos );
        while( speciesIdx < speciesLeaves.size() 
               && speciesLeaves[speciesIdx]->hasFather()
               && speciesName != speciesLeaves[speciesIdx]->getName() )
        {
            speciesIdx++;
        }
        if( speciesIdx >= speciesLeaves.size() ) {
            cout << leafName << " with species " << speciesName << endl;
            throw bpp::Exception( 
                    "NetAlg::LCA: Could not find species for gene" );
        }
        recon[geneLeaf->getId()] = speciesLeaves[speciesIdx];

        if( prevSpecies == NULL )
            prevSpecies = speciesLeaves[speciesIdx];
        else if( allSame && prevSpecies != speciesLeaves[speciesIdx] )
            allSame = false;
    }

    double cost = 0;
    /*
    if( geneTree->getNumberOfNodes() == 2 ) {
        // in two node trees, root is marked as a leaf
        int geneRootId = geneTree->getRootNode()->getId();
        int geneLeafId = geneTree->getRootNode()->getSon(0)->getId();
        recon[geneRootId] = recon[geneLeafId]->getFather();
        while( recon[geneRootId]->hasFather() 
               && recon[geneRootId]->getNumberOfSons() != 2 )
        {
            recon[geneRootId] = recon[geneRootId]->getFather();
        }
        if( recon[geneRootId]->getNumberOfSons() == 2 ) {
            numLosses++;
            cost += mLossCost;
        }
    */

    if( allSame ) {
        vector<MyGeneNode*> geneNodes = geneTree->getNodes();
        BOOST_FOREACH( MyGeneNode *geneNode, geneNodes) {
            recon[geneNode->getId()] = prevSpecies;
            if( geneNode->getNumberOfSons() == 2 ) {
                numDupli++;
                cost += mDupliCost;
            }
        }
    } else {
        MyGeneNode *root = geneTree->getRootNode();
        if( root->isLeaf() ) {
            geneTree->printMe();
            throw bpp::Exception( "NetAlg::LCA: Root is a leaf." );
        }
        
        // find LCA recursively
        cost = LCAaux( geneTree->getRootNode(),
                       numLosses, numDupli, numTransfers, recon );
    }

    return cost;
}



/**
 * Calculate the minimum switching by checking the cost of each switching.
 *
 * @return minimum cost
 */
double NetAlg::runBrute( 
    int &numLosses, ///< number of losses in the minimum switching
    int &numDupli,  ///< number of duplications in the minimum switching
    int &numTransfers ///< number of transfers in the minimum switching
    ) 
{
    // create swtiching using all nodes (biComp=0 for all nodes)
    vector<int> biCompMap( mSpeciesNetwork->getNumberOfNodes(), 0 );
    vector<MyNetwork*> switchings = mSpeciesNetwork->allSwitchings(
                                    mSpeciesNetwork->getRootNode(), 0, 
                                    biCompMap );
            
            

        
//if( gDebug ) 
    cout << switchings.size() << " switchings" << endl;

    // check the cost of each switching
    vector<MySpeciesNode*> recon;
    double lowestCost = numeric_limits<double>::max();
    BOOST_FOREACH( MyNetwork*switching, switchings ) {
        
        
        int switchNumLosses =0 ;
        int switchnumDupli=0 ;
        int switchnumTransfers=0 ;
        vector<MySpeciesNode*> switchRecon( mGeneTree->getNumberOfNodes() );
        double cost = LCA( switching, mGeneTree, switchNumLosses, 
                           switchnumDupli, switchnumTransfers, switchRecon );
//if( gDebug ) 
    cout << "LCA cost=" << cost << " losses=" << switchNumLosses 
         << " dupli=" << switchnumDupli << " transfers=" << switchnumTransfers << endl;

        if( cost < lowestCost ) {
            numLosses = switchNumLosses;
            numDupli = switchnumDupli;
            numTransfers =switchnumTransfers;
            lowestCost = cost;
            recon = switchRecon;
        }
    }

    // error check
    for( size_t i=0; i<recon.size(); i++ ) {
        if( recon[i] == NULL ) {
            cerr << "gene " << i << " not assigned" << endl;
            throw bpp::Exception(
                    "NetAlg::runBrute: reconciliation not complete" );
        }
if( gDebug )
    cout << "gene " << i << " -> " << recon[i]->getId() << endl;
    }

    BOOST_FOREACH( MyNetwork *tree, switchings ) 
        delete tree;
        

    return lowestCost;
}


/**
 * Recursively make a tree of biconnected components from original tree.
 *
 * Internal node ids are the biconnected component number (biNum). 
 * Leafs have the same name as the orignal tree.
 *
 * @return Created node
 */
MySpeciesNode *NetAlg::makeBicompTree(
    MySpeciesNode *node,           ///< node of original tree
    MySpeciesNode *parent,         ///< parent in new tree
    vector<bool> &seenIt,   ///< nodes visited
    vector<int> &biCompMap ) ///< given map of species tree to bicomponents
{
    int id = node->getId();
    if( seenIt[id] )
        return NULL;

    seenIt[id] = true;
    int biNum = biCompMap[id];

    // create a new node the first time each biNum is seen
    MySpeciesNode *newNode = NULL;
    if( parent == NULL || biNum != parent->getId() ) {
        newNode = new MySpeciesNode();
        if( node->isLeaf() ) {
            newNode->setName( node->getName() );
            newNode->setId( -1 );
        } else 
            newNode->setId( biNum );
        parent = newNode;
    }

    // recursion
    for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
        MySpeciesNode *son = makeBicompTree( node->getSon(i), parent, 
                            seenIt, biCompMap );
        if( son != NULL )
            parent->addSon( son );
    }

    return newNode;
}


/**
 * Recursive part of makeGn
 *
 * @return New node corresponding to the given node.
 */
MyGeneNode *NetAlg::makeGnAux(
    MyGeneNode *node,           ///< gene tree node
    vector<MySpeciesNode*> &recon, ///< reconciliation with bicomp tree
    vector<int> &mu,        ///< mapping of gene nodes to bicomp id
    int &cnt,               ///< counter for ids
    MySpeciesNode *parentBiCompNode ) ///< parent bicomp node
{
    MyGeneNode *newNode = new MyGeneNode();
    MyGeneNode *topNode = newNode; // node to return

    // add artifical nodes if necessary
    MySpeciesNode *biCompNode = recon[node->getId()];
    int biCompId = biCompNode->getId();
    if( biCompNode->isLeaf() )
        biCompId = biCompNode->getFather()->getId();

    newNode->setId( cnt++ );
    mu.push_back( biCompId );
if( gDebug ) {
cout << newNode->getId() << " mu=" << biCompId << " from " << node->getId() << endl;
if( parentBiCompNode != NULL )
cout << "gN bicomp=" << biCompId << " parent: " << parentBiCompNode->getId() << endl;
}
    if( node->hasFather() && biCompNode != parentBiCompNode ) {
        while( biCompNode->getFather() != parentBiCompNode )
        {
            biCompNode = biCompNode->getFather();
            biCompId = biCompNode->getId(); // can't be a leaf
            MyGeneNode *artificialNode = new MyGeneNode();
            artificialNode->setId( cnt++ );
            mu.push_back( biCompId );

            // add above
            artificialNode->addSon( topNode );
            topNode = artificialNode;
if( gDebug ) 
cout << topNode->getId() << " mu=" << biCompId << " from " << node->getId() << " ***" << endl;
        }
    }

    if( node->isLeaf() ) {
        newNode->setName( node->getName() );
        mu[newNode->getId()] = -1;
if( gDebug ) 
cout << "changing " << (cnt-1) << " to -1" << endl;
    }

    // recursion
    for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
        newNode->addSon( makeGnAux( node->getSon(i), recon, mu, cnt,
                                    recon[node->getId()] ) );
    }

    return topNode;
}


/**
 * Create a gene tree corresponding to the network species tree.
 *  
 *  Map lowest Bi whose leaves encompass gene subtree
 *  Note: each Bi might map to several (possibly disconnected) nodes
 *
 * @return the new gene tree
 */
MyGeneTree *NetAlg::makeGn(
    vector<int> &mu, ///< mapping of gene nodes to bicomp id
    vector<int> &biCompMap ) ///< map of species tree to bicomponents
{
    // 1. make a tree of the biconnected componets
    vector<bool> seenIt( mSpeciesNetwork->getNumberOfNodes(), false );
    MySpeciesNode *root = makeBicompTree( mSpeciesNetwork->getRootNode(), NULL, 
                              seenIt, biCompMap );
    MyNetwork *biCompTree = new MyNetwork( *root );
if( gDebug ) {
cout << "========biCompTree======================" << endl;
biCompTree->printMe();
}
    // 2. use LCA to map gene tree to bicomp tree
    int notUsed1;
    int notUsed2;
    int notUsed3;
    vector<MySpeciesNode*> recon( mGeneTree->getNumberOfNodes() );
//cout << "calling LCA" << endl;
    LCA( biCompTree, mGeneTree, notUsed1, notUsed2, notUsed3, recon );
//cout << "out calling LCA" << endl;

if( gDebug ) {
cout << "=== given gene tree ====" << endl;
mGeneTree->printMe();
cout << "=== bicomp recon===" << endl;
}
    // confirm that all genes are there
    vector<MyGeneNode*> nodes = mGeneTree->getNodes();
    BOOST_FOREACH( MyGeneNode *node, nodes ) {
        int id = node->getId();
        if( recon[id] == NULL ) {
            cerr << "gene " << id << " not assigned" << endl;
            throw bpp::Exception(
                    "NetAlg::makeGn: reconciliation not complete" );
        }
        if( gDebug )
            cout << "gene " << id << " -> " << recon[id]->getId() << endl;
    }

    
    // 3. create Gn, adding missing bicomponents
    int cnt = 0;
    MyGeneNode *gnRoot = makeGnAux( mGeneTree->getRootNode(), recon, mu, cnt );
    delete biCompTree;

    return new MyGeneTree( *gnRoot );
}

/**
 * Calculate the minimum cost of reconciliation of the bicomponent, 
 * consisting of the given switchings, and the given gene gene trees.
 *
 * The cost of each switching is the sum of reconciliations cost with
 * each gene tree. The best cost is returned.
 *
 * @return Best cost
 */
double NetAlg::biCompCost(
        vector<MyNetwork*> &switchings, ///< switchings of a bicomponent
        vector<MyGeneTree*> &biCompGeneTrees, ///< corresponding gene trees
        bool isRootBiComp,      ///< true if bicomponent is the root
        int &bestLosses,        ///< number of losses for the minium cost
        int &bestDupes,         ///< number duplications for the min cost
        int &bestTransfers,         ///< number transfers for the min cost
        int &bestSwitching) ///< number of best switching 
{
    // for each switching
    int cnt = 0;
    double bestCost = -1;    
    
    int bs=0;
    BOOST_FOREACH( MyNetwork *switching, switchings ) {
    	
		if( gDebug ) 
			cout << "== switching === " << cnt++ << endl;
		//switching->printMe();
       // sum costs (DL) of each gene tree with switching 
       int sumLosses = 0;
       int sumDupes = 0;
       int sumTransfers = 0;
       double sumCost = 0;
       BOOST_FOREACH( MyGeneTree* g, biCompGeneTrees ) {
            // LCA algorithm
            int lcaLosses=0;
            int lcaDupes=0;
            int lcaTransfers=0;
            vector<MySpeciesNode*> recon( g->getNumberOfNodes() );
            double cost = LCA( switching, g, lcaLosses, lcaDupes, lcaTransfers, recon );

            // add losses (=distance) and transfers to switch root, if not bicomp root
            if( !isRootBiComp ) {
                MySpeciesNode *node = recon[g->getRootNode()->getId()];
                while( node != switching->getRootNode() ) {
if( gDebug ) 
cout << "**root up " << node->getId() << endl;
                    
                    if( (node->getInfos().secondaryFather!=NULL)&& (node->getInfos().secondaryFather->getId()== node->getFather()->getId())){
                    	lcaTransfers++; 
                        cost += mTransferCost;
                    }
                    
                    node = node->getFather();
                    if( node->getNumberOfSons() > 1 ) {
                        lcaLosses++; 
                        cost += mLossCost;
                    }

                }
            }

if( gDebug ) 
cout << "LCA cost=" << cost << " losses=" << lcaLosses << " dupli=" << lcaDupes << " transfers=" << lcaTransfers << endl;
for( size_t i=0; i<recon.size(); i++ )
if( gDebug && recon[i] != NULL )
cout << "   gene " << i << " ->sp " << recon[i]->getId() << endl;
            sumLosses += lcaLosses;
            sumDupes += lcaDupes;
            sumTransfers += lcaTransfers;
            sumCost += cost;
       }
       if( bestCost == -1 || sumCost < bestCost ) {
           bestLosses = sumLosses;
           bestDupes = sumDupes;
           bestTransfers = sumTransfers;
           bestCost = sumCost;
           bestSwitching= bs;
       }
if( gDebug ) 
cout << "SWITCH cost=" << bestCost << " losses=" << bestLosses << " dupli=" << bestDupes << " transfers="   <<  bestTransfers << endl << endl;
    bs++;
    }

    return bestCost;
}


/**
 * Find the switching of mSpeciesNetwork whose reconciliation with mGeneTree
 * has the least cost.
 *
 * @return The least cost.
 */
double NetAlg::runMinSwitch( 
    int &numLosses, ///< number of losses for the minimum switching
    int &numDupli, ///< number of duplications for the minimum switching
    int &numTransfers,  ///< number of transfers  for the minimum switching
     vector<std::pair <int,int> > & edgesBestSwitchings)  ///< bestSwitching to store
{
    // calculate biconnected components B and create map speciesId->biNum
    vector<MySpeciesNode *> biRoots;
    vector<int> biCompMap = mSpeciesNetwork->getBiconnectedCompMap( biRoots );

    // map B to gene tree (mu), adding artifical nodes (GN)
    vector<int> mu; // mapping of gN id to bi comp set
    MyGeneTree *gN = makeGn( mu, biCompMap );

    // get highest node of Gn associated with each bicomp
    vector< vector<MyGeneNode*> > gNroots( biRoots.size() );
    size_t rootBiSetId = 0;
    BOOST_FOREACH( MyGeneNode *node, gN->getNodes() ) {
        int setId = mu[node->getId()];
        if( !node->hasFather() ) {
            //add root
if( gDebug )
cout << node->getId() << " root for " << setId << endl;
            gNroots[setId].push_back( node );
            rootBiSetId = mu[node->getId()];
        } else if( !node->isLeaf() ) {
            int fatherSetId = mu[node->getFather()->getId()];
            if( setId != fatherSetId ) {
                gNroots[setId].push_back( node );
if( gDebug )
cout << node->getId() << " root for " << setId << endl;
            }
        }
    }

    // Calculate and sum cost for each bicomponent.
    double totalCost = 0;
    for( size_t i=0; i<biRoots.size(); i++ ) {
        if( biRoots[i] == NULL ) 
            throw bpp::Exception( "runMinSwitch: bicomp set with no root" );

        vector<MyNetwork*> switchings = 
            mSpeciesNetwork->allSwitchings( biRoots[i], i, biCompMap );
if( gDebug ) 
cout << "========= proc biComp set " << i << " has " << switchings.size()
    << " switchings =====" << endl;

        // get gene trees associated with this bicomp
        vector<MyGeneTree*> biCompGeneTrees;
        BOOST_FOREACH( MyGeneNode *root, gNroots[i] ) {
            MyGeneTree* tree = gN->switching( root, i, mu );
            tree->resetNodesId();
            biCompGeneTrees.push_back( tree );
if( gDebug ) {
cout << " == bicomp Gn tree for root " << root->getId() << endl;
tree->printMe();
}
        }

        // calculate the cost for the bicomponent
        int bestLosses = 0;
        int bestDupes = 0;
        int bestTransfers =0;
        int bestSwitching =0;

        totalCost += biCompCost( switchings, biCompGeneTrees, 
                                      rootBiSetId==i, bestLosses, bestDupes,bestTransfers, bestSwitching );
           
		
		vector<MySpeciesNode*> allNodesTemp = switchings[bestSwitching]->getNodes();
		

		BOOST_FOREACH( MySpeciesNode *node, allNodesTemp ){
			if (node->hasFather() && node->getInfos().secondaryFather!=NULL){ // hybrid eges that are kept
				//cout << node->getFather()->getId() << " " << node->getId() << endl;
				std::pair <int,int> edge = std::make_pair(node->getFather()->getId(),node->getId());
				edgesBestSwitchings.push_back(edge);
			}	
		}
         
                              
        numLosses += bestLosses;
        numDupli += bestDupes;
        numTransfers +=  bestTransfers;
        
        BOOST_FOREACH( MyNetwork *tree, switchings )
            delete tree;
        BOOST_FOREACH( MyGeneTree *tree, biCompGeneTrees )
            delete tree;
    }
    


    delete gN;

    return totalCost;
}


/**
 * Fill the ancestors vector with ancestors of node.
 */
void markParents(
        MySpeciesNode *node,       ///< a node
        vector<MySpeciesNode*> &otherParents, ///< secondary parents
        vector<MySpeciesNode*> &ancestors ) ///< vector to fill
{
    if( ancestors[node->getId()] == NULL) {
        ancestors[node->getId()] = node;
        if( node->hasFather() )
            markParents( node->getFather(), otherParents, ancestors );
        if( otherParents[node->getId()] != NULL ) 
            markParents( otherParents[node->getId()], otherParents, ancestors );
    }
}


/**
 * Calculate the number of paths from node1 to node2, counting the
 * number of times each nodes occurs in a path.
 *
 * @return number of paths from node1
 */
int allPaths( 
        MySpeciesNode *node1, ///< start node
        MySpeciesNode *node2,  ///< end node
        vector<MySpeciesNode*> &ancestors, ///< ancestors of node1
        vector<int> &nodeCounts ) ///< number of times each node occurs
{
    size_t paths = 0;

    if( node1 == node2 ) {
        paths = 1;
    } else {
        for( size_t i = 0; i<node1->getNumberOfSons(); i++ ) {
            MySpeciesNode *son = node1->getSon( i );
            if( ancestors[son->getId()] != NULL )
                paths += allPaths( son, node2, ancestors, nodeCounts );
        }
    }

    nodeCounts[node1->getId()] += paths;
    return paths;
}


/**
 * Brute force calculation of computeMin. If the results 
 * don't match those of computeMin (in toCheck) an exception is thrown
 */
void NetAlg::checkComputeMin( 
        MySpeciesNode *y,  ///< one node
        MySpeciesNode *z,  ///< a second node
        vector<MySpeciesNode*> &otherParents, ///< secondary parents
        vector<MySpeciesNode*> &toCheck ) ///< those found by computeMin 
{
    vector<bool> inToCheck( mSpeciesNetwork->getNumberOfNodes(), false );
    BOOST_FOREACH( MySpeciesNode *node, toCheck ) 
        inToCheck[node->getId()] = true;

    vector<MySpeciesNode*> yAncestors( 
                            mSpeciesNetwork->getNumberOfNodes() );
    vector<MySpeciesNode*> zAncestors( 
                            mSpeciesNetwork->getNumberOfNodes() );

    // 1. Find all shared ancestors
    markParents( y, otherParents, yAncestors );
    markParents( z, otherParents, zAncestors );
    
    // 2. Find those with separate paths (MuN),
    //    e.g. x in MuN if one child of x leads to y and the other to z,
    //    i.e. one child is ancestor of y and the other of z.
    vector<MySpeciesNode*> MuN;
    BOOST_FOREACH( MySpeciesNode *node, yAncestors ) {
        if( node != NULL && node->getNumberOfSons() == 2 ) {
            int id0 = node->getSon(0)->getId();
            int id1 = node->getSon(1)->getId();
            if( (yAncestors[id0] != NULL && zAncestors[id1] != NULL)
                || (yAncestors[id1] != NULL && zAncestors[id0] != NULL) )
            {
                MuN.push_back( node );
            }
        }
    }

    // 3. Find x in MuN that do not share nodes on the path.
    //    For each path to y, check if there is a path to z
    //    that doesn't share nodes.
    size_t checked = 0;
    BOOST_FOREACH( MySpeciesNode *mu, MuN ) {
        vector<int> pathCountsY( mSpeciesNetwork->getNumberOfNodes(), 0);
        int numPathsY = allPaths( mu, y, yAncestors, pathCountsY );

        vector<int> pathCountsZ( mSpeciesNetwork->getNumberOfNodes(), 0);
        int numPathsZ = allPaths( mu, z, zAncestors, pathCountsZ );

        // check if there are any shared nodes in all paths
        bool isMin = true;
        int allPathsId = -1;
        BOOST_FOREACH( MySpeciesNode *m, MuN ) { // check just MUs
            int id = m->getId();
            if( mu->getId() != id
                && pathCountsY[id] == numPathsY 
                && pathCountsZ[id] == numPathsZ ) 
            {
                allPathsId = id;
                isMin = false;
            }
        }
        if( inToCheck[mu->getId()] ) {
            if( !isMin ) {
                cout << "MUs= ";
                BOOST_FOREACH( MySpeciesNode *m, MuN ) {
                    if( inToCheck[mu->getId()] )
                        cout << " " << m->getId() << "*";
                    else
                        cout << " " << m->getId();
                }
                cout << endl;

                cout << "mu=" << mu->getId() << " is not a Min" << endl;
                cout << allPathsId << " in all paths" << endl;
                throw bpp::Exception( "checkComputeMin: Min not right" );
            }
            checked++;
        } else if ( isMin ) {
            throw bpp::Exception( "checkComputeMin: missing Min" );
        }

    }
    if( checked != toCheck.size() ) 
        throw "checkComputeMin: not all Min in MuN";
}


/**
 * Find the distances to every ancestor of node.
 */
void ancestorDistances(
        MySpeciesNode *node,   ///< a node
        vector<MySpeciesNode*> &otherParents, ///< secondary parents
        vector<int> &dists, ///< distances to each ancestor by id
        int dist=0 )    ///< current distance in the recursion
{
    int id = node->getId();
    if( dists[id] != -1 && dist >= dists[id] ) 
        return; // already set and less than current

    dists[id] = dist;

    if( node->hasFather() ) {
        int nextDist = dist;
        if( node->getFather()->getNumberOfSons() == 2 )
            nextDist++;
        ancestorDistances( node->getFather(), otherParents, dists, nextDist );
    }
    if( otherParents[id] != NULL ) {
        int nextDist = dist;
        if( otherParents[id]->getNumberOfSons() == 2 )
            nextDist++;
        ancestorDistances( otherParents[id], otherParents, dists, nextDist );
    }
}


/**
 * Recursion to find mins for y and z by propogating labels, given
 * the distances to the ancestors of y and z.
 *
 * The label is the id of the deepest common descendant of y and z
 * or -1 if there is no such node.
 *
 * @return label for node
 */
int NetAlg::computeMinAux(
        MySpeciesNode *node,   ///< current node in the recursion
        vector<int> &yDists, ///< given distances to the ancestors of y
        vector<int> &zDists, ///< given distances to the ancestors of z
        vector<int> &labels, ///< labels
        vector<MySpeciesNode*> &min ) ///< list of mins to return
{
    int id = node->getId();
    if( yDists[id] <= 0 || zDists[id] <= 0 ) 
        // not an ancestor
        return -1;

    if( labels[id] > -1 )
        return labels[id]; // seen it

    int labelSon0 = computeMinAux( node->getSon(0), yDists, zDists, 
                                   labels, min );

    // one son, propagate label for hybrid nodes
    if( node->getNumberOfSons() == 1 ) {
        labels[id] = labelSon0;
if( gDebug )
cout << " propagate1 " << labels[id] << " to " << id << endl;
        return labelSon0;
    }

    int labelSon1 = computeMinAux( node->getSon(1), yDists, zDists, 
                                   labels, min );

    // check if in MuN (each son has a path to y and z separately)
    int id0 = node->getSon(0)->getId();
    int id1 = node->getSon(1)->getId();
    if( (yDists[id0] != -1 && zDists[id1] != -1)
        || (yDists[id1] != -1 && zDists[id0] != -1) )
    {
        if( labelSon0 != -1 && labelSon0 == labelSon1 ) {
            labels[id] = labelSon0; // MuN, but not a min
if( gDebug )
cout << " propagate2 " << labels[id] << " to " << id << endl;
        } else {
            labels[id] = id;
            min.push_back( node ); // Found a min.
if( gDebug )
cout << " MIN " << id << endl;
        }
    } else {
        // one child is not an ancestor, propagate the label
        if( yDists[id0] == -1 && zDists[id0] == -1 )
            labels[id] = labelSon1;
        else if( yDists[id1] == -1 && zDists[id1] == -1 )
            labels[id] = labelSon0;
        else
            throw bpp::Exception( "computeMin: both labels -1" );
if( gDebug )
cout << " propagate3 " << labels[id] << " to " << id << endl;
    }

    return labels[id];
}


/**
 * Find the minimal set of shared ancestors (with separate paths to y and z).
 *
 * @return minimal set
 */
vector<MySpeciesNode*> NetAlg::computeMin( 
        MySpeciesNode *y,          ///< a node 
        MySpeciesNode *z,          ///< a second node
        vector<MySpeciesNode*> &otherParents, ///< secondary parents
        vector<int> &yDists,    ///< ancestor distances to y
        vector<int> &zDists )   ///< ancestor distances to z
{
if( gDebug ) 
cout << "compute Min " << y->getId() << " and " << z->getId() << endl;

    // get distances to ancestors
    yDists.resize( mSpeciesNetwork->getNumberOfNodes(), -1 );
    zDists.resize( mSpeciesNetwork->getNumberOfNodes(), -1 );
    ancestorDistances( y, otherParents, yDists );
    ancestorDistances( z, otherParents, zDists );
   
    // search for mins recursively
    vector<MySpeciesNode*> min;
    vector<int> labels( yDists.size(), -1 );
    computeMinAux( mSpeciesNetwork->getRootNode(), 
                   yDists, zDists, labels, min );

    return min;
}


/**
 * Add newPair to cList if itis not there.
 */
void NetAlg::merge(
        vector<NodeVal> &cList, ///< list of node/values
        NodeVal newPair )   ///< a node/value
{
    for( size_t i=0; i<cList.size(); i++ ) {
        if( cList[i].node == newPair.node ) {
            // found it
            if( newPair.val < cList[i].val ) {
                cList[i].val = newPair.val;
                cList[i].dupe = newPair.dupe;
                cList[i].losses = newPair.losses;
                cList[i].c1nodeId = newPair.c1nodeId;
                cList[i].c2nodeId = newPair.c2nodeId;
            }
            return;
        }
    }

    // not found, add it
    cList.push_back( newPair );
}


/**
 * Calculates cost of combining two node/cost children.
 */
void NetAlg::processInternalNode(
        NodeVal c1, ///< node/value pair
        NodeVal c2, ///< second node/value pair
        vector<MySpeciesNode*> &otherParents, ///< secondary parents
        vector<NodeVal> &cList ) ///<  cost list
{
    int c1nodeId = c1.node->getId();
    if( c1.node->isLeaf() )
        c1nodeId = 0;
    int c2nodeId = c2.node->getId();
    if( c2.node->isLeaf() )
        c2nodeId = 0;

    vector<int> dists1;
    vector<int> dists2;
    vector<MySpeciesNode*> minList = computeMin( c1.node, 
            c2.node, otherParents, dists1, dists2 );
if( gDebug )
checkComputeMin( c1.node, c2.node, otherParents, minList );

    BOOST_FOREACH( MySpeciesNode *minNode, minList ) {
        NodeVal c;
        int minSon1id = minNode->getSon(0)->getId();
        int minSon2id = minNode->getSon(1)->getId();
if( gDebug ) 
cout << " " << minNode->getId() << " is min of " <<
c1.node->getId() << " and " << c2.node->getId() 
<< " with sons " << minSon1id << " and "
<< minSon2id << endl;
        double sumDist1 = -1;
        double sumDist2 = -1;
        if( dists1[minSon1id] != -1 && dists2[minSon2id] != -1 ) 
            sumDist1 = dists1[minSon1id] + dists2[minSon2id];
if( gDebug ) 
cout << " sumDist1 = " << sumDist1 << " = " << dists1[minSon1id] << " + " << dists2[minSon2id] << endl;
        
        if( dists1[minSon2id] != -1 && dists2[minSon1id] != -1 )
            sumDist2 = dists1[minSon2id] + dists2[minSon1id];
if( gDebug ) 
cout << " sumDist2 = " << sumDist2 << " = " << dists1[minSon2id] << " + " << dists2[minSon1id] << endl;

        if( sumDist1 == -1 && sumDist2 == -1 )
            throw bpp::Exception( "runMinReconAux: both dists -1" );

        if( sumDist2 == -1 || (sumDist1 != -1 && sumDist1 < sumDist2) ) {
            c.val = c1.val + c2.val + mLossCost*(sumDist1);
            c.losses = sumDist1;
        } else {
            c.val = c1.val + c2.val + mLossCost*(sumDist2);
            c.losses = sumDist2;
        }
        c.c1nodeId = c1nodeId;
        c.c2nodeId = c2nodeId;
        c.dupe = false;
        c.node = minNode;
if( gDebug ) 
cout << "  -> " << minNode->getId() << "/" << c.val << endl;
        merge( cList, c );
    }


    if( dists2[c1.node->getId()] != -1 ) {
        // c2.node is desc of c1.node
        NodeVal c;
        c.node = c1.node;
        c.c1nodeId = c1nodeId;
        c.c2nodeId = c2nodeId;
        c.dupe = true;
        c.losses = dists2[c1.node->getId()];
        c.val = c1.val + c2.val + mDupliCost 
              + mLossCost*dists2[c1.node->getId()];
if( gDebug ) 
cout << " duplication at " << c1.node->getId() << " to " 
<< c2.node->getId() << " dists=" << dists2[c1.node->getId()] 
<< "  -> " << c1.node->getId() << "/" << c.val << endl;
        if( dists1[c2.node->getId()] != -1 ) 
            throw bpp::Exception ("runMinReconAux: mutual descendants" );
        merge( cList, c );
    } else if( dists1[c2.node->getId()] != -1 ) {
        // c1.node is desc of c2.node
        NodeVal c;
        c.node = c2.node;
        c.c1nodeId = c1nodeId;
        c.c2nodeId = c2nodeId;
        c.dupe = true;
        c.losses = dists1[c2.node->getId()];
        c.val = c1.val + c2.val + mDupliCost 
              + mLossCost*dists1[c2.node->getId()];
if( gDebug ) 
cout << " duplication at " << c2.node->getId() << " to " 
<< c1.node->getId() << " dists=" << dists1[c2.node->getId()] 
<< "  -> " << c2.node->getId() << "/" << c.val << endl;
        merge( cList, c );
    }
}

/**
 * DFS traversal of network to calculate minimum reconciliation.
 */
void NetAlg::runMinReconAux( 
    MyGeneNode* geneNode, ///< current node
    vector<MySpeciesNode*> &otherParents, ///< secondary parents
    vector<MySpeciesNode*> &leafMapping ) 
            ///< mapping of gene leave to specie leaves
{
    vector<NodeVal> cList;
    if( geneNode->isLeaf() ) {
        // gene leaf
        NodeVal c;
        c.dupe = false;
        c.losses = 0;
        c.val = 0;
        c.node = leafMapping[geneNode->getId()];
        cList.push_back( c );
    } else {
                    
        if( geneNode->getNumberOfSons() != 2 )
            throw bpp::Exception( "runMinReconAux: non-binary gene node" );

        // DFS
        MyGeneNode *son1 = geneNode->getSon(0);
        MyGeneNode *son2 = geneNode->getSon(1);
        runMinReconAux( son1, otherParents, leafMapping );
        runMinReconAux( son2, otherParents, leafMapping );

        if( gDebug ) 
            cout << "processing gene " << geneNode->getId()
                 << ": " << son1->getId() << "," << son2->getId() << endl;
        // combine cost lists of children
        vector<NodeVal> cList1 = mCLists[son1->getId()];
        vector<NodeVal> cList2 = mCLists[son2->getId()];
        BOOST_FOREACH( NodeVal c1, cList1 ) {
            int c1nodeId = c1.node->getId();
            if( c1.node->isLeaf() )
                c1nodeId = 0;
            BOOST_FOREACH( NodeVal c2, cList2 ) {
                int c2nodeId = c2.node->getId();
                if( c2.node->isLeaf() )
                    c2nodeId = 0;
                if( c1.node == c2.node ) {
                    NodeVal c;
                    c.dupe = true;
                    c.losses = 0;
                    c.c1nodeId = c1nodeId;
                    c.c2nodeId = c2nodeId;
                    c.node = c1.node;
                    c.val = c1.val + c2.val + mDupliCost;
                    merge( cList, c );
if( gDebug ) 
cout << " duplication at " << c2.node->getId() << " to " 
<< c1.node->getId() << " dists=0  -> " << c2.node->getId() 
<< "/" << c.val << endl;
                } else 
                    processInternalNode( c1, c2, otherParents, cList );
            }
        }
    }

    if( gDebug ) {
        cout << "cList " << geneNode->getId() << ":";
        BOOST_FOREACH( NodeVal c, cList ) {
            cout << " " << c.node->getId() << "/" << c.val;
            if( !geneNode->isLeaf() )
                 cout << " spSons="  << c.c1nodeId << "," << c.c2nodeId
                 << " d=" << c.dupe << " l=" << c.losses;
            cout << endl;
        }
    }

    mCLists[geneNode->getId()] = cList;
}


/**
 * Find the minimum network reconciliation between mSpeciesNetwork
 * and mGeneTree
 *
 * @return cost of the minimal reconciliation
 */
double NetAlg::runMinRecon()
{
    // get non-father parents
    vector<MySpeciesNode*> otherParents( mSpeciesNetwork->getNumberOfNodes() );
    vector<MySpeciesNode*> nodes = mSpeciesNetwork->getNodes();
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
        for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
            MySpeciesNode *son = node->getSon(i);
            if( son->getFather()->getId() != node->getId() )
                otherParents[son->getId()] = node;
        }
    }

    vector<MySpeciesNode*> leafMapping( mGeneTree->getNumberOfNodes() );
    vector<MySpeciesNode*> speciesLeaves = mSpeciesNetwork->getSortedLeaves();
    size_t speciesIdx = 0;
    BOOST_FOREACH( MyGeneNode* geneLeaf, mGeneTree->getSortedLeaves() ) {
        // gene leaf
        string leafName = geneLeaf->getName();
        size_t pos = leafName.find( mCharSep );
        string speciesName = leafName.substr( 0, pos );

        while( speciesIdx < speciesLeaves.size() 
               && speciesName != speciesLeaves[speciesIdx]->getName() )
        {
            speciesIdx++;
        }
        if( speciesIdx == speciesLeaves.size() )
            throw bpp::Exception( 
                    "runMinRecon: no species found for a gene leaf" );
        leafMapping[geneLeaf->getId()] = speciesLeaves[speciesIdx];
    }

    mCLists.resize( mGeneTree->getNumberOfNodes() );
    runMinReconAux( mGeneTree->getRootNode(), otherParents, leafMapping );
// What is the recon?
// What about losses and dupli?


    // find the root

    double minCost = -1;
    NodeVal minNode;
    BOOST_FOREACH( NodeVal c, mCLists[mGeneTree->getRootNode()->getId()] ) {
        if( minCost == -1 || c.val < minCost ) {
           minCost = c.val; 
           minNode = c;
        }
    }

    return minCost;
}

void NetAlg::printRecon(  
    MyGeneNode* geneNode, ///< current node
    int cListIdx )
{
    if( geneNode == NULL ) {
        cout << "-----Species network-----------" << endl;
        mSpeciesNetwork->printNetwork();
        cout << endl;
        geneNode = mGeneTree->getRootNode();
        double minCost = -1;
        cListIdx = 0;
        int idx = 0;
        BOOST_FOREACH( NodeVal c, mCLists[geneNode->getId()] ) {
            if( minCost == -1 || c.val < minCost ) {
               minCost = c.val; 
               cListIdx = idx;
            }
            idx++;
        }
    } else if( geneNode->isLeaf() ) {
        MySpeciesNode *spNode = mCLists[geneNode->getId()][0].node;
        cout << geneNode->getId() << "(" << geneNode->getName() 
             << "): " << spNode->getName() << endl;
        return;
    } 

    NodeVal minNode = mCLists[geneNode->getId()][cListIdx];
    string minExtId = bpp::TextTools::toString(
                            minNode.node->getInfos().realPostOrder );
    if( minNode.node->isLeaf() )
        minExtId = minNode.node->getName();
    cout << geneNode->getId() << ": " << minExtId //minNode.node->getId()
         << " cost=" << minNode.val;
    if( minNode.dupe )
        cout << " duplication,";
    cout << " " << minNode.losses << " losses: "
         << geneNode->getSon(0)->getId() << ","
         << geneNode->getSon(1)->getId() << endl;

    vector<NodeVal> cList1 = mCLists[geneNode->getSon(0)->getId()];
    vector<NodeVal> cList2 = mCLists[geneNode->getSon(1)->getId()];
    for( size_t idx1=0; idx1<cList1.size(); idx1++ ) {
        NodeVal c1 = cList1[idx1];
        int c1nodeId = c1.node->getId();
        if( c1.node->isLeaf() )
            c1nodeId = 0;
        if( minNode.c1nodeId != c1nodeId ) 
            continue; // species must match
        for( size_t idx2=0; idx2<cList2.size(); idx2++ ) {
            NodeVal c2 = cList2[idx2];
            int c2nodeId = c2.node->getId();
            if( c2.node->isLeaf() )
                c2nodeId = 0;
            if( minNode.c2nodeId != c2nodeId )
                continue; // species must match
            double cost = c1.val + c2.val + mLossCost*minNode.losses;
            if( minNode.dupe )
                cost += mDupliCost;
            if( (minNode.val - cost) <= COST_DIFF ) {
                printRecon( geneNode->getSon(0), idx1 );
                printRecon( geneNode->getSon(1), idx2 );
                return;
            }
        }
    }
    throw Exception ("NetAlg::printRecon: failed to backtrack" );
}

