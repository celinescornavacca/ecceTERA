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
Implementation of MyTreeTemplate for species networks.
*/


#include "MyNetwork.h"

#include <iostream>
#include <queue>
#include <boost/foreach.hpp>


/** 
 * Constructor from a string in Newick format.
 */
MyNetwork::MyNetwork( 
        string description, ///< Newick string 
        string &errString ) ///< error description
    : MySpeciesTree(), mEnewick(true)
{
    readTree( description, errString, true, false );
}


/**
 * Deep copy construtor.
 */
MyNetwork::MyNetwork( MyNetwork &origNetwork ) 
{
    vector<MySpeciesNode*> origNodes = origNetwork.getNodes();

 	mCorrespondence.resize( origNodes.size() );

    // create new nodes
    map<MySpeciesNode*,MySpeciesNode*> nodeMap;
    BOOST_FOREACH( MySpeciesNode *origNode, origNodes ) {
        MySpeciesNode *newNode = new MySpeciesNode();
        NodeInfos info( origNode->getInfos() );
        newNode->setInfos( info );

        int nodeId = origNode->getId();
        newNode->setId( nodeId );
        mCorrespondence[nodeId] = newNode;
        
        if( origNode->hasName() ) 
            newNode->setName( origNode->getName() );

        // Must copy all properties too:
        vector<std::string> names = origNetwork.getNodePropertyNames(nodeId);
        for (size_t i = 0; i < names.size(); i++) {
            newNode->setNodeProperty( names[i], 
                    *origNetwork.getNodeProperty( nodeId, names[i]) );
        }
        names = origNetwork.getBranchPropertyNames(nodeId);
        for (size_t i = 0; i < names.size(); i++) {
            newNode->setBranchProperty(names[i], 
                        *origNetwork.getBranchProperty(nodeId, names[i]));
        }

        nodeMap[origNode] = newNode;
    }

    // connect new nodes
    map<MySpeciesNode*,MySpeciesNode*>::iterator iter;
    BOOST_FOREACH( MySpeciesNode *origNode, origNodes ) {
        iter = nodeMap.find( origNode );
        if( iter == nodeMap.end() ) 
            throw bpp::Exception( "MyNetwork::MyNetwork:"
                        " could not find original node in map" );
        MySpeciesNode *newNode = iter->second;

        if( origNode->hasFather() ) {
            iter = nodeMap.find( origNode->getFather() );
            if( iter == nodeMap.end() ) 
                throw bpp::Exception( "MyNetwork::MyNetwork:"
                            " could not find father node in map" );
            newNode->setFather( iter->second );
        } else {
            this->setRootNode( newNode );
        }

        if( origNode->getInfos().primaryFather != NULL ) {
            iter = nodeMap.find( origNode->getInfos().primaryFather );
            if( iter == nodeMap.end() ) 
                throw bpp::Exception( "MyNetwork::MyNetwork:"
                            " could not find primary father node in map" );
            newNode->getInfos().primaryFather = iter->second;
        }
        else
        	newNode->getInfos().primaryFather = NULL;

        if( origNode->getInfos().secondaryFather != NULL ) {
            iter = nodeMap.find( origNode->getInfos().secondaryFather );
            if( iter == nodeMap.end() ) 
                throw bpp::Exception( "MyNetwork::MyNetwork:"
                            " could not find secondary father node in map" );
            newNode->getInfos().secondaryFather = iter->second;
        }
        else
        	newNode->getInfos().secondaryFather = NULL;
    }


}



void MyNetwork::setCorrespondance(vector <MySpeciesNode *> allNodes){
	    mCorrespondence.resize( allNodes.size() );
	    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
	    		mCorrespondence[node->getId()] = node;	
	    }
} 


    /**
     * Overrides MyTreeTemplate function to copy parents, for networks 
     */

void MyNetwork::copyPrimSecParentsInfo( 
        	MySpeciesNode *node, 
        	MySpeciesNode *nodeToCopy )
{
         	node->getInfos().primaryFather 
                    = nodeToCopy->getInfos().primaryFather;
         	node->getInfos().secondaryFather 
                    = nodeToCopy->getInfos().secondaryFather;
         } 


  /**
     * Make the tree binary.
     Note that we cannot have root node with outdegree 1 in networks/switchings because the newick format (extended and not) does not permit it
     */
    void MyNetwork::makeBinarySwitching( 
            MySpeciesNode *node, bool onlyDeleteUnlabeledLeaves ) ///< a tree node
    { 
      bool atRoot = false;
        if( node == NULL ) {
            atRoot = true;
            node = this->getRootNode();
        }

        int sonCount = node->getNumberOfSons();
            	
        //cout << "node->id " << node->getId() << endl;
        
        if( sonCount==0 ){
                //cout << "leaf" << endl;
            if( node->hasName() && node->getName().find('#')==std::string::npos)    
            	return;
			else if(!( node->hasName()) || node->getName().find('#')!=std::string::npos ){
            	MySpeciesNode *father = node->getFather();
				father->removeSon(node);
				delete(node);
				//cout << "delete unlabeled\n";	
				makeBinarySwitching( father,onlyDeleteUnlabeledLeaves );          		
            	return;
            }		
         }   


        if( sonCount == 1 ) {
            // remove this node
            if(!atRoot){
                MySpeciesNode *son = node->getSon(0);
				if(!onlyDeleteUnlabeledLeaves){
            		//cout << "here collapsing "  << node->getId() << " " << node->getSon(0)->getId() << "\n";
            		MySpeciesNode *father = node->getFather();
            		father->addSon(son);
            		father->removeSon(node);
            		delete(node);
				}
            	makeBinarySwitching( son,onlyDeleteUnlabeledLeaves );
            }
        }

        if( sonCount > 1 ) {
            //cout << "here recursing "  << node->getId() << " " << node->getSon(0)->getId() << " " << node->getSon(1)->getId() << "\n";
			MySpeciesNode * son1 = node->getSon(0);
			MySpeciesNode *son2 =  node->getSon(1);
            makeBinarySwitching( son1,onlyDeleteUnlabeledLeaves );
            makeBinarySwitching( son2 ,onlyDeleteUnlabeledLeaves);
        }

		if( atRoot) {// to avoid to have outdegree-1 root
		    //cout << "root cleaning\n";
		     if(!onlyDeleteUnlabeledLeaves){
				//cout << "root cleaning\n";
				while(sonCount == 1){
		        	node = this->getRootNode();
		    		MySpeciesNode *son = node->getSon(0);
				this->setRootNode(son);
				node->removeSon(son);
				son->removeFather();
            		delete(node);	
            		sonCount = son->getNumberOfSons();

            	}
			}
            makeBinarySwitching( this->getRootNode(),onlyDeleteUnlabeledLeaves);

		}		
		
    	if( atRoot ) {
            if( !isBinary() )
                throw bpp::Exception("makeBinarySwitching NOT WORKING" );
            this->resetNodesId();
        }
    }

/**
 * Overrides MyTreeTemplate function to process to store the information that the arc that we are going to add  is a secondary one. 
 *  
*/

    void  MyNetwork::handleSecondaryArcs( 
            MySpeciesNode *node    ///< node 

            ) 
    { 
    		if(node->getNumberOfSons()==1 && node->getInfos().primaryFather==NULL){
    			//cout << node->getFather() << endl;
				node->getInfos().primaryFather = node->getFather();
				//cout << node->getInfos().primaryFather << endl;

			}	
    }


/**
 * Overrides MyTreeTemplate function to process eNewick tags
 */
MySpeciesNode *MyNetwork::handleNewickTag( 
        bool isLeaf,
        MySpeciesNode *node, 
        string tag, 
        string &errString ) 
{
    map<string,MySpeciesNode*>::iterator tagIter = tags.find( tag );
    //cout << "tag " << tag << endl;
	tags.insert( make_pair(tag, node) );
    return node;
}


/**
 * Delete all but root.
 */
void MyNetwork::deleteNetwork()
{
    vector<MySpeciesNode*> allNodes = getNodes();
    getRootNode()->removeSons();
    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
        if( node->hasFather() ) {
            delete node;
        }
    }
}


/**
 * Read newick tree from a file.
 *
 * @return a tree
 */
MyNetwork* MyNetwork::readMyNetwork( 
    const char *treePathChar, ///< path to file with newick trees
    string &errString ) ///< error description
{
    vector<string> treeStrings = readNewickStrings( treePathChar );
    if( treeStrings.size() == 0 )
        return NULL;

    return new MyNetwork( treeStrings[0], errString );
}



/** 
 * Get sorted leaves for a network.
 * @return sorted leaves
 */
vector<MySpeciesNode*> MyNetwork::getSortedLeaves() 
{
    vector<MySpeciesNode*> leaves;
    if( getRootNode()->getNumberOfSons() == 0 )
        return leaves;

    vector<MySpeciesNode*> allNodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
        if( node->getNumberOfSons() == 0 ) 
            leaves.push_back( node );
    }

    sort( leaves.begin(), leaves.end(), leaf_sort() );

    return leaves;
}

/**
 * Assign ids to a network using a breadth-first postorder.
 *
 * Also, set the mCorrespondence vector mapping ids to nodes.
 */
void MyNetwork::assignNetworkPostOrderIds()
{
    if( getNumberOfNodes() == 0 )
        return ;

    // get all nodes
    vector<MySpeciesNode*> allNodes;
    vector<MySpeciesNode*> curNodes;
    MySpeciesNode *root = getRootNode();
    curNodes.push_back( root );
    allNodes.push_back( root );
    root->setId( -999 );
    while( curNodes.size() > 0 ) {
        vector<MySpeciesNode*> sons;
        BOOST_FOREACH( MySpeciesNode *node, curNodes ) {
            for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
                MySpeciesNode *son = node->getSon(i);
                if( son->getId() != -999 ) {
                    son->setId( -999 );
                    allNodes.push_back( son );
                    sons.push_back( son );
                }
            }
        }
        curNodes = sons;
    }

    mCorrespondence.resize( allNodes.size() );

    int id = 0;
    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
        mCorrespondence[id] = node;
        node->setId( id++ );
     }
     
/* Breadth first
    int nodeCount = allNodes.size();
    vector<bool> seenIt( nodeCount, false );

    // get non-father parents
    vector<MySpeciesNode*> otherParents( nodeCount ); 
    vector<MySpeciesNode*> nodes; // leaves
    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
        if( node->getNumberOfSons() == 0 )
            nodes.push_back( node ); 
        for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
            MySpeciesNode *son = node->getSon(i);
            if( son->getFather() != node )
                otherParents[son->getId()] = node;
        }
    }
        
    int pOrd = 0;
    // bottom-up breadth first method
    while( nodes.size() > 0 ) {
        vector<MySpeciesNode*> parentNodes; // nodes for next level up
        BOOST_FOREACH( MySpeciesNode *node, nodes ) {
            if( !seenIt[node->getId()] ) { 
                seenIt[node->getId()] = true;
                mCorrespondence[pOrd] = node;
                pOrd++;
                if( node->hasFather() ) {
                    parentNodes.push_back( node->getFather() );
                    if( otherParents[node->getId()] != NULL 
                        && otherParents[node->getId()] != node->getFather() ) 
                    {
                        parentNodes.push_back( otherParents[node->getId()] );
                    }
                }
            }
        }
        nodes = parentNodes;
    }

    for( size_t i=0; i<mCorrespondence.size(); i++ ) 
        mCorrespondence[i]->setId( i );
*/

}



/**
 * Assign ids to a network using a breadth-first postorder.
 *
 * Also, set the mCorrespondence vector mapping ids to nodes.
 */
int MyNetwork::assignNetworkPostOrderIdsDFS()
{
    if( getNumberOfNodes() == 0 )
        return 0;

    // get all nodes
    vector<MySpeciesNode*> allNodes =getNodes();
    vector<MySpeciesNode*> orderedNodes ;

    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
		node->getInfos().seenIt = false;
	}
	
    			
    vector<MySpeciesNode*> stack ;

    mCorrespondence.resize( allNodes.size() );
    mCorrespondenceTS.resize( allNodes.size() );
    	
	stack.push_back(getRootNode() );

    while(stack.size()>0){
    	MySpeciesNode* node= stack.front(); 
		stack.erase(stack.begin());
    	
    	if(!node->getInfos().seenIt){
			//cout << "id node entering " << node->getId() << " " << node->getNumberOfSons() << "\n";
			//node->getInfos().seenIt = true;
			if(node->getNumberOfSons() == 0 ){
				if(!node->getInfos().seenIt ){
					orderedNodes.push_back( node ); 
					node->getInfos().seenIt = true;
					//cout << "id node " << node->getId() << "\n";
				}

			}
			else{
				size_t seenSons=0;
				for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
			
					MySpeciesNode *son = node->getSon(i); 

					if( ! son->getInfos().seenIt  )   {
						stack.push_back( son );
					}
					else
						seenSons++;
				
				}
				
				if( node->getNumberOfSons() == seenSons){
					if(node->getSon(0)->getInfos().seenIt){//  && node->getInfos().secondaryFather->getInfos().seenIt && node->getInfos().primaryFather->getInfos().seenIt ) { 
						node->getInfos().seenIt =true;
						orderedNodes.push_back( node ); 
						//cout << "id node int " << node->getId() << "\n";
					} 
				}

				if (!node->getInfos().seenIt){
					stack.push_back( node );
				}
			}
				
               
    	}
    }	
    
    for( size_t i=0; i<orderedNodes.size(); i++ ){ 
    	//cout << orderedNodes[i]->getId( ) << " ";
    	orderedNodes[i]->setId( i );
    	//cout << orderedNodes[i]->getId( ) << endl;
       	mCorrespondence[i]=orderedNodes[i];
       	mCorrespondenceTS[i].push_back(i);
    }
	
	return allNodes.size();

}

void MyNetwork::compute_RealPostOrder() 
{
    const vector<MySpeciesNode*>& nodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
    	node->getInfos().realPostOrder = node->getId();  
    }
}

/**
 * Recursion for printNetwork.
 */
void MyNetwork::printNetworkAux( 
        MySpeciesNode *node, ///< current node in recursion
        vector<bool> &seenIt, ///< list of nodes visited
        string indent ) ///< indent for current line
{
    cout << indent;

    int id = node->getId();
    if( seenIt[id] ) 
        cout << "*";
    cout << id << " (rpo=" << node->getInfos().realPostOrder << ")";
    if( node->getNumberOfSons()==0 ) 
        cout << " leaf " << node->getName();

    cout << endl;

    if( !seenIt[id] ) {
        seenIt[id] = true;
        int sonCount = node->getNumberOfSons();
        for( int i=0; i<sonCount; i++ ) 
            printNetworkAux( node->getSon(i), seenIt, indent+"  " );
    }
}

/**
 * Print this network.
 */
void MyNetwork::printNetwork() {
    if( mCorrespondence.size() == 0 )
        throw bpp::Exception("MyNetwork::printNetwork: postOrder not set");

    int nodeCount = getNumberOfNodes(); // this count is greater or equal
    vector<bool> seenIt( nodeCount, false );

    printNetworkAux( getRootNode(), seenIt );
}


/**
 * Get all of the nodes in the network.
 *
 * The normal getNodes function doesn't work because it adds subtrees
 * with multiple fathers more than once.
 *
 * @return network nodes.
 */
vector<MySpeciesNode*> MyNetwork::getNodes()
{
    vector<MySpeciesNode*> nodes;
    vector<MySpeciesNode*> rawNodes = 
        MyTreeTemplate<MySpeciesNode>::getNodes(); // some nodes twice
    BOOST_FOREACH( MySpeciesNode *node, rawNodes )
        node->getInfos().seenIt = false;
    BOOST_FOREACH( MySpeciesNode *node, rawNodes )
        if( !node->getInfos().seenIt ) {
            node->getInfos().seenIt = true;
            nodes.push_back( node );
            //cout << node->getId() << endl;
        }

    return nodes;
}


/**
 * Check for a valid network. 
 *
 * Internal nodes have in-degree/out-degree 1/2 or 2/1
 *
 * Insert a node for hybrid nodes with two sons.
 *
 * @return true if network valid
 */
bool MyNetwork::checkNetwork( 
        string &errStr ) ///< return description of any errors
{
   if( mCorrespondence.size() == 0 )
        throw bpp::Exception("MyNetwork::checkNetwork: postOrder not set");

    vector<MySpeciesNode*> allNodes = getNodes();
    vector<MySpeciesNode*> otherParent( allNodes.size() );
    boost::unordered_map<string,int> leafMap;
    int maxId = 0;
    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {

        if( maxId < node->getId() )
            maxId = node->getId();

        int sonCount = node->getNumberOfSons();
        if( sonCount > 2 ) {
            errStr = "A network node has more than two children";
            return false;
        }
        for( int i=0; i<sonCount; i++ ) {
            MySpeciesNode *son = node->getSon( i );
            if( node->getId() != son->getFather()->getId() ) {
                if( otherParent[son->getId()] != NULL ) {
                    errStr = "A hybrid node has more than two parents";
                    return false;
                }
                otherParent[son->getId()] = node;
                //cout << "adding " << node->getId() << " as other parent of " << son->getId()<< endl; 
            }
        }
	
        // check for unique leaves
        if( node->isLeaf() ) {
            boost::unordered_map<string,int>::iterator iter 
                = leafMap.find( node->getName() );
            if ( iter == leafMap.end() ) 
                leafMap.insert(std::pair<string,int>(node->getName(),0));
            else  if(node->getName()!="") {
                errStr = "The species network does not have unique leaves.";
                return false;
            }
        }

        
    }


    // check in-degree and out-degree
    bool changed = false;
    maxId++;
    BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
    	node->getInfos().seenIt = false;
        int sonCount = node->getNumberOfSons();
        if( !node->hasFather() ) {
            if( sonCount != 2 ) {
                errStr = "The species root must have two children";
                return false;
            }
        } else if( node->getInfos().secondaryFather!=NULL && sonCount != 1 ){//otherParent[node->getId()] != NULL && sonCount != 1 ) { 
            


            // add hybrid node if there is another parent and son count
            // is not one
            MySpeciesNode *newNode = new MySpeciesNode();
            newNode->setId( maxId++ );
            if( sonCount == 0 ) { // leaf
                string name = node->getName();
                size_t pos = name.find( "#" ); 
                string leafName = name.substr(0,pos);
                string hybridName = name.substr(pos);
                newNode->setName( leafName );
                node->setName( hybridName );
                node->deleteName();
                
            } else {
                MySpeciesNode *son0 = node->getSon( 0 );
                MySpeciesNode *son1 = node->getSon( 1 );
                node->removeSons();
                newNode->addSon( son0 );
                newNode->addSon( son1 );
            }
            node->addSon( newNode );
            changed = true;
        } else if(node->getInfos().secondaryFather==NULL && sonCount == 1){ //if( otherParent[node->getId()] == NULL && sonCount == 1 ) {
            // remove nodes with one parent and one son
            MySpeciesNode *father = node->getFather();
            cout << father->getId() << " " <<  node->getId() << " " << node->getSon(0)->getId() << " " << endl;

            
            father->removeSon( node );
            //node->removeFather();
            father->addSon( node->getSon(0) );
            node->getSon(0)->getInfos().primaryFather  = node->getInfos().primaryFather;
         	node->getSon(0)->getInfos().secondaryFather  = node->getInfos().secondaryFather;
            //cout << node->getSon(0)->getNumberOfSons() << endl;
            //cout << node->getSon(0)->getInfos().secondaryFather << endl;

            delete node;

            changed = true;
        }
    }
	
	
    vector<MySpeciesNode*> allNodesNew;
    vector<MySpeciesNode*> curNodes;
    MySpeciesNode *root = getRootNode();
    root->getInfos().seenIt=true;
    curNodes.push_back( root );
    while( curNodes.size() > 0 ) {
        vector<MySpeciesNode*> sons;
        BOOST_FOREACH( MySpeciesNode *node, curNodes ) {
            for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
                MySpeciesNode *son = node->getSon(i);
                if( son->getInfos().seenIt == false ) {
                    son->getInfos().seenIt = true;
                    allNodes.push_back( son );
                    sons.push_back( son );
                }
            }
        }
        curNodes = sons;
    }


		
    if( changed ) {
        assignNetworkPostOrderIds();
    	cout << "id changed";    
    }
	
    return true;
}


/**
 * Construct a subtree for the given id, i.e.,
 * only nodes where biCompMap[nodeId] == biCompId are included.
 *
 * If there are hybrid nodes, a tree is constructed for each possible
 * switching.
 * 
 * @return switching trees for the given biCompId
 */
vector<MyNetwork*> MyNetwork::allSwitchings(
        MySpeciesNode *netRoot, ///< root of network (speciation node)
        const int biCompId, ///< biconnected set to use 
        const vector<int> &biCompMap ) ///< mapping of node ids to set
{
    // find all hybrid nodes and record which nodes (of all) are in network
    vector<int> hybridIndices( getNumberOfNodes() ); 
    int hybridCount = 0;
    vector<MySpeciesNode*> nodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
        if( biCompMap[node->getId()] == biCompId 
            && node->getNumberOfSons() == 1 ) 
        {
            hybridIndices[node->getId()] = hybridCount;
            hybridCount++;
        }
    }

    // create swtichings
    int numSwitchings = pow( (float) 2, (int) hybridCount );
    vector<MyNetwork*> switches; 
    for( int switchNum=0; switchNum<numSwitchings; switchNum++ ) {
        // maps id->new
        vector<MySpeciesNode*> newNodes( getNumberOfNodes() ); 
        vector<bool> deletedNodes( getNumberOfNodes(), false ); 
        allSwitchingsAux( netRoot, switchNum, biCompId, biCompMap,
                        hybridIndices, newNodes, deletedNodes );
        MySpeciesNode *newRoot = newNodes[netRoot->getId()];
        MyNetwork *switching = new MyNetwork( *newRoot );
        switches.push_back( switching );
    }

    return switches;
}



/**
* Assign biNum to all nodes in subtree that are not already assigned.
*/
void MyNetwork::assignBiconnectedComp(
        MySpeciesNode *node, ///< current node in the recursion
        const int biNum, ///< bicomp number
        vector<int> &map ) ///< map of node-ids to bicomp number
{
    if( !node->isLeaf() && map[node->getId()] < 0 )
        map[node->getId()] = biNum;
    else
        return;

    for( size_t i=0; i<node->getNumberOfSons(); i++ ) 
        assignBiconnectedComp( node->getSon(i), biNum, map );
}


/**
* Find all biconnected components and save results to map.
*
* @return low - least id in recursion
*/
int MyNetwork::getBiconnectedCompMapAux( 
        MySpeciesNode *node, ///< current node in the recursion
        int &counter, ///< pre-numbering counter
        vector<int> &map, ///< map of node-ids to bicomp number
        vector<int> &pre, ///< pre-number for each node
        vector<MySpeciesNode*> &otherParents, ///< second parent of hybrids
        MySpeciesNode *parent ) ///< calling node 
{
    if( node->isLeaf() ) 
        return counter;

    int id = node->getId();
    if( pre[id] == -1 ) // not visited
        pre[id] = counter++;
    else 
        return pre[id]; // seen it

    vector<MySpeciesNode*> neighbors;
    for( size_t i=0; i<node->getNumberOfSons(); i++ ) 
        neighbors.push_back( node->getSon(i) );
    if( node->hasFather() && node->getFather() != parent )
        neighbors.push_back( node->getFather() );
    if( otherParents[id] != NULL && otherParents[id] != parent ) 
        neighbors.push_back( otherParents[id] );

    int low = pre[id];
    BOOST_FOREACH( MySpeciesNode *n, neighbors ) {
        int nLow = getBiconnectedCompMapAux( n, counter, map, pre,
                                             otherParents, node );
        if( nLow < low ) 
            low = nLow;
    }
    if( low == pre[id] ) 
        assignBiconnectedComp( node, node->getId(), map );

    return low;
}


/**
* Find all biconnected components.
*
* @return map of node ids to biconnected component number
*/
vector<int> MyNetwork::getBiconnectedCompMap( 
    vector<MySpeciesNode*> &biRoots ) ///< root node of each connected component
{
    int nodeCount = getNumberOfNodes();

    // get non-father parents
    vector<MySpeciesNode*> otherParents( nodeCount );
    vector<MySpeciesNode*> nodes = getNodes();
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
        for( size_t i=0; i<node->getNumberOfSons(); i++ ) {
            MySpeciesNode *son = node->getSon(i);
            if( son->getFather()->getId() != node->getId() )
                otherParents[son->getId()] = node;
        }
    }

    // create map (set id is post order id of root of component)
    vector<int> map( nodeCount, -1 );
    vector<int> pre( nodeCount, -1 );
    int counter = 0;
    getBiconnectedCompMapAux( getRootNode(), counter, map, pre, otherParents );

    // renumber consectutively
    int setIdx = 0;
    vector<int> setIdMap( nodeCount, -1 );
    for( int i=0; i<nodeCount; i++ ) {
        int origId = map[i];
        if( origId != -1 && setIdMap[origId] == -1 ) {
            // not set
            setIdMap[origId] = setIdx++;
        }
    }

    // redo map and get root nodes 
    biRoots.resize( setIdx );
    BOOST_FOREACH( MySpeciesNode *node, nodes ) {
        int origId = map[node->getId()];
        if( origId != -1 ) {
            int newId = setIdMap[origId];
            if( origId == node->getId() ) {
                biRoots[newId] = node; // sets named after id of root
            }
            map[node->getId()] = newId;
         }
    }

    return map;
}


