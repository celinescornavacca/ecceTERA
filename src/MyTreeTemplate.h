#ifndef MYTREETEMPLATE_H_
#define MYTREETEMPLATE_H_

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

#include <vector>

#include <Bpp/Phyl/NodeTemplate.h>
#include <Bpp/Phyl/TreeTemplate.h>

#include <boost/unordered_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

using namespace std;

typedef boost::tokenizer<boost::char_separator<char> > Tokenizer; 

/**
 * Adds basic functions to bpp tree template.
 */
template<class N>
class MyTreeTemplate: public bpp::TreeTemplate<N> {

protected:
    vector<N*> mCorrespondence;   ///< mapping from id to node

public:
    /**
     * Empty constructor.
     */
 	MyTreeTemplate() : bpp::TreeTemplate<N>() {}



    /**
     * Constructor that takes root node.
     */
 	MyTreeTemplate(
            const MyTreeTemplate<N>& root ) ///< root node
        : bpp::TreeTemplate<N>(root)
    {}  



    /**
     * Constructor that takes root node pointer.
     */
 	MyTreeTemplate(
            N* root ) ///< root node
        : bpp::TreeTemplate<N>(root)
    {}  

protected:
    

    /**
     * Reads a newick string into this tree.
     */
    void readTree( 
            string description, ///< Newick string 
            string &errString, ///< error description
            bool bootstrap ) ///< names are bootstrap values
    {
        // tokenize tree input
        boost::char_separator<char> sep("", "(),:;");
        Tokenizer tok( description, sep );
        Tokenizer::iterator iter=tok.begin();

        // create tree
        N *root = readMyNode( tok, iter, errString, bootstrap );
        if( errString != "" ) 
            return;

        this->setRootNode( root );

        // check for correct format
        if( iter==tok.end() || (*iter).compare(";") ) {
            errString = "Bad format: no semi-colon found.";
            return;
        }
    
        this->resetNodesId(); // assign ids
    }



    /**
     * Read a string tokenizer for a node name and properties.
     *
     * Assign the node name if the node is a leaf, else assign
     * the bootstrap value.
     *
     * @return eNewick tag if found (so always empty because NwetAlgo is deprecated in this version)
     */
    static string readNodeName( 
            Tokenizer &tok, ///< string tokenizer
            Tokenizer::iterator &iter, ///< string tokenizer iterator
            N *node, ///< node
            bool isLeaf,  ///< true if node is a leaf
            string &errString, ///< error description
            bool bootstrap ) ///< names are bootstrap values
    {
        string tag = "";

        // check for node name and branch
        bool afterColon = false;
        while( iter != tok.end() ) {
            if( !iter->compare(")") || !iter->compare(",") 
                || !iter->compare(";") ) 
            {
                break; // processed in next level up
            }

            if( !iter->compare("(") ) {
               errString = "Bad format: unexpected '('.";
               return tag;
            }

            if( !iter->compare(":") )
                afterColon = true;
            else if( afterColon ) 
                // branch length, (toDouble throws nan exception)
                node->setDistanceToFather( boost::lexical_cast<double>(*iter) );
            else {
                // node name
                string val = 
                        bpp::TextTools::removeSurroundingWhiteSpaces(*iter);
                if( isLeaf || !bootstrap ) {
                    node->setName(val);

                } else {
                    // internal node bootstrap
                    try {
                        double bs = bpp::TextTools::toDouble( val );
                        if( bs < 0 ) {
                            errString = "Negative bootstrap value.";
                            return tag;
                        }
                        if( bs == 0 ) {
                            errString = "Internal node has a bootstrap of 0.";
                            return tag;
                        }
                        node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, 
                                                bpp::Number<double>(bs));
                    } catch( bpp::Exception e ) {
                        //cout << "bootstrap = " << val << endl;
                        errString = "Non-numeric bootstrap value.";
                        return tag;
                    }
                }
             }
             iter++;
        }
        return tag;
    }
   




    /**
     * Create a node from a tokenized string.
     *
     * The function is called recursively to read all of the nodes
     * and build a tree.
     *
     * @return The new node.
     */
    N* readMyNode( 
            Tokenizer &tok, ///< string tokenizer
            Tokenizer::iterator &iter, ///< string tokenizer iterator
            string &errString, ///< error description
            bool bootstrap ) ///< names are bootstrap values
    {
        // initialize basic information
        N *node = new N();

        // States:
        // ( -> first son, descend
        // ) -> return after getting node info
        // , -> next son, descend
        // else node info <name><:#>

        // process any children
        bool inNode = false;
        bool isLeaf = false;
        while( iter!=tok.end() ) {
            if( !(*iter).compare("(") || !(*iter).compare(",") ) {
                // start a child node, don't advance iter
                if( !(*iter).compare("(") ) 
                    inNode = true;
                iter++;
                N *son = readMyNode( tok, iter, errString, 
                        bootstrap );
                if( errString != "" ) 
                    return NULL;
                node->addSon( son );
            } else if( !(*iter).compare(")") ) {
                // finished node
                iter++;
                inNode = false;
                break;
            } else if( !(*iter).compare(";") ) {
                errString = "Bad format: unexpected ';' found.";
                return NULL;
            } else {
                // leaf
                isLeaf = true;
                break;
            }
        }

        if( inNode ) {
            errString = "Bad format: missing ')'.";
            return NULL;
        }

        string tag = readNodeName ( tok, iter, node, isLeaf, errString,
                                bootstrap );
                                     
                                 
        if( errString != "" ) 
            return NULL;

        return node;
    }


    /**
     * Delete given node and all descendants.
     *
     * @return Number of nodes delete.
     */
    int deleteSubtree( 
            N *node ) ///< subtree root to delete
    {
        int deleteCount = 1;

        // get all sons before deleting them
        vector<N*> sons;
        int sonCount = node->getNumberOfSons();
        for( int i=0; i<sonCount; i++ ) 
            sons.push_back( node->getSon(i) );

        // delete son's trees
        BOOST_FOREACH( N *son, sons )
            deleteCount += deleteSubtree( son );

        delete node;
        return deleteCount;
    }



    /** 
     * Remove a node and merge the sibling with a cousin.
     *
     * Repeats with father if the father has only a single son. The 
     * distance to father and bootstrap values are added if transferProperties
     * is true.
     *
     * @return False if tree is empty.
     */
    bool collapseEdge( 
            N* node, ///< node to collapse
            bool transferProperties = false, 
                     ///< propagate distance and bootstrap
            bool keepBootstrap = false  ) ///< don't change the bootstrap value
    {
        bool hasLeaves = true;

        if( transferProperties ) {
            if( node->hasFather() && node->getFather()->hasFather() ) {
                N *father = node->getFather();
                if( !keepBootstrap 
                    && node->hasBranchProperty(bpp::TreeTools::BOOTSTRAP))
                {
                    double bs = dynamic_cast<const bpp::Number<double> *> 
                            (node->getBranchProperty(bpp::TreeTools::BOOTSTRAP))
                            ->getValue();
                    if( father->hasBranchProperty( bpp::TreeTools::BOOTSTRAP ) ) 
                        bs += dynamic_cast<const bpp::Number<double> *> 
                          (father->getBranchProperty(bpp::TreeTools::BOOTSTRAP))
                              ->getValue();
                    father->setBranchProperty( bpp::TreeTools::BOOTSTRAP, 
                                               bpp::Number<double>(bs) );
                }
                if( node->hasDistanceToFather() ) { 
                    unsigned sonCount = node->getNumberOfSons();
                    double dist = node->getDistanceToFather();
                    for (unsigned i=0; i<sonCount; i++) {
                        if( node->getSon(i)->hasDistanceToFather() ) 
                            node->getSon(i)->setDistanceToFather( node->getSon(i)->getDistanceToFather() +dist);
                    }
                }
            } // else drop values at root
        }

        if( !node->hasFather() ) {
            //NOTE: isLeaf() will return true for a root with one son
            if( node->getNumberOfSons() > 0 ) {
                N *newRoot = node->getSon(0);
                newRoot->removeFather();
                this->setRootNode( newRoot );
            } else {
                // deleting last node
                hasLeaves = false;
            }
            delete node;
        } else {
            N *father= node->getFather();
            unsigned sonCount = node->getNumberOfSons();	
            for (unsigned i=0; i<sonCount; i++) {
                // we take always the first son...it's always a different one
                N *son = node->getSon(0);
                node->removeSon( son );
                father->addSon( son );
            } 
            father->removeSon( node );
            delete node;

            // collapse again if necessary
            if( father->getNumberOfSons()==1 )
                hasLeaves = collapseEdge( father );
        }

        return hasLeaves;
    }


    /**
     * Read newick strings from a file.
     */
    static vector<string> readNewickStrings(
        const char *treePathChar ) ///< path to file with newick trees
    {
        ifstream treePath( treePathChar, ios::in );  
        if( !treePath ) { throw bpp::IOException ( "MyTreeTemplate"
                "::readNewickStrings: failed to read from stream"); }

        string description;
        vector<string> treeStrings;
        while( !treePath.eof() ) {
            // Copy current line in temporary string
            string temp;
            getline(treePath, temp, '\n');  

            // skip empty lines and empty trees
            if(temp != "" && temp != "()" ) {
                string::size_type index = temp.find(";");
                if(index != string::npos) {
                    description += temp.substr(0, index + 1);
                    treeStrings.push_back( description );
                    description = "";
                }
                else description += temp;
            }
        }

        return treeStrings;
    }



    struct leaf_sort {
        /** leaf sorting operator */
        inline bool operator() (const N* a, const N* b ) {
            return b->getName() > a->getName();
        }
    };

public:
    /**
     * Get tree leaves, sorted by name.
     *
     * @return sorted leaves
     */
    vector<N*> getSortedLeaves() 
    {
        N *root = this->getRootNode();

        vector<N*> leaves;
        if( root->getNumberOfSons() == 0 )
            return leaves;

        if( root->getNumberOfSons() < 2 ) {
            // 2 nodes - root is technically a leaf, but don't return it
            leaves = bpp::TreeTemplateTools::getLeaves( *(root->getSon(0)) );
        } else {
            leaves = this->getLeaves();
        }
        sort( leaves.begin(), leaves.end(), leaf_sort() );
        return leaves;
    }


    /** 
     * Check if the tree leaves are unique.
     *
     * @return True if the leaves are unique.
     */
    bool uniqueLeaves( string &dupName ) 
    {
        boost::unordered_map<string,int> leafMap;
        vector<N*> leaves = this->getLeaves();		
        N *leaf;
        BOOST_FOREACH( leaf, leaves ) {
            boost::unordered_map<string,int>::iterator iter 
                = leafMap.find( leaf->getName() );
            if ( iter == leafMap.end() ) 
                leafMap.insert(std::pair<string,int>(leaf->getName(),0));
            else {
                dupName = leaf->getName();
                return false;
            }
        }

        return true;
    }


    /** 
     * Check if the tree is binary.
     *
     * @return true if the tree is binary
     */
    bool isBinary() 
    {
        int count = getBinaryCount();
        return count == 0;
    }


    /** 
     * Check if the tree is binary.
     *
     * @return Number of non-binary nodes.
     */
    int getBinaryCount() 
    {
        int nonBinCount = 0;
        vector<N*> nodes = this->getNodes();		
        for( size_t i=0; i<nodes.size(); i++ ) {
            if( nodes[i]->getNumberOfSons() != 0 && nodes[i]->getNumberOfSons() != 2 )  
                nonBinCount++;
        }

        return nonBinCount;
    }


    /**
     * Return the node with the given postorder id, as calculated with
     * assignPostOrderIds. If the node is artficial, return the the
     * real child descendant.
     *
     * @return the node with the given post order
     */
    N* getNodeById( 
            int id ) ///< id of a node.
    {
        if( mCorrespondence.size() == 0 )
            throw bpp::Exception("MyTreeTemplate::getNodeById: not set");

        if( (size_t) id >= mCorrespondence.size() || id < 0 ) {
            cout << "demand for " << id 
                << ", array size " << mCorrespondence.size() << endl;
            throw bpp::Exception("MyTreeTemplate::getNodeById:"
                            " id out of range");
        }

        return mCorrespondence[id];
    }	
 


    /**
     * useful debugging function
     */
    void printMe( 
            N *node=NULL, ///< current node
            int depth=0 ) ///< node depth
    {
        if( node==NULL )
            node = this->getRootNode();

        for( int i=0; i<depth; i++ )
            cout << "  ";
        int sonCount = node->getNumberOfSons();
        cout << node->getId() << " (" << sonCount << ")";
        if( !node->hasFather() )
            cout << "*";
        if( node->hasName() )
            cout << " : " << node->getName();
        else if( node->isLeaf() ) 
            cout << " : L";
        cout << endl;
        for( int i=0; i<sonCount; i++ )
            printMe( node->getSon(i), depth+1 );
    }


    /**
     * Make the tree binary.
     */
    void makeBinary( 
            N *node = NULL ) ///< a tree node
    {
        bool atRoot = false;
        if( node == NULL ) {
            atRoot = true;
            node = this->getRootNode();
        }

        if( node->isLeaf() )
            return;

        int sonCount = node->getNumberOfSons();
        if( sonCount > 3 ) {
            // split son list in two 
            N *firstNewNode = new N();
            N *secondNewNode = new N();
            N *newNode = firstNewNode;
            for( int i=0; i<sonCount; i++ ) {
                if( i == int(sonCount/2) ) 
                    newNode = secondNewNode;
                N *son = node->getSon(0);
                node->removeSon( son  );
                newNode->addSon( son );
            }
            node->addSon( firstNewNode );
            node->addSon( secondNewNode );
        } else if( sonCount == 3 ) {
            // add one node
            N *newNode = new N();
            for( int i=0; i<2; i++ ) {
                N *son = node->getSon(1); // get the second son twice
                node->removeSon( son  );
                newNode->addSon( son );
            }
            node->addSon( newNode );
        } else if( sonCount == 1 ) {
        	cout << "here collapsing\n";
            // remove this node
            N *son = node->getSon(0);
            collapseEdge( node, true );
            makeBinary( son );
        }

        if( sonCount > 1 ) {
            makeBinary( node->getSon(0) );
            makeBinary( node->getSon(1) );
        }

        if( atRoot ) {
            if( !isBinary() )
                throw bpp::Exception("makeBinary NOT WORKING" );
            this->resetNodesId();
        }
    }


    /**
     * Print the tree root at node to the file.
     */
    void printNewick( string fileName, bool append=false ) 
    {
        N *node = this->getRootNode();
        string newickStr = bpp::TreeTemplateTools::nodeToParenthesis( *node );
        ios_base::openmode mode = ios::out | ios::binary;
        if( append )
            mode |= ios::app; // append
        ofstream out( fileName.c_str(), mode );
        out << newickStr << ";" << endl;
        out.close();
       
    }

};	
#endif /*MYTREETEMPLATE_H_*/
