#ifndef MYGENETREE_H_
#define MYGENETREE_H_
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

struct GeneNodeInfos  {
};
typedef bpp::NodeTemplate<GeneNodeInfos> MyGeneNode;


// Implementation of BPP tree
class MyGeneTree: public MyTreeTemplate<MyGeneNode> {
private:

    void postorder( int &pOrd, MyGeneNode *node );

public:

    MyGeneTree(){} //Empty constructor for DeCo

    /**
     * Constructor that takes root node.
     */
 	MyGeneTree(
            MyGeneNode &root ) ///< root node
        : MyTreeTemplate<MyGeneNode>(& root)
    {}  

    /**
     * Constructor that reads a tree from a Newick string.
     */
    MyGeneTree( 
        string description, ///< Newick string 
        string &errString, ///< error description
        bool readBootstrap=false ) ///< read bootstrap values
        : MyTreeTemplate<MyGeneNode>()
    {
        readTree( description, errString, readBootstrap );
    }

    void assignPostOrderIds();
    int collapseOnTree( float threshold, int mode); 
    bool rootTree( MyGeneNode *rootSon=NULL );
    bool reroot( vector<string> &outGroup, int proportion );
    bool rootWithThreshold( float threshold, int mode );

    bool restrictTreeToASetOfTaxa( boost::unordered_map<string,int> &taxaNames , boost::unordered_map<string, string> &mapNames,
                                   char charSep='x', bool verbose=false );

    static vector<MyGeneTree*> readMyGeneTrees( const char *treePathChar,
                                string &errString, bool readBoostrap=false );

};	

#endif /*MYGENETREE_H_*/
