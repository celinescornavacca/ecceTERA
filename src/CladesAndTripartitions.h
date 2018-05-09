#ifndef CLADESANDTRIPARTITIONS_H_
#define CLADESANDTRIPARTITIONS_H_
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




//#define DEBUG

// Define BITS to use vector< bitset<BLOCK_SIZE> > to represent clades
// rather than vector<bool>
#define BITS //celineS


#include <bitset>

#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <MyGeneTree.h>


typedef std::vector<int> VectorInt; 
    ///< typedef for split type, which can be used as type in maps
	



/**
 * Clade collection
 */
class Clades {

    friend class CladesAndTripartitions;

    private:
    static const int BLOCK_SIZE = 32;
        ///< Size of block for bit representation of clades
    bool mNodesOnly; ///< clades are just tree nodes
    int mLabelCount; ///< number of leaves; same for all gene trees
    int mCladeCount; ///< Current number of clades.

#ifdef BITS
    typedef vector< bitset<BLOCK_SIZE> > CladeType; ///< clade representation
    typedef vector<int> CladeHashType;  ///< type of a hashed clade

    // Methods and members for bitset representation of clades.
    int mBitBlocks;    ///< number of blocks (=labelCount/BLOCK_SIZE+1)
    bitset<BLOCK_SIZE> mMask;  ///< mask of used bits for highest block 

    CladeHashType hashClade( CladeType cladeV, bool flip );
#else
    typedef vector<bool> CladeType;
    typedef vector<bool> CladeHashType;
#endif

    // Clade information
    boost::unordered_map<CladeHashType,int> mCladeToInt;
        ///< mapping of clade array to number
    vector<string> mSpeciesNameMapping; 
        ///< maps clade number to first part gene name (species name)
    vector<string> mLeafNameMapping; 
        ///< maps clade number to leaf name
	boost::unordered_map<std::string,int> mLeafClades;
        ///< maps gene names to leaf clade numbers


    public:

    // Clade info used during calculation of clades
    struct Clade {
        int bits;       ///< number of leaves in clade 
        int cladeInt;   ///< clades integer representaion
        int antiCladeInt;  ///< integer representation of anti-clade
        CladeType cladeV;  ///< array representation of clade
    };

    private:
    void printClade( CladeType cladeV ); // for debugging
    int insertCladeInt( CladeType cladeV, bool flip );
    int getCladeInt( CladeType cladeV, bool flip );
    bool fillCladeInts( Clade &clade, bool useAntiClades );
    int createSimpleCladesAux( MyGeneNode *node, 
            vector< vector< pair<int,int> > > &cladeSplits,
            vector<double> &branchLengths );

    // creation functions that don't work with nodesOnly
    void precomputeLeaves( bool useAntiClades );
    void precomputeLeafSplits( bool useAntiClades, int treeCount, 
  	        boost::unordered_map<VectorInt,double> 
                    &splitsOccurrences );
    Clade createClade( string leafName, bool useAntiClades );
    bool mergeClades( vector<Clade> &clades, Clade &clade, bool useAntiClades );
    vector<int> deleteClades( vector<int> &cladesToDelete );

    // nodesOnly clade creator
    vector< vector< pair<int,int> > > createSimpleClades( MyGeneTree &geneTree,
            vector<double> &branchLengths );

    vector<string> getCladesAsStrings();
    vector< vector<int> > getNoncompatible();

    Clades(): mNodesOnly(false), mLabelCount(0), mCladeCount(0) {}
    void init( char charSep, MyGeneTree &geneTree, bool nodesOnly = false );
    void init( char charSep, vector<string> &leafNames, bool nodesOnly );

    public:

    void mapSpeciesNames( vector<string> &speciesNames,
                          vector<string> &geneNames, string &errStr );


    /******* Get clade information *************/
   
    /** 
     * Return true if clade identified by cladeNum is a leaf.  
     *
     * @return True if the clade is a leaf.
     */
    bool isLeaf( int cladeNum ) ///< a clade identifier
    { 
        if( cladeNum == 0 || cladeNum > mLabelCount ) 
            return false;
        else 
            return true;
    }
    
    /** 
     * Returns the species name associated with a clade number.
     *
     * This is the first part of the gene (leaf) name separated by charSep.
     *
     * @return gene name
     */
    string getSpeciesName( int cladeNum ) ///< a clade identifier
    { 
        if( cladeNum == 0 || cladeNum > mLabelCount ) {
            cout << cladeNum << endl;
            throw bpp::Exception("Clades.getSpeciesame given a non-leaf clade");
        }
        return mSpeciesNameMapping[cladeNum]; 
    }

    /**
     * Returns the leaf name associated with a clade number.
     *
     * @return leaf name
     */
    string getLeafName( int cladeNum ) ///< a clade identifier
    {
        if( cladeNum == 0 || cladeNum > mLabelCount ) {
            cout << cladeNum << endl;
            throw bpp::Exception("Clades.getLeafName given a non-leaf clade");
        }
        return mLeafNameMapping[cladeNum]; 
    }

    /**
     * Returns the number of clades.
     *
     * @return number of clades
     */
    int getCladeCount() { return mCladeCount; }

    /**
     * Returns the id of the root clade.
     *
     * @return the id of the root clade 
     */
    int getRootClade() { return 0; } // all clade number 

    vector<int> getSortedClades();

    /**
     * Returns the number of leaf clades.
     *
     * @return number of leaf clades
     */
    int getLeafCount() { return mLabelCount; };


    /**
    Takes:
        - vector<string> LeafNames: list of leaf names
    Returns:
        int: id of the corresponding clade
    **/
    int getCladeIntFromLeafList(vector<string> LeafNames);

}; // end Clades class






/** 
 * Class to create clades and tripartitions (splits).
 * The clades are stored in the Clades class.
 * Tripartitions exists as a list of tripartitions for each clade
 * as well as the associated number of occurences and splits ratio.
 */
class CladesAndTripartitions {

    protected: //switched to protected in order to inherit
    bool mRooted;   ///< is the gene tree rooted
    int mGeneTreeCount;
    string mNewickGeneTree; ///< save newick string for ALE output

    boost::unordered_map<VectorInt,double> mSplitsOccurrencesMap;
        ///< Occurrences per split divided by occurrences of clade.
        ///< Used to calculate mSplitsRatio and mCladeOccurrences, then cleared
    vector< vector< pair<int,int> > > mCladeSplits; 
        ///< All splits for each clade (by clade id).
    vector< vector<double> > mSplitsRatio; 
        ///< Occurence of each split as a vector, divided by clade occurrences.
    vector<double> mCladeOccurrences;
        ///< Occurence of each clade, indexed by clade number.
    vector<double> mBranchLengths; ///< branch length for each clade
    vector< pair<int,int> > mPolytomySizes; 
        ///< size of polytomy a clade is in <cladeNum,size>
    vector< pair<int,int> > mPolytomyCladeSizes; 
        ///< size of polytomy for a clade <cladeNum,size>

    // The two occurrences arrays are both ordered
    // by clade. Furthermore, the secondary vectors of splitsOccurrences
    // are ordered by cladeSplits.

    void countSplit( VectorInt split, double occurrences );
    void calculateSplits( double occurrences, 
            const MyGeneNode *node, Clades::Clade clade, 
            Clades::Clade firstChild, Clades::Clade secondChild,
            bool notRootSon=false );
    static unsigned long long doubleFactorial( int n, bool &overflow );
    void polytomy( const MyGeneNode *node, const int polytomySize, 
            vector<Clades::Clade> &childClades,
            Clades::Clade allClade, boost::unordered_map<int,bool> &seenIt,
            bool &overflow );
    Clades::Clade computeSplits( const MyGeneNode *node, 
                    bool polytomic, bool &overflow, int level=0 );
    int getCladesBySizeAux( int curClade, vector<int> &cladeSizes ); 
    void postOrderMappingAux( vector<int> &mapping, int &pOrd, int id_u );
    void getRootMpAmAux( int id_u, vector<double> &vals );
    MyGeneNode *getTreeAux( int &pOrd, int idU );

    bool toInt( string in, int &out );
    bool toDouble( string in, double &out );
    string loadALE( string aleFileName, char charSep, int &splitCount );
    vector<int> redoTripartitions( vector<int> &deletedClades,
                                   vector<int> &cladeIdMapping );

    // bootstrap split ratio
    int checkBootstrapValues( MyGeneTree* geneTreeForBootstrap, 
                              string &errStr );
    void bootstrapSplitRatio( MyGeneTree* geneTreeForBootstrap,
                        int modifier, string &errStr, bool printCladeInfo,
                        bool useBootstrapValues=true );
    Clades::Clade mapCladesToNodesAux( MyGeneNode* node,
            boost::unordered_map<string,int> &nameToInt,
            vector<MyGeneNode*> &cladeToNode, 
            bool useBootstrapValues, string &errStr );
    vector<MyGeneNode*> 
        mapCladesToNodes( MyGeneTree* geneTree, bool useBootstrapValues,
                string &errStr );

    public:
    Clades mClades;  ///< Clades class containing clades information.


    CladesAndTripartitions( char charSep, vector<MyGeneTree*> &geneTrees, 
                            bool verbose, bool &overflow, string &errStr,
                            bool polytomic = false, 
                            MyGeneTree *bootstrapGeneTree = NULL,
                            bool polytomicUnrooted = false,
                            bool printCladeInfo = false );


    // rooted version
    CladesAndTripartitions( char charSep, MyGeneTree &geneTree );

    // Ale reader
    CladesAndTripartitions( char charSep, string aleFileName, string &errStr,
                                bool verbose );

    // debugging
    void printMe();



    /**
     * Get number splits for the given clade.
     *
     * @return number of for the clade.
     */
    int getSplitCount( int cladeNum ) ///< a clade identifier
    {
        return mCladeSplits[cladeNum].size();
    }

    /** 
     * Get the specified clade split.
     *
     * @return clade split
     */
    pair<int,int> getCladeSplit(
            int cladeNum, ///< a clade identifier
            int i ) ///< split index
    {
        return mCladeSplits[cladeNum][i];
    }

    /** 
     * Return occurrences of split divided by total for clade.
     *
     * @return split ratio.
     */
    double getSplitRatio(
            int cladeNum, ///< a clade identifier
            int i ) ///< split index
    {
        return mSplitsRatio[cladeNum][i];
    }

    /** 
     * Return occurrences of clade
     *
     * @return clade occurrences
     */
    unsigned long long
    getCladeOccurrences( int cladeNum ) ///< a clade identifier
    {
        return mCladeOccurrences[cladeNum];
    }
   

    double getRootMpAm();
    vector<int> getPostOrderMapping();
    vector< vector<int> > getCladesBySize();
    vector< vector<int> > getCladeParents();
    MyGeneTree *getTree();
    int excludeClades( float percentToExclude );
    static unsigned long long findPolytomies( 
                const MyGeneNode *node, bool &overflow );


    void saveALE( string aleFileName );
};	



#endif 
