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

The CladesAndTripartition class takes a gene tree or a set of
gene trees and calculates the clades and splits (tripartitions).
This calculation is done in the constructors. The remaining functions
return information associated with clades and splits.

The nodesOnly constructor taking a single gene tree creates a single clade
per node in the gene tree. The other constructor calculates all possible
clades and splits for the one or more gene trees and calculates statistics
on the frequency of the clades and splits in the collection of trees.


Clades can be defined as a set of leaves. As such, they are represented
by a set of leaf ids.

Leaf clade number are assigned first, as well as anti-clade (all other leaves). 
Thus, they can be calculated, saving a map find. 0 is used for the root.
Leaf cladeNums become, leaf index + 1. Anti-clade numbers are the 
clade number + leafCount.

lc = number of leaves (labels)
Clade numbering:
    0 - all clades (all 1's)
    1 to lc: leafs
    lc+1 to 2*lc: anti-leafs
    >2lc: use next available number, inserting anti-clade at same time
      insert one with the least bits first (equal bits don't care about order).

Clades can be stored in two ways. 
   1. vector<bool>
   2. vector< bitset<BLOCK_SIZE> > 
   Define BITS to use the second option. The second option is more 
   complicated, but more efficient, espcially with larger trees.
   
Note: clades are not explicitly saved. They exists as clade ids with
associated information (leaf names, associated splits).

**/

//#define DEBUG

#include <iostream>
#include <limits>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <Bpp/Text/TextTools.h>
#include "CladesAndTripartitions.h"


#ifdef BITS
/**
 * Create and return a hash code used for the clade (set of leaves)
 * given as a CladeType.
 *
 * @return hash code
 */
Clades::CladeHashType Clades::hashClade( 
        CladeType cladeV, ///< the clade to hash
        bool flip ) ///< flip bits before hashing
{
    CladeHashType hash;
    for( int i=0; i<mBitBlocks; i++ ) {
        int h;
        if( flip ) {
            bitset<BLOCK_SIZE> flipped( cladeV[i] );
            flipped.flip();
            if( i==mBitBlocks-1 ) 
                flipped &= mMask; // undo unused bits
            h = int( flipped.to_ulong() );
        } else {
            h = int( cladeV[i].to_ulong() );
        }
        hash.push_back( h ); 
    }

    return hash;
}
#endif

/**
 * Output clade to stdout.
 */
void Clades::printClade( 
        CladeType cladeV ) ///< clade as a set of leaves
{
#ifdef BITS
    for( int k=mBitBlocks-1; k>=0; k-- ) {
        cout << " " << cladeV[k].to_string();
    }
    cout << endl;
#else
    BOOST_FOREACH( bool t, cladeV )
        cout << t << " ";
    cout << endl;
#endif
}

/**
 * Assign clade number to the given clade represented as an array. 
 *
 * NOTE: This function does not
 * check if the clade has already been assigned a number.
 *
 * @return clade number assigned to clade 
 */
int Clades::insertCladeInt( 
        CladeType cladeV, ///< clade as a set of leaves
        bool flip ) ///< flip bits before inserting (inserts anti-clade) 
{
#ifdef BITS
    CladeHashType hash = hashClade( cladeV, flip ); 
    mCladeToInt.insert(make_pair( hash, mCladeCount++ ));
#else
    if( flip ) cladeV.flip();
    mCladeToInt.insert(make_pair( cladeV, mCladeCount++ ));
    if( flip ) cladeV.flip();
#endif

#ifdef DEBUG
    // print out assigned number and array
    cout << "insert clade " << mCladeCount-1;
    #ifdef BITS
    for( int k=mBitBlocks-1; k>=0; k-- ) 
        cout << " " << hash[k] <<" "  << cladeV[k];
    #else
    BOOST_FOREACH( bool t, cladeV )
        cout << t << " ";
    cout << endl;
    #endif
    cout << endl;
#endif

    return mCladeCount-1;
}


/**
 * Return the clade number associated with the given clade represented
 * as CladeType. If flip is true, return the clade number of the 
 * anti-clade (all leaves not in cladeV).
 * 
 * @return Clade number or -1 if it doesn't exits.
 */
int Clades::getCladeInt( 
        CladeType cladeV, ///< clade as a set of leaves
        bool flip ) ///< find anti-clade if true (i.e. flip all bits)
{
    boost::unordered_map<CladeHashType,int>::iterator iterCladesToInt;
#ifdef BITS
    CladeHashType hash = hashClade( cladeV, flip );
    iterCladesToInt = mCladeToInt.find( hash );
#else
    if( flip ) cladeV.flip();
    iterCladesToInt = mCladeToInt.find( cladeV );
    if( flip ) cladeV.flip();
#endif

    int cladeInt = -1; // -1 for not found
    if( iterCladesToInt != mCladeToInt.end() ) 
        cladeInt = iterCladesToInt->second;

    return cladeInt;
}


/**
 * Find cladeInt and anticladeInt for the given clade 
 * and put them in cladeInt and antiCladeInt of the Clade structure. 
 *
 * This function
 * assumes that the clade array has already been set.
 * Given one clade number, the anti-clade can be calculated except
 * in a special instance. See notes on clade numbering in the Clades
 * class description.
 * If the clade number doesn't exist, it is added, as well as the 
 * anti-clade. The anti-clade number is always the one with more 
 * bits.
 *
 * @return true if a new clade
 */
bool Clades::fillCladeInts( 
    Clade &clade,  ///< clade to fill
    bool useAntiClades ) ///< set anti-clade with clade if true
{
    if( clade.bits == mLabelCount ) {
        // root
        clade.cladeInt = 0;
        return false;
    }

    clade.cladeInt = getCladeInt( clade.cladeV, false );
    bool isNew = false;
    if( clade.bits == 1 ) {
        // leaf (calculate anti-clade)
        if( clade.cladeInt == -1 )
            throw bpp::Exception( 
                    "Clades::fillCladeInts: Leaf clade not found" );
        clade.antiCladeInt = clade.cladeInt+mLabelCount;
    } else if( useAntiClades && clade.bits == mLabelCount-1 ) {
        // anti-leaf (calculate clade)
        if( clade.cladeInt == -1 ) 
            throw bpp::Exception(
                    "Clades::fillCladeInts: Anti-leaf clade not found");
        clade.antiCladeInt = clade.cladeInt - mLabelCount;
    } else if( clade.cladeInt == -1 ) {
        // Doesn't exist, so add it. Insert clade and anti-clade together so
        // that second one is consecutively numbered. Insert clade or anti 
        // with least bits first.
        if( useAntiClades ) {
            if( clade.bits <= mLabelCount/2 ) {
                clade.cladeInt = insertCladeInt( clade.cladeV, false );
                clade.antiCladeInt = insertCladeInt( clade.cladeV, true );
            } else {
                clade.antiCladeInt = insertCladeInt( clade.cladeV, true );
                clade.cladeInt = insertCladeInt( clade.cladeV, false );
            }
        } else
            clade.cladeInt = insertCladeInt( clade.cladeV, false );
        isNew = true;
    } else if( useAntiClades ) {
        // found clade number, calculate anticlade number
        if( clade.bits == mLabelCount/2 ) {
            // This could be calculated, too, but rare so not worth it.
            clade.antiCladeInt = getCladeInt( clade.cladeV, true );
            if( clade.antiCladeInt == -1 ) 
                throw bpp::Exception( 
                        "Clades::fillCladeInts: Anti-clade not found");
        } else if( clade.bits < mLabelCount/2 ) {
            clade.antiCladeInt = clade.cladeInt+1;
        } else {
            clade.antiCladeInt = clade.cladeInt-1;
        }
    } 

    return isNew;
}


/**
 * Clades initialization. 
 *
 * Collects leaf names from gene tree and calls other init function.
 */
void Clades::init( 
        char charSep, ///< Seperator used in gene names.
        MyGeneTree &geneTree, ///< A gene tree 
        bool nodesOnly ) 
            ///< Mimic gene tree by have only nodes instead of clades.
{
	vector<MyGeneNode *> allLeafNodes = geneTree.getLeaves();	
    vector<string> leafNames;
    BOOST_FOREACH( MyGeneNode *leaf, allLeafNodes ) {
        // nodes need names
        if( !leaf->hasName() ) 
            throw bpp::Exception("Clades::init: Leaf node missing name");
        leafNames.push_back( leaf->getName() );
    }
    init( charSep, leafNames, nodesOnly );
}

/**
 * Clades initialization. 
 *
 * Create clade numbers for root, leaves, and anti-leaves.
 *
 * The gene tree is used to create a mapping between leaf
 * names and gene names, and vice versa. The mapping from leaf name
 * to clade number uses the full leaf name while the mapping from
 * number to name only uses the first part of the name, before charSep.
 */
void Clades::init( 
        char charSep, ///< Seperator used in gene names.
        vector<string> &leafNames, ///< list of leaf leaf names
        bool nodesOnly ) 
            ///< Mimic gene tree by have only nodes instead of clades.
{
    mCladeCount = 0;
    mNodesOnly = nodesOnly;
    mLabelCount = leafNames.size();

    if( !nodesOnly ) {
#ifdef BITS
        // basic info required by bits implementation
        mBitBlocks = int(mLabelCount/BLOCK_SIZE);
        if( mLabelCount%BLOCK_SIZE != 0 ) 
            mBitBlocks++;

        // create mask to hide unused portion of highest block
        int hideBits = mLabelCount%BLOCK_SIZE;
        if( mLabelCount%BLOCK_SIZE == 0 ) // highest block is full
            hideBits = BLOCK_SIZE;
        for( int i=0; i<hideBits; i++ ) 
            mMask[i] = 1;
        #ifdef DEBUG
        cout << mBitBlocks << " blocks, mask=" << mMask.to_string() << endl;
        #endif
#endif
    }

    // maps clade number to first part of gene name
    mSpeciesNameMapping.push_back( "__root" ); // clade 0 is root
    mLeafNameMapping.push_back( "__root" ); // clade 0 is root
    for( int i=0; i<mLabelCount; i++ ) {
        string leafName = leafNames[i];

        // check for duplicates
        boost::unordered_map<string,int>::const_iterator iter = 
                                            mLeafClades.find( leafName ); 
        if( iter != mLeafClades.end() ) {
            cout << "duplicate leaf name: " << leafName << endl;
            throw bpp::Exception("Clades::init : duplicate leaf name");
        }

        // leaf clades are numbered as leaf number + 1
	    mLeafClades.insert( make_pair( leafName, i+1 ) );

        mLeafNameMapping.push_back( leafName );
        size_t pos = leafName.find( charSep );
        string speciesName = leafName.substr( 0, pos );
        mSpeciesNameMapping.push_back( speciesName );
    }
}

/**
 * Create a correspondance between species and gene names.
 *
 * The two given arrays must be parallel, e.g. speciesNames[i]
 * is the species name for gene geneNames[i].
 */
void Clades::mapSpeciesNames(
    vector<string> &speciesNames,
    vector<string> &geneNames,
    string &errStr )
{
    if( speciesNames.size() != geneNames.size() )
        throw bpp::Exception( "Clades::mapSpeciesNames:"
                " input arrays are not the same size" );

    // clear any current values
    for( int i=1; i<mLabelCount; i++ ) 
        mSpeciesNameMapping[i] = "";

    int addedCount = 0;
    for( size_t i=0; i<geneNames.size(); i++ ) {
        string leafName = geneNames[i];

        // check for duplicates
        boost::unordered_map<string,int>::const_iterator iter = 
                                            mLeafClades.find( leafName ); 
        if( iter == mLeafClades.end() ) {
            // found it
            int cladeId = iter->second;
            if( mSpeciesNameMapping[cladeId] == "" ) {
                mSpeciesNameMapping[cladeId] = speciesNames[i];
                addedCount++;
            } else if( mSpeciesNameMapping[cladeId] != speciesNames[i] ) {
                errStr = "Inconsistent names";
                return;
            } 
        } 
    }

    if( addedCount != mLabelCount )
        errStr = "Missing gene name mapppings";
}


/**
 * Creates clades for leaves and anti-leaves.
 */
void Clades::precomputeLeaves( 
        bool useAntiClades ) ///< precompute anticlades of leaves
{
    // this function does not apply with nodes only
    if( mNodesOnly ) 
       throw bpp::Exception ("Clades::precomputeLeaves: called with nodesOnly");
    
    if( mLabelCount == 0) 
       throw bpp::Exception ("Clades::precomputeLeaves: called with no leaves");

#ifdef BITS
    CladeType cladeV( mBitBlocks ); // 0;
#else
    CladeType cladeV( mLabelCount, false ); // all false 
#endif

    // assign leaves and root as clade 
    // insert root 
    insertCladeInt( cladeV, true );

    // insert leaves (first loop), and anti leaves (second loop)
    bool flip = false;
    int loops = 1;
    if( useAntiClades )
        loops = 2;
    for( int j=0; j<loops; j++ ) {
        for( int idx=0; idx<mLabelCount; idx++ ) {
#ifdef BITS
            int blockIdx = int(idx/BLOCK_SIZE);
            cladeV[blockIdx][idx%BLOCK_SIZE] = 1;
            insertCladeInt( cladeV, flip );
            cladeV[blockIdx][idx%BLOCK_SIZE] = 0;
#else
            cladeV[idx] = true;
            insertCladeInt( cladeV, flip );
            cladeV[idx] = false;
#endif
        }
        flip = true; // flip bits to insert anti-leaves
    }
}

/**
 * Insert leaf split occurrences for leaves.
 *
 * The number of gene trees is used to precompute some splitsOccurrences.
 */
void Clades::precomputeLeafSplits( 
        bool useAntiClades, ///< precompute anti-clade splits
        int treeCount, ///< number of gene trees
  	    boost::unordered_map<VectorInt,double> &splitsOccurrencesMap ) 
                ///< map of number of times each split occurs
{
    // this function does not apply with nodes only
    if( mNodesOnly ) 
        throw bpp::Exception ("Clades::precomputeSplits called with nodesOnly");
    
    if( mLabelCount == 0) 
        throw bpp::Exception ("Clades::precomputeSplits called with no leaves");

    // Precompute split occurences for leaves. 
    VectorInt split(3);
    // add leaf occurences
    for( int i=0; i<mLabelCount; i++ ) {
        split[0] = i+1;
        split[1] = split[2] = -1;
        splitsOccurrencesMap.insert(
                std::pair<VectorInt,double>(split,treeCount));
    }

    if( useAntiClades ) {
        // add all split on leaf/anti-leaf occurences
        for( int i=0; i<mLabelCount; i++ ) {
            split[0] = 0;
            split[1] = i+1;
            split[2] = i+1+mLabelCount;
            splitsOccurrencesMap.insert(
                    std::pair<VectorInt,double>(split,treeCount));
        }
    }
}



/**
 * Initialize a clade structure for the given leaf.
 *
 * @return new calde
 */
Clades::Clade Clades::createClade( 
    string leafName, ///< leaf name of the clade to create
    bool useAntiClades ) ///< set anti-clade with clade if true
{
    if( mNodesOnly ) 
        throw bpp::Exception( "Clades::createClade: called for nodesOnly");

    boost::unordered_map<string,int>::const_iterator iter = 
                                mLeafClades.find( leafName ); 
    if( iter == mLeafClades.end() ) 
        throw bpp::Exception(
                "Clades::createClade: clade number not found for leaf");
    int leafCladeNum = iter->second;

    // Leaves are assigned clade numbers 1-numLeaves
    // The leaf number is the position of the leaf in the bit
    // array, hence base 0 and 1 minus the clade number.
    int leafNum = leafCladeNum -1;
    if( leafNum >= mLabelCount ) 
        throw bpp::Exception("Clades::createClade: leafNum >= leafCladeNum");

    // set bit for leaf in clade array
    Clade clade;
#ifdef BITS
    clade.cladeV.resize( mBitBlocks );
    int blockIdx = int(leafNum/BLOCK_SIZE);
    clade.cladeV[blockIdx][leafNum%BLOCK_SIZE] = 1;
#else
    clade.cladeV.resize( mLabelCount );
    clade.cladeV[leafNum] = true;
#endif

    // get/create cladeInts
    clade.bits = 1;
    clade.cladeInt = leafNum + 1;
    if( useAntiClades )
        clade.antiCladeInt = leafNum + mLabelCount + 1;

    return clade;
}



/**
 * Combines clades to create a new clade.
 *
 * The new clade is the union of all the given leaves.
 *
 * @return true if the clade is new
*/
bool Clades::mergeClades( 
    vector<Clades::Clade> &clades, ///< set of clades to merge
    Clade &clade, ///< merged clade
    bool useAntiClades ) ///< set anti-clade with clade if true
{
    if( mNodesOnly ) 
        throw bpp::Exception( "Clades::mergeClades: called for nodesOnly");

    if( clades.size() == 0 )
        throw bpp::Exception( "Clades::mergeClades: given no clades");

    if( clades.size() == 1 ) {
        clade = clades[0];
        return false;
    }

    Clade newClade;
    newClade.bits = 0;
    // This assumes the clades are disjoint, which they should be.
    BOOST_FOREACH( Clades::Clade cl, clades ) 
        newClade.bits += cl.bits;

#ifdef BITS
    newClade.cladeV.resize( mBitBlocks );
    for( int i=0; i<mBitBlocks; i++ )  {
        BOOST_FOREACH( Clades::Clade cl, clades )
            newClade.cladeV[i] |= cl.cladeV[i];
    }
#else
    newClade.cladeV.resize( mLabelCount );
    for( size_t i=0; i<newClade.cladeV.size(); i++ ) {
        BOOST_FOREACH( Clades::Clade cl, clades )
            if( cl.cladeV[i] )
                newClade.cladeV[i] = true;
    }
#endif

    // get clade and anti-clade numbers
    bool isNew = fillCladeInts( newClade, useAntiClades ); 

    clade = newClade;
    return isNew;
}


#ifdef BITS
// special compartor for sorting bit clade arrays
struct less_than_key {
    inline bool operator() (const vector<int>& a, 
                            const vector<int>& b ) {
        int aSize = a.size();
        for( int i=aSize-1; i>= 0; i-- ) {
            if( a[i] < 0 && b[i] >= 0 ) 
                return false; 
            if( a[i] >= 0 && b[i] < 0 ) 
                return true;
            if( a[i] != b[i] || i== 0) 
                return a[i] < b[i];
        }
        // Never get here
        throw bpp::Exception( "CladesAndTripartitions::less_than_key: "
                " should never be here" );
        return false;
    }
};
#endif

/**
 * Return sorted list of clade ids.
 *
 * Sort order is numerical (based on conversion of clade array to number).
 * The requirement is to see each clade before it's merged clades (parent).
 *
 * @return sorted clades
 */
vector<int> Clades::getSortedClades() 
{
    vector<int> v;
    if( mNodesOnly ) {
        // 0 = nodes, 1 to mLabelCount = leaves,
        // rest assigned in post order 
        for( int i=1; i<mCladeCount; i++ ) 
            v.push_back( i );
        v.push_back( 0 );
    } else {
        // Create a vector of clade arrays or hashes
        pair<CladeHashType,int> me; // what a map is made of
        vector< CladeHashType > C; // vector of clade arrays or hashes
        BOOST_FOREACH(me, mCladeToInt) 
            C.push_back(me.first);

        // sort it
        #ifdef BITS
            // bits requires a special comparator
            sort(C.begin(),C.end(), less_than_key());
        #else
            sort(C.begin(),C.end());
        #endif

        // retrieve clades in order, get the clade number and add to v
        CladeHashType cladeHash;
        BOOST_FOREACH(cladeHash, C) {
            boost::unordered_map<CladeHashType,int>::iterator iterCladesToInt;
            iterCladesToInt = mCladeToInt.find( cladeHash );
            if( iterCladesToInt == mCladeToInt.end() ) 
                throw bpp::Exception( "Clade not found  in Clades.sorted" );
            v.push_back( iterCladesToInt->second );
        }
    }

    return v;
}

/**
 * Return a list of strings representing each clade,
 * indexed by clade number.
 *
 * @return clade strings
 */
vector<string> Clades::getCladesAsStrings()
{
    vector<string> cladeStrings( mCladeCount );
    if( mNodesOnly ) 
        return cladeStrings;

    pair<CladeHashType,int> me; // what a map is made of
    BOOST_FOREACH(me, mCladeToInt) {
        CladeHashType c = me.first;
        string bits;
#ifdef BITS
        int block = 0;
        int bit = 1;
        for( int i=0; i<mLabelCount; i++ ) {
            // vector<int>
            if( c[block] & bit ) 
                bits += "1";
            else
                bits += "0";
            bit <<= 1;
            if( bit == 0 ) {
                bit = 1;
                block++;
            }
        }
#else
        for( int i=0; i<mLabelCount; i++ ) {
            // vector<bool>
            if( c[i] )
                bits += "1";
            else
                bits += "0";
        }
#endif
        reverse( bits.begin(), bits.end() );

        int cladeNum = me.second;
        cladeStrings[cladeNum] = bits;
    }

    return cladeStrings;
}


/**
 * Return the list of non-compatible clades for each clade.
 *
 * Non-compatible clades are those that interesect but are not
 * strict subsets.
 *
 * @return non-compatiblity lists
 */
vector< vector<int> > Clades::getNoncompatible()
{
    if( mNodesOnly ) // doesn't apply (nothing non-compatible)
        throw bpp::Exception( "Clades:getNoncompabile called with nodes only" );

    // get a list of clades and a list of clade numbers
    vector<CladeHashType> allClades;
    vector<int> allCladeNums;
    pair<CladeHashType,int> me; // what a map is made of
    BOOST_FOREACH(me, mCladeToInt) {
        allClades.push_back( me.first );
        allCladeNums.push_back( me.second );
    }

    vector< vector<int> > noncompatibleLists( allClades.size() );

    for( size_t i=0; i<allClades.size(); i++ ) {
        int cladeNum = allCladeNums[i];
        CladeHashType cladeI = allClades[i];
        for( size_t j=i+1; j<allClades.size(); j++ ) {
            CladeHashType cladeJ = allClades[j];
            // check for intersection and subset
            bool hasCommon = false;
            bool iNotSubset = false;
            bool jNotSubset = false;
#ifdef BITS
            for( int b=0; b<mBitBlocks; b++ ) {
                int intersection = cladeI[b] & cladeJ[b];
                if( intersection != 0 )
                    hasCommon = true;
                if( cladeI[b] != intersection )
                    iNotSubset = true;
                if( cladeJ[b] != intersection )
                    jNotSubset = true;
                if( hasCommon && iNotSubset && jNotSubset )
                    break;
                    
            }
#else
            for( size_t b=0; b<cladeI.size(); b++ ) {
                if( cladeI[b] && cladeJ[b] )
                    hasCommon = true;
                else if( cladeI[b] && !cladeJ[b] )
                    iNotSubset = true;
                else if( !cladeI[b] && cladeJ[b] )
                    jNotSubset = true;
                if( hasCommon && iNotSubset && jNotSubset )
                    break;
            }
#endif
            if( hasCommon && iNotSubset && jNotSubset ) {
                noncompatibleLists[cladeNum].push_back( allCladeNums[j] );
                noncompatibleLists[allCladeNums[j]].push_back( cladeNum );
            }
        }
    }

    return noncompatibleLists;
}

/**
 * Recursive function for createSimpleClades
 *
 * Automatically assign root to clade 0 and 
 * number leaves consecutively.
 *
 * @return clade number assigned to node
 */
int Clades::createSimpleCladesAux( 
    MyGeneNode *node, ///< gene tree node
    vector< vector< pair<int,int> > > &cladeSplits,
        ///< list of child pairs for each node
    vector<double> &branchLengths ) ///< branch length of each node
{
    int nodeCladeNum;
    if( node->isLeaf() ) {
        string leafName = node->getName();
        boost::unordered_map<string,int>::const_iterator iter = 
                                mLeafClades.find( leafName ); 
        if( iter == mLeafClades.end() ) {
            cerr << leafName << " not found" << endl;
            throw bpp::Exception( "createimpleCladesAux: "
                    "clade number not found for leaf");
        }
        nodeCladeNum = iter->second;
        cladeSplits[nodeCladeNum].push_back( make_pair(-1, -1) );

        if( node->hasDistanceToFather() ) 
            branchLengths[nodeCladeNum] = node->getDistanceToFather();
    } else {

        if( node->getNumberOfSons() != 2 ) 
            throw bpp::Exception("Non-binary tree node in simpleClades");
    
        int childClade0 = createSimpleCladesAux( node->getSon(0), cladeSplits,
                                                 branchLengths );
        int childClade1 = createSimpleCladesAux( node->getSon(1), cladeSplits,
                                                 branchLengths );


        if( node->hasFather() ) {
            nodeCladeNum = mCladeCount++; // post-order

            // branch length
            double dist = -1;
            if( node->hasDistanceToFather() ) 
                dist = node->getDistanceToFather();
            branchLengths.push_back( dist );
            if( branchLengths.size() != (size_t) mCladeCount ) 
                throw bpp::Exception ("NOPE: branchLengths" );
        } else {
            nodeCladeNum = 0; // root
        }

        cladeSplits[nodeCladeNum].push_back( 
                make_pair(childClade0, childClade1) );
    }

    if( nodeCladeNum < 0 || (size_t) nodeCladeNum >= branchLengths.size() ) {
        cout << "BAD " << nodeCladeNum << " vs size=" 
             << branchLengths.size() << endl;
        throw bpp::Exception( "FAILURE in branchLengths" );
    }
    
    return nodeCladeNum;
}

/**
 * Assign nodes to clades for a rooted gene tree.
 *
 * @return  A list of splits for each clade.
 */
vector< vector< pair<int,int> > > 
Clades::createSimpleClades(
    MyGeneTree &geneTree, ///< a gene tree
    vector<double> &branchLengths ) ///< branch length of each node
{
    if( !mNodesOnly ) 
        throw bpp::Exception( "createSimpleClades called on not nodesOnly" );

    int nodeCount = geneTree.getNodes().size();
    vector< vector< pair<int,int> > > cladeSplits( nodeCount ); 

    mCladeCount = geneTree.getLeaves().size() + 1;
    branchLengths.resize( mCladeCount, -1 );
    createSimpleCladesAux( geneTree.getRootNode(), cladeSplits, branchLengths );

    return cladeSplits;
}


/**
 * Delete the give clades according to the given index.
 *
 * @return A mapping from the old id to the new id.
 */
vector<int> Clades::deleteClades( 
        vector<int> &cladesToDelete ) ///< list of clades to delete
{
    // convert deletion list to a map
    vector<bool> deleteV( mCladeCount, false );
    BOOST_FOREACH( int idU, cladesToDelete )
        deleteV[idU] = true;

    vector<int> idMapping( mCladeCount );

    mCladeCount = mLabelCount+1; 
    boost::unordered_map<CladeHashType,int>::iterator iter;
    for( iter = mCladeToInt.begin(); iter != mCladeToInt.end(); iter++ ) {
        int id = iter->second;
        if( deleteV[id] ) {
            if( id <= mLabelCount )
                throw bpp::Exception( "Clades:deleteClades: request to "
                                 "remove leaf or root" );
            mCladeToInt.erase( iter );
        } else if( id <= mLabelCount ) {
            // don't change leaves and root
            idMapping[id] = id;
        } else {
            // reassign id
            idMapping[id] = mCladeCount;
            iter->second = mCladeCount++;
        }
    }

    //CHECK IT
    if( (int) mCladeToInt.size() != mCladeCount )
        throw bpp::Exception("Clades::deleteClades: map size wrong" );

    return idMapping;
}


/**
Takes:
    - vector<string> LeafNames: list of leaf names
Returns:
    int: id of the corresponding clade
**/
int Clades::getCladeIntFromLeafList(vector<string> LeafNames)
{
    if( mNodesOnly ) 
        throw bpp::Exception( "MyClades::getCladeIntFromLeafList: called for nodesOnly");

    vector <int> LeafNums;

    BOOST_FOREACH(string leafName, LeafNames)
    {
        boost::unordered_map<string,int>::const_iterator iter = 
                                    mLeafClades.find( leafName ); 
        if( iter == mLeafClades.end() ) 
            throw bpp::Exception("MyClades::getCladeIntFromLeafList: clade number not found for leaf");
        int leafCladeNum = iter->second;

        // Leaves are assigned clade numbers 1-numLeaves
        // The leaf number is the position of the leaf in the bit
        // array, hence base 0 and 1 minus the clade number.
        int leafNum = leafCladeNum -1;
        if( leafNum >= mLabelCount ) 
            throw bpp::Exception("MyClades::getCladeIntFromLeafList: leafNum >= leafCladeNum");
        LeafNums.push_back(leafNum);
    }
    // set bit for leaf in clade array
    Clade clade;


#ifdef BITS
    BOOST_FOREACH(int leafNum , LeafNums)
    {
    clade.cladeV.resize( mBitBlocks );
        int blockIdx = int(leafNum/BLOCK_SIZE);
        clade.cladeV[blockIdx][leafNum%BLOCK_SIZE] = 1;
    }
#else
    BOOST_FOREACH(int leafNum , LeafNums)
    {
        clade.cladeV.resize( mLabelCount );
        clade.cladeV[leafNum] = true;
    }
#endif

    //The cladeV instance is ready to be fed to getCladeInt
    return getCladeInt(clade.cladeV,false); 
}




/**
 * If the split (triparition) exists, add to it's occurrences.  
 * If not, add  it.
 */
void CladesAndTripartitions::countSplit( 
        VectorInt split, ///< split, given by three clade numbers
        double occurrences ) ///< precomputed occurrences (polytomy)
{
    // make sure split is ordered correctly - no need to swap then
    if( split[1] > split[2] ) {
        int tmp = split[1];
        split[1] = split[2];
        split[2] = tmp;
    }
 
    // add or update split occurrences
 	boost::unordered_map<VectorInt,double>::iterator iterSplits 
        = mSplitsOccurrencesMap.find( split );
    if (iterSplits != mSplitsOccurrencesMap.end() ) {
        iterSplits->second += occurrences;
        #ifdef DEBUG
        cout << " [" << split[0] << "/" << split[1] << "," << split[2] << "]"
             << " +" << occurrences << endl;
        #endif
    } else {
        // not there, insert it
        mSplitsOccurrencesMap.insert(std::pair<VectorInt, 
                double>(split,occurrences)); 	
        #ifdef DEBUG
        cout << " [" << split[0] << "/" << split[1] << "," << split[2] << "]"
             << " " << occurrences << endl;
        #endif
    }
}

/**
 * Calculate possible splits of this node and add to split counts.
 */
void CladesAndTripartitions::calculateSplits(
    double occurrences, 
    ///< occurrences of split for polytomic rooted trees
    const MyGeneNode *node, ///< a gene tree node in recursion
    Clades::Clade clade, ///< clade for which to count splits
    Clades::Clade firstChild, ///< child of clade
    Clades::Clade secondChild, ///< other child of clade
    bool notRootSon ) ///< force adding root split if unrooted (polytomy)
{
    // Calculate splits for the node.
    VectorInt split(3);

    if( !mRooted ) {
        // if not root, split each anti-child with sibling 
        if( node->hasFather() ) {
            //for each child, split anti-child on sibling and (all - this clade)
            split[0] = firstChild.antiCladeInt; // antiChild
            split[1] = secondChild.cladeInt; // sibling clade int
            split[2] = clade.antiCladeInt; // antiChild - sibling clade
            countSplit( split, occurrences );

            split[0] = secondChild.antiCladeInt; // antiChild
            split[1] = firstChild.cladeInt; // sibling clade int
            split[2] = clade.antiCladeInt; // antiChild - sibling clade
            countSplit( split, occurrences );
        }
    }

    // For splitting clade with children.
    // Except if the tree is not polytomic (mRooted=false)
    // and this is the root node and a child is a leaf, which
    // is handled by the Clades inititzation.
    if( mRooted || clade.cladeInt != 0 
            || (firstChild.bits > 1 && secondChild.bits > 1) ) 
    {
        split[0] = clade.cladeInt;
        split[1] = firstChild.cladeInt;
        split[2] = secondChild.cladeInt;
        countSplit( split, occurrences );
    }


    if( !mRooted ) {
        // split all on this clade with anti-clade if not root or root-son
        if( node->hasFather() 
            && (notRootSon || node->getFather()->hasFather() ) )
        {
            split[0] = 0;
            split[1] = clade.cladeInt;
            split[2] = clade.antiCladeInt;
            countSplit( split, occurrences );
        }
    }
}


/**
 * Calculate double factorial.
 *
 * @returns n!!
 */
unsigned long long 
CladesAndTripartitions::doubleFactorial(
        int n, ///< input number
        bool &overflow ) ///< true if overflowed
{
    // first several
    static const unsigned long long values[24]
        = { 1, 1, 
            2, 3, 
            8, 15, 
            48, 105, 
            384, 945, 
            3840, 10395, 
            46080, 135135, 
            645120 ,2027025, 
            10321920, 34459425, 
            185794560, 654729075, 
            3715891200, 13749310575, 
            81749606400, 316234143225 };

    if( n == -1 )
        return 1;
    if( n < -1 ) {
        cerr << "doubleFactorial given " << n << endl;
        throw bpp::Exception( "CladesAndTripartitions::doubleFactorial:"
                " this implmentation doesn't handle values less than -1" );
    }
    if( n < 24 ) 
        return values[n];
   
    unsigned long long val = 0; 
    try {
        double df = boost::math::double_factorial<double>(n); 
        if( df > ULLONG_MAX ) 
            overflow = true;
        else
            val = (unsigned long long) df;
    } catch( bpp::Exception e ) {
        // throws overflow exception
        overflow = true;
    }



    return val;
}

/**
 * Clade creation with polytomies.
 *
 * Creates splits and computes their occurrences for all possible
 * combinations.
 *
 * Note: this is not the true number of occurrences, but rather
 * #occurrences/allPossibleTrees.
 *
 */
void
CladesAndTripartitions::polytomy(
    const MyGeneNode *node, ///< current gene tree node in recursion
    const int polytomySize, ///< size of polytomy
    vector<Clades::Clade> &childClades, ///< clades to split
    Clades::Clade allClade, ///< merge of child clades
    boost::unordered_map<int,bool> &seenIt, ///< clades already visited
    bool &overflow )  ///< overflow check
{
    unsigned long childCount = (unsigned long) childClades.size();

// Is this still true?
    if( childCount > log2(ULONG_MAX) ) {
        cerr << "Polytomy currently limited to " << log2(ULONG_MAX)
            << " children" << endl;
        overflow = true;
        return;
    }

    seenIt.insert( make_pair( allClade.cladeInt, true ) );

    mPolytomySizes.push_back( make_pair( allClade.cladeInt, polytomySize ) );
    if( (size_t) polytomySize != childClades.size() )  {
        // not needed for root of polytomy, add would overlap previous
        mPolytomyCladeSizes.push_back( 
                make_pair( allClade.cladeInt, childClades.size() ) );
        // Not needed for anti-clades in support calculation.
    }

    // part of occurrences equation (only calculate once)
    double df = doubleFactorial(2*polytomySize-3, overflow);
    if( overflow )
        return;

    // create all possible divisions of children into two parts
    vector<Clades::Clade> firstPartClades;
    Clades::Clade firstPart;
    float limit = pow( (float) 2, int (childCount-1) );
    for( unsigned long partition=1; partition<limit; partition++ ) {

        // assign each child to one of two partitions
        firstPartClades.clear();
        vector<Clades::Clade> secondPartClades;
        // each bit position corresponds to a child
        unsigned long bitRep = 1;  
        for( unsigned long i=0; i<childCount; i++ ) {
            // If child's bit is 1, include in first partition.
            // Last child is always in second partition.
            if( partition & bitRep ) 
                firstPartClades.push_back( childClades[i] );
            else 
                secondPartClades.push_back( childClades[i] );
            bitRep = bitRep << 1; 
        }

        // merge two partitions
        mClades.mergeClades( firstPartClades, firstPart, !mRooted );
        Clades::Clade secondPart;
        mClades.mergeClades( secondPartClades, secondPart, !mRooted );

        // calculate split occurrences 
        double occurrences = 
                doubleFactorial(2*firstPartClades.size()-3, overflow)
                * doubleFactorial( 2*secondPartClades.size()-3, overflow )
                * doubleFactorial( 
                        2*(polytomySize-childCount+1)-3, overflow)
                / df;
// should never happen since possibleTreeCounts would overflow first
        if( overflow ) 
            return;

        calculateSplits( occurrences, node, allClade, firstPart, secondPart,
                         true );

        boost::unordered_map<int,bool>::const_iterator iter;
        iter = seenIt.find( firstPart.cladeInt ); 
        if( firstPartClades.size() > 1 && iter == seenIt.end() ) {
            polytomy( node, polytomySize, firstPartClades,
                      firstPart, seenIt, overflow );
            if( overflow )
                return;
        }

        iter = seenIt.find( secondPart.cladeInt ); 
        if( secondPartClades.size() > 1 && iter == seenIt.end() ) {
            polytomy( node, polytomySize, secondPartClades,
                      secondPart, seenIt, overflow );
            if( overflow )
                return;
        }
    }
}


/**
 * Given the root node of a tree, calculate clades and the occurences
 * of every possible clade split. 
 *
 * This function recursively traverses the the tree.
 *
 * If there are more than two children of the node (polytomic), all pairwise
 * combinations of the children are used to create child clades.
 *
 * @return clade structure created for this node
 */
Clades::Clade
CladesAndTripartitions::computeSplits(
    const MyGeneNode *node, ///< current gene tree node in recursion
    bool polytomic, ///< gene tree can be polytomic
    bool &overflow,  ///< true if overflow occurrences
    int level )   ///< current recursion depth
{
    #ifdef DEBUG 
    for( int i=0; i<level; i++ ) cout << "  ";
    cout << level << " ";
    #endif

    Clades::Clade newClade; // clade to return
    if( node->isLeaf() ) {
        newClade = mClades.createClade( node->getName(), !mRooted );
        // NOTE: splits for leaves are precomputed by Clades
        #ifdef DEBUG 
        cout << "leaf: " << node->getName() << "->" 
            << newClade.cladeInt;
        if( !mRooted )
            cout << "," << newClade.antiCladeInt;
        cout << endl;
        #endif
    } else {
        // node
        int sonCount = node->getNumberOfSons();
	    if( !polytomic && sonCount != 2 ) {
            cout << "level=" << level << " sons=" 
                 << node->getNumberOfSons() << endl;
            throw bpp::Exception("CladesAndTripartitions::computeSplits:"
                    " non-binary tree node");
        }

        vector<Clades::Clade> childClades( sonCount );
        for( int idx=0; idx<sonCount; idx++ ) {
            childClades[idx] = computeSplits( node->getSon(idx), 
                                    polytomic, overflow, level+1 );    
            if( overflow )
                return newClade;
        }

        if( sonCount > 2 ) {
            mClades.mergeClades( childClades, newClade, !mRooted );
            // create polytomic subclades and calculate their splits
            boost::unordered_map<int,bool> seenIt;
            polytomy( node, sonCount, childClades, newClade, seenIt, overflow );
        } else {
            int occurrences = 1;
            mClades.mergeClades( childClades, newClade, !mRooted );
            calculateSplits( occurrences, node, newClade, childClades[0],
                             childClades[1] );
        }

        #ifdef DEBUG 
        cout << "node " << newClade.cladeInt;
        if( !mRooted )
            cout << "/" << newClade.antiCladeInt;
        cout << endl;
        #endif
    }

    if( node->hasDistanceToFather() ) {
        if( (size_t) newClade.cladeInt >= mBranchLengths.size() )
            mBranchLengths.resize( newClade.cladeInt + 1, -1 );
        mBranchLengths[newClade.cladeInt] += node->getDistanceToFather();
        if( !mRooted && newClade.cladeInt != mClades.getRootClade() ) {
            if( (size_t) newClade.antiCladeInt >= mBranchLengths.size() )
                mBranchLengths.resize( newClade.antiCladeInt + 1, -1 );
            mBranchLengths[newClade.antiCladeInt]+=node->getDistanceToFather();
        }
    }

    return newClade;
}



/**
 * Find polytomies and count the possible number of trees
 * that can be constructed from them.
 *
 * @return the possible number of trees
 */
unsigned long long
CladesAndTripartitions::findPolytomies(
    const MyGeneNode *node, ///< current gene tree node in recursion
    bool &overflow ) ///< true if number of trees overlowed 
{
    overflow = false;

    unsigned long long possibleTreeCount = 1;
    if( !node->isLeaf() ) {
        int sonCount = node->getNumberOfSons();
        if( sonCount > 2 ) {
            possibleTreeCount = doubleFactorial( 2*sonCount-3, overflow );
            if( overflow )
                return 0;
        }

        for( int idx=0; idx<sonCount; idx++ ) {
            unsigned long long childTrees
                = findPolytomies( node->getSon(idx), overflow );
            if( overflow )
                return 0;

            // Check for overflow from multiplication.
            // Not exact, but close.
            if( log2( childTrees ) + log2( possibleTreeCount )
                    > log2( ULLONG_MAX ) ) 
            {
                overflow = true;
                return 0;
            }

            possibleTreeCount *= childTrees;
        }
    }

    return possibleTreeCount;
}


/**
 * Constructor - Compute all clades and tripartitions for the given 
 * gene trees. 
 */
CladesAndTripartitions::CladesAndTripartitions( 
    char charSep, ///< separator of gene names
    vector<MyGeneTree*> &geneTrees, ///< a set of gene trees
    bool verbose, ///< print timing and stats 
    bool &overflow, ///< overflowed number of trees
    string &errStr, ///< returns any erros
    bool polytomic,   ///< trees can be polytomic (and are rooted) 
    MyGeneTree *bootstrapGeneTree,
            ///< gene tree with values for splitRatio calculation
    bool polytomicUnrooted, ///< polytomic tree is unrooted (default=false)
    bool printCladeInfo ) ///< output clade info
{
    clock_t start = clock();

    if( geneTrees.size() == 0 ) 
        throw bpp::Exception( "CladesAndTripartitions::constructor: "
                " empty geneTrees vector" );

    mGeneTreeCount = geneTrees.size();
    if( bootstrapGeneTree != NULL &&  mGeneTreeCount != 1 ) {
        errStr = "Bootstrap weighting cannot be applied to amalgamations.";
        return;
    }
    
    overflow = false;
    errStr = "";
    mRooted = false; // amalgamated trees are considered unrooted

    // Initialize clades and create a map between leaf name and 
    // its clade number
    mClades.init( charSep, *(geneTrees[0]) );

    mNewickGeneTree = bpp::TreeTemplateTools::nodeToParenthesis( 
                            *(geneTrees[0]->getRootNode()) ); 

    if( polytomic ) {
        if( mGeneTreeCount != 1 ) 
            throw bpp::Exception( "CladesAndTripartitions:: constructor "
                        "Polytomy only applies to a single gene tree." );

        if( !polytomicUnrooted )
            mRooted = true; // rooted

        // check for polytomies 
        if( verbose ) {
            bool overflowCount = false;
            unsigned long long allPossibleTrees = 
                findPolytomies( geneTrees[0]->getRootNode(), overflowCount );
            if( overflowCount || allPossibleTrees > 1 ) {
                // find max degree
                vector<MyGeneNode*> nodes = geneTrees[0]->getNodes();
                size_t maxDegree = 0;
                BOOST_FOREACH( MyGeneNode *node, nodes ) 
                    if( node->getNumberOfSons() > maxDegree )
                        maxDegree = node->getNumberOfSons();

                cout << "gene tree is polytomic: " ;
                if( overflowCount )
                    cout << "overflowed number of trees";
                else
                    cout << allPossibleTrees << " possible trees";
                cout << ", maximum out-degree = " << maxDegree << endl;
            }
        }
    }

    // compute splits
    mClades.precomputeLeaves( !mRooted );
    mClades.precomputeLeafSplits( !mRooted, mGeneTreeCount, 
                                  mSplitsOccurrencesMap );
	for( int g=0; g < mGeneTreeCount; g++) {
        computeSplits( geneTrees[g]->getRootNode(), polytomic, overflow );
        if( overflow )
            return;
    }
	
   	int cladeCount = mClades.getCladeCount(); 
    if( verbose ) 
   	    cout << "clades and splits:  " << cladeCount << " " 
            << mSplitsOccurrencesMap.size() << endl;
    

    // Convert mSplitsOccurencesMap into more efficient structures.
    // loop over all splits to find:
    //    1. occurrences: sum of split occurrencnces for each clade
    //    2. mCladeSplits: list of splits for each clade 
    //    3. splitRatio: occurence of each split as a vector, divided
    //      by clade occurrences
    //       ordered the same as mCladeSplits 
    //  NOTE: polytomic occurrences are the fraction of allPossibleTrees,
    //      i.e. the true occurrences = mSplitsRatio * allPossibleTrees
 	boost::unordered_map<VectorInt,double>::iterator iterSplits;
    mCladeSplits.resize( cladeCount );
    mSplitsRatio.resize( cladeCount );
    mCladeOccurrences.resize( cladeCount );
    mBranchLengths.resize( cladeCount, -1 );
  	for(iterSplits = mSplitsOccurrencesMap.begin(); 
            iterSplits != mSplitsOccurrencesMap.end(); iterSplits++ )
    {
        vector<int> split = iterSplits->first;
        int cladeNum = split[0];
        mCladeSplits[cladeNum].push_back( make_pair(split[1], split[2]) );
        mSplitsRatio[cladeNum].push_back( iterSplits->second );
        mCladeOccurrences[cladeNum] += iterSplits->second;
   	}

    // finish split ratio calculation with clade totals
   	if( bootstrapGeneTree != NULL ) {	  
        if( printCladeInfo ) {
            cout << "======== collapsed gene tree ==============" << endl;
            geneTrees[0]->printMe();
        }
        int modifier = checkBootstrapValues( bootstrapGeneTree, errStr );
        if( errStr != "" )
            return;
		bootstrapSplitRatio( bootstrapGeneTree, modifier, errStr,
                                printCladeInfo );
        if( errStr != "" )
            return;
    } else if( polytomic ) {
		bootstrapSplitRatio( geneTrees[0], 1, errStr, printCladeInfo, false );
        if( errStr != "" )
            return;
    } else {
		for( size_t cladeNum=0; cladeNum<mSplitsRatio.size(); cladeNum++ ) 
		    for( size_t i=0; i<mSplitsRatio[cladeNum].size(); i++ ) {
		        mSplitsRatio[cladeNum][i] /= mCladeOccurrences[cladeNum];
            }
    }

    for( size_t cladeNum=0; cladeNum<mSplitsRatio.size(); cladeNum++ ) {

        // average branch lengths
        if( cladeNum != (size_t) mClades.getRootClade() 
            && mBranchLengths[cladeNum] != -1 ) 
        {
            mBranchLengths[cladeNum] /= mCladeOccurrences[cladeNum];
        }
    }

    mSplitsOccurrencesMap.clear(); // done with it, free the space

    if( printCladeInfo ) {
        cout << "============= Clades and Tripartitions =============" << endl;
        printMe();
    }

    if( verbose ) {
        clock_t end = clock();
        double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
        cout << "split time = " << time << endl;
    }

}


/**
 * Recursively create a (cladeToNode) map of clade numbers to nodes.
 *
 * @return Clade for the given node
 */
Clades::Clade CladesAndTripartitions::mapCladesToNodesAux(
        MyGeneNode* node, ///< Current node in bootstrap gene tree.
        boost::unordered_map<string,int> &nameToInt, 
            ///< a map of leave names to clade numbers
        vector<MyGeneNode*> &cladeToNode,
            ///< mapping of clade numbers to nodes (returned result)
        bool useBootstrapValues, ///< use bootstrap values in geneTree
        string &errStr ) ///< description of any errors
{
    Clades::Clade clade;
    if( node->isLeaf() ) {
        boost::unordered_map<string,int>::iterator iter
            = nameToInt.find( node->getName() );
        if( iter == nameToInt.end() ) {
            errStr = "bootstrap gene tree doesn't have the same leaves"
                     " as the original tree";
        } else
            clade = mClades.createClade( node->getName(), false );
    } else if( useBootstrapValues && node->getNumberOfSons() != 2  ) {
        errStr = "bootstrap gene tree is not binary";
        return clade;
    } else {
        // merge clades to get parent clade
        vector<Clades::Clade> cladesToMerge;
        size_t sonCount = node->getNumberOfSons();
        for( size_t i=0; i<sonCount; i++ ) {
            Clades::Clade cladeSon = mapCladesToNodesAux( node->getSon(i), 
                                            nameToInt, cladeToNode, 
                                            useBootstrapValues, errStr );
            if( errStr != "" ) return clade;
            cladesToMerge.push_back( cladeSon );
        }

        bool isNew = mClades.mergeClades( cladesToMerge, clade, false );
        if( isNew ) {
            // clade didn't exist if it is new
            errStr = "bootstrap gene tree does not correpsond to gene tree";
            return clade;
        }

    }

    cladeToNode[clade.cladeInt] = node;

    return clade;
}


/**
 * Create map of clade numbers to nodes of the given tree.
 *
 * @return map 
 */
vector<MyGeneNode*> CladesAndTripartitions::mapCladesToNodes(
        MyGeneTree* geneTree, ///< gene tree with nodes to map
        bool useBootstrapValues, ///< use bootstrap values in geneTree
        string &errStr ) ///< description of any errors
{
    // map leaf names to clade numbers
    boost::unordered_map<string,int> nameToInt;
    for( size_t cladeNum=0; cladeNum<mSplitsRatio.size(); cladeNum++ ) {
        if( mClades.isLeaf( cladeNum ) ) {
            string leafName = mClades.getLeafName( cladeNum );
            nameToInt.insert(
                std::pair<string,int>(leafName,cladeNum));
        }
    }

    // create map of clades to bootstrap nodes
    vector<MyGeneNode*> cladeToNode( mClades.getCladeCount(), NULL );
    mapCladesToNodesAux( geneTree->getRootNode(), 
                      nameToInt, cladeToNode, useBootstrapValues, errStr );

    return cladeToNode;
}


/**
 * Validate bootstrap values and return modifier.
 *
 * @return false modifier
 */
int CladesAndTripartitions::checkBootstrapValues(
        MyGeneTree* geneTreeForBootstrap, ///< gene tree with bootstrap values
        string &errStr )  ///< return any errors
{
    // check if bootstrap values are 0-1 or 0-100
    vector<MyGeneNode*> nodes = geneTreeForBootstrap->getNodes();
    int modifier = 1;
    BOOST_FOREACH( MyGeneNode *node, nodes ) {
        if( !node->isLeaf() 
            && node->hasBranchProperty( bpp::TreeTools::BOOTSTRAP) )
        {
            double bootstrapValue = dynamic_cast< const bpp::Number<double>*> 
                 (node->getBranchProperty(bpp::TreeTools::BOOTSTRAP))
                 ->getValue(); 
            if( bootstrapValue < 0 ) {
                errStr = "Invalid bootstrap value: < 0";
            } else if( bootstrapValue > 100 ) {
                errStr = "Invalid bootstrap value: > 100";
            }  else if( bootstrapValue > 1 ) 
                modifier = 100;
        }
    }

    return modifier;
}



/**
 * Calculate split ratios using the bootstrap value from the given
 * tree. Splits not in the tree are normalized.
 */
void CladesAndTripartitions::bootstrapSplitRatio(
    MyGeneTree* geneTreeForBootstrap, ///< gene tree with bootstrap values
    int modifier, ///< adjustment to bootstrap value
    string &errStr,  ///< return any errors
    bool printCladeInfo, ///< print split ratio calculations
    bool useBootstrapValues ) ///< use bootstrap values in geneTreeForBootstrap
{
    size_t numClades = mSplitsRatio.size();

    // clear any previous values
    for( size_t C1=0; C1<numClades; C1++ ) {
        int numberSplits = mSplitsRatio[C1].size();
        for( int i=0; i<numberSplits; i++ ) 
            mSplitsRatio[C1][i] = -1;
    }

    vector< vector<int> > noncompatibleLists = mClades.getNoncompatible();
    vector<MyGeneNode*> cladeToNode 
        = mapCladesToNodes( geneTreeForBootstrap, useBootstrapValues, errStr );
    if( errStr != "" )
        return;

    // anti clade set for unrooted = those without node or polytomySize
    vector<bool> isAnti( numClades, true ); 
    if( !mRooted ) {
        for( size_t i=0; i<numClades; i++ ) 
            if( cladeToNode[i] != NULL ) 
                isAnti[i] = false;
    }

    // organize polytomy sizes (size of polytomy of which clade is embedded)
    vector<int> polytomySizes( numClades, 1 );
    for( size_t i=0; i<mPolytomySizes.size(); i++ ) {
        int cladeNum = mPolytomySizes[i].first;
        polytomySizes[cladeNum] = mPolytomySizes[i].second;

        // only clades created by polytomy have a size (not their anti-clades)
        if( !mRooted ) 
            isAnti[cladeNum] = false;
    }


    vector<int> polyCladeSizes( numClades, 1 );
    for( size_t i=0; i<mPolytomyCladeSizes.size(); i++ ) 
        polyCladeSizes[mPolytomyCladeSizes[i].first] 
                                = mPolytomyCladeSizes[i].second;

    vector<double> fC( numClades, -1 ); // clade support

    // set fC to bootstrap if clade is in gene tree
    for( size_t cladeNum=0; cladeNum<numClades; cladeNum++ ) {
        MyGeneNode *node = cladeToNode[cladeNum];
        if( node != NULL ) {
            if( node->hasBranchProperty( bpp::TreeTools::BOOTSTRAP) )
            { 
                fC[cladeNum] = dynamic_cast< const bpp::Number<double>*> 
                     (node->getBranchProperty(bpp::TreeTools::BOOTSTRAP))
                     ->getValue(); 
                fC[cladeNum] /= modifier;
            } else {
                // no bootstrap value defaults to 1
                fC[cladeNum] = 1;
            }
        }
    }

    bool overflow = false;
    // set fC to adjusted value if not in gene tree
    for( size_t cladeNum=0; cladeNum<numClades; cladeNum++ ) {
        // if not seen and not an anticlade 
        if( fC[cladeNum] == -1 && (mRooted || !isAnti[cladeNum] ) ) { 
            // find maximum non-comptabile f()
            int Cmax = -1;
            double maxVal = 1;
            for( size_t i=0; i<noncompatibleLists[cladeNum].size(); i++ ) {
                int C = noncompatibleLists[cladeNum][i];
                if( cladeToNode[C] != NULL 
                    && (Cmax == -1 || fC[C] > maxVal) ) 
                {
                    Cmax = C;
                    maxVal = fC[C];
                }
            }
            if( Cmax == -1 ) {
                if( useBootstrapValues ) 
                    throw bpp::Exception("CladesAndTripartitions::"
                        "bootstrapSplitRatio: no FC() max found" );
                Cmax = 0; // if no bootstrap values given
            }

            int polytomySize = polytomySizes[cladeNum];
            int nC1 = polyCladeSizes[cladeNum];
            int secondPart1 = polytomySize - nC1 + 1;
            double fbarC1 = doubleFactorial(2*nC1-3, overflow)
                          * doubleFactorial( 2*secondPart1-3, overflow );

            int nC = polyCladeSizes[Cmax];
            int secondPart = polytomySize - nC + 1;
            double fbarC = doubleFactorial(2*nC-3, overflow)
                          * doubleFactorial( 2*secondPart-3, overflow );
            fC[cladeNum] = (1-maxVal) * fbarC1;
            fC[cladeNum] /=  (doubleFactorial( 2*polytomySize-3, overflow )
                               - fbarC );
            if( overflow ) { 
                errStr = "bootstrapSplitRatio calculation overflowed";
                return;
            }
        }
    }

    vector< vector<double> > f( numClades );

    // split ratio calculations (C1,{C2,C3})
    for( size_t C1=0; C1<numClades; C1++ ) {
        if( !mRooted && isAnti[C1] ) // handle these later
            continue; 

        MyGeneNode *node = cladeToNode[C1];
        int numberSplits = mSplitsRatio[C1].size();
        int polytomySize = polytomySizes[C1];

        int foundSplit = -1;
        if( node != NULL && node->isLeaf() )
            foundSplit = 0;
        else if( node != NULL && node->getNumberOfSons() == 2 ) {
            MyGeneNode *child0 = node->getSon(0);
            MyGeneNode *child1 = node->getSon(1);
            // find corresponding split
            for( int i=0; i<numberSplits; i++ ) {
                int C2 = mCladeSplits[C1][i].first;
                int C3 = mCladeSplits[C1][i].second;
                MyGeneNode *splitNode0 = cladeToNode[C2];
                MyGeneNode *splitNode1 = cladeToNode[C3];
                // pointer comparision
                if( (splitNode0 == child0 && splitNode1 == child1) 
                   || (splitNode0 == child1 && splitNode1 == child0) ) 
                {
                    foundSplit = i;
                    break;
                }
            }
            // It must there. Something's wrong if not found.
            if( foundSplit == -1 ) 
                throw bpp::Exception( "CladesAndTripartitions"
                        "::bootstrapSplitRatio: could not find split" );
        }

        // do the found split first
        double gHatPrime = 0;
        double foundF = 0;
        if( foundSplit != -1 ) {
            double min = fC[C1];
            if( !node->isLeaf() ) {
                int C2 = mCladeSplits[C1][foundSplit].first;
                int C3 = mCladeSplits[C1][foundSplit].second;
                if( fC[C2] < min )
                    min = fC[C2];
                if( fC[C3] < min )
                    min = fC[C3];
                if( polytomySize > 1 ) {
                    int nC2 = polyCladeSizes[C2];
                    int nC3 = polyCladeSizes[C3];
                    gHatPrime = doubleFactorial(2*nC2-3, overflow)
                              * doubleFactorial(2*nC3-3, overflow)
                              * doubleFactorial(
                                      2*(polytomySize-nC2-nC3+1)-3, overflow);
                    gHatPrime /= doubleFactorial(2*polytomySize-3, overflow);
                }
            }
            foundF = min;
            mSplitsRatio[C1][foundSplit] = min / fC[C1];

        }

        // do other root splits later
        if( !mRooted && C1 == (size_t) mClades.getRootClade() )
            continue;
    
        for( int i=0; i<numberSplits; i++ ) {
            if( foundSplit == i ) {
                if( printCladeInfo )
                    f[C1].push_back( foundF );
                continue;
            }


            int C2 = mCladeSplits[C1][i].first;
            int C3 = mCladeSplits[C1][i].second;
            int nC1 = polyCladeSizes[C1];
            int nC2 = polyCladeSizes[C2];
            int nC3 = polyCladeSizes[C3];
            if( fC[C2] == -1 || fC[C3] == -1 )
                throw bpp::Exception( "CladesAndTripartitions"
                        "::bootstrapSplitRatio: fC not set for a clade" );

            double gC1 = doubleFactorial(2*nC1-3, overflow)
                       * doubleFactorial(2*(polytomySize-nC1+1)-3, overflow);
            gC1 /= doubleFactorial(2*polytomySize-3, overflow);
                        
            double gHat = doubleFactorial(2*nC2-3, overflow)
                        * doubleFactorial(2*nC3-3, overflow)
                        * doubleFactorial(
                                2*(polytomySize-nC2-nC3+1)-3, overflow );
            gHat /= doubleFactorial(2*polytomySize-3, overflow);
            if( node != NULL ) {
                // C1 in binary tree
                mSplitsRatio[C1][i] = ( 1 - foundF/fC[C1] )
                                    * gHat / (gC1 - gHatPrime);
                                    // Bug fix 13/6/16 - EHJ
                                    //* gHat / ( 1 - foundF ); 
                                    // old equation
                                    //* gHat / ( 1 - gHatPrime );
            } else {
                mSplitsRatio[C1][i] = gHat/gC1;
                /** first equation change on 13/16/16 
                fbar = doubleFactorial(2*nC2-3, overflow)
                     * doubleFactorial(2*nC3-3, overflow);
                fBar /= doubleFactorial(2*(nC2+nC3)-3, overflow);
                mSplitsRatio[C1][i] = fBar;
                */
                // Old equation 
                //mSplitsRatio[C1][i] = gHat;
            }

            if( printCladeInfo )
                f[C1].push_back( mSplitsRatio[C1][i] * fC[C1] );

            if( overflow ) { 
                errStr = "bootstrapSplitRatio calculation overflowed";
                return;
            }
        }
    }

    if( !mRooted ) {
        size_t rootClade = mClades.getRootClade();

        // unrooted, find root splits and set anti-clade to the same,
        // for clades in Gb (node != NULL)
        int numberRootSplits = mSplitsRatio[rootClade].size();
        vector<int> antiCladeIds( numClades );
        for( int i=0; i<numberRootSplits; i++ ) {
            int C2 = mCladeSplits[rootClade][i].first;
            int C3 = mCladeSplits[rootClade][i].second;
            antiCladeIds[C2] = C3;
            antiCladeIds[C3] = C2;
            if( fC[C2] == -1 && fC[C3] == -1 ) 
                throw bpp::Exception("CladesAndTripartitions::"
                    "bootstrapSplitRatio: unrooted split has no fC" );
            double fVal = -1;
            if( fC[C3] == -1 && fC[C2] != -1 ) {
                fC[C3] = fC[C2];
                fVal = fC[C2];
            }
            if( fC[C2] == -1 && fC[C3] != -1 ) {
                fC[C2] = fC[C3];
                fVal = fC[C3];
            }

            mSplitsRatio[rootClade][i] = 1;
            if( printCladeInfo )
                f[rootClade].push_back( fVal );
        }

  
        // All anti-clade splits are rotated versions of the clade
        // split and have the same support as the unrotated version.
        for( size_t C1=0; C1<numClades; C1++ ) {
            if( !isAnti[C1] )
                continue;

            int numberSplits = mSplitsRatio[C1].size();
            for( int i=0; i<numberSplits; i++ ) {
                int C2 = mCladeSplits[C1][i].first;
                int C3 = mCladeSplits[C1][i].second;

                // One child will be the anti-clade of the C1 of the
                // split needed. The other child will be a child
                // in the unrotated split.
                int anti = -1;
                int other = -1;
                if( isAnti[C2] ) {
                    anti = antiCladeIds[C2];
                    other = C3;
                } else if( isAnti[C3] ) {
                    anti = antiCladeIds[C3];
                    other = C2;
                } else if( !isAnti[antiCladeIds[C2]] ) {
                    anti = antiCladeIds[C2];
                    other = C3;
                } else if( !isAnti[antiCladeIds[C3]] ) {
                    anti = antiCladeIds[C3];
                    other = C2;
                } else
                    throw bpp::Exception( "CladesAndTripartitions"
                        "::bootstrapSplitRatio: no anti-clade found" );

                if( isAnti[other] )
                    throw bpp::Exception( "CladesAndTripartitions"
                        "::bootstrapSplitRatio: other clade not right" );

                // search unrotated split
                int n = mCladeSplits[anti].size();
                bool found = false;
                for( int j=0; j<n; j++ ) {
                    if( mCladeSplits[anti][j].first == other 
                        || mCladeSplits[anti][j].second == other ) 
                    {
                        found = true;
                        mSplitsRatio[C1][i] = mSplitsRatio[anti][j];
                        if( printCladeInfo )
                            f[C1].push_back( f[anti][j] );
                        break;
                    }
                }
                if( !found )
                    throw bpp::Exception( "CladesAndTripartitions"
                        "::bootstrapSplitRatio: rotated split not found" );
            }
        }
    }

    // check that all split ratios are filled
    for( size_t C1=0; C1<numClades; C1++ ) {
        int numberSplits = mSplitsRatio[C1].size();
        for( int i=0; i<numberSplits; i++ ) 
            if( mSplitsRatio[C1][i] == -1 ) {
                throw bpp::Exception( "CladesAndTripartitions"
                        "::bootstrapSplitRatio: splits ratio not set" );
            }
    }

    if( printCladeInfo ) {
        cout << "===========BOOTSTRAP WEIGHTING=======" << endl;
        for( int cladeNum=0; cladeNum<mClades.getCladeCount(); cladeNum++ ) {
            cout << "clade id=" << cladeNum;
            MyGeneNode *node = cladeToNode[cladeNum];
            if( mClades.isLeaf( cladeNum ) )
                cout << " LEAF:" << mClades.getLeafName( cladeNum );
            cout << endl;
            if( node != NULL ) { 
                cout << "  gene tree id=" << node->getId();
                if( node->hasBranchProperty( bpp::TreeTools::BOOTSTRAP) )
                { 
                    double bs = dynamic_cast< const bpp::Number<double>*> 
                      (node->getBranchProperty(bpp::TreeTools::BOOTSTRAP))
                         ->getValue(); 
                    cout << " bootstrap=" << bs;
                }
                cout << endl;
            }
            cout << "  f(C) = " << fC[cladeNum] << endl;
            cout << "  splits:" << endl;
            for( size_t i=0; i<mCladeSplits[cladeNum].size(); i++ ) {
                pair<int,int> split = getCladeSplit( cladeNum, i );
                cout << "    " << split.first << ", " << split.second
                     << "  f=" << f[cladeNum][i] 
                     << "   cond prob=" << mSplitsRatio[cladeNum][i] << endl;
            }
        }
    }
}


/**
 * Constructor - Create a clades and tripartitions structure for the given 
 * gene tree, but instead of clades, there are just nodes,
 * i.e. use rooted tree.
 *
 */
CladesAndTripartitions::CladesAndTripartitions( 
        char charSep, ///< separator of gene names
        MyGeneTree &geneTree ) ///< a gene tree
{ 
    mGeneTreeCount = 1;

    // assign internal node clade numbers
    mClades.init( charSep, geneTree, true );
    mCladeSplits = mClades.createSimpleClades( geneTree, mBranchLengths );

    //  set all splitsRatios to 1 as
    //  each split is a gene node and is seen once
   	int cladeCount = mClades.getCladeCount(); 
    mSplitsRatio.resize( cladeCount );
    for( int i=0; i<cladeCount; i++ ) 
        mSplitsRatio[i].push_back( 1 );
}


/**
 * ALE constructor
 */
CladesAndTripartitions::CladesAndTripartitions( 
    char charSep, ///< separator of gene names
    string aleFileName,  ///< path to an ALE file
    string &errStr, ///< return any errors in reading the file
    bool verbose ) ///< print timing and stats 
{
    int splitCount = 0;
    string resultErr = loadALE( aleFileName, charSep, splitCount );
    if( resultErr != "" )
        errStr = resultErr;

    if( verbose ) 
   	    cout << "clades and splits:  "
   	        << mClades.getCladeCount() << " " << splitCount << endl;
}

/**
 * Print clades and splits.
 */
void CladesAndTripartitions::printMe() 
{
    vector<string> cladeStrings = mClades.getCladesAsStrings();

    for( int cladeNum=0; cladeNum<mClades.getCladeCount(); cladeNum++ ) {
        cout << cladeNum << " " << cladeStrings[cladeNum];
        if( mClades.isLeaf( cladeNum ) )
            cout << " LEAF:" << mClades.getLeafName( cladeNum );
        cout << endl;
        for( size_t i=0; i<mCladeSplits[cladeNum].size(); i++ ) {
            pair<int,int> split = getCladeSplit( cladeNum, i );
            double ratio = getSplitRatio( cladeNum, i );
            cout << "  " << split.first << ", " << split.second
                 << "  " << ratio << endl;
        }
    }
}

/**
 * Recursive helper function for getRootMpAm.
 *
 * To avoid recalculation, vals stores the values are calculated
 * and must be intialized to -1 with cladeCount size.
 */
void CladesAndTripartitions::getRootMpAmAux( 
    int id_u,   ///< current clade in the recursion
    vector<double> &vals ) ///< values already calculated
{
    int splitCount = getSplitCount( id_u );
    double bestCost = std::numeric_limits<double>::max();

    // find split with the best cost
    for( int i=0; i<splitCount; i++ ) {
        pair<int,int> cladeSplit = getCladeSplit( id_u, i );
        if( cladeSplit.first == -1 ) {
            bestCost = 0;
            break;
        }

        // fill children values if not there
        if( vals[cladeSplit.first] == -1 ) 
            getRootMpAmAux( cladeSplit.first, vals );
        if( vals[cladeSplit.second] == -1 ) 
            getRootMpAmAux( cladeSplit.second, vals );

        // calculate cost
        double cost = vals[cladeSplit.first] + vals[cladeSplit.second]
                    - log10( getSplitRatio( id_u, i ) );

        if( cost < bestCost )
            bestCost = cost;
    } 
    vals[id_u] = bestCost;
}


/**
 * MpAm is a value used to recalculate costs.
 *
 * @return root MpAm
 */
double CladesAndTripartitions::getRootMpAm() 
{
    int rootClade = mClades.getRootClade();
    vector<double> vals( mClades.getCladeCount(), -1 );
    getRootMpAmAux( rootClade, vals );

    return vals[rootClade];
}



/**
 * Helper for postorderMapping
 */
void CladesAndTripartitions::postOrderMappingAux(
    vector<int> &mapping, ///< map of clade ids to DF post-order numbering
    int &pOrd, ///< current depth-first post-order number
    int id_u ) ///< current clade in the recursion
{
    pair<int,int> cladeSplit = getCladeSplit( id_u, 0 );
    if( cladeSplit.first != -1 ) {
        postOrderMappingAux( mapping, pOrd, cladeSplit.first );
        postOrderMappingAux( mapping, pOrd, cladeSplit.second);
    }
    mapping[id_u] = pOrd++;
}


/**
 * Create a map from clade ids to a depth-first post-order numbering
 * of clades.
 *
 * Only returns the first split, of which there is only one in
 * a simple tree (non-amalgamated).
 *
 * @return A vector mapping clade ids to a post order traversal.
 */
vector<int> CladesAndTripartitions::getPostOrderMapping() 
{
    vector<int> mapping( mClades.getCladeCount() ); // map to fill
    int pOrd = 0; // counter for recursion
    postOrderMappingAux( mapping, pOrd, mClades.getRootClade() );
    return mapping;
}



/**
 * Finds size of clade (number of leaves in subtree rooted at clade ).
 *
 * @return size of clade
 */
int CladesAndTripartitions::getCladesBySizeAux( 
    int cladeNum, ///< current clade in the recursion
    vector<int> &cladeSizes ) ///< size list to fill
{
    if( cladeSizes[cladeNum] != -1 )
        return cladeSizes[cladeNum];

    // get any clade split for this clade (all pairs have the same sum)
    int size = 1; // leaf size, default
    pair<int,int> cladeSplit = mCladeSplits[cladeNum][0];
    if( cladeSplit.first != -1 ) {
        size = getCladesBySizeAux( cladeSplit.first, cladeSizes )
             + getCladesBySizeAux( cladeSplit.second, cladeSizes );
    }
    cladeSizes[cladeNum] = size;

    return size;
}



/**
 * Get a list of pairs <size of clade, clade id> and sort
 * by size, smallest first.
 *
 * Fill sizePairs size of clade and id pairs. The size of the
 * clade is the number of leaves in the subtree. Sort the list
 * by size, smallest first.
 *
 * @return sorted list of <clade size, clade id> pairs.
 */
vector< vector<int> > CladesAndTripartitions::getCladesBySize()
{
    // store sizes to avoid recalculation during recursion
    vector<int> cladeSizes( mCladeSplits.size(), -1 );
    int maxSize = getCladesBySizeAux( mClades.getRootClade(), cladeSizes );
    vector< vector<int> > cladesBySize( maxSize+1 );
    for( size_t cladeNum = 0; cladeNum < mCladeSplits.size(); cladeNum++ ) {
        int size = getCladesBySizeAux( cladeNum, cladeSizes );
        if( size > 0 ) 
            cladesBySize[size].push_back( cladeNum );
    }

    return cladesBySize;
}

/**
 * Create list of clade parents, indexed by clade.
 *
 * Accessing the ith list element gives the list of clade ids 
 * for the ith clade's parentx.
 *
 * @return A list of parent ids for each clade.
 */
vector< vector<int> > CladesAndTripartitions::getCladeParents() {
    vector< vector<int> > cladeParents;
    cladeParents.resize( mCladeSplits.size() );
    pair<int,int> split;
    for( size_t u = 0; u<mCladeSplits.size(); u++ ) {
        int splitCount = getSplitCount( u );
        for( int i=0; i<splitCount; i++ ) {
            pair<int,int> split = getCladeSplit( u, i );
            if( split.first != -1 ) {
                cladeParents[split.first].push_back( u );
                cladeParents[split.second].push_back( u );
            }
        }
    }
    return cladeParents;
}

/**
 * Create a gene tree with post order ids in the bootstrap field.
 *
 * @return root of the new tree
 */
MyGeneTree *CladesAndTripartitions::getTree() {
    int pOrd = 0;
    MyGeneNode *node = getTreeAux( pOrd, mClades.getRootClade() );
    return new MyGeneTree( *node );
}

/**
 * Auxilary function for getTree.
 *
 * @return newly created node
 */
MyGeneNode *CladesAndTripartitions::getTreeAux(
    int &pOrd, ///< current depth-first post-order number
    int idU ) ///< current clade in the recursion
{
    MyGeneNode *node = new MyGeneNode();
    pair<int,int> cladeSplit = getCladeSplit( idU, 0 );
    if( cladeSplit.first != -1 ) {
        node->addSon( getTreeAux( pOrd, cladeSplit.first ) );
        node->addSon( getTreeAux( pOrd, cladeSplit.second ) );
        node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, 
                                bpp::Number<double>(pOrd));
    } else 
        node->setName( mClades.getLeafName( idU ) );

    if( mBranchLengths[idU] != -1 )
        node->setDistanceToFather( mBranchLengths[idU] );

    pOrd++;

    return node;
}

/**
 * Parse a double from a string.
 *
 * @return True if an double is found.
 */
bool CladesAndTripartitions::toDouble( 
        string in,  ///< string to parse
        double &out )  ///< parsed double 
{
    try {
        out = bpp::TextTools::toDouble( in.c_str() );
    } catch( bpp::Exception e ) {
        return false;
    }
    return true;
}

/**
 * Parse an integer from a string.
 *
 * @return True if an integer is found.
 */
bool CladesAndTripartitions::toInt( 
        string in,  ///< string to parse
        int &out )  ///< parsed integer
{
    try {
        out = bpp::TextTools::toInt( in.c_str() );
    } catch( bpp::Exception e ) {
        return false;
    }
    return true;
}


/**
 * Read clades and tripartition in ALE format from a file.
 *
 * @return A description of an error or the empty string if successful.
 */
string CladesAndTripartitions::loadALE(
        string aleFileName,   ///< file path of ALE
        char charSep, ///< separator of gene names
        int &splitCount ) ///< number of splits read from file
{
    ifstream fileStream (aleFileName.c_str());
    if( !fileStream.is_open() ) 
        return "Could not open ALE file.";
    
    vector<string> tokens;
    string reading="#nothing";
    int cladeCount;
    vector<string> bipStrings;
    vector<string> blsStrings;
    vector<string> dipStrings;
    vector<string> leafStrings;
    vector<string> setStrings;
    while( !fileStream.eof() ) {
        string line;
        getline( fileStream, line );
        boost::trim( line );

        if( line == "" ) {
            // ignore blank lines
        } else if( boost::find_first(line, "#") ) {
            reading = line;
        } else if( reading == "#constructor_string" ) {
            // use tree string to initialize clades
            // NOT USED
            //MyGeneTree geneTree = MyGeneTree::readMyTreeFromString( line );
            //mClades.init( charSep, geneTree, true );
            mNewickGeneTree = line;
            reading = "#nothing";
        } else if( reading == "#observations" ) {
            // number of trees used to build ALE
            if( !toInt( line, mGeneTreeCount ) )
                return "Observations is not an integer."; 
            reading = "#nothing";
        } else if( reading == "#Bip_counts" ) {
            // clade counts,  two ints 
            bipStrings.push_back( line );
        } else if( reading == "#Bip_bls" ) {
            // int float
            // Sum of the branch lengths associated to the bipartitions.
            blsStrings.push_back( line );
        } else if( reading == "#Dip_counts" ) {
            // Split counts (four ints)
            // Contains the frequency of triplets: mother clade and 
            // its two daughter clades. Map between the bipartition id 
            // of the mother clade and another map containing a set of 
            // bipartition ids (couldn't it be just a pair, in the case 
            // of bifurcating trees?) and the frequency of the 
            // associated triplet of bipartitions. 
            dipStrings.push_back( line );
        } else if( reading == "#last_leafset_id" ) {
            // clade count, excluding root (int)
            if( !toInt( line, cladeCount ) )
                return "last_leafset_id is not an integer."; 
            cladeCount++;
            reading = "#nothing";
        } else if( reading == "#leaf-id" ) {
            // leaf ids, string int
            leafStrings.push_back( line );
        } else if( reading == "#set-id" ) {
            // clade id ->  leafs in clade, int : int list
            setStrings.push_back( line );
        } else if( reading == "#END" ) {
            break; // done
        } else {
            return "Unknown state - bad file format.";
        }
    }
    fileStream.close();

    // map leaf ids, which are numbered from 1
    vector<string> leafNames( leafStrings.size() );
    BOOST_FOREACH( string line, leafStrings ) {
        boost::split( tokens, line, boost::is_any_of("\t "),
                      boost::token_compress_on );
        if( tokens.size() != 2 )
            return "A leaf-id does not have two tokens.";
        string leafName = tokens[0];
        int leafId;
        if( !toInt( tokens[1], leafId ) )
            return "A leaf_id is not an integer.";
        leafNames[leafId-1] = leafName;
    }

    int rootNum = mClades.getRootClade();
    mCladeOccurrences.resize( cladeCount );
    mCladeSplits.resize( cladeCount );
    mBranchLengths.resize( cladeCount, -1 );
    mSplitsRatio.resize( cladeCount );
    mCladeOccurrences[rootNum] = 0; // root count not in ALE 

    // initialize clades and precompute leaf ids
    mClades.init( charSep, leafNames, false );
    mClades.precomputeLeaves( true );

    // create a clade for each leaf
    vector<Clades::Clade> cladeMap( cladeCount ); // ale # -> cladeNum
    vector<Clades::Clade> leafClades( leafNames.size() );
    for( size_t i=0; i<leafNames.size(); i++ ) {
        leafClades[i] = mClades.createClade( leafNames[i], true );

        // stats for leaf clades
        int cladeNum = leafClades[i].cladeInt;
        mCladeOccurrences[cladeNum] = mGeneTreeCount;
        mCladeSplits[cladeNum].push_back( make_pair(-1, -1) );
        mSplitsRatio[cladeNum].push_back( mGeneTreeCount );
        int antiCladeNum = leafClades[i].antiCladeInt;
        mSplitsRatio[rootNum].push_back( mGeneTreeCount );
        mCladeOccurrences[rootNum] += mGeneTreeCount;
        if( cladeNum < antiCladeNum )
            mCladeSplits[rootNum].push_back( make_pair(cladeNum, antiCladeNum));
        else
            mCladeSplits[rootNum].push_back( make_pair(antiCladeNum, cladeNum));
    }

    BOOST_FOREACH( string line, setStrings ) {
        boost::split( tokens, line, boost::is_any_of("\t "),
                      boost::token_compress_on );
        if( tokens.size() < 3 )
            return "A set-id has too few tokens.";
        if( tokens[1] != ":" )
            return "A set-id is missing colon."; 
        int setId; 
        if( !toInt( tokens[0], setId ) )
            return "A set-id is not an integer.";
        if( setId >= cladeCount )
            return "Bad set-id.";

        vector<Clades::Clade> cladesToMerge;
        for( size_t i = 2; i<tokens.size(); i++ ) {
            int leafId;
            if( !toInt( tokens[i], leafId ) )
                return "A set-id contains a non-integer.";
            if( leafId > (int) leafClades.size() )
                return "A set-id has an integer that does not correspond to"
                       " a leaf.";
            cladesToMerge.push_back( leafClades[leafId-1] );
        }
        if( cladesToMerge.size() == 1 ) {
            cladeMap[setId] = cladesToMerge[0];
        } else {
            Clades::Clade clade;
            mClades.mergeClades( cladesToMerge, clade, true );
            cladeMap[setId] = clade;
        }
    }

    // create root
    Clades::Clade rootClade;
    mClades.mergeClades( leafClades, rootClade, true );

    if( mClades.getCladeCount() != cladeCount )
        return "Created clade count does not match last_leafset_id";
    // bip counts to clade counts
    BOOST_FOREACH( string line, bipStrings ) {
        boost::split( tokens, line, boost::is_any_of("\t "),
                      boost::token_compress_on );
        if( tokens.size() != 2 )
            return "A bip count does not have two tokens.";
        int setId;
        double val;
        if( !toInt( tokens[0], setId) || !toDouble( tokens[1], val ) )
            return "A bip count contains a invalid number.";
        if( setId > cladeCount ) 
            return "A bip count has an invalid set number.";
        int cladeNum = cladeMap[setId].cladeInt;
        mCladeOccurrences[cladeNum] = (unsigned long long) val;

        int antiCladeNum = cladeMap[setId].antiCladeInt;
        if( cladeNum < antiCladeNum ) { // don't add split twice
            mCladeSplits[rootNum].push_back( make_pair(cladeNum, antiCladeNum));
            mSplitsRatio[rootNum].push_back( val );
            mCladeOccurrences[rootNum] += val;
        }
    }

    // branch lengths
    BOOST_FOREACH( string line, blsStrings ) {
        boost::split( tokens, line, boost::is_any_of("\t "),
                      boost::token_compress_on );
        if( tokens.size() != 2 )
            return "A bip bls count does not have two tokens.";
        int setId;
        double val;
        if( !toInt( tokens[0], setId) || !toDouble( tokens[1], val ) )
            return "A bip bls count contains a invalid number.";
        if( setId > cladeCount ) 
            return "A bip bls count has an invalid set number.";
        int cladeNum = cladeMap[setId].cladeInt;
        mBranchLengths[cladeNum] = val/mCladeOccurrences[cladeNum];
    }

    // dip counts to splits
    BOOST_FOREACH( string line, dipStrings ) {
        boost::split( tokens, line, boost::is_any_of("\t "),
                      boost::token_compress_on );
        if( tokens.size() != 4 )
            return "A dip count does not have four tokens.";
        int setId, first, second;
        double val;
        if( !toInt( tokens[0], setId ) || !toInt( tokens[1], first ) 
            || !toInt( tokens[2], second ) || !toDouble( tokens[3], val ) )
            return "A dip count contains a non-integer.";
        if( setId > cladeCount || first > cladeCount || second > cladeCount )
            return "A dip count has an invalid set number.";
        int cladeNum = cladeMap[setId].cladeInt;
        mCladeSplits[cladeNum].push_back( 
                            make_pair( cladeMap[first].cladeInt, 
                                       cladeMap[second].cladeInt) );
        mSplitsRatio[cladeNum].push_back( val );
    }


    // finish split ratio calculation with clade totals
    for( size_t cladeNum=0; cladeNum<mSplitsRatio.size(); cladeNum++ ) {
        unsigned long long sum = 0;
        for( size_t i=0; i<mSplitsRatio[cladeNum].size(); i++ ) {
            splitCount++;
            sum += mSplitsRatio[cladeNum][i];
            mSplitsRatio[cladeNum][i] /= mCladeOccurrences[cladeNum];
        }
        if( sum != mCladeOccurrences[cladeNum] ) {
            cout << cladeNum << ": sum=" << sum << " != " 
                 << mCladeOccurrences[cladeNum] << endl;
            return "Dip and bip counts do not correspond.";
        }
    }

    return "";
}


/**
 * Save clade and tripartition information to a file in ALE format.
 *
 */
void CladesAndTripartitions::saveALE(
        string aleFileName )   ///< file path for output of ALE
{
    ofstream fout( aleFileName.c_str() );

   	int cladeCount = mClades.getCladeCount(); 

    //must be first!
    fout << "#constructor_string" << endl;
    fout << mNewickGeneTree << ";" << endl;

    fout << "#observations" <<endl; // number of trees
    fout << mCladeOccurrences[1] << endl; // leaves are in all trees

    fout << "#Bip_counts" << endl;
    for( size_t i=0; i<mCladeOccurrences.size(); i++ ) 
        if( !mClades.isLeaf( i ) && i != (size_t) mClades.getRootClade() )
            fout << i << "\t" << mCladeOccurrences[i] << endl;

    fout << "#Bip_bls" << endl;
    for( size_t i=0; i<mCladeOccurrences.size(); i++ ) 
        if( i != (size_t) mClades.getRootClade() ) 
            fout << i << "\t" << mBranchLengths[i]*mCladeOccurrences[i] << endl;

    fout << "#Dip_counts" << endl;
    for( size_t i=0; i<mCladeSplits.size(); i++ ) {
        if( !mClades.isLeaf( i ) && i != (size_t) mClades.getRootClade() ) {
            vector< pair<int,int> > splits = mCladeSplits[i];
            for( size_t j=0; j<splits.size(); j++ ) {
                int occ = floor(mSplitsRatio[i][j] * mCladeOccurrences[i]+0.5);
                fout << i << "\t" << mCladeSplits[i][j].first << "\t" 
                    << mCladeSplits[i][j].second << "\t" << occ << endl;
            }
        }
    }

    fout << "#last_leafset_id" << endl;
    fout << cladeCount-1 << endl;

    fout << "#leaf-id" << endl;
    for( int i=0; i<cladeCount; i++ ) {
        if( mClades.isLeaf( i ) )
            fout << mClades.getLeafName( i ) << "\t" << i << endl;
    }

    fout << "#set-id" << endl;
    vector<int> sortedClades = mClades.getSortedClades();
    vector< vector<int> > leafLists( cladeCount );
    BOOST_FOREACH( int u, sortedClades ) {
        if( u == mClades.getRootClade() )
            continue;
        if( mClades.isLeaf( u ) ) {
            leafLists[u].push_back( u );
        } else if( mCladeSplits[u].size() > 0 ) {
            int firstSplit = mCladeSplits[u][0].first;
            int secondSplit = mCladeSplits[u][0].second;
            leafLists[u].insert( leafLists[u].end(),
                leafLists[firstSplit].begin(), leafLists[firstSplit].end() );
            leafLists[u].insert( leafLists[u].end(),
                leafLists[secondSplit].begin(), leafLists[secondSplit].end() );
            sort( leafLists[u].begin(), leafLists[u].end() );
        } else {
            throw bpp::Exception ("CladesAndTripartitions::saveALE: "
                     " clade split has no elements" );
        }
        fout << u << "\t:";
        BOOST_FOREACH( int id, leafLists[u] ) 
            fout << "\t" << id;
        fout << endl;
    }

    fout << "#END" << endl;
    fout.close();
}

/**
 * After deleting clades, update mCladeSplits and mSplitsRatio, 
 * removing splits with deleted clades. Also update mCladeOccurrences.
 *
 * @return List of clades with no remaining splits.
 */
vector<int> CladesAndTripartitions::redoTripartitions(
        vector<int> &deletedClades, ///< list of deleted clades
        vector<int> &cladeIdMapping ) ///< new clade ids
{
    int newCount = mClades.getCladeCount();
    int oldCount = mCladeSplits.size();

    // convert deletion list to a map
    vector<bool> deletedV( oldCount, false );
    BOOST_FOREACH( int idU, deletedClades )
        deletedV[idU] = true;

    vector<int> emptyClades; // list of empty splits, if any
    vector<bool> inSplits( newCount, false );

    vector< vector< pair<int,int> > > newCladeSplits( newCount );
    vector< vector<double> > newSplitsRatio( newCount );
    vector<double> newCladeOccurrences( newCount );
    for( int idU=0; idU<oldCount; idU++ ) {
        if( !deletedV[idU] ) {
            int newId = cladeIdMapping[idU];
            for( size_t i=0; i<mCladeSplits[idU].size(); i++ ) {
                int firstId = mCladeSplits[idU][i].first;
                int secondId = mCladeSplits[idU][i].second;
                if( mCladeSplits[idU][i].first == -1 // leaf
                    || (!deletedV[firstId] && !deletedV[secondId] ) )
                {
                    // note clades appearing in splits
                    inSplits[cladeIdMapping[firstId]] = true;
                    inSplits[cladeIdMapping[secondId]] = true;

                    // neither excluded - keep it
                    pair<int,int> split;
                    if( mCladeSplits[idU][i].first == -1 )
                        // leaves don't change
                        split = mCladeSplits[idU][i];
                    else 
                        split = make_pair( cladeIdMapping[firstId],
                                           cladeIdMapping[secondId] );
                        newCladeSplits[newId].push_back( split );
                    newSplitsRatio[newId].push_back( mSplitsRatio[idU][i] );
                }
            }
            newCladeOccurrences[newId] = mCladeOccurrences[idU];
            if( newCladeSplits[newId].size() == 0 ) {
                inSplits[newId] = true; // not there, but marked for deletion
                emptyClades.push_back( newId );
            }
        }
    }

    for( int i=1; i<newCount; i++ ) {
        if( !inSplits[i] ) 
            emptyClades.push_back( i );
    }

    mCladeSplits = newCladeSplits;
    mSplitsRatio = newSplitsRatio;
    mCladeOccurrences = newCladeOccurrences;

    return emptyClades;
}

/**
 * Remove clades with occurrences less than the given amount.
 *
 * Additionally, clades that become empty (no splits) or clades
 * not in any split are also removed. If this occurs to a leaf
 * or the root, -1 is returned as this is an invalid state.
 *
 * @return The number of clades removed or -1 if too many are removed.
 */
int CladesAndTripartitions::excludeClades( float percentToExclude ) 
{
    // identify clades to delete
    int oldCladeCount = mClades.getCladeCount();
    vector<int> cladesToDelete;
    for( int idU=0; idU<oldCladeCount; idU++ ) {
        double frac = (double) mCladeOccurrences[idU]/mGeneTreeCount*100;
        if( frac <= percentToExclude ) 
            cladesToDelete.push_back( idU );
    }

    // Loop to propagate deletions if a clade has no more splits.
    while( cladesToDelete.size() > 0 ) {
        // check if a leaf or root is marked for deletion
        int rootIdU = mClades.getRootClade();
        BOOST_FOREACH( int idU, cladesToDelete ) {
            if( idU == rootIdU || mClades.isLeaf( idU ) ) 
                return -1;
        }

        // delete clades
//cout << "deleting " << cladesToDelete.size() << " clades"  << endl;
        vector<int> cladeIdMapping = mClades.deleteClades( cladesToDelete );

        // returns list of clades with no splits
        cladesToDelete = redoTripartitions( cladesToDelete, cladeIdMapping );
    }

    return oldCladeCount - mClades.getCladeCount();
}

