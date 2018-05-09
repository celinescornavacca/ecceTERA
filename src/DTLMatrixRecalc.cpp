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

@section DESCRIPTION
Functions to compute amalgamation matrix wrapped in a class in order
to have state variables accesible from sub-functions.  
The constructors are used to set the input and parameters before
calling the function calculateMatrix (or recalculateMatrix).

To calculate the matrix, clades are visited from leaves to root.
For each time slice, costs for each clade species node combination
are calculated, initially without transferLoss, which is calculated
after the eestReciever (lowest cost species node) is found.

The printCSVmatrix outputs the matrix for subsequent recalculations.
The matrix and bestReceiver data can be written to a file and
read for a fast recalculation for changes in the species tree.

If alpha is added to calculate transfer from the dead, it is 
assumed that alpha was added as the last species node in the time
slice.

The species trees and gene trees are assumed to be binary.

The backtrack function returns the consensus gene tree. To calculate
it, all of the bestSplits must be known, which isn't the case
for recalculated matrices.

**/


#include "DTLMatrixRecalc.h"
#include <boost/foreach.hpp>


// recalculate events used to create number which indicates the
// events to calculate
static const int N_EVENT = 1;  ///< null 
static const int SL_EVENT = 2; ///< split/loss
static const int D_EVENT = 4;  ///< duplication
static const int S_EVENT = 8;  ///< split
static const int T_EVENT = 16;  ///< transfer
static const int ALL_EVENT = 31; ///< all of the above 


/**
 * Compute costs for the current split.
 */
void DTLMatrixRecalc::computeOptimaForCladeSplit(
        int toCompute,  ///< bit index with costs to compute
        int splitIdx,   ///< current split index
        pair<int,int> cladeSplit, ///< current split
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{

    // if not a leaf
    if( cladeSplit.first != -1 && 
       (toCompute & T_EVENT || toCompute & D_EVENT || toCompute & S_EVENT ))
    {
        int idUl = cladeSplit.first;
        int idUr = cladeSplit.second;
  
#ifdef SPEED
        double costThisSplit = 0;
#else
        double splitRatio = mCladesTrips->getSplitRatio( state.idU, splitIdx );
        double costThisSplit = -(mSplitWeight * log10( splitRatio ));
#endif

        if( state.timeSlice != mMaxTS && toCompute & T_EVENT ) 
            computeTransferCost(idUl, idUr, costThisSplit, state, 
                    optCost, bestSplit );
//cout << "  T: " << optCost << endl;
 
        if( toCompute & D_EVENT ) 
            computeDuplicationCost( idUl, idUr, costThisSplit,
                    state, optCost, bestSplit );
//cout << "  D: " << optCost << endl;

        vector< pair<int,int> > speciesSplits
            = mSpeciesTree->getSplits( state.idX );
        if( speciesSplits.size() > 1 ) 
            throw bpp::Exception( "DTLMatrixRecalc: ils not supported" );
        if( speciesSplits.size() == 1 && speciesSplits[0].second != -1
            && toCompute & S_EVENT ) 
            computeSpeciationCost( idUl, idUr, speciesSplits[0].first, 
                        speciesSplits[0].second,
                        costThisSplit, state, optCost, bestSplit );
//cout << "  S: " << optCost << endl;
    }
}

/**
 * Compute costs for the current split.
 */
void DTLMatrixRecalc::computeOptimaForSpeciesSplit(
        int toCompute,  ///< bit index with costs to compute
        pair<int,int> cladeSplit, ///< current split
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    vector< pair<int,int> > speciesSplits
        = mSpeciesTree->getSplits( state.idX );
    if( speciesSplits.size() > 1 ) 
        throw bpp::Exception( "DTLMatrixRecalc: ils not supported" );

    if( speciesSplits.size() == 1 && speciesSplits[0].second != -1 ) {
        if( toCompute & SL_EVENT ) 
            computeSpeciationPlusLossCost( speciesSplits[0].first,
                        speciesSplits[0].second,
                        state, optCost, bestSplit );
//cout << " SL:" << optCost << endl;
    } else if( speciesSplits.size() == 1 && speciesSplits[0].second == -1 ) {
        if( toCompute & N_EVENT ) 
            computeNullCost( speciesSplits[0].first,
                    state, optCost, bestSplit );
//cout << "  N: " << optCost << endl;
    } else if( cladeSplit.first == -1 ) { 
        // no sons 
        // if both are leaf, check if they are the "same"
        MySpeciesNode *xNode = mSpeciesTree->getNodeById( state.idX );
        if( mCladesTrips->mClades.getSpeciesName(state.idU) 
                == xNode->getName() ) {
#ifndef SPEED
            if( mSubOpt ) 
                state.costList.push_back( 0 );
            else if( mTriplets ) { 
                addZeroCost();
            }
#endif
            optCost = 0;
            bestSplit.event = 'f';
//cout << "  L: 0" << endl;
        } 
//else 
//cout << "  L: max" << endl;
    }
//cout << "  O: " << optCost << " event=" << bestSplit.event << endl;
}


/**
 * Recomputes some costs without TL.
 *
 * Recomputed if state.idX,state.idU is in event queue, recomputeAll
 * is true, or this is alpha and other nodes have changed.
 *
 * Uses recal event queue, matrix, and mBestReceivers/Cost
 * Modifies speciesNodeIdsTS, mCalculatedSpeciesIds.
 *
 * @return true if costs changed or this is the best receiver
 */
bool DTLMatrixRecalc::recomputeOptima( 
    bool &changed, ///< True if any values have changed.
    bool recomputeAll, ///< Force recomputation regardless of queue.
    double &cost,   ///< calculated cost
    DTLMatrixState &state ) ///< various current info
{
    bool recompute = false; 
    int events = ALL_EVENT; // default to recalc all events
    if( recomputeAll ) {
        recompute = true;
    } else if ( mSpeciesTree->hasAlpha() && mSpeciesTree->isAlpha( state.idX ) 
               && changed ) 
    {
        // add alpha if this node is alpha and other nodes have
        // changed (alpha cost before TL is needed for TL cost)
        recompute = true;
    } else {
        // search map for event 
 	    boost::unordered_map<pair<int,int>,int>::iterator 
            iterRecal = mRecalQ.find( pair<int,int>(state.idU,state.idX) );
        if( iterRecal != mRecalQ.end() ) {
            recompute = true;
            events = iterRecal->second;
        }
    }

    // recompute
    bool add = false;
    if( recompute ) {
        add = true;
        cost = computeOptimaForCell( state.idX, events, state );
        double oldCost = mMatrix.getValueSure(state.idU, state.idX);
        if( !COST_EQUAL(cost,oldCost) && events != ALL_EVENT ) 
        {
            // cost went up, check other costs for less
            double otherCost = computeOptimaForCell( 
                                    state.idX, ~events, state );
            if( otherCost < cost ) 
                cost = otherCost;
        }
        if( !COST_EQUAL(cost,oldCost) )
            changed = true;

    } else if( state.idX == mBestReceiver[state.timeSlice][state.idU] ) {
        add = true;
        // needed for checking if best receiver changed
        cost = mBestReceiverCost[state.timeSlice][state.idU];
    }

    return add;
}


/**
 * For the current nodes, calculate which parent nodes
 * need to change and add them to the recal queues.
 *
 * Changes are for the current node changing and the
 * best receiver changing.
 */
void DTLMatrixRecalc::updateEventQueue( 
    bool changed, ///< true if the current node changed
    vector<int> &parentsU, ///< parents of state.idU
    DTLMatrixState &state ) ///< various current info
{
    boost::unordered_map<pair<int,int>,int>::iterator iterRecal;

    // add p(x),NULL,SL to u,t+1 
    bool speciesChange = false;
    int parentX = mSpeciesTree->getArtificialFatherId( state.idX );
    if( changed && parentX != -1 ) {
        pair<int,int> pxPair(state.idU,parentX);
        iterRecal = mRecalNextQ.find( pxPair );
        int events = N_EVENT | SL_EVENT;
        if( iterRecal == mRecalNextQ.end() ) {
            mRecalNextQ.insert( make_pair(pxPair,events) ); 
        } else {
            iterRecal->second |= events;
        }
        speciesChange = true;
    }

    // If not root clade add x,D,T to p(u),t. 
    // only T if BR1 changed.
    bool cladeChange = false;
    if( state.idU != 0 ) {
        BOOST_FOREACH( int parentU, parentsU ) {
            pair<int,int> puPair(parentU,state.idX);
            iterRecal = mRecalQ.find( puPair );
            int events = T_EVENT;
            if( changed ) 
                events |= D_EVENT;
            if( iterRecal == mRecalQ.end() ) 
                mRecalQ.insert( make_pair(puPair,events) ); 
            else 
                iterRecal->second |= events;
            cladeChange = true;
        }
    }

    // add p(x),S to p(u),t+1 
    if( changed && speciesChange && cladeChange ) {
        BOOST_FOREACH( int parentU, parentsU ) {
            pair<int,int> puxPair(parentU,parentX);
            iterRecal = mRecalNextQ.find( puxPair );
            if( iterRecal == mRecalNextQ.end() ) 
                mRecalNextQ.insert( make_pair(puxPair,S_EVENT) ); 
            else
                iterRecal->second |= S_EVENT;
        }
    }
}

/**
 * Parse a variable line from the matrix file.
 */
void DTLMatrixRecalc::readMatrixVariable(
    string line ) ///< a line from the matrix file
{
    // read variable
    string::size_type index = line.find("=");
    if(index != string::npos) {
        string var = line.substr(1, index-1 );
        string val = line.substr(index + 1);
        try { 
            if( var == "fixedCosts" )  {
                if( val == "false" )
                    mFixedCosts = false;
            } else if( var == "computeTL" )  {
                if( val == "false" )
                    mComputeTL = false;
            } else if( var == "dupliCost" )  
                mDupliCost = bpp::TextTools::toDouble(val.c_str());
            else if( var == "hgtCost" )  
                mHGTCost = bpp::TextTools::toDouble(val.c_str());
            else if( var == "lossCost" ) 
                mLossCost = bpp::TextTools::toDouble(val.c_str());
            else if( var == "WGDCost" )  
                mWGDCost = bpp::TextTools::toDouble(val.c_str());
            else if( var == "weight" )  
                mSplitWeight = bpp::TextTools::toDouble(val.c_str());
            else if( var == "transferDead" ) { 
                if( val == "true" ) {
                    if( !mSpeciesTree->hasAlpha() )
                        throw bpp::Exception( "transferDead=true"
                            " in matrix file but alpha not added" );
                } else if( val == "false" ) {
                    if( mSpeciesTree->hasAlpha() )
                        throw bpp::Exception( "transferDead=false"
                            " in matrix file but alpha added" );
                } else {
                    throw bpp::Exception( "bad value for transferDead"
                            " in matrix file" );
                }
            }
        } catch( bpp::Exception e ) {
            throw bpp::Exception(
                "Number conversion error reading matrix parameters" );
        }
        //cout << "<" << var << "/" << val << ">" << endl;
    }
}


/**
 * Read line of matrix, best receivers or best splits.
 */
void DTLMatrixRecalc::readMatrixLine(
    string line, ///< a line from the matrix file
    int spNodeCnt,   ///< current species node count
    int &brTs,  ///< best receiver time slice
    int &state )    ///< reading state
{
    string::size_type index = -1;
    int cladeNum = -1;
    int prevIndex= -1;
    string val;
    while( true ) {
        index = line.find(",", index+1 );
        if( index == string::npos ) 
            val = line.substr(prevIndex+1); // last token 
        else 
            val = line.substr(prevIndex+1, index-prevIndex-1);

        if( cladeNum < 0 ) {
            // first column
            if( val == "bestReceivers" )  
                state = READING_BEST_REC;
            else if( val == "bestReceiversCost" ) {
                state = READING_BEST_REC_COST;
                brTs = -1;
            }
        } else {
            // insert values
            if( state == READING_MATRIX )  {
                try {
                    double matrixVal = bpp::TextTools::toDouble(val.c_str());
                    mMatrix.setValue(cladeNum, spNodeCnt, matrixVal ); 
                } catch( bpp::Exception e ) {
                    throw bpp::Exception("Number conversion error reading matrix");
                }
            } else if (state == READING_BEST_REC_COST) {
                size_t index2 = val.find(";");
                if( index2 == string::npos ) 
                    throw bpp::Exception(
                            "Best receivers or splits missing ;");
                double val1, val2;
                try {
                    val1 = bpp::TextTools::toDouble(
                            (val.substr(0,index2)).c_str());
                    val2 = bpp::TextTools::toDouble(
                            (val.substr(index2+1)).c_str());
                } catch( bpp::Exception e ) {
                    throw bpp::Exception( 
                            "Not a number in Best receivers cost");
                }
                mBestReceiverCost[brTs][cladeNum] = val1; 
                mSecondBestReceiverCost[brTs][cladeNum] = val2;
            } else {
                size_t index2 = val.find(";");
                if( index2 == string::npos ) 
                    throw bpp::Exception( "Best receivers or splits missing ;");
                int val1, val2;
                try {
                    val1 = bpp::TextTools::toInt( 
                            (val.substr(0,index2)).c_str() );
                    val2 = bpp::TextTools::toInt( 
                            (val.substr(index2+1)).c_str() );
                } catch( bpp::Exception e ) {
                    throw bpp::Exception(
                            "Not a number in Best receivers or splits");
                }

                // set values
                if( state == READING_BEST_REC ) {
                    mBestReceiver[brTs][cladeNum] = val1; 
                    mSecondBestReceiver[brTs][cladeNum] = val2;
                } 
            } 
        }
        if( index == string::npos ) 
            break;
        cladeNum++;
        prevIndex = index;
    }
}


/**
 * Constructor to initialize variables from existing matrix.
 *
 * Reads matrix file.
 */
DTLMatrixRecalc::DTLMatrixRecalc( 
    MySpeciesTree *speciesTree, ///< species trees
    CladesAndTripartitions *cat, ///< clades and tripartitions
    const char *matrixFileName, 
        ///< file containing matrix and data from previous run
    int maxTS ) ///< maximum time slice 
{

    mMaxTS = maxTS;
    mSpeciesTree = speciesTree;
    mCladesTrips = cat;
	mSTnodes = (mSpeciesTree->getNodes()).size();

    ifstream csvfile( matrixFileName );
    if (!csvfile.is_open()) {
        cout << "Unable to open matrix file " << matrixFileName << endl;
        throw bpp::Exception("Bad matrix file");
    }

    mFixedCosts = true;
    mComputeTL = true;
    mDupliCost = -std::numeric_limits<double>::max();
    mHGTCost = -std::numeric_limits<double>::max();
    mLossCost = -std::numeric_limits<double>::max();
    mWGDCost = -std::numeric_limits<double>::max();
    mSplitWeight = -std::numeric_limits<double>::max();
    int readClades = 0;
    int spNodeCnt = 0;
    int brTs = -1; // best receiver time slice
    int state = READING_MATRIX;
    string line;
    while( getline( csvfile, line ) ) {
        if( line.empty() ) 
            continue;

        if( line[0] == '#' ) {
            readMatrixVariable( line );
        } else if( !readClades ) {
            // finished reading variables; initialize structures
	        mSTnodes = (mSpeciesTree->getNodes()).size();
            init();

            // check clade count
            string::size_type index = -1;
            int cnt=0;
            while( (index = line.find(",",index+1)) != string::npos ) 
                cnt++;
            if( cnt != mCladeCount ) {
                cerr << "Clade count from geneTreeFile is " << mCladeCount
                    << " but " << cnt << " in matrix file";
                throw bpp::Exception("Clade count in matrix file doesn't match");
            }
            readClades++;
        } else {
            readMatrixLine( line, spNodeCnt, brTs, state );
            if( state == READING_MATRIX ) 
                spNodeCnt++;
            if( state == READING_BEST_REC )  
                brTs++;
            if( state == READING_BEST_REC_COST )  
                brTs++;
        }
    }
    csvfile.close();

    // check the values
    if( mDupliCost == -std::numeric_limits<double>::max() ) 
        throw bpp::Exception("dupliCost parameter not found in matrix file");
    if( mLossCost == -std::numeric_limits<double>::max() ) 
        throw bpp::Exception("lossCost parameter not found in matrix file");
    if( mHGTCost == -std::numeric_limits<double>::max() ) 
        throw bpp::Exception("hgtCost parameter not found in matrix file");
    if( mWGDCost == -std::numeric_limits<double>::max() ) 
        throw bpp::Exception("WGDCost parameter not found in matrix file");
    if( mSplitWeight == -std::numeric_limits<double>::max() ) 
        throw bpp::Exception("weight parameter not found in matrix file");
    if( spNodeCnt != mSTnodes ) {
        cerr << "Species node count from speciesTreeFile is " << mSTnodes
                    << " but " << spNodeCnt << " in matrix file" << endl;
        throw bpp::Exception("Species node count in matrix file doesn't match");
    }
    if( brTs != mMaxTS+1 ) {
        cerr << "max time slice is " << mMaxTS 
             << ", but read " << brTs << " time slices for best receivers"
             << " in the matrix file" << endl;
        throw bpp::Exception("Time slice count in matrix file doesn't match");
    }

    mMatrixCalculated = true;
}


/**
 * Read parameters from a matrix file.
 *
 * static function to read specific parameters
 * from a matrix file, though reduced to just
 * transferDead at the moment
 *
 */
void DTLMatrixRecalc::getCSVparams( 
   const char *fileName,    ///< matrix file
   bool &transferDead ) ///< parameter to find
{
    ifstream csvfile( fileName );
    if (!csvfile.is_open()) {
        cout << "Unable to open matrix file " << fileName << endl;
        throw bpp::Exception("Bad matrix file");
    }
    string line;
    while( getline( csvfile, line ) ) {
        if( line.empty() ) 
            continue;
        if( line[0] != '#' ) // no more parameters
            break;

        // read variable
        string::size_type index = line.find("=");
        if(index != string::npos) {
            string var = line.substr(1, index-1);
            string val = line.substr(index+1);
            if( var == "transferDead" ) {
                if( val == "true" ) 
                    transferDead = true;
                else if( val == "false" )  
                    transferDead = false;
                else
                    throw bpp::Exception(
                        "getCSVparams, bad value for transferDead" );
            }
         }
    }
    csvfile.close();
}


/**
 * Reompute optimal costs for all nodes in the time slice, excluding
 * the transfer loss cost.
 *
 * @return True if any value changed.
 */
bool DTLMatrixRecalc::recomputeOptimaForTS( 
    vector<int> &speciesNodeIdsTS, ///< all specie nodes in the time slice
    vector<double> &optVectorTS,    ///< calculated costs for the nodes
    boost::unordered_map<int,double> &calculatedSpeciesIds,
        ///< map of nodes already calculated
    double &alphaCost, ///< cost of alpha for this time slice (if it exists)
    bool recomputeAll, ///< recompute everything if true
    DTLMatrixState &state ) ///< various current info
{
    // calculate suboptimal for those in queue
    vector<int> allSpeciesNodeIdsTS 
        = mSpeciesTree->getVectorWithTS(state.timeSlice);

    bool changed = false;
    for ( size_t i=0; i<allSpeciesNodeIdsTS.size(); i++ ) {
        // state variables in class
        state.idX = allSpeciesNodeIdsTS[i];	
        if( !mFixedCosts ) {
            MySpeciesNode *xNode = mSpeciesTree->getNodeById( state.idX );
            state.dupliCost = xNode->getInfos().duplicationCost;
            state.hgtCost = xNode->getInfos().hgtCost;
            state.lossCost = xNode->getInfos().lossCost;
        }

        double cost;
        if( recomputeOptima( changed, recomputeAll, cost, state ) ) {
            speciesNodeIdsTS.push_back( state.idX );
            optVectorTS.push_back( cost ); 

            calculatedSpeciesIds.insert( make_pair(state.idX, cost) );
            if ( mSpeciesTree->hasAlpha() && mSpeciesTree->isAlpha( state.idX ))
                alphaCost = cost;
        }
    } // end suboptimal loop

    return changed;
}



/**
 * Fill optVector with recalculated costs using those
 * costs already calculated where available.
 */
void DTLMatrixRecalc::redoOptVector( 
    boost::unordered_map<int,double> &calculatedSpeciesIds,
        ///< map of recalculated species ids and their costs
    vector<double> &optVectorTS,  ///< calculated cost for speciesNodeIdsTS
    DTLMatrixState &state ) ///< various current info
{
    boost::unordered_map<int,double>::iterator seenIt; 
    seenIt = calculatedSpeciesIds.find( state.idX );
    if( seenIt == calculatedSpeciesIds.end() ) {
        double optAllOverTriples = 
            computeOptimaForCell( state.idX, ALL_EVENT, state );
        optVectorTS.push_back(optAllOverTriples); 
    } else {
        double tmpOpt = seenIt->second; // already calculated
        // compute events not already calculated
 	    boost::unordered_map<pair<int,int>,int>::iterator 
            iterRecal = mRecalQ.find( pair<int,int>(state.idU,state.idX) );
        int events;
        if( iterRecal != mRecalQ.end() ) {
            events = iterRecal->second;
            double optAllOverTriples = computeOptimaForCell( 
                                           state.idX, ~events, state );
            if( tmpOpt > optAllOverTriples )
                tmpOpt = optAllOverTriples;
        } // else alpha was recalculated 
        // or bestReceivers are old versions
        optVectorTS.push_back( tmpOpt );
    }
}


/**
 * Recalculate values if best receivers have changed.
 */
void DTLMatrixRecalc::recomputeBestReceivers( 
    int brChanged, ///< -1 if new best receiver is less, 1 if more
    vector<int> &speciesNodeIdsTS, ///< all specie nodes in the time slice
    vector<double> &optVectorTS,  ///< calculated cost for speciesNodeIdsTS
    boost::unordered_map<int,double> &calculatedSpeciesIds,
        ///< map of recalculated species ids and their costs
    DTLMatrixState &state ) ///< various current info
{
    optVectorTS.clear();
    speciesNodeIdsTS = mSpeciesTree->getVectorWithTS( state.timeSlice );
    for ( size_t i=0; i<speciesNodeIdsTS.size(); i++) {
        state.idX = speciesNodeIdsTS[i];
        if( brChanged == -1 ) {
            // BR is less, fill with old values
            boost::unordered_map<int,double>::iterator seenIt; 
            seenIt = calculatedSpeciesIds.find( state.idX );
            if( seenIt == calculatedSpeciesIds.end() ) 
                optVectorTS.push_back(
                        mMatrix.getValueSure(state.idU,state.idX));
            else 
                optVectorTS.push_back( seenIt->second ); 
        } else { 
            // BR is more, so all costs must be redone
            redoOptVector( calculatedSpeciesIds, optVectorTS, state );
        }
    }
    // update bestReceiver
    computeBestReceivers( state.timeSlice, state.idU, speciesNodeIdsTS, 
                          optVectorTS );
}


/**
 * Realculate and set a matrix cell values for the given clade
 * and nodes in the given time slice.
 *
 */
void DTLMatrixRecalc::recalculateMatrixTS(
    int idU,    ///< clade id
    int timeSlice, ///< current time slice
    vector<int> &idUparents, ///< parents of idU
    bool recomputeAll ) ///< recompute everything if true
{
    DTLMatrixState state;
    state.idU = idU;
    state.timeSlice = timeSlice;
    if( mFixedCosts ) {
        state.dupliCost = mDupliCost;
        state.hgtCost = mHGTCost;
        state.lossCost = mLossCost;
    }

    vector<int> speciesNodeIdsTS; // list of species nodes
    vector <double> optVectorTS; // costs for speciesNodeIdsTS
    boost::unordered_map<int,double> calculatedSpeciesIds; 
        // for saving calculated nodes
    double alphaCost = -1;
    bool changed = recomputeOptimaForTS( speciesNodeIdsTS, optVectorTS,
                calculatedSpeciesIds, alphaCost, recomputeAll, state );
    if( !changed && !recomputeAll ) // recomputeAll can have cost changes
        return; // nothing more to recalculate

    // Recalculate best receivers.  -1 if new is less, 1 if more
    int brChanged = computeBestReceivers( timeSlice, idU, speciesNodeIdsTS, 
                                          optVectorTS );

    // Recalculate other nodes if BR changes.
    if( !recomputeAll && brChanged != 0 )  
        recomputeBestReceivers( brChanged, speciesNodeIdsTS, optVectorTS, 
                                  calculatedSpeciesIds, state );

    // check if alpha has changed
    bool alphaUpdated = false;
    if( mSpeciesTree->hasAlpha() && speciesNodeIdsTS.size() > 0 
        && timeSlice != mMaxTS ) 
    {
        if( alphaCost == -1 ) 
            throw bpp::Exception("DTLMatrix::recalculateMatrixTS: "
                    " alphaCost not set" );
        int alphaId =  mSpeciesTree->getAlphaIdForTS( timeSlice );
        if( !COST_EQUAL( alphaCost, mMatrix.getValueSure(idU,alphaId) ) )
            alphaUpdated = true;
    }

    // Set matrix values. Calculate TL costs if necessary.
    for( size_t i=0; i<speciesNodeIdsTS.size(); i++ ) {
        // costs without considering transfer loss events
        double opt = optVectorTS[i]; 
        state.idX = speciesNodeIdsTS[i];	
        if( !mFixedCosts ) {
            MySpeciesNode *xNode = mSpeciesTree->getNodeById( state.idX );
            state.dupliCost = xNode->getInfos().duplicationCost;
            state.hgtCost = xNode->getInfos().hgtCost;
            state.lossCost = xNode->getInfos().lossCost;
        }

        bool changed = false;
        if( !COST_EQUAL( opt, mMatrix.getValueSure(idU,state.idX) ) )
            changed = true;
        if( brChanged != 0 || changed || alphaUpdated || recomputeAll ) {
            if( mComputeTL && opt != 0 && speciesNodeIdsTS.size() > 1 ) 
                if( computeTransferLossCost( alphaCost, opt, state) )
                    changed = true;
            mMatrix.setValue(idU,state.idX,opt); 
            updateEventQueue( changed, idUparents, state );
        }
    } 
}


/**
 * Recalculate the matrix for the parameters given in the the constructor. 
 *
 * All time slices between the minimum and maximum in tsChangeList
 * are recalculated. Matrix cells in time slices greater than the
 * maximum (up to the root) are also recalculated if their
 * dependancies have changed.
 *
 */
void DTLMatrixRecalc::recalculateMatrix( 
        vector<int> tsChangeList, ///< List of time slices that changed.
        bool verbose)  ///< print timing and stats
{
    clock_t start = clock();

    if( tsChangeList.size() == 0 )  
        throw bpp::Exception( "DTLMatrix::recalculateMatrix: "
                 " input has no changed time slice" );

    mRecal = true;
    sort( tsChangeList.begin(), tsChangeList.end() );
    int minChangedTS = tsChangeList[0];

    // get sorted clade ids (sorted so that smaller clades are seen first)
    vector<int> sortedCladeNums = mCladesTrips->mClades.getSortedClades();
    // create parent list for clades, used to find cells needing recalculation
    vector< vector<int> > cladeParents = mCladesTrips->getCladeParents();

    // main loop over time slices
    size_t changeIdx = 0;
    for ( int timeSlice=minChangedTS; timeSlice<=mMaxTS; timeSlice++ ) {

        // move to recal list computed in the last time slice
        mRecalQ.clear(); 
        mRecalQ = mRecalNextQ;
        mRecalNextQ.clear();

        bool recomputeAll = false;
        if( changeIdx < tsChangeList.size() 
            && timeSlice == tsChangeList[changeIdx] ) 
        {
            recomputeAll = true;
            changeIdx++;
        } else if( changeIdx == tsChangeList.size() && mRecalQ.size() == 0 ) {
            break; //done
        }

        // clade loop
        BOOST_FOREACH( int idU, sortedCladeNums ) 
            recalculateMatrixTS( idU, timeSlice, cladeParents[idU],
                                 recomputeAll );
    }
		

    if( verbose ) {
        double time = (double) (clock()-start) / CLOCKS_PER_SEC * 1000.0;
        cout << "recal matrix time = " << time << endl;
    }
}
