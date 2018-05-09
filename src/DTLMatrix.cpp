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

**/


#include <iostream>
#include <limits>

#include <Bpp/Text/TextTools.h>
#include <boost/foreach.hpp>

#include "DTLMatrix.h"

#ifdef SPEED
    #define UPDATE_OPT(optCost, bestSplit, cost, optimumSub, idX, c) if( cost < optCost) optCost = cost;
#else
    #define UPDATE_OPT(optCost, bestSplit, cost, optimumSub, idX, c); updateOpt(optCost, bestSplit, cost, optimumSub, idX, c);
#endif


const double DTLMatrix::COST_DIFF = 0.00001;

/**
 * Cost comparison accounting for rounding.
 *
 * @return true if a is equal to b
 */
bool DTLMatrix::COST_EQUAL( double a, double b ) {
    if( abs(a-b) <= COST_DIFF )
        return true;
    return false;
}
/**
 * Cost comparison accounting for rounding.
 *
 * @return true if a is greater than b
 */
bool DTLMatrix::COST_GREATER( double a, double b ) {
    if( a > (b+COST_DIFF) ) 
        return true;
    return false;
}
        


/**
 * Given a pair of clades and species, a and b, calculate costs
 * from the matrix.
 *
 * If mSubOpt is true
 * find costsA and costsB, from the vector matrix and make all combinations
 * of additions, add otherCost, and push onto costList.
 *
 * @return lowest cost
 */
double DTLMatrix::computeCostList( 
        double otherCost, ///< additional costs to add
        int idUa, ///< clade a
        int idXa, ///< species a
        int idUb, ///< clade b 
        int idXb, ///< species b
        EventTriplet eventTrip, ///< events to add
        DTLMatrixState &state ) ///< various current info
{
    double lowestCost = otherCost
                     + mMatrix.getValueSure(idUa,idXa)
                     + mMatrix.getValueSure(idUb,idXb);
#ifndef SPEED
    if( mSubOpt ) { 
        // add other costs to costList
        vector<double> costsA = mMatrixV->getValueSure(idUa,idXa); 
        vector<double> costsB = mMatrixV->getValueSure(idUb,idXb); 
        BOOST_FOREACH( double costA, costsA ) {
            BOOST_FOREACH( double costB, costsB ) {
                double cost = costA + costB + otherCost; 
                if( COST_GREATER( cost, state.curBestEpsCost ) )
                    continue;
 
                state.costList.push_back( cost );
               
                if( COST_GREATER( lowestCost, cost ) ) {
                    cout << "lowestCost=" << lowestCost << endl;
                    cout << "cost=" << cost << "=" << costA 
                        << "+" << costB << "+" << otherCost << endl; 
                    throw bpp::Exception ("DTLMatrix::computeCostList:"
                            " impossible cost");
                }
            }
        }
    } else if( mTriplets ) {
        computeList( idUa, idXa, idUb, idXb, 
                lowestCost, state.curBestEpsCost, otherCost, eventTrip );
    }
#endif

    return lowestCost;
}

/**
 * Given a clade and specie, calculate costs from the matrix.
 *
 * If mSubOpt is true
 * find costsA, from the vector matrix and make all combinations
 * of additions, add otherCost, and push onto costList.
 *
 * @return lowest cost
 */
double DTLMatrix::computeCostList( 
        double otherCost, ///< additional costs to add
        int idUa, ///< clade
        int idXa, ///< species
        EventTriplet eventTrip, ///< events to add
        DTLMatrixState &state ) ///< various current info
{
    double lowestCost = otherCost + mMatrix.getValueSure(idUa,idXa);
#ifndef SPEED
    if( mSubOpt ) {
        // calculate all costs
        vector<double> costsA = mMatrixV->getValueSure(idUa,idXa); 
        BOOST_FOREACH( double costA, costsA ) {
            double cost = costA + otherCost; 
            if( COST_GREATER( cost, state.curBestEpsCost ) )
                continue;

            state.costList.push_back( cost );
            
            if( COST_GREATER( lowestCost, cost ) ) {
                cout << "lowestCost = " << lowestCost << endl;
                cout << "cost=" << cost << "=" << costA 
                        << "+" << otherCost << endl; 
                throw bpp::Exception ("computeCostList: impossible cost");
            }
        }
    } else if( mTriplets ) {
        computeList( idUa, idXa, 
                lowestCost, state.curBestEpsCost, otherCost, eventTrip );
    }
#endif

    return lowestCost;
}

/**
 * Fill the best split structure if the given cost is less than
 * optCost.
 */
void DTLMatrix::updateOpt( 
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit,  ///< current best split
        double cost, ///< cost to check
        int firstSplit, ///< species split
        int secondSplit, ///< species split
        char event ) ///< event associated with the given cost
{
    if( cost < optCost) 
    {
        optCost = cost;
        if( mUseBestSplits ) {
            bestSplit.speciesSplit.first = firstSplit;
            bestSplit.speciesSplit.second = secondSplit;
            bestSplit.event = event;
        }
    }
}



/**
 * Compute transfer cost for the current split.
 */
void DTLMatrix::computeTransferCost( 
        int idUl,     ///< clade child
        int idUr,     ///< clade child
        double costThisSplit,  ///< current split cost
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    int idUsub = idUl; 
    int idUother = idUr;
    EventTriplet eventTrip;
    for( int i=0; i<2; i++ ) { // Calculate transfer in each direction 
        // optimumSub can be zero for second best receivers with 
        // no optimal cost (in leaves), skip it
        int optimumSub = mBestReceiver[state.timeSlice][idUsub];  
        if (optimumSub == -1)
            continue;
#ifndef SPEED
        if( mSubOpt ) { // for all VT, do a cost combination
            if( optimumSub == state.idX ) {
                BOOST_FOREACH( double c, mVT2[state.timeSlice][idUsub] ) 
                    computeCostList( costThisSplit + state.hgtCost + c,
                                     idUother, state.idX, eventTrip, state ); 
            }
			// transfer all optimal to this node (except self)
            pair<double,int> p;
            BOOST_FOREACH( p, mVT1[state.timeSlice][idUsub]) {
                // don't add if it is this species node's cost
                if( p.second != state.idX ) {
                    computeCostList( costThisSplit + state.hgtCost + p.first, 
                                     idUother, state.idX, eventTrip, state ); 
                }
            }
        } else if( mTriplets ) { // for all VT, do a cost combination
            computeTransferCostSub( optimumSub, idUsub, idUother, 
                                 costThisSplit, state );
        }
#endif
		// if the bestReceiver is x, we take the second one
		if( optimumSub == state.idX ) 
		    optimumSub = mSecondBestReceiver[state.timeSlice][idUsub]; 
	
		// compute the cost, even for mSubOpt
		if (optimumSub != -1) {
		    // if state.idX is in Alpha or not it stays the same
		    double cost = costThisSplit 
						+ mMatrix.getValueSure(idUother,state.idX) 
		        		+ mMatrix.getValueSure(idUsub,optimumSub) 
                        + state.hgtCost; 
			if( i == 0 ) {
               UPDATE_OPT(optCost, bestSplit, cost, optimumSub, state.idX, 't');
            } else {
               UPDATE_OPT(optCost, bestSplit, cost, state.idX, optimumSub, 't');
            }
		}

        // switch to other child for second iteration
        idUsub = idUr;
        idUother = idUl;
    }

    if( state.timeSlice != mMaxTS && mSpeciesTree->hasAlpha() 
        && !mSpeciesTree->isAlpha( state.idX ) ) 
    {
        int alphaForThisTS = mSpeciesTree->getAlphaIdForTS(state.timeSlice);

        //Transfers to alpha cost 0
        //c4	Transfer to alpha
        double cost = computeCostList( costThisSplit, idUl, state.idX, 
                                idUr, alphaForThisTS, eventTrip, state ); 
        UPDATE_OPT( optCost, bestSplit, cost, state.idX, alphaForThisTS, 'a' );

        cost = computeCostList( costThisSplit, idUr, state.idX, 
                                idUl, alphaForThisTS, eventTrip, state ); 
        UPDATE_OPT( optCost, bestSplit, cost, alphaForThisTS, state.idX, 'a' );
    }
}



/**
 * Compute duplication cost for the current split.
 */
void DTLMatrix::computeDuplicationCost( 
        int idUl,     ///< clade child
        int idUr,     ///< clade child
        double costThisSplit,  ///< current split cost
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    double otherCosts = costThisSplit;
    EventTriplet eventTrip;

//    if( state.xNode->getInfos().WGD ) 
//        otherCosts += mWGDCost;	// WG Duplication
//    else 
    if( !mSpeciesTree->hasAlpha() || !mSpeciesTree->isAlpha( state.idX )) {
        otherCosts += state.dupliCost;	
        eventTrip.d = 1;
    }
    //else duplication in alpha costs 0
    
    double cost = computeCostList( otherCosts, idUl, state.idX, 
                                   idUr, state.idX, eventTrip, state );
    if( mSpeciesTree->hasAlpha() && mSpeciesTree->isAlpha( state.idX ) ) {
        UPDATE_OPT( optCost, bestSplit, cost, state.idX, state.idX, 'a' );
    } else {
        UPDATE_OPT( optCost, bestSplit, cost, state.idX, state.idX, 'd' );
    }
}



/**
 * Compute the speciation cost for the current split.
 */
void DTLMatrix::computeSpeciationCost( 
        int idUl,     ///< clade child
        int idUr,     ///< clade child
        int idXl,     ///< species child
        int idXr,     ///< species child
        double costThisSplit,  ///< current split cost
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    EventTriplet eventTrip;

    // Cospeciation impossible in Alpha, c2 Cospeciation
    double cost = computeCostList( costThisSplit,
                    idUl, idXl, idUr, idXr, eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXl, idXr, 's' );

    cost = computeCostList( costThisSplit,
                    idUl ,idXr, idUr, idXl, eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXr, idXl, 's' );
}


/**
 * Compute the ils cost for the current split.
 */
void DTLMatrix::computeIlsCost( 
        int idUl,     ///< clade child
        int idUr,     ///< clade child
        int idXl,     ///< species child
        int idXr,     ///< species child
        double costThisSplit,  ///< current split cost
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    EventTriplet eventTrip;

    double otherCosts = costThisSplit + mIlsCost;

    // Cospeciation impossible in Alpha, c2 Cospeciation
    double cost = computeCostList( otherCosts,
                    idUl, idXl, idUr, idXr, eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXl, idXr, 'i' );

    cost = computeCostList( otherCosts,
                    idUl ,idXr, idUr, idXl, eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXr, idXl, 'i' );
}
 

/**
 * Calculate speciation plus loss cost.
 */ 
void DTLMatrix::computeSpeciationPlusLossCost( 
        int idXl,     ///< species child
        int idXr,     ///< species child
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{ 
    EventTriplet eventTrip;
    eventTrip.l = 1;

    double otherCosts = state.lossCost;

    // Cospeciation  + Loss impossible in Alpha
    //c3	Cospeciation + Loss
    double cost = computeCostList( otherCosts, state.idU, idXl, 
                                    eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXl, -1, 'l' );

    cost = computeCostList( otherCosts, state.idU, idXr, eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXr, -1, 'l' );
}


/**
 * Calculate ILS plus loss cost.
 */ 
void DTLMatrix::computeIlsPlusLossCost( 
        int idXl,     ///< species child
        int idXr,     ///< species child
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{ 
    EventTriplet eventTrip;
    eventTrip.l = 1;

    double otherCosts = state.lossCost + mIlsCost;

    // Cospeciation  + Loss impossible in Alpha
    //c3	Cospeciation + Loss
    double cost = computeCostList( otherCosts, state.idU, idXl, 
                                    eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXl, -1, 'j' );
/*
if( cost < 1000 ) {
cout << state.idU << "," << state.idX << " from " << state.idU << "," << idXl << " IL 1 cost=" << cost << ":";
BOOST_FOREACH( double d, state.costList )
    cout << " " << d;
cout << endl;
}
*/

    cost = computeCostList( otherCosts, state.idU, idXr, eventTrip, state );
    UPDATE_OPT( optCost, bestSplit, cost, idXr, -1, 'j' );
/*
if( cost < 1000 ) {
cout << state.idU << "," << state.idX << " from " << state.idU << "," << idXr << " IL 2 cost=" << cost << ":";
BOOST_FOREACH( double d, state.costList )
    cout << " " << d;
cout << endl;
}*/
}


/**
 * Calculate null cost for current matrix cell.
 */ 
void DTLMatrix::computeNullCost( 
        int idXl,     ///< species child
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    EventTriplet eventTrip;

    double otherCosts = 0;
//    if( state.xNode->getInfos().WGD ) {
        //c3 No Event in the case of a WGD, we pay for a loss 
        //and mWGDCost because the duplication happened
//        otherCosts = mWGDCost + state.lossCost;		
//        eventTrip.l = 1;
//    }
    //else //c3 without loss i.e. No Event
    double cost = computeCostList( otherCosts, state.idU, idXl, 
                                    eventTrip, state );		

//    if( state.xNode->getInfos().WGD ) {
//        UPDATE_OPT( optCost, bestSplit, cost, idXl, -1, 'w' );
//    } else {
        UPDATE_OPT( optCost, bestSplit, cost, idXl, -1, 'n' );
//    }
}


//#define PRINT_COSTS 1;
/**
 * Compute costs for the current split.
 */
void DTLMatrix::computeOptimaForCladeSplit(
        int toCompute,  ///< bit index with costs to compute (not used
                        ///<  in this class)
        int splitIdx,   ///< current split index
        pair<int,int> cladeSplit, ///< current split
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    int idUl = cladeSplit.first;
    int idUr = cladeSplit.second;

#ifdef SPEED
    double costThisSplit = 0;
#else
    double splitRatio = mCladesTrips->getSplitRatio( state.idU, splitIdx );
    double costThisSplit = -(mSplitWeight * log10( splitRatio ));
#endif
    // maxTS has one node (therefore no transfer) normally, but with
    // ILS, there can be multiple idX at maxTS
    if( mComputeT && ( state.timeSlice != mMaxTS || mUseILS )) 
        computeTransferCost(idUl, idUr, costThisSplit, state, 
                optCost, bestSplit );

#ifdef PRINT_COSTS
cout << "  T: " << optCost << endl;
#endif

    computeDuplicationCost( idUl, idUr, costThisSplit,
            state, optCost, bestSplit );
#ifdef PRINT_COSTS
cout << "  D: " << optCost << endl;
#endif

    // speciations
    vector< pair<int,int> > speciesSplits 
                                = mSpeciesTree->getSplits( state.idX );
    for( size_t i=0; i<speciesSplits.size(); i++ ) {
        if( speciesSplits[i].second == -1 ) // null, ignore 
            continue;

        computeSpeciationCost( idUl, idUr, speciesSplits[i].first,
            speciesSplits[i].second, costThisSplit, 
            state, optCost, bestSplit );
#ifdef PRINT_COSTS
cout << "  S: " << optCost << endl;
#endif
    }

    if( mUseILS ) {
        // ILS 
        vector< pair<int,int> > ilsSplits 
                                = mSpeciesTree->getIlsSplits( state.idX );
        for( size_t i=0; i<ilsSplits.size(); i++ ) {
// Is this possible
            if( ilsSplits[i].second == -1 )
                continue;
            computeIlsCost( idUl, idUr, ilsSplits[i].first, ilsSplits[i].second, 
                            costThisSplit, state, optCost, bestSplit );
#ifdef PRINT_COSTS
cout << "  I: " << optCost << endl;
#endif
        }
    }
}


/**
 * Compute costs for the current split.
 */
void DTLMatrix::computeOptimaForSpeciesSplit(
        int toCompute,  ///< bit index with costs to compute (not used
                        ///<  in this class)
        pair<int,int> cladeSplit, ///< current split
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    vector< pair<int,int> > speciesSplits 
                                = mSpeciesTree->getSplits( state.idX );
    for( size_t i=0; i<speciesSplits.size(); i++ ) {
        if( speciesSplits[i].second == -1 ) { 
            // null 
            computeNullCost( speciesSplits[i].first, state, optCost, bestSplit);
#ifdef PRINT_COSTS
cout << "  N: " << optCost << endl;
#endif
        } else {
            computeSpeciationPlusLossCost( speciesSplits[i].first, 
                        speciesSplits[i].second, state, optCost, bestSplit );
#ifdef PRINT_COSTS
cout << "  SL: " << optCost << endl;
#endif
        }
    }

    // ILS loss
    if( mUseILS ) {
        vector< pair<int,int> > ilsSplits 
                                = mSpeciesTree->getIlsSplits( state.idX );
        for( size_t i=0; i<ilsSplits.size(); i++ ) {
            if( ilsSplits[i].second == -1 ) {
                computeNullCost( ilsSplits[i].first, state, optCost, bestSplit);
#ifdef PRINT_COSTS
cout << "  IN: " << optCost << endl;
#endif
            } else {
                computeIlsPlusLossCost( ilsSplits[i].first, ilsSplits[i].second,
                                        state, optCost, bestSplit );
#ifdef PRINT_COSTS
cout << "  IL: " << optCost << endl;
#endif
            }
        }
    }

    if( speciesSplits.size() == 0 && cladeSplit.first == -1 ) { 
        // both are leaves
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
#ifdef PRINT_COSTS
cout << "  L: 0" << endl;
#endif
        } 
#ifdef PRINT_COSTS
else 
cout << "  L: max" << endl;
#endif
    }
}




/**
 * Compute the temp optima for the matrix cell c(state.idX,state.idU) 
 * for the given species node, which is the optimal costs 
 * without considering the transfer loss cost.
 *
 * toCompute is used to only calculate some costs for the
 * recalcaulteMatrix function
 *
 * If best splits are needed, the split and event associated
 * with the best cost is saved.
 *
 * @return optimal cost
 */
double DTLMatrix::computeOptimaForCell(
        int idX, ///< species node id
        int toCompute,  ///< bit index with costs to compute
        DTLMatrixState &state ) ///< various current info
{
    // state variables in class
    state.idX = idX;
    //cout << state.idX <<" " << state.idU << endl;	
    if( mSubOpt ) 
        state.costList.clear();
    if( !mFixedCosts ) {
        MySpeciesNode *xNode = mSpeciesTree->getNodeById( idX );
        state.dupliCost = xNode->getInfos().duplicationCost;
        state.hgtCost = xNode->getInfos().hgtCost;
        state.lossCost = xNode->getInfos().lossCost;
    }

    double optCost = std::numeric_limits<double>::max();
    BestSplit bestSplit;
    bestSplit.event = 'u'; // non-existent event
    double lastCost = optCost;

    // Loop over all splits of the current clade from
    // cladeSplits. splitsRatio is a parallel array with cladeSplits
    int splitCount = mCladesTrips->getSplitCount( state.idU );
    if( splitCount < 1 ) 
        throw bpp::Exception("DTLMatrix::computeOptimaForCell:"
                    " no splits for a node" );

    pair<int,int> cladeSplit;
    for( int splitIdx=0; splitIdx<splitCount; splitIdx++ ) {
        cladeSplit = mCladesTrips->getCladeSplit( state.idU, splitIdx );

        if( cladeSplit.first != -1 ){ // if not a leaf
            //cout << "computeOptimaForCladeSplit\n";
            computeOptimaForCladeSplit( toCompute, splitIdx, cladeSplit, state,
                                   optCost, bestSplit );
         	                          
        }                           

        if( splitIdx == 0 ){ // only call once
            //cout << "computeOptimaForSpeciesSplit\n";
            computeOptimaForSpeciesSplit( toCompute, cladeSplit, state,
                               optCost, bestSplit );
         	                      
        }                       

        // if cost was updated, save this clade split
        if( optCost < lastCost ) 
            bestSplit.cladeSplit = cladeSplit;
        lastCost = optCost;
    
    }

    if( mUseBestSplits ) 
        mBestSplits[state.idU][state.idX] = bestSplit;

    if( mSubOpt ) {
        // sort and remove duplicates 
        sort( state.costList.begin(), state.costList.end() );
        size_t prev = 0;
        int eraseCount = 0;
        size_t listSize = state.costList.size();
        for( size_t i=1; i<listSize; i++ ) {
            if( !COST_EQUAL( state.costList[i], state.costList[prev] )) {
                prev++;
                if( prev != i )
                    state.costList[prev] = state.costList[i];
            } else
                eraseCount++;
        }
        if( eraseCount > 0 )
            state.costList.erase( state.costList.end()-eraseCount,
                                  state.costList.end() );
    } else if( mTriplets ) {
        pushList( state );
    }

    return optCost;
} 



/**
 * Save the given best receivers.
 *
 * @return -1 if new best receiver is less, 1 for greater, or 0 for no change
 *          (used by recalc code)
 */ 
int DTLMatrix::saveBestReceivers( 
    int timeSlice,      ///< time slice to save
    int idU,            ///< gene id
    int bestReceiver,   ///< best receiver
    int secondBestReceiver, ///< second best receiver
    double bestCost,    ///< cost of best receiver
    double secondBestCost,  ///< cost of best receiver
    vector<int> &otherBests ) ///< other second best receivers
{

    // All other best receivers are needed for graph calculation, but
    // not the matix calculation. Therefore, just add all others
    // to a giant list and put it into a more usable format 
    // as needed later in the getAllBestReceivers() function.
#ifndef PARALLEL
#ifndef SPEED 
    // not thread safe
    BOOST_FOREACH( int idX, otherBests ) {
        if( idX != -1 ) {
            mAllOtherBestReceivers.push_back( timeSlice );
            mAllOtherBestReceivers.push_back( idU );
            mAllOtherBestReceivers.push_back( idX );
        }
    }
#endif
#endif

    // save the result
    // reading and writing can change the rounding, so COST_DIFF 
    // is needed to check for equality when  doing recalculations
    // Note, the actual bestReceiver can change as long as 
    // the cost remains the same.
    // -1 values occur for the root and leaves. These are not used.
    int bestRecChanged = 0; 
    if( bestReceiver != -1 ) {
        if( mRecal && 
            !COST_EQUAL( 
                mBestReceiverCost[timeSlice][idU], bestCost ))
        {
            // check if the best receiver changed,
            // don't care about second best in the recalculation algorithm
            if( bestCost <= mBestReceiverCost[timeSlice][idU] )
                bestRecChanged = -1; // new score is less
            else
                bestRecChanged = 1; 
        }
        mBestReceiver[timeSlice][idU] = bestReceiver; 
        mBestReceiverCost[timeSlice][idU] = bestCost;
    } 

    if( secondBestReceiver != -1 ) {
        mSecondBestReceiver[timeSlice][idU] = secondBestReceiver; 
        mSecondBestReceiverCost[timeSlice][idU] = secondBestCost;
    }

    return bestRecChanged;
} 




/**
 * Compute and save bestReceiver and second best recievers.
 *
 * Find species nodes in this time slice (or just the nodes in
 * speciesNodeIdsTS) associated with least cost 
 * (bestReceiver) and the second least (or equal) cost.
 * The results are stored in member variables by calling
 * saveBestReceivers.
 *
 * @return -1 if new best receiver is less, 1 for greater, or 0 for no change
 */
int DTLMatrix::computeBestReceivers( 
    int timeSlice,      ///< time slice to save
    int idU,            ///< gene id
    vector<int> &speciesNodeIdsTS, ///< nodes to considier
    vector<double> &optVectorTS )  ///< costs for speciesNodeIdsTS
{

    int bestReceiver = -1;
    int secondBestReceiver = -1;
    double bestCost = std::numeric_limits<double>::max();
    double secondBestCost = std::numeric_limits<double>::max();

    // Save other bests for graph construction.
    // if best==secondBest, save all nodes equal to best cost
    // if best!=secondBest, save all nodes equal to second best cost
    vector<int> otherBests;
    bool first = true;
    for ( size_t z=0; z<speciesNodeIdsTS.size(); z++ ) {	

        int idXf = speciesNodeIdsTS[z];

        // ILS nodes cannot be the source of a transfer
        if( mSpeciesTree->hasILS() && mSpeciesTree->isILS( idXf ) ) 
            continue;

        double spNodeOpt = optVectorTS[z];  // current optimal
        if( first ) {
            first = false;
            bestReceiver = idXf;
            bestCost = spNodeOpt;  
        } else if( COST_GREATER( bestCost, spNodeOpt ) ) {

            // move best to second best
            if( timeSlice > 0 &&
                !COST_EQUAL( bestCost, numeric_limits<double>::max() ) )  
            {
#ifndef SPEED
                if(  COST_EQUAL( bestCost, secondBestCost ) ) 
                { 
                    // other bests were all equal to best receiver, 
                    // now keep same list, but now it is of second bests
                    otherBests.push_back( secondBestReceiver );
                } else {
                    // only a single second best
                    otherBests.clear();
                }
#endif
                // put old best as second best id and save old z
                // if not timeSlice=0 || not MAX
                secondBestReceiver = bestReceiver;
                secondBestCost = bestCost;
            }

            // new best
            bestReceiver = idXf;
            bestCost = spNodeOpt;  

        } else if( timeSlice == 0 
                && COST_EQUAL( spNodeOpt, numeric_limits<double>::max() ) )
        {
            //  don't add max cost to second best in time slice 0 
        } else if( COST_GREATER( secondBestCost, spNodeOpt ) ) {
            // save best receiver for best recevier
            secondBestReceiver = idXf;
            secondBestCost = spNodeOpt;

#ifndef SPEED
            // only a single second best, now
            otherBests.clear();
        } else if( COST_EQUAL( spNodeOpt, secondBestCost ) ) {
            otherBests.push_back( idXf );
#endif
        } 

    }

    if( first )
        throw bpp::Exception ("DTLMatrix::computeBestReceivers:"
                          "best receiver not set" );

    return saveBestReceivers( timeSlice, idU, bestReceiver, 
                    secondBestReceiver, bestCost, secondBestCost, otherBests );
}




/**
 * Calculate all suboptimal best and second beceivers (VT1 and VT2) 
 * for the current clade and time slice.
 *
 * The list of species ids for the current time slice and their costs are
 * given.
 *
 * The results are stored in the following member variables:
 * mVT1 - all costs within epsilon of the optimal for current time slice, paired
 *  with an id if the cost is unique for that id, else paired
 *  with -1 (meaning multiple ids for this cost)
 * mVT2 - if mVT1 contains a single id, then mVT2 contains all costs
 *        (> BRcost <= BR2 cost)
 * mVT1all, mVT2all - list of ALL ids with costs in VT1 (VT2) for
 *  generating graph.
 *  
 */
void DTLMatrix::computeVTs( 
    bool noSub, ///< no subdivision if true
    bool recompute, ///< replace existing values
    int timeSlice, ///< time slice (not used if noSub is true)
    int idU, ///< gene id
    vector<int> &speciesNodeIdsTS, ///< nodes to considier
    vector< vector<double> > &allCostsTS )
        ///< all suboptimal costs associated with speciesNodeIdsTS
{
    // create cost/id pairs from the optVector
    vector< pair<double,int> > costPairs;
    vector<int> idXs; // id lists ordered the same as speciesNodeIdsTS
    for( size_t nodeIndex=0; nodeIndex<speciesNodeIdsTS.size(); nodeIndex++ ) {
        int idX = speciesNodeIdsTS[nodeIndex];

        // ILS nodes cannot be the source of a transfer
        // DOES THIS COVER ILS CHANGES FOR SUBOPT?
        if( mSpeciesTree->hasILS() && mSpeciesTree->isILS( idX ) ) 
            continue;

        idXs.push_back( idX );
        if( recompute ) {
            vector<double> costs = mMatrixV->getValueSure(idU,idX); 
            BOOST_FOREACH( double cost, costs )
                costPairs.push_back( make_pair( cost, idX ) );
        } else {
            BOOST_FOREACH( double cost, allCostsTS[nodeIndex] )
                costPairs.push_back( make_pair( cost, idX ) );
        }
    }

    if( !noSub && recompute ) {
        mVT1[timeSlice][idU].clear();
        mVT2[timeSlice][idU].clear();
        mVT1all[timeSlice][idU].clear();
        mVT2all[timeSlice][idU].clear();
    }

    // sort 
    sort( costPairs.begin(), costPairs.end() );

    int br = -1;
    int br2 = -1;
    double brCost = 0;
    double br2cost = 0;
    if( !noSub ) {
        br = mBestReceiver[timeSlice][idU];
        brCost = mBestReceiverCost[timeSlice][idU];
        br2 = mSecondBestReceiver[timeSlice][idU];
        br2cost = mSecondBestReceiverCost[timeSlice][idU]; 
    }
    size_t loops = 1;
    if( noSub )
        loops = speciesNodeIdsTS.size();
    for( size_t i=0; i<loops; i++ ) {
        int ts = timeSlice;
        if( noSub ) {
            ts  = mSpeciesTree->getTimeSlice( speciesNodeIdsTS[i] );
            brCost = mBestReceiverCost[ts][idU];
            if( recompute ) {
                mVT1[ts][idU].clear();
                mVT2[ts][idU].clear();
                mVT1all[ts][idU].clear();
                mVT2all[ts][idU].clear();
            }
        }
        bool firstCost = true;
        double curCost = -1;
        double lastBestCost = -1;
        int bestId = -1;
        int costCount = 0;
        bool costGreater = false;
        for( size_t j=0; j<costPairs.size(); j++ ) {
            if( noSub && !mComparable[idXs[i]][costPairs[j].second] ) 
                continue;
            if( COST_GREATER( brCost, costPairs[j].first ) ) 
                throw bpp::Exception("DTLMatrix::computeVTs:"
                        " vt cost less than best receiver" );
            if( COST_GREATER( costPairs[j].first, brCost + mEpsilon ) ) {
                if( noSub ) 
                    break; // seen everything
               
                if( COST_GREATER( costPairs[j].first, br2cost + mEpsilon ))
                {
                    costGreater = true;
                    break; // second best finished, too
                } else if( br2 != -1  // second best is not optimal
                           && costPairs[j].second != br  )
                {
                    // second best receiver

                    // add cost to second optimums list (mVT2) 
                    // if within mEpsilon
                    if( !COST_EQUAL( costPairs[j].first, curCost ) ) {
                        mVT2[ts][idU].push_back( costPairs[j].first );
                        curCost = costPairs[j].first;
                    }
                    // all ids with a cost within epsilon for second best
                    // (except best receiver)
                    if( recompute ) 
                        mVT2all[ts][idU].push_back( costPairs[j].second );
                }
            } else { 
                // best receiver (keep one pair for each cost)
                if( firstCost || !COST_EQUAL( costPairs[j].first, curCost ) )
                {
                    if( firstCost )
                        firstCost = false;
                    else {
                        // cost is different, add the id for previous cost
                        mVT1[ts][idU].push_back( make_pair( curCost, bestId ) );
                    }
                    bestId = costPairs[j].second;
                    curCost = costPairs[j].first;
                    lastBestCost = curCost;
                    costCount = 1;
                } else {
                    // An id of -1 indicates that there are multiple species
                    // nodes with this cost.
                    costCount++;
                    bestId = -1;
                }
                // all ids with a cost within epsilon
                if( recompute ) 
                    mVT1all[ts][idU].push_back( costPairs[j].second );
            }
        }
        if( costCount > 0 ) 
            mVT1[ts][idU].push_back( make_pair( lastBestCost, bestId ) );

        // remove duplicates
        if( recompute ) {
            removeDuplicates( mVT1all[ts][idU] );
            if( !noSub ) {
                removeDuplicates( mVT2all[ts][idU] );
                if( costGreater )
                    break; // done all costs greater then br
           }
        }
    }
}



/**
 * Sort a list and remove duplicates.
 */
void DTLMatrix::removeDuplicates( 
        vector<int> &v ) ///< the list
{
    sort( v.begin(), v.end() );
    v.erase( unique( v.begin(), v.end() ), v.end() );
}






/**
 * Compute transfer cost with loss.
 *
 * @return true if computed cost is lower than opt
*/
bool DTLMatrix::computeTransferLossCost(
    double alphaCost, ///< cost for alpha of this time slice
    double &opt, ///< current optimal cost
    DTLMatrixState &state ) ///< various current info
{
    bool changed = false;
    double bestRecCost = mBestReceiverCost[state.timeSlice][state.idU]; 
    int bestReceiverId = mBestReceiver[state.timeSlice][state.idU];
    if( bestReceiverId == state.idX ) {
        //if the bestReceiver is x, we take the second one ..
        bestRecCost = mSecondBestReceiverCost[state.timeSlice][state.idU]; 
        bestReceiverId = mSecondBestReceiver[state.timeSlice][state.idU];
    }

    // Transfer from optimal node. Add loss cost if this is not alpha.
    char event = 'u';
    if( bestReceiverId != -1 ) {	
        double cost = bestRecCost + state.hgtCost;		
        event = 'y';  // TLFD
        if( !mSpeciesTree->hasAlpha() || !mSpeciesTree->isAlpha( state.idX ) )
        {
            cost += state.lossCost;  // don't add this for alpha
            event = 'x';
        }
        if( cost < opt ) {
            opt = cost;
            changed = true;
        }
    }

    // Transfer from alpha TLTD.
    if( mSpeciesTree->hasAlpha() && !mSpeciesTree->isAlpha( state.idX ) ) {
        double cost = alphaCost + state.lossCost; 
        if( cost < opt ) {
            opt = cost;
            event = 'z';
            changed = true;
        }
    }

    if( changed && mUseBestSplits ) {
       mBestSplits[state.idU][state.idX].event = event;
       mBestSplits[state.idU][state.idX].speciesSplit.first = bestReceiverId;
       mBestSplits[state.idU][state.idX].speciesSplit.second = -1;
    }

    return changed;
}


/**
 * Add all possible transfer costs within suboptimal epsilon to
 * allCosts.
 */
void DTLMatrix::computeTransferLossCostVT(
    vector<double> &alphaCosts, ///< alpha costs
    vector<double> &allCosts, ///< all suboptimal costs
    DTLMatrixState &state ) ///< various current info
{
    if( mBestReceiver[state.timeSlice][state.idU] == -1 ) 
        return; // no best receiver

    double otherCosts = state.hgtCost;

    //  transfer to alpha (TLTD) has no extra cost
    if( !mSpeciesTree->hasAlpha() || !mSpeciesTree->isAlpha( state.idX ) ) 
        otherCosts += state.lossCost;
   
    // for all VT, do a cost combination
    if( state.idX == mBestReceiver[state.timeSlice][state.idU] ) {
        // transfering to optimal species node
        BOOST_FOREACH( double c, mVT2[state.timeSlice][state.idU] ) {
            double cost = otherCosts + c;
            if( !COST_GREATER( cost, state.curBestEpsCost) ) 
                allCosts.push_back( cost );
        }
    } 

    pair<double,int> p;
    BOOST_FOREACH( p, mVT1[state.timeSlice][state.idU] ) {
        // don't add if it is this species node's cost
        if( p.second != state.idX ) {
            double cost = otherCosts + p.first;
            if( !COST_GREATER( cost, state.curBestEpsCost) ) 
                allCosts.push_back( cost );
        }
    }

    // Transfer from alpha.
    if( mSpeciesTree->hasAlpha() && !mSpeciesTree->isAlpha( state.idX ) ) { 
        BOOST_FOREACH( double alphaCost, alphaCosts ) {
            double cost = alphaCost + state.lossCost;
            if( !COST_GREATER( cost, state.curBestEpsCost) ) 
                allCosts.push_back( cost );
        }
    }
}



/**
 * Initialize matrix and best receivers.
 *
 */
void DTLMatrix::init() 
{
    mCladeCount = mCladesTrips->mClades.getCladeCount();

    // create matrix
    mMatrixV = NULL;
 	try{	
        // we can have at most 3GTnodes different clades all 
        // over all possible rootings of G
		mMatrix.setDim(mCladeCount,mSTnodes); 
	}
	catch (exception& genericException) {
     	throw bpp::Exception("----Low memory!!! The dataset is too big-------");
    }
		
    // First two best receivers and cost for each clade and timeslice.
    mBestReceiver = new int*[mMaxTS+1]; 
    mBestReceiverCost = new double*[mMaxTS+1]; 
    mSecondBestReceiver = new int*[mMaxTS+1]; 
    mSecondBestReceiverCost = new double*[mMaxTS+1]; 
	for( int timeSlice=0; timeSlice<=mMaxTS; timeSlice++ ) {
	    mBestReceiver[timeSlice] = new int[mCladeCount];
	    mBestReceiverCost[timeSlice] = new double[mCladeCount];
	    mSecondBestReceiver[timeSlice] = new int[mCladeCount];
	    mSecondBestReceiverCost[timeSlice] = new double[mCladeCount];
    }
    clearBestReceivers(); // initializes all values to -1
    

    // not used in all cases (e.g. for recalculation)
    mBestSplits = NULL;
    mMatrixCalculated = false;
    mEpsilon = 0;
    mSubOpt = false;
    mTriplets = false;
    mRecal = false;
    mUseILS = mSpeciesTree->hasILS();
}


/**
 * Reset best receiver values
 */
void DTLMatrix::clearBestReceivers() 
{
	for( int timeSlice=0; timeSlice<=mMaxTS; timeSlice++ ) 
	    for( int u=0; u<mCladeCount; u++ ) {
            mBestReceiver[timeSlice][u] = -1;
            mBestReceiverCost[timeSlice][u] = -1;
            mSecondBestReceiver[timeSlice][u] = -1;
            mSecondBestReceiverCost[timeSlice][u] = -1;
        }
}



/**
 * Basic constructor.
 *
 */
DTLMatrix::DTLMatrix( 
    MySpeciesTree *speciesTree, ///< species trees
    CladesAndTripartitions *cat, ///< clades and tripartitions
    double WGDCost,     ///< cost of whole genome duplication
    bool fixedCosts,    ///< use given costs rather than node specific
    bool computeT,      ///< compute transfer costs
    bool computeTL,     ///< compute TL costs
    double dupliCost,   ///< fixed duplication cost
    double hgtCost,     ///< fixed transfer cost
    double lossCost,    ///< fixed loss cost
    int maxTS,          ///< maximum time slice
    double weight,      ///< split weight
    bool useBestSplits, ///< record splits
    double ilsCost )     ///< cost of an ils event 
:
    mUseBestSplits(useBestSplits),
    mMaxTS(maxTS), 
    mSpeciesTree(speciesTree), 
    mCladesTrips(cat),
    mFixedCosts(fixedCosts), 
    mComputeT(computeT),
    mComputeTL(computeTL),
    mDupliCost(dupliCost), 
    mHGTCost(hgtCost), mLossCost(lossCost), mSplitWeight(weight), 
    mWGDCost(WGDCost),
    mIlsCost(ilsCost)
{	
    if( !mComputeT )
        mComputeTL = false;

	mSTnodes = mSpeciesTree->getNumberOfIds();

    init();

    // Create with same dimensions as matrix to store
    // best split for each cell. This is a pointer intto
    // cladeSplits structure. Note that this is not done when
    // reading the matrix. BestSplit is not stored and therefore
    // is incomplete.
    if( mUseBestSplits ) {
        mBestSplits = new BestSplit*[mCladeCount];
        for( int u=0; u<mCladeCount; u++ ) 
            mBestSplits[u] = new BestSplit[mSTnodes];
    }
}


/**
 * Destructor
 */
DTLMatrix::~DTLMatrix() 
{
    // clean up best receivers
	for( int ts=0; ts<=mMaxTS; ts++) {
        delete [] mBestReceiver[ts];
        delete [] mBestReceiverCost[ts];
        delete [] mSecondBestReceiver[ts];
        delete [] mSecondBestReceiverCost[ts];
        if( mSubOpt ) {
            delete [] mVT1[ts];
            delete [] mVT2[ts];
            delete [] mVT1all[ts];
            delete [] mVT2all[ts];
        }
    }
    delete [] mBestReceiver;
    delete [] mBestReceiverCost;
    delete [] mSecondBestReceiver;
    delete [] mSecondBestReceiverCost;

    if( mSubOpt ) {
        delete [] mVT1;
        delete [] mVT2;
        delete [] mVT1all;
        delete [] mVT2all;
    }

    if( mBestSplits != NULL ) {
        for( int u=0; u<mCladeCount; u++ ) 
            delete [] mBestSplits[u];
        delete [] mBestSplits;
    }

    if( mMatrixV != NULL ) 
        delete mMatrixV;

}



/**
 * Initialize data structures for suboptimal calculations.
 */
void DTLMatrix::initSubOpt() 
{
#ifdef SPEED
    throw bpp::Exception( "DTLMatrix: subopt called with SPEED defined" );
#endif
    if( mTriplets ) 
        throw bpp::Exception( "DTLMatrix: subopt called after triplets" );

    if( mSubOpt ) {
        // clean up previous
        for( int ts=0; ts<=mMaxTS; ts++) {
            delete [] mVT1[ts];
            delete [] mVT2[ts];
            delete [] mVT1all[ts];
            delete [] mVT2all[ts];
        }
        delete [] mVT1;
        delete [] mVT2;
        delete [] mVT1all;
        delete [] mVT2all;

        delete mMatrixV;
    }


    mSubOpt = true;
    mUseBestSplits = false;
    mMatrixV = new MyMatrixV();
    try{	
        // we can have at most 3GTnodes different clades all 
        // over all possible rootings of G
        mMatrixV->setDim(mCladeCount,mSTnodes); 
    }
    catch (exception& genericException){
        throw bpp::Exception("----Low memory!!! The dataset is too big-------");
    }

    mVT1 = new vector< pair<double,int> > *[mMaxTS+1];
    mVT2 = new vector<double> *[mMaxTS+1];
    mVT1all = new vector<int> *[mMaxTS+1];
    mVT2all = new vector<int> *[mMaxTS+1];
    for( int timeSlice=0; timeSlice<=mMaxTS; timeSlice++ ) {
        mVT1[timeSlice] = new vector< pair<double,int> > [mCladeCount];
        mVT2[timeSlice] = new vector<double> [mCladeCount];
        mVT1all[timeSlice] = new vector<int> [mCladeCount];
        mVT2all[timeSlice] = new vector<int> [mCladeCount];
    }
}

void DTLMatrix::printAllCostsError(
    double optVal, ///< optimum
    vector<double> &allCosts, ///< all suboptimal costs
    DTLMatrixState &state ) ///< various current info
{
    cout << "ts=" << state.timeSlice << " " << state.idU << "," << state.idX
     << " allCosts size=" << allCosts.size() << ": ";
    if( allCosts.size() < 1 ) 
        cout << "no values" << endl;
    else {
        cout << allCosts[0] << " != optVal=" << optVal << endl;
        BOOST_FOREACH( double d, allCosts)
            cout << "  " << d << endl;
    }
}


/**
 * Remove dupilcates and set matrix cost list.
 */
void DTLMatrix::suboptSetFinalCosts( 
    double optVal, ///< optimum
    vector<double> &allCosts, ///< all suboptimal costs
    DTLMatrixState &state, ///< various current info
    bool checkOptVal )
{

    sort( allCosts.begin(), allCosts.end() );

    // check matrix vector first value
    if( checkOptVal && optVal < numeric_limits<double>::max()
        && (allCosts.size() < 1 || !COST_EQUAL( allCosts[0], optVal) )) 
    {
        printAllCostsError( optVal, allCosts, state );
        throw bpp::Exception("DTLMatrix::suboptSetFinalCosts:"
                " matrix vector opt not right");
    }

    // remove duplicates
    vector<double> finalList;
    double prevValue = 0;
    bool firstLoop = true;
    BOOST_FOREACH( double c, allCosts ) 
    {
        // another check
        if( checkOptVal && optVal < numeric_limits<double>::max()
            && COST_GREATER( optVal, c ) ) 
        {
            printAllCostsError( optVal, allCosts, state );
            throw bpp::Exception("adding an impossible cost to the matrix");
        }

        if( COST_GREATER( c, optVal+mEpsilon ) )
            break;

        // always add the first one, but don't add values
        // close to previous values
        if( firstLoop || !COST_EQUAL( prevValue, c ) ) {
            finalList.push_back( c );
            prevValue = c;
        }

        firstLoop = false;
    }

    mMatrixV->setValues( state.idU, state.idX, finalList );
//cout << state.idU << "," << state.idX << " final:";
//BOOST_FOREACH( double d, finalList)
//    cout << "  " << d;
//cout << endl;
}

/**
 * Calculate and set a matrix cell values for the given clade
 * and nodes in the given time slice.
 */
void DTLMatrix::calculateMatrixTS( 
    int idU,    ///< clade id
    int timeSlice ) ///< current time slice
{
    // cell information to be passed to the cost calculation functions
    DTLMatrixState state;
    state.idU = idU;
    state.timeSlice = timeSlice;
    if( mFixedCosts ) {
        state.dupliCost = mDupliCost;
        state.hgtCost = mHGTCost;
        state.lossCost = mLossCost;
    }

    // nodes for this time slice
    vector<int> speciesNodeIdsTS = mSpeciesTree->getVectorWithTS( timeSlice );

    // Find maximum cost to keep for this time slice.
    if( mSubOpt || mTriplets ) {
        for ( size_t i=0; i<speciesNodeIdsTS.size(); i++ ) {
            state.idX = speciesNodeIdsTS[i];	
            if( i==0 || state.curBestEpsCost 
                    < mMatrix.getValueSure(state.idU,state.idX) )
            {
               state.curBestEpsCost = mMatrix.getValue(state.idU,state.idX);
            }
        }
        state.curBestEpsCost += mEpsilon + state.hgtCost + state.lossCost;
    }

    // compute costs 
    // score clade and store results in optVectorTS,
    // which contains tempOptima for each node in timeslice
	initTS(); // this does nothing in this class
    vector<double> optVectorTS; // costs of speciesNodeIdsTS
    vector< vector<double> > allCostsTS; // subopt costs
    vector<double> alphaCosts; // suboptimal alpha costs
    double alphaCost = 0;
    for ( size_t i=0; i<speciesNodeIdsTS.size(); i++ ) {
//cout << idU << "," << speciesNodeIdsTS[i] << " ts=" << timeSlice << endl;
        double optAllOverTriples = computeOptimaForCell( 
                                        speciesNodeIdsTS[i], 0, state );
//cout << idU << "," << speciesNodeIdsTS[i] << " ts=" << timeSlice << " opt=" << optAllOverTriples << endl;
        optVectorTS.push_back( optAllOverTriples ); 

        if( mSpeciesTree->isAlpha( speciesNodeIdsTS[i]) ) { 
            alphaCost = optAllOverTriples;
            int alphaIdx = (int) optVectorTS.size() - 1;
            if( mTriplets )
                setAlpha( alphaIdx );
            else if( mSubOpt ) 
                alphaCosts = state.costList;
        }
        if( mSubOpt ) 
            allCostsTS.push_back( state.costList );

        // Needed for ILS loss
        if( mUseILS ) {
            mMatrix.setValue( state.idU, speciesNodeIdsTS[i], 
                              optAllOverTriples );
            if( mSubOpt ) 
                suboptSetFinalCosts( optAllOverTriples, state.costList, state, 
                                     false );
            else if( mTriplets ) 
                setFinalCosts( i, state );
        }

    }

    // calculate and set best receivers (first and second least cost nodes)
    if( speciesNodeIdsTS.size() > 1 ) {
        if( mSubOpt ) 
            computeVTs( false, false, timeSlice, idU, speciesNodeIdsTS, 
                        allCostsTS ); 
        else if( mTriplets ) 
            computeVTs( false, false, timeSlice, idU, speciesNodeIdsTS ); 
        else
            computeBestReceivers(timeSlice, idU, speciesNodeIdsTS, optVectorTS);
    }


    // Calculate transfer loss costs and set matrix cell.
    for( size_t i=0; i<speciesNodeIdsTS.size(); i++ ) {
        double opt = optVectorTS[i]; // cost before considering TL event
        state.idX = speciesNodeIdsTS[i];	
        if( !mFixedCosts ) {
            MySpeciesNode *xNode = mSpeciesTree->getNodeById( state.idX );
            state.dupliCost = xNode->getInfos().duplicationCost;
            state.hgtCost = xNode->getInfos().hgtCost;
            state.lossCost = xNode->getInfos().lossCost;
        }

        // check for better opt by transfer
        if( mComputeTL && opt != 0 && speciesNodeIdsTS.size() > 1 ) {
            if( (!mSubOpt && !mTriplets) || mUseILS )
                computeTransferLossCost( alphaCost, opt, state );

            if( mSubOpt ) 
                computeTransferLossCostVT( alphaCosts, allCostsTS[i], state );
             else if( mTriplets ) 
                computeTransferLossCostVT( i, state );
        }

        if( (!mSubOpt && !mTriplets) || mUseILS )
            mMatrix.setValue( idU, state.idX, opt ); 

        // set matrix value
        if( mSubOpt ) 
            suboptSetFinalCosts( mMatrix.getValueSure( idU, state.idX ),
                                allCostsTS[i], state );
        else if( mTriplets ) 
            setFinalCosts( i, state );
    }

    // update VTs with costs from transfer losses
    if( speciesNodeIdsTS.size() > 1 ) {
        if( mSubOpt ) 
            computeVTs( false, true, timeSlice, idU, speciesNodeIdsTS, 
                        allCostsTS ); 
        else if( mTriplets ) 
            computeVTs( false, true, timeSlice, idU, speciesNodeIdsTS ); 
    }
}



/**
 * Calculate the matrix.
 *
 * If maxIterations is greater than 1, than the parameters
 * are re-estimated and the matrix re-calculated until
 * maxIterations or the cost converges.
 *
 * For suboptimal (epsilon != -1), calculateMatrix must be called
 * first without suboptimal (default epsilon = -1).
 *
 * If epsilon is zero, than all optimal costs nodes are calculated.
 *
 * The matrix can be partially recalculated using tsStart and tsEnd, which 
 * default to tsStart=0 (leaves) and tsEnd=-1 (mMaxTS=root). The matrix must
 * be fully calculated first.
 *
 */
void DTLMatrix::calculateMatrix( 
    bool verbose,  ///< print timing and stats
    int maxIterations,
        ///< Number of times to calculate matrix with recomputed parameters.
    bool updateWeightOnly, 
        ///< only update the weight costs, not DTL costs when doing iterations
    int dated, ///< 0=undated, 1=partially, 2=fully
    double epsilon, ///< suboptimal epsilon (ignored if -1)
    int tsStart,    ///< time slice from which to start calculation (default=0)
    int tsEnd )     ///< last time slice to calculate (default is -1 = mMaxTS).
{
    clock_t start = clock();

    if( epsilon != -1 ) {
        mEpsilon = epsilon;
        initCalculation(); // initializes structures for subopt calculations
        if( verbose ) {
            if( mEpsilon == std::numeric_limits<double>::max() ) 
                cout << "EPSILON unlimited" << endl;
            else
                cout << "EPSILON initialized " << mEpsilon << endl;
        }

        if( !mMatrixCalculated )
            throw bpp::Exception("DTLMatrix::calculateMatrix: non-suboptimal"
                " matrix calculation must be called before suboptimal call" );
        if( maxIterations != 1 ) 
            throw bpp::Exception("DTLMatrix::calculateMatrix: suboptimal"
               " calculations only can do one iteration." );
    }

    if( tsEnd == -1 )
        tsEnd = mMaxTS;
    else if( tsStart != 0 && !mMatrixCalculated ) 
        throw bpp::Exception("DTLMatrix::calculateMatrix: the matrix must be"
            " calculated first before using partial time-slices." );

    if( verbose )
		cout << mCladeCount << " x " << mSTnodes << " matrix" << endl; 

    // precompute comparable nodes for best recevier calculation
    if( dated == 1 ) // partially dated
        initComparable( true );
    else if( dated == 0 ) // undated
        initComparable( false );

    int iteration = 1;
    while( iteration <= maxIterations ) {

        // collect other best receivers needed for graph construction
        mAllOtherBestReceivers.clear();

//double t1= (double) (clock()-start) / CLOCKS_PER_SEC * 1000.0;
//cout << "matrix start time = " << t1 << endl;
        // loop over clades of same size
        vector< vector<int> > cladesBySize = mCladesTrips->getCladesBySize();
        BOOST_FOREACH( vector<int> &sameSizeClades, cladesBySize ) {
#ifdef PARALLEL
            #pragma omp parallel for // for OpenMP parallelization
#endif
            for( size_t p=0; p<sameSizeClades.size(); p++ ) {
                int idU = sameSizeClades[p];
                //cout << "idU " << idU << endl;
                if( dated != 2 )
                    calculateMatrixNoSub( idU );
                else
                    // loop over all time slices
                    for( int ts=tsStart; ts<=tsEnd; ts++ ) 
                        calculateMatrixTS( idU, ts );
            }
        }
//double t2= (double) (clock()-start) / CLOCKS_PER_SEC * 1000.0;
//cout << "matrix end time = " << t2 << endl;

        /* // Testing loop
        vector<int> sortedCladeNums = mCladesTrips->mClades.getSortedClades(); 
        for( int ts=tsStart; ts<=tsEnd; ts++ ) 
            BOOST_FOREACH( int idU, sortedCladeNums ) 
                    calculateMatrixTS( idU, ts );
        */

        iteration++;
        if( iteration <= maxIterations ) {
            if( !updateCosts( updateWeightOnly, verbose ) ) 
                break; // update costs; break loop if no changes
        }
    } // end iterations
		
    if( verbose ) {
        double time = (double) (clock()-start) / CLOCKS_PER_SEC * 1000.0;
        cout << "iterations = " << iteration-1 << endl;
        cout << "matrix time = " << time << endl;
    }

    mMatrixCalculated = true;
}



/**
 * Calculate matrix for species trees with no subdivisions.
 *
 * Species Tree has no subdivision and time slice = idX 
 * (set this here or before?).
 */
void DTLMatrix::calculateMatrixNoSub( 
    int idU )    ///< clade id
{

    DTLMatrixState state;
    state.idU = idU;
    if( mFixedCosts ) {
        state.dupliCost = mDupliCost;
        state.hgtCost = mHGTCost;
        state.lossCost = mLossCost;
    }

    if( mSubOpt || mTriplets ) {
        for( int ts=0; ts<=mMaxTS; ts++ ) {
            state.timeSlice = ts;
            vector<int> singleNode = mSpeciesTree->getVectorWithTS( ts );
            state.idX = singleNode[0];
            if( ts==0 || state.curBestEpsCost 
                    < mMatrix.getValueSure(state.idU,state.idX) )
            {
               state.curBestEpsCost = mMatrix.getValue(state.idU,state.idX);
            }
        }
        state.curBestEpsCost += mEpsilon + state.hgtCost + state.lossCost;
    }

    vector<double> optVector; // costs of speciesNodeIdsTS
    vector< vector<double> > allCosts; // subopt costs for each TS
    vector<int> speciesNodeIds;
    vector<double> curBestCosts; // for mSubOpt
    double alphaCost = 0;
    vector<double> alphaCosts;
	initTS(); // clear triplets
    for( int ts=0; ts<=mMaxTS; ts++ ) {
        state.timeSlice = ts;

        // node for this time slice (time slices are depth first ordering
        // of species nodes)
        vector<int> singleNode = mSpeciesTree->getVectorWithTS( ts );
        if( singleNode.size() != 1 )
            throw bpp::Exception( "DTLMatrix::calculateMatrixNoSub: "
                    " a time slice does not have a single node." );

        // compute costs (it appends state.costList)
        double optAllOverTriples = computeOptimaForCell( 
                                        singleNode[0], 0, state );
        optVector.push_back( optAllOverTriples ); 

        if( mSubOpt ) 
            allCosts.push_back( state.costList );

        if( mSpeciesTree->isAlpha( state.idX ) ) { 
            alphaCost = optVector[0];
            if( mTriplets ) {
                int alphaIdx = (int) optVector.size() - 1;
                setAlpha( alphaIdx );
            } else if( mSubOpt ) 
                alphaCosts = state.costList;
        }

        // set the matrix so that SL has something to use
        if( mSubOpt || mTriplets ) {

            if( mSubOpt && optVector[ts] !=  numeric_limits<double>::max() 
                && allCosts[ts].size() == 0 ) 
            {
                throw bpp::Exception( "DTLMatrix::calculateMatrixNoSub:" 
                        " allCosts is zero" );
            } 

            // This can be zero if best is max or if nothing
            // is better than epsilon (only after TL)
            if( mSubOpt && allCosts[ts].size() > 0 )
                suboptSetFinalCosts( optVector.back(), allCosts.back(), state,
                                     false );
            else if( mTriplets ) 
                setFinalCosts( ts, state );

            // Save value for setting cost list.
            // - pre TL value need for speciation calculations
            curBestCosts.push_back( mMatrix.getValueSure( idU, state.idX ) );
        }

        mMatrix.setValue( idU, state.idX, optVector.back() ); 
        speciesNodeIds.push_back( singleNode[0] );
    }

    // set best receiever and cost for each node
    computeBestReceiversNoSub( idU, false, speciesNodeIds, optVector );
    if( mSubOpt ) 
        computeVTs( true, false, -1, idU, speciesNodeIds, allCosts ); 
    else if( mTriplets ) 
        computeVTs( true, false, -1, idU, speciesNodeIds ); 



    if( !mComputeTL && !mSubOpt && !mTriplets )
        return;

    // Calculate transfer loss costs and set matrix cell.
    for( int ts=0; ts<=mMaxTS; ts++ ) {
        double opt = optVector[ts]; // cost before considering TL event
        state.idX = speciesNodeIds[ts];	

        // check for better opt by transfer
        if( mComputeTL && opt != 0 )
        {
            state.timeSlice = mSpeciesTree->getTimeSlice( state.idX );
            if( !mFixedCosts ) {
                MySpeciesNode *xNode = mSpeciesTree->getNodeById( state.idX );
                state.dupliCost = xNode->getInfos().duplicationCost;
                state.hgtCost = xNode->getInfos().hgtCost;
                state.lossCost = xNode->getInfos().lossCost;
            }

            if( mSubOpt ) {
                computeTransferLossCostVT( alphaCosts, allCosts[ts], state );
                suboptSetFinalCosts( curBestCosts[ts], allCosts[ts], state );
            } else if( mTriplets ) {
                computeTransferLossCostVT( ts, state );
                setFinalCosts( ts, state );
            } else {
                computeTransferLossCost( alphaCost, opt, state );
                mMatrix.setValue( idU, state.idX, opt ); 
            }
        }

        // make sure values are reset, even if TL is not computed
        if( mSubOpt || mTriplets ) 
            mMatrix.setValue( idU, state.idX, curBestCosts[ts] ); 
    }

    computeBestReceiversNoSub( idU, true, speciesNodeIds, optVector );
    if( mSubOpt ) 
        computeVTs( true, true, -1, idU, speciesNodeIds, allCosts ); 
    else if( mTriplets ) 
        computeVTs( true, true, -1, idU, speciesNodeIds ); 
}

/**
 * Fill mComparable. mComparable[xi][xj] is true if species
 * nodes xi and xj are comparable.
 * 
 * Comparable nodes are siblings of ancestors and their descendants.
 */
void DTLMatrix::initComparable(
        bool partialDates ) ///< species tree is partially dated 
{ 
    vector<MySpeciesNode*> speciesNodes = mSpeciesTree->getNodes();

    // make a matrix of comparable nodes
    mComparable.resize( speciesNodes.size() );
    vector<double> bootstraps;
    vector<double> parentBootstraps;
    for( size_t i=0; i<speciesNodes.size(); i++ ) {
        mComparable[i].resize( speciesNodes.size(), true );


        if( partialDates ) {
            MySpeciesNode *node = speciesNodes[i];
            double bs = -1;
            if( node->hasBranchProperty(bpp::TreeTools::BOOTSTRAP) ) 
                bs = dynamic_cast<const bpp::Number<double> *> 
               (node->getBranchProperty(bpp::TreeTools::BOOTSTRAP))->getValue();
            bootstraps.push_back( bs );

            double parentBs = -1;
            if( node->hasFather() ) {
                MySpeciesNode *parent = node->getFather();
                if( parent->hasBranchProperty(bpp::TreeTools::BOOTSTRAP) ) {
                    parentBs = dynamic_cast<const bpp::Number<double> *> 
                        (parent->getBranchProperty(bpp::TreeTools::BOOTSTRAP))
                            ->getValue();
                }
            }
            parentBootstraps.push_back( parentBs );
        }
    }

    for( size_t i=0; i<speciesNodes.size(); i++ ) {
        MySpeciesNode *node = speciesNodes[i];
        int idX = node->getId();
        mComparable[idX][idX] = false;
        MySpeciesNode *fatherNode = node;
        while( fatherNode->hasFather() ) {
            fatherNode = fatherNode->getFather();
            int parentX = fatherNode->getId();
            mComparable[idX][parentX] = false;
            mComparable[parentX][idX] = false;
        }

        if( partialDates && bootstraps[i] != -1 ) {
            // compare bootstrap values, valid if x>y and x<P(y)
            for( size_t j=0; j<speciesNodes.size(); j++ ) {
                if( bootstraps[j] == -1 || parentBootstraps[j] == -1 )
                    continue;

                MySpeciesNode *y = speciesNodes[j];
                int idY = y->getId();
                if( bootstraps[i] > parentBootstraps[j] 
                    || bootstraps[i] > parentBootstraps[j] ) 
                {
                    mComparable[idX][idY] = false;
                    mComparable[idY][idX] = false;
                }
            }
        }
    }
}


/**
 * Find best receivers for each species node.
 */
void DTLMatrix::computeBestReceiversNoSub( 
        int idU, ///< gene id
        bool recompute, ///< use matrix value instead of optVector
        vector<int> &speciesNodeIds, ///< nodes to considier
        vector<double> &optVector )  ///< costs for speciesNodes
{
    // create cost/id pairs from the optVector
    vector< pair<double,int> > costPairs;
    int aMaxId = -1; // any node id whose cost is numeric limit max
    vector<int> idXs; // id lists ordered the same as speciesNodeIds
    for( size_t nodeIndex=0; nodeIndex<speciesNodeIds.size(); nodeIndex++ ) {
        int idX = speciesNodeIds[nodeIndex];
        idXs.push_back( idX );
        double cost;
        if( recompute ) {
            cost = mMatrix.getValueSure( idU, idX );
        } else {
            cost = optVector[nodeIndex];
        }
#ifndef SPEED 
        // don't include all max for speed (only need one)
        if( cost < std::numeric_limits<double>::max() ) 
            costPairs.push_back( make_pair( cost, idX ) );
        else if (aMaxId == -1 ) 
            aMaxId = idX;
#else
        costPairs.push_back( make_pair( cost, idX ) );
#endif
    }

    sort( costPairs.begin(), costPairs.end() );

    // find best receivers for each node
    for( size_t i=0; i<speciesNodeIds.size(); i++ ) {
        int bestReceiver = aMaxId;
        double bestCost = std::numeric_limits<double>::max();
        bool first = true;
        vector<int> otherBests;
        vector<bool> &comp = mComparable[idXs[i]]; 
        int secondBestReceiver = -1;
        double secondBestCost = std::numeric_limits<double>::max();
        for( size_t j=0; j<costPairs.size(); j++ ) {
            if( !comp[costPairs[j].second] ) 
                continue;

            if( first ) {
                first = false;
                bestCost = costPairs[j].first;
                bestReceiver = costPairs[j].second;
#ifndef SPEED 
            // second and otherBests only needed for graph
            } else if( COST_EQUAL( costPairs[j].first, bestCost)) {
                if( secondBestReceiver == -1 ) {
                    secondBestReceiver = costPairs[j].second;
                    secondBestCost = bestCost;
                } else if( recompute )
                    otherBests.push_back( costPairs[j].second );
            } else {
#endif
                break; // all costs greater
            }
        }

        int ts = mSpeciesTree->getTimeSlice( speciesNodeIds[i] );
        saveBestReceivers( ts, idU,
                bestReceiver, secondBestReceiver, bestCost, 
                secondBestCost, otherBests );
    }
}



/**
 * Print the matrix in CSV format.
 */
void DTLMatrix::printMatrixCSV( 
    const char *fileName,       ///< output file name
    const char *speciesTreeFile, ///< species tree file
    const char *genesTreeFile ) ///< gene tree file
{
    ofstream csvfile( fileName );
    if (!csvfile.is_open()) {
        cout << "Unable to open file " << fileName << endl;
        return;
    }

    csvfile << "#speciesTreeFile=" << speciesTreeFile << endl;
    csvfile << "#genesTreeFile=" << genesTreeFile << endl;
    if( mFixedCosts )
        csvfile << "#fixedCosts=true" << endl;
    else
        csvfile << "#fixedCosts=false" << endl;
    if( mComputeTL )
        csvfile << "#computeTL=true" << endl;
    else
        csvfile << "#computeTL=false" << endl;
    csvfile << "#dupliCost=" << mDupliCost << endl;
//    csvfile << "#WGDCost=" << mWGDCost << endl;
    csvfile << "#hgtCost=" << mHGTCost << endl;
    csvfile << "#lossCost=" << mLossCost << endl;
    csvfile << "#weight=" << mSplitWeight << endl;
    if( mSpeciesTree->hasAlpha() ) 
        csvfile << "#transferDead=true" << endl;
    else 
        csvfile << "#transferDead=false" << endl;

    // clade numbers
    csvfile << "node/clade";
    for( int i=0; i<mCladeCount; i++ ) 
        csvfile << "," << i;
    csvfile << endl;

    // print matrix
    csvfile.precision( PRINT_MATRIX_PRECISION );
    for( int j=0; j<mSTnodes; j++ ) {
        csvfile << j;
        for( int i=0; i<mCladeCount; i++ ) 
            csvfile << "," << mMatrix.getValueSure(i,j);
        csvfile << endl;
    }

    //print best receivers
    csvfile << "bestReceivers" << endl;
    for( int t=0; t<=mMaxTS; t++ ) {
        csvfile << t;
        for( int u=0; u<mCladeCount; u++ ) 
            csvfile << "," << mBestReceiver[t][u] 
                    << ";" << mSecondBestReceiver[t][u];
        csvfile << endl;
    }

    //print best receivers cost
    csvfile << "bestReceiversCost" << endl;
    for( int t=0; t<=mMaxTS; t++ ) {
        csvfile << t;
        for( int u=0; u<mCladeCount; u++ ) 
            csvfile << "," << mBestReceiverCost[t][u] 
                    << ";" << mSecondBestReceiverCost[t][u];
        csvfile << endl;
    }

    csvfile.close();
}




/**
 * Finds lowest cost species node for gene tree root (all clade).
 *
 * @return Best cost most parsimonous reconciliation.
 */
double DTLMatrix::getBestCost( 
    int &rootClade, ///< root clade
    int &x ) ///< species id for best cost
{
    rootClade = mCladesTrips->mClades.getRootClade();
    x = 0;
    double bestCost = mMatrix.getValue( rootClade, 0 );
    for ( int j=1; j<mSTnodes; j++) {	
        double cost = mMatrix.getValue( rootClade, j );	
        if( cost < bestCost ) { 
            bestCost = cost;
            x = j;
        }
    }

    return bestCost;
}


/**
 * Finds lowest cost species node for gene tree root (all clade).
 *
 * @return Best cost most parsimonous reconciliation.
 */
double DTLMatrix::getBestCost( 
        bool keepRoot ) ///< use species root
{
    if( keepRoot ) {
        int rootClade = mCladesTrips->mClades.getRootClade();
        int x = mSpeciesTree->getRootNode()->getId();
        return mMatrix.getValue( rootClade, x );
    }

    int u;
    int x;
    return getBestCost( u, x );
}



/**
 * Create a consensus gene tree through backtracking of best
 * splits.
 *
 * @return return current node (root node returned finally)
 */
MyGeneNode *DTLMatrix::backtrack( 
    bool printTree, ///< print tree to cout if true
    bool keepRoot, ///< use species root
    int idU, ///< current clade id in recursion
    int idX, ///< current species id in recursion
    int level ) ///< recursion depth
{
    // Need best splits, but not calculated for recalculation or suboptimal.
    if( mBestSplits == NULL ) 
        throw bpp::Exception( 
                "DTLMatrix::backtrack: Backtracking not available." );

    if( idU < 0 ) { // root, get ids for best cost
        if( keepRoot ) {
            idU = mCladesTrips->mClades.getRootClade();
            idX = mSpeciesTree->getRootNode()->getId();
        } else
            getBestCost( idU, idX );
    }

    BestSplit bestSplit = mBestSplits[idU][idX];
    char bse = bestSplit.event;
    MyGeneNode *node = NULL;
    string leafName = "";
    if( bse == 'f' ) { // leaf
        node = new MyGeneNode();
        leafName = mCladesTrips->mClades.getLeafName( idU );
        node->setName( leafName );
        leafName = " ===LEAF " + leafName;
    }

    // print node info
    if( printTree && bse != 'n' ) {
        for( int i=0; i<level; i++ ) cout << "  ";
        cout << idU << "/" << idX << "  "  << mMatrix.getValue( idU, idX ) 
             << "(" << bse << ") -> u="
             << bestSplit.cladeSplit.first << "/"<< bestSplit.cladeSplit.second
             << " x=" << bestSplit.speciesSplit.first << "/"
             << bestSplit.speciesSplit.second << leafName << endl;
    }

    if( bse == 't' || bse == 'd' || bse == 's' || bse == 'a' || bse == 'i' ) {
        // clade split
        node = new MyGeneNode();
        MyGeneNode *son = backtrack( printTree, keepRoot, 
                                bestSplit.cladeSplit.first, 
                                 bestSplit.speciesSplit.first, level+1 );
        node->addSon( son );
        son = backtrack( printTree, keepRoot, bestSplit.cladeSplit.second, 
                         bestSplit.speciesSplit.second, level+1 );
        node->addSon( son );
    } else if( bse == 'l' || bse == 'n' || bse == 'w' || bse == 'x' 
            || bse == 'y' || bse == 'z' || bse == 'j'  )
    {
        // species event
        node = backtrack( printTree, keepRoot, idU, 
                          bestSplit.speciesSplit.first, level );
    }
    
    if( node == NULL ) {
        cout << idU << "," << idX
             << " unknown bestSplit.event = " << bse << endl;
        throw bpp::Exception( "DTLMatrix::backtrack: unknown bestSplit event" );
    }

    return node;
}


/**
 * Calculate new costs based on event counts.
 *
 * @return true if any costs changed
 */
bool DTLMatrix::updateCosts( 
    bool weightOnly, ///< only update split weight
    bool verbose  )  ///< print timing and stats
{
    // Node specific costs not yet added to this function.
    if( !mFixedCosts ) 
        throw bpp::Exception("DTLMatrix::updateCosts: fixed costs only" );

    int duplications = 0;
    int transfers = 0;
    int losses = 0;
    int ils = 0;
    double N_am = countEvents( duplications, transfers, losses, ils );
    
    double alpha = 1; // pseudocounts
    double Nf_am = abs(N_am - mCladesTrips->getRootMpAm() );
    double eventCount = duplications + transfers + losses + ils + Nf_am;
    double newWeight = -log10( (double) (Nf_am+alpha)/(eventCount + 4*alpha) ); 

    if( verbose ) {
        cout << "nAm=" << N_am << ", duplications=" << duplications 
             << ", transfers=" << transfers << ", losses=" << losses;
        if( mSpeciesTree->hasILS() )
            cout << ", ils=" << ils;
        cout << endl;
    }

    bool changed = false;
    if( !COST_EQUAL( mSplitWeight, newWeight ) ) {
        cout << "  weight=" << mSplitWeight << " -> " << newWeight << endl;
        mSplitWeight = newWeight;
        changed = true;
    }

    if( !weightOnly ) {
        double denom = eventCount + 4*alpha;
        double newDupliCost = -log10( (double) (duplications+alpha) / denom);
        if( !COST_EQUAL( mDupliCost, newDupliCost ) ) {
            if( newDupliCost == 0 )
                newDupliCost = 0.001;
            cout << "  duplication cost = " << mDupliCost 
                << " -> " << newDupliCost<< endl;
            mDupliCost = newDupliCost;
            changed = true;
        }
        double newHGTCost = -log10( (double) (transfers+alpha) / denom );
        if( !COST_EQUAL( mHGTCost, newHGTCost )) {
            if( newHGTCost == 0 )
                newHGTCost = 0.001;
            cout << "  HGT cost = " << mHGTCost << " -> " << newHGTCost << endl;
            mHGTCost = newHGTCost;
            changed = true;
        }
        double newLossCost = -log10( (double) (losses+alpha) / denom );
        if( !COST_EQUAL( mLossCost, newLossCost ) ) {
            if( newLossCost == 0 )
                newLossCost = 0.001;
            cout << "  loss cost = " << mLossCost << " -> " 
                << newLossCost << endl;
            mLossCost = newLossCost;
            changed = true;
        }
        double newIlsCost = -log10( (double) (ils+alpha) / denom );
        if( !COST_EQUAL( mIlsCost, newIlsCost ) ) {
            if( newIlsCost == 0 )
                newIlsCost = 0.001;
            cout << "  ils cost = " << mIlsCost << " -> " 
                << newIlsCost << endl;
            mIlsCost = newIlsCost;
            changed = true;
        }
    }

    return changed;
}



/**
 * Count events from the best consensus gene tree 
 * through backtracking of mBestSplits.
 *
 * @return nAm - sum of split costs 
 */
double DTLMatrix::countEvents( 
    int &duplications, ///< number of duplications
    int &transfers, ///< number of duplications
    int &losses, ///< number of duplications
    int &ils )  ///< number of ils events
{
    // Need best splits, but not calculated for recalculation or suboptimal.
    if( !mUseBestSplits ) 
        throw bpp::Exception( "DTLMatrix::countEvents:"
                " events not enabled (useBestSplits)." );

    int idU;
    int idX;
    getBestCost( idU, idX );

    double cost; // not returned
    return countEventsAux( idU, idX, duplications, transfers, losses, 
                            ils, cost );
}


/**
 * Recursive part of countEvents.
 *
 * Traverses the tree stored in mBestSplits.
 *
 * @return nAm - sum of split costs 
 */
double DTLMatrix::countEventsAux( 
    int idU, ///< current clade number in recursion
    int idX, ///< current species id in recursion
    int &duplications, ///< calculated number of duplications
    int &transfers, ///< calculated number of duplications
    int &losses, ///< calculated number of duplications
    int &ils, ///< calculated number of ils events
    double &cost ) ///< cost of current node
{
    BestSplit bestSplit = mBestSplits[idU][idX];
    cost = mMatrix.getValueSure(idU,idX);
    if( cost == std::numeric_limits<double>::max() ) 
        return 0; // OVERFLOW

    char bse = bestSplit.event;
//cout << idU << "," << idX << " " << bse << endl;

    // get event counts and costs
    double eventCost = 0;
    if( bse == 'd' ) {      // duplication 
        duplications++;
        eventCost = mDupliCost;
// Add this when WGD is implemented.
//        if( state.xNode->getInfos().WGD ) 
//            cost += mWGDCost;	// WG Duplication
    } else if( bse == 't' || bse == 'y' ) {  // transfer or tl from alpha
        transfers++;
        eventCost = mHGTCost;
    } else if( bse == 'l'  || bse == 'z' || bse == 'w') 
    {  // sl, tl to alpha, null: wgd + loss 
        losses++;
        eventCost = mLossCost;
    } else if( bse == 'i' ) { // ils
        ils++;
        eventCost = mIlsCost;
    } else if( bse == 'j' ) { // ils + loss
        ils++;
        losses++;
        eventCost = mIlsCost + mLossCost;
    } else if( bse == 'x' ) { // transfer loss (tl)
        transfers++;
        losses++;
        eventCost = mLossCost+mHGTCost;
    } else if( bse == 'f' || bse == 's' || bse == 'n' || bse == 'a' ) { 
        // leaf, speciation, null, alpha dupl or transfer  -- no costs 
    } else { 
        cout << "event: <" << bestSplit.event << ">" << endl;
        throw bpp::Exception("DTLMatrix::countEvents: undefined event type" );
    }

    // recursion 
    double nAm = 0;
    double calcCost = -1;
    if( bse == 't' || bse == 'd' || bse == 's' || bse == 'a' || bse == 'i' ) {
        // clade split
        double firstCost, secondCost;
        nAm += countEventsAux( 
                 bestSplit.cladeSplit.first, bestSplit.speciesSplit.first,
                 duplications, transfers, losses, ils, firstCost );
        nAm += countEventsAux( 
                 bestSplit.cladeSplit.second, bestSplit.speciesSplit.second,
                 duplications, transfers, losses, ils, secondCost );
        calcCost = firstCost+secondCost+eventCost;

        // find the ordering for this split to find split occurrences ratio
        int splitCount = mCladesTrips->getSplitCount( idU );
        int i=0;
        for( ; i<splitCount; i++ ) {
            pair<int,int> cladeSplit = mCladesTrips->getCladeSplit( idU, i );
            if( cladeSplit == bestSplit.cladeSplit ) 
                break;
        }
        if( i == splitCount ) 
            throw bpp::Exception("DTLMatrix::countEvents:"
                    " unable to recover split" );
        
        nAm += -log10( mCladesTrips->getSplitRatio( idU, i ) );

    } else if( bse == 'l' || bse == 'n' || bse == 'w' || bse == 'x' 
               || bse == 'y' || bse == 'z' || bse == 'j' )
    {
        // species split
        int childId = bestSplit.speciesSplit.first;

        if( bse == 'z' ) {
            // get alphaId for the current time slice
            int ts = mSpeciesTree->getTimeSlice( idX );
            childId = mSpeciesTree->getAlphaIdForTS( ts );
        }

        double childCost;
        nAm = countEventsAux( idU, childId, duplications, transfers, 
                              losses, ils, childCost );
        calcCost = childCost+eventCost;
    }

    if( calcCost != -1 && mSplitWeight == 0 && !COST_EQUAL( cost, calcCost) ) {
        cout << idU << "/" << idX << ":" << bse << "=" << cost << " != " 
             << calcCost << endl;
        throw bpp::Exception("DTLMatrix::countEvents: event cost wrong" );
    }

    return nAm;
}




/**
 * Add graph vertices if computed costs match best costs.
 *
 * Create an event vertex off of the given vertex z and create
 * pair vertex from the given ids. The cost of this event is
 * computed. If the cost matches the current cost, add it to
 * the graph.
 */
void DTLMatrix::addVertices( 
    DTLGraph &graph, ///< the graph
    DTLGraph::MyGraph::Vertex z, ///< vertex to attach to 
    vector<DTLGraph::MyGraph::Vertex> &qList, ///< queue new vertices
    int u,     ///< clade id 
    int x,     ///<  species id 
    const char *event, ///< event type to add
    double curCost,   ///< cost to compare 
    double otherCost, ///< costs for this event
    EventTriplet zTrip ) ///< event triplet
{
	vector<double> values;
    if( mSubOpt ) 
		values = mMatrixV->getValue(u,x);
    else if( mTriplets ) {
throw bpp::Exception ("SHOULD NEVER SEE THIS");
    } else
		values.push_back( mMatrix.getValueSure(u,x) );

    for( size_t i=0; i<values.size(); i++ ) {
        // check if this cost matches the matrix cost
	    double cost = values[i] + otherCost; 
	    if( !COST_EQUAL( cost, curCost) ) 
	        continue;

	    DTLGraph::MyGraph::Vertex eventVertex 
	                = graph.addEventVertex( z, event );
        DTLGraph::MyGraph::Vertex pairVertex;
        bool isNew = graph.addPairVertex( eventVertex, u, x, values[i], 
                                         pairVertex );
        if( isNew )
	        qList.push_back(pairVertex); // add to queue if new vertex
	}

}


/**
 * Add graph vertices for a split if computed costs match best costs.
 *
 * Create an event vertex with off of the given vertex z and create
 * pair vertex from the given ids for each cost combination.
 */
void DTLMatrix::addVertices( 
    DTLGraph &graph, ///< the graph
    DTLGraph::MyGraph::Vertex z, ///< vertex to attach to 
    vector<DTLGraph::MyGraph::Vertex> &qList, ///< queue new vertices
    int u1,     ///< first clade id 
    int x1,     ///< first species id 
    int u2,     ///< second clade id 
    int x2,     ///< second species id
    const char*event, ///< event type to add
    double curCost,   ///< cost to compare 
    double otherCost, ///< costs for this event
    EventTriplet zTrip ) ///< event triplet
{
	vector<double> values1;
	vector<double> values2;
    if( mSubOpt ) { 
		values1 = mMatrixV->getValue(u1,x1);
		values2 = mMatrixV->getValue(u2,x2);
    } else if( mTriplets ) {
throw bpp::Exception ("SHOULD NEVER SEE THIS2");
	} else {
		values1.push_back( mMatrix.getValueSure(u1,x1) );
		values2.push_back( mMatrix.getValueSure(u2,x2) );
    }
    // create an event vertex and 2 pair vertices for each combination of costs
    for( size_t i=0; i<values1.size(); i++ ) {
        for( size_t j=0; j<values2.size(); j++ ) {
            // check if this cost matches the matrix cost
            double cost = values1[i] + values2[j]+ otherCost; 
            if( !COST_EQUAL( cost, curCost) ) 
                continue;


            DTLGraph::MyGraph::Vertex eventVertex 
                        = graph.addEventVertex( z, event );
            DTLGraph::MyGraph::Vertex pairVertex;

			// ordering to make output comparisions easier,
            // order insertions the same way
            bool isNew;
            if( (u1 > u2) || (u1==u2 && x1 > x2)
                || (u1==u2 && x1==x2 && values1[i]>values2[j]) )
            {
                // add second first
                isNew = graph.addPairVertex( eventVertex, u2, x2, 
                                            values2[j], pairVertex );
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
                isNew = graph.addPairVertex( eventVertex, u1, x1, 
                                            values1[i], pairVertex ); 
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
			} else {
                isNew = graph.addPairVertex( eventVertex, u1, x1, 
                                        values1[i], pairVertex ); 
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
                isNew = graph.addPairVertex( eventVertex, u2, x2, 
                                        values2[j], pairVertex ); 
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
            }
		}
    }
}



/**
 * Creates a structure for all best receivers.
 *
 * Matrix calculation needs only best and second best receiver.
 * This function creates allBestReceivers using mAllOtherBestReceivers 
 * as well as best and second best.
 *
 * @return allBestReceivers
 */
vector<int> ***DTLMatrix::getAllBestReceivers() 
{
    vector<int> ***allBestReceivers = new vector<int>**[mMaxTS+1];

    // Find alphas and all best receivers.
    for( int ts=0; ts<=mMaxTS; ts++ ) { 

        allBestReceivers[ts] = new vector<int>*[mCladeCount];
        for( int idU=0; idU<mCladeCount; idU++ ) {
            allBestReceivers[ts][idU] = new vector<int>[2];

            // put best in first list if it exists
            if( mBestReceiver[ts][idU] == -1 ) 
                continue;
            allBestReceivers[ts][idU][0].push_back( mBestReceiver[ts][idU] );

            // second best in first list if the same cost (don't add it)
            if( COST_EQUAL( mBestReceiverCost[ts][idU],
                    mSecondBestReceiverCost[ts][idU] ) )
            {
                allBestReceivers[ts][idU][0].push_back( 
                        mSecondBestReceiver[ts][idU] );
            } else if( mSecondBestReceiver[ts][idU] != -1 ) {
                allBestReceivers[ts][idU][1].push_back( 
                        mSecondBestReceiver[ts][idU] );
            }
        }
    }

    // add other best receivers beyond second best
    size_t i = 0;
    while( i < mAllOtherBestReceivers.size() ) {
        int ts = mAllOtherBestReceivers[i++];
        int idU = mAllOtherBestReceivers[i++];
        int idX = mAllOtherBestReceivers[i++];
        if( COST_EQUAL( mBestReceiverCost[ts][idU], 
                         mSecondBestReceiverCost[ts][idU] ))
        {
            allBestReceivers[ts][idU][0].push_back( idX );
        } else 
            allBestReceivers[ts][idU][1].push_back( idX );
    }

    return allBestReceivers;
}


/**
 * Get the root nodes of the reconciliation graph.
 *
 * The root nodes are all of the species nodes with minimal cost
 * at the clade root.
 *
 * @return root nodes
 */
vector<DTLGraph::MyGraph::Vertex> DTLMatrix::getRootNodes(
    DTLGraph &graph, ///< the graph
    bool keepRoot ) ///< use species root only
{
    vector<DTLGraph::MyGraph::Vertex> rootList;
	double bestCost = getBestCost();
    if( mSubOpt )
        bestCost += mEpsilon;
    int rootClade = mCladesTrips->mClades.getRootClade();
    int j=0;
    if( keepRoot )
        j = mSpeciesTree->getRootNode()->getId();
    for( ; j<mSTnodes; j++ ) {	
		vector<double> values;
		if( mSubOpt ) {
			values = mMatrixV->getValue(rootClade,j);
		} else {
			values.push_back( mMatrix.getValue(rootClade,j) );
		}
        for( size_t i=0; i<values.size(); i++ ) {
            if( !COST_GREATER( values[i], bestCost ) ) {
                DTLGraph::MyGraph::Vertex vertex
                    = graph.addRoot( rootClade, j, values[i]);
                rootList.push_back( vertex );
        	}
        }
        if( keepRoot ) 
            break;
    }

    return rootList;
}


/**
 * Add all possible vertices to the given vertex.
 *
 * Check all possible events and add vertices for events
 * that have the same cost as the vertex.
 */
void DTLMatrix::createVertices(
    DTLGraph &graph, ///< the graph
    DTLGraph::MyGraph::Vertex pairVertex, ///< vertex to attach to 
    vector<DTLGraph::MyGraph::Vertex> &qList,  ///< queue new vertices
    vector<int> *** allBestReceivers ) ///< all best receivers
{
    int idU, idX;
    double cost;
    EventTriplet trip;
    if( mTriplets )
        graph.getVertexIdentfiers( pairVertex, idU, idX, cost, 
                                        trip.d, trip.t, trip.l );
    else
        graph.getVertexIdentfiers( pairVertex, idU, idX, cost );

    int ts = mSpeciesTree->getTimeSlice( idX );
    pair<int,int> split;
    vector< pair<int,int> > speciesSplits = mSpeciesTree->getSplits( idX );

    int alphaIdForTS = -1;
    if( mSpeciesTree->hasAlpha() )
        alphaIdForTS = mSpeciesTree->getAlphaIdForTS( ts );


    // Only graphs on a single tree (no amalgamations).
    // So, there is only a single split.
    int splitCount = mCladesTrips->getSplitCount( idU );
    if( splitCount != 1 )
        throw bpp::Exception( "DTLMatrix::createVertices: splitCount not one" );
    pair<int,int> cladeSplit = mCladesTrips->getCladeSplit( idU, 0 );

    // leaf
    MySpeciesNode *xNode = mSpeciesTree->getNodeById( idX );
    if( speciesSplits.size() == 0 && cladeSplit.first == -1
        && mCladesTrips->mClades.getSpeciesName(idU) == xNode->getName() )
    {
        graph.addEventVertex( pairVertex, "C" ); 
    }

    // if not a leaf
    if( cladeSplit.first != -1 ) {
        int idUl = cladeSplit.first;
        int idUr = cladeSplit.second;

        double costIdU = -(mSplitWeight * 
                 log10( mCladesTrips->getSplitRatio( idU, 0 )));
        if( costIdU != 0 ) 
            throw bpp::Exception ("DTLMatrix::createVertices: bad costIdU");

        // Speciation costs
        BOOST_FOREACH( split, speciesSplits ) {
            if( split.second != -1 ) {
                addVertices( graph, pairVertex, qList, idUl, split.first, 
                             idUr, split.second, "S", cost, 0, trip );
                addVertices( graph, pairVertex, qList, idUl, split.second, 
                             idUr, split.first, "S", cost, 0, trip );
            }
        }

        // Duplication costs
        addVertices( graph, pairVertex, qList, idUl, idX, 
                     idUr, idX, "D", cost, mDupliCost, trip ); 
        // Duplication not in Alpha cost 0
        if( mSpeciesTree->hasAlpha() && idX == alphaIdForTS ) 
            addVertices( graph, pairVertex, qList, 
                         idUl, idX, idUr, idX, "DD", cost, 0,  trip ); 
        //if( xNode->getInfos().WGD ) // WGD
        //    addVertices( graph, pairVertex, qList, idUl, 
        //                 idX, idUr, idX, "WGD", cost, mWGDCost, trip ); 
       
        // Transfer costs
        int idUsub = idUl; 
        int idUother = idUr;
        if( mComputeT ) {
            for( int i=0; i<2; i++ ) {
                const char *event = "T";
                if( mSpeciesTree->hasAlpha() && alphaIdForTS == idX ) 
                    event = "TFD";
                if( mSubOpt || mTriplets ) {
                    // if idX is best, include all secondary
                    if( idX == mBestReceiver[ts][idUsub] ) {
                        BOOST_FOREACH( int id, mVT2all[ts][idUsub] ) 
                            addVertices( graph, pairVertex, qList, idUsub,
                                id, idUother, idX, "T", cost, mHGTCost, trip ); 
                    }
                    // include all best receivers not equal to idX
                    BOOST_FOREACH( int brId, mVT1all[ts][idUsub] ) 
                        if( brId != idX ) 
                            addVertices( graph, pairVertex, qList, idUsub, brId,
                                   idUother, idX, event, cost, mHGTCost, trip );
                } else {
                    vector<int> allBRs = allBestReceivers[ts][idUsub][0];
                    if( allBRs.size()==1 && allBRs.at(0)==idX ) 
                       allBRs = allBestReceivers[ts][idUsub][1];
                    BOOST_FOREACH( int brId, allBRs ) {
                        if( brId != idX ) 
                            addVertices( graph, pairVertex, qList, idUsub, brId,
                                  idUother, idX, event, cost, mHGTCost, trip );
                    }
                }
                 
                //transfers to alpha cost 0
                //transfers from alpha cost hgtCost so nothing to change
                if( mSpeciesTree->hasAlpha() && alphaIdForTS != idX 
                    && alphaIdForTS != -1 )
                {
                    addVertices( graph, pairVertex, qList, idUsub, 
                       alphaIdForTS, idUother, idX, "TTD", cost, 0, trip );
                }

                // switch to other child for second iteration
                idUsub = idUr;
                idUother = idUl;
            }
        }
    }
   
    // Speciation plus loss costs
    BOOST_FOREACH( split, speciesSplits ) {
        if( split.second != -1 ) {
            addVertices( graph, pairVertex, qList, idU, split.first, "SL", 
                         cost, mLossCost, trip );
            addVertices( graph, pairVertex, qList, idU, split.second, "SL", 
                         cost, mLossCost, trip );
        } else {
            // Null costs
                //if( xNode->getInfos().WGD ) {
                    //null event in a WGD node cost L + WGDCost
                    // the null events are always the last ones, 
                    // I call it O and not WGDL because 
                    // I need it for canonicals
                //    addVertices( graph, pairVertex, qList, idU,
                //          split.first, "0", cost, mLossCost+mWGDCost, trip );
                //} else 
            // the null events are always the last ones
            addVertices( graph, pairVertex, qList, idU, split.first, "O", 
                    cost, 0, trip );
        }
    }

    // Transfer plus loss costs
    if( mComputeTL ) {
        if( mSubOpt || mTriplets ) {
            if( idX == mBestReceiver[ts][idU] ) {
                BOOST_FOREACH( int id, mVT2all[ts][idU] ) 
                    addVertices( graph, pairVertex, qList, 
                            idU, id, "TL", cost, mLossCost+mHGTCost, trip ); 
            }
            BOOST_FOREACH( int id, mVT1all[ts][idU] ) 
                if( id != idX ) {
                    addVertices( graph, pairVertex, qList,
                            idU, id, "TL", cost, mLossCost+mHGTCost, trip ); 
                    if( mSpeciesTree->hasAlpha() && alphaIdForTS == idX ) 
                        addVertices( graph, pairVertex, qList,
                            idU, id, "TLFD", cost, mHGTCost, trip ); 
                }
        } else {        
            vector<int> allBRs = allBestReceivers[ts][idU][0];
            if( allBRs.size()==1 && allBRs.at(0)==idX ) 
                allBRs = allBestReceivers[ts][idU][1];
            BOOST_FOREACH( int brId, allBRs ) 
                if( brId != idX ) {
                    addVertices( graph, pairVertex, qList, idU,  
                                 brId, "TL", cost, mLossCost+mHGTCost, trip );
                    if( mSpeciesTree->hasAlpha() && alphaIdForTS == idX )
                        // transfer + loss  from alpha, they cost hgtCost
                        addVertices( graph, pairVertex, qList, 
                                    idU, brId, "TLFD", cost, mHGTCost, trip );
                }
        }
    }

    //transfers to alpha plus loss cost lossCost
    if( mSpeciesTree->hasAlpha() && alphaIdForTS!=idX && alphaIdForTS!=-1 ) 
        addVertices( graph, pairVertex, qList, idU, 
                alphaIdForTS, "TLTD", cost, mLossCost, trip );



    if( mUseILS ) {
        // ILS - there can be multiple paths through the species, unlike
        // the gene tree, where there must be one tree, i.e. no best
        // split is chosen for the species tree.
        // In thre graph, we keep the ambiguity. For the reconciliation,
        // fake ILS nodes collapse to the real parent.
        vector< pair<int,int> > ilsSplits 
                                = mSpeciesTree->getIlsSplits( idX );
        for( size_t i=0; i<ilsSplits.size(); i++ ) {
            if( ilsSplits[i].second == -1 ) {
                addVertices( graph, pairVertex, qList, idU, ilsSplits[i].first,
                    "O", cost, 0, trip );
            } else { 
                // ILS
                if( cladeSplit.first != -1 ) {
                    addVertices( graph, pairVertex, qList, cladeSplit.first, 
                            ilsSplits[i].first, cladeSplit.second, 
                            ilsSplits[i].second, "I", cost, mIlsCost, trip );
                    addVertices( graph, pairVertex, qList, cladeSplit.first, 
                            ilsSplits[i].second, cladeSplit.second,
                            ilsSplits[i].first, "I", cost, mIlsCost, trip );
                }

                // ILS + loss
                addVertices( graph, pairVertex, qList, idU, 
                            ilsSplits[i].first, "IL", cost, 
                            mLossCost+mIlsCost, trip );
                addVertices( graph, pairVertex, qList, idU, 
                             ilsSplits[i].second, "IL", cost, 
                             mLossCost+mIlsCost, trip );
            }
        }
    }
}


/**
 * Construct a DTLGraph from the matrix.
 *
 * Construct events graph by finding all best costs paths 
 * through matrix, considering all events with lowest cost, 
 * which is the matrix cell value.
 *
 * @return DTLGraph
 */
DTLGraph DTLMatrix::constructGraph( 
    bool verbose, ///< print timing information and stats
    bool keepRoot ) ///< use species root
{
    clock_t start = clock();

    // Node specific costs not yet added to this function.
    if( !mFixedCosts ) 
        throw bpp::Exception("DTLMatrix::constructGraph: fixed costs only" );
#ifdef PARALLEL
    // mAllOtherBestReceivers is not thread safe
    throw bpp::Exception("DTLMatrix::constructGraph: "
            "not implemented with parallel matrix computation" );
#endif
  
    bool subopt = false;
    if( mSubOpt || mTriplets )
        subopt = true;
    // create graph 
    DTLGraph graph( subopt, mTriplets, mSpeciesTree, mCladesTrips ); 
    vector<int> *** allBestReceivers = getAllBestReceivers(); 
   
    // initialize the queue, qList, with the root nodes
    vector<DTLGraph::MyGraph::Vertex> rootList 
                                        = getRootNodes( graph, keepRoot );
//    if( verbose ) 
//        cout << "root count: " << rootList.size() << endl;
    vector<DTLGraph::MyGraph::Vertex> qList( rootList ); 

    // main loop - qList is pairs of u,x,cost which are also vertex names
    long vertexCount = 0;
    while( !qList.empty() ) {
        // pop queue
        DTLGraph::MyGraph::Vertex pairVertex = qList.back();
        qList.pop_back();
        vertexCount++;
        createVertices( graph, pairVertex, qList, allBestReceivers );
    } // end qList

    // clean up best receivers
    for( int ts=0; ts<=mMaxTS; ts++ ) { 
        for( int idU=0; idU<mCladeCount; idU++ ) 
            delete [] allBestReceivers[ts][idU];
        delete [] allBestReceivers[ts]; 
    }
    delete [] allBestReceivers;

    if( verbose ) {
        double time = (double) (clock()-start) / CLOCKS_PER_SEC * 1000.0;
        cout << "constructGraph: " << time << " ms, " << vertexCount 
             << " mapping vertices" << endl;
    }

    return graph;
}






///////////////////////////////////////////////////////////////////
//////////////////// DEPRECATED //////////////////////////////////


/**
 * Confirm validity of a solution.  DEPRECATED
 */
void DTLMatrix::checkSolution( 
    int splitNum, int id_u, int z,
    vector< vector<int> > &recAlphaVector,
    double &costRec, bool &canonical, bool &isRec )
{
    // Node specific costs not yet added to this function.
    if( !mFixedCosts ) 
        throw bpp::Exception("DTLMatrix::checkSolution: fixed costs only" );

    int id_x= recAlphaVector[id_u][z];
    vector<int> xSonIds; // = mSpeciesTree->getSplit( id_x );
    vector< pair<int,int> > splits = mSpeciesTree->getSplits( id_x );
    if( splits.size() > 1 ) 
        throw bpp::Exception( "DTLGraph::checkSolution:"
                " multiple species splits not supported" );
    if( splits.size() == 1 ) {
        xSonIds.push_back( splits[0].first );
        if( splits[0].second != -1 )
            xSonIds.push_back( splits[0].second );
    }

    pair<int,int> cladeSplit = mCladesTrips->getCladeSplit( id_u, splitNum );
    int id_u_l = cladeSplit.first;
    int id_u_r = cladeSplit.second;


    int size= recAlphaVector[id_u].size();

    // if x is in Alpha particular event costs 
    // will be applied 
    bool inAlpha=false;
    int xTimeSlice = mSpeciesTree->getTimeSlice( id_x );
    int alphaIdForTS = -1;
    if( mSpeciesTree->hasAlpha() ) {
        alphaIdForTS = mSpeciesTree->getAlphaIdForTS( xTimeSlice );
        if( id_x == alphaIdForTS )
            inAlpha = true; 
    }



    if( z==(size-1) ) {
        if( xSonIds.size() == 0 && id_u_l == -1 ) {
            //C
            MySpeciesNode *xNode = mSpeciesTree->getNodeById(id_x); 
            if( mCladesTrips->mClades.getSpeciesName(id_u) == xNode->getName() )
                isRec=true;
        } else {
            int id_alpha1_ul= recAlphaVector[id_u_l][0];
            int id_alpha1_ur= recAlphaVector[id_u_r][0];
    
            if( xSonIds.size() == 2  ) {
                //S
                if( (id_alpha1_ul==xSonIds[0] && id_alpha1_ur==xSonIds[1]) 
                     || (id_alpha1_ur==xSonIds[0] && id_alpha1_ul==xSonIds[1]))
                    isRec=true;
            }

            bool dupliOrHgt=false;
            if(id_alpha1_ul==id_x && id_alpha1_ur==id_x) {
                isRec=true;
                //if( xNode->getInfos().WGD )
                //    costRec += mWGDCost;
                //else if(!inAlpha)
                    //dupli in alpha cost 0
                    costRec +=(mDupliCost );
                dupliOrHgt=true; //D
            }
        
            if( (id_alpha1_ul==id_x  
                && id_alpha1_ur!=id_x 
                && mSpeciesTree->getTimeSlice(id_alpha1_ur) == xTimeSlice ) 
                || (id_alpha1_ur==id_x  
                    && id_alpha1_ul !=id_x 
                    && mSpeciesTree->getTimeSlice(id_alpha1_ul) == xTimeSlice))
            {
                bool inAlphaXF=false;
                // same time slice as id_x
                // x is in Alpha so particulat event costs 
                // will be applied 
                if( id_alpha1_ur == alphaIdForTS
                    || id_alpha1_ul == alphaIdForTS )
                    inAlphaXF =true; 
                //transfer to alpha cost 0
                if(!(!inAlpha && inAlphaXF) ) 
                    costRec +=(mHGTCost );
                isRec=true;
                dupliOrHgt=true; //T
            }
        
            if( dupliOrHgt
                && recAlphaVector[id_u_l].size() > 1 
                && recAlphaVector[id_u_r].size() > 1
                && xSonIds.size() == 1 
                && recAlphaVector[id_u_l][1] == xSonIds[0] 
                && recAlphaVector[id_u_r][1] == xSonIds[0] )
            {
                //non canonical D or T event
                canonical= false; 
            }
                                                    
        }

    }
    else
    {
        int id_alphaiplus1_u= recAlphaVector[id_u][z+1];
        if( xSonIds.size() == 1  ) {
            if(id_alphaiplus1_u==xSonIds[0]){
                isRec=true;
            }
            int rootClade = mCladesTrips->mClades.getRootClade();
            if( id_u == rootClade && z==0){
                //it starts with a null event
                canonical= false; 
            }
            //if( xNode->getInfos().WGD ) {
                // null event costs loss 
                // + WGDCost in a WGD node
            //    costRec +=(mLossCost+mWGDCost); 
            //}	
        } else if( xSonIds.size() == 2  ) {
            if(id_alphaiplus1_u==xSonIds[0]|| id_alphaiplus1_u==xSonIds[1]){
                isRec=true;
                costRec += mLossCost;
            }
        }

        if( mSpeciesTree->getTimeSlice( id_alphaiplus1_u ) == xTimeSlice ) {
            isRec=true;
            bool inAlphaXF=false;
            // x is in Alpha so particulat event costs 
            // will be applied 
            if( id_alphaiplus1_u == alphaIdForTS )
                inAlphaXF =true; 

            //transfer loss to alpha cost l
            if(!inAlpha && inAlphaXF) 
                costRec += mLossCost ; 
            else if (inAlpha && !inAlphaXF)
                //transfer loss from alpha cost h
                costRec += mHGTCost; 
            else
                costRec += (mHGTCost + mLossCost);	

            if( z+2 < size && xSonIds.size() == 1  ) {
                //vector<int> aSonIds = mSpeciesTree->getSplit( 
                vector< pair<int,int> > aSonIds = mSpeciesTree->getSplits( 
                                            recAlphaVector[id_u][z+1] );
                if( aSonIds.size() == 1 && aSonIds[0].second == -1
                    && recAlphaVector[id_u][z+2] == aSonIds[0].first ) 
                {
                    //non canonical TL event 
                    canonical= false; 
                }
            }
        }
    }
}


// Useful function for comparing two vectors
// Return true if all a vectors are in b.
// Requires that the input vectors are sorted.
template<class T>
bool isIncludedOrEqual(const vector<T> & a, const vector<T> &  b){
    
    if( a.size() > b.size() ) return false;

    size_t cont = 0;
    for( size_t i=0; i<a.size(); i++ ) {
        size_t j = 0;
        while( j < b.size() && b[j] <=a [i] ) {
            if( b[j] == a[i] ) cont++;
            j++;
        }
    }

    if( cont == a.size() )
        return true;
    else
        return false;
}




