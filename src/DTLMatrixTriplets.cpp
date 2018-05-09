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

**/


#include "DTLMatrixTriplets.h"
#include <boost/foreach.hpp>


void DTLMatrixTriplets::computeList( 
    int idUa,
    int idXa,
    int idUb,
    int idXb,
    double lowestCost,
    double curBestEpsCost, 
    double otherCost,
    EventTriplet eventTrip )
{
    // merge lists
    vector<EventTriplet> tripletsA = mMatrixT->getValueSure(idUa,idXa); 
    vector<EventTriplet> tripletsB = mMatrixT->getValueSure(idUb,idXb); 
    for( size_t a=0; a<tripletsA.size(); a++ ) {
        for( size_t b=0; b<tripletsB.size(); b++ ) {
            double cost = otherCost + tripletsA[a].cost 
                        + tripletsB[b].cost;
            if( COST_GREATER( cost, curBestEpsCost ) )
                continue;

            EventTriplet newTrip( tripletsA[a].d + tripletsB[b].d 
                    + eventTrip.d,
                    tripletsA[a].t + tripletsB[b].t + eventTrip.t,
                    tripletsA[a].l + tripletsB[b].l + eventTrip.l, cost );
            mTripletList.push_back( newTrip );
            
            if( COST_GREATER( lowestCost, cost ) ) {
                cout << idUa << "," << idXa << " and " 
                     << idUb << "," << idXb
                     << " lowestCost=" << lowestCost << " > " 
                     << "cost=" << cost << endl;
                throw bpp::Exception ("DTLMatrixTriplets::computeList: "
                        "impossible cost");
            }
        }
    }
}

void DTLMatrixTriplets::computeList( 
    int idUa,
    int idXa,
    double lowestCost,
    double curBestEpsCost, 
    double otherCost,
    EventTriplet eventTrip )
{
    // check costs of new list
    vector<EventTriplet> tripletsA = mMatrixT->getValueSure(idUa,idXa); 
    for( size_t a=0; a<tripletsA.size(); a++ ) {
        double cost = otherCost + tripletsA[a].cost;
        if( COST_GREATER( cost, curBestEpsCost ) )
            continue;

        EventTriplet newTrip( tripletsA[a].d + eventTrip.d,
            tripletsA[a].t + eventTrip.t, tripletsA[a].l + eventTrip.l, cost );
        mTripletList.push_back( newTrip );
        
        if( COST_GREATER( lowestCost, cost ) ) {
            cout << idUa << "," << idXa
                 << " lowestCost=" << lowestCost << " > " 
                 << "cost=" << cost << endl;
            throw bpp::Exception ("DTLMatrixTriplets::computeList: "
                    "impossible cost");
        }
    }
}

void DTLMatrixTriplets::computeTransferCostSub( 
        int optimumSub, ///< BR
        int idUsub,     ///< clade child
        int idUother,     ///< clade child
        double costThisSplit,  ///< current split cost
        DTLMatrixState &state ) ///< various current info
{
    EventTriplet eventTrip;
    if( optimumSub == state.idX ) {
        // transfer all sub-optimal to optimal species node
        BOOST_FOREACH( EventTriplet t, mVTt2[state.timeSlice][idUsub] )
        {
            eventTrip.d = t.d;
            eventTrip.t = t.t + 1;
            eventTrip.l = t.l;
            computeCostList( costThisSplit + state.hgtCost + t.cost, 
                    idUother, state.idX, eventTrip, state ); 
        }
    }

    // transfer all optimal to this node (except self)
    pair<EventTriplet,int> p;
    BOOST_FOREACH( p, mVTt1[state.timeSlice][idUsub]) {
        // don't add if it is this species node's cost
        if( p.second != state.idX ) {
            eventTrip.d = p.first.d;
            eventTrip.t = p.first.t + 1;
            eventTrip.l = p.first.l;
            computeCostList( costThisSplit + state.hgtCost + 
                p.first.cost, idUother, state.idX, eventTrip, state ); 
        }
    }
}

/*
void DTLMatrixTriplets::computeVTs( 
    vector<MySpeciesNode*> &speciesNodesTS, ///< nodes to considier
    DTLMatrixState &state ) ///< various current info
{
    int br = mBestReceiver[state.timeSlice][state.idU];
    double brCost = mBestReceiverCost[state.timeSlice][state.idU];
    int br2 = mSecondBestReceiver[state.timeSlice][state.idU];
    double br2cost = mSecondBestReceiverCost[state.timeSlice][state.idU]; 


    vector< pair<EventTriplet,int> > &VTt1 = mVTt1[state.timeSlice][state.idU]; 
    for( size_t z = 0; z < speciesNodesTS.size(); z++ ) {
        int idX = speciesNodesTS[z]->getId();
        bool inVTt1all = false;
        bool inVTt2all = false;
        for( size_t t=0; t<mAllTripletsTS[z].size(); t++ ) {
            EventTriplet trip = mAllTripletsTS[z][t];

            // add to list if within mEpsilon
            if( !COST_GREATER( trip.cost, brCost + mEpsilon )  ) {
            //if( trip.cost <= brCost + mEpsilon + COST_DIFF ) {
                // check if in the list
                pair<EventTriplet,int> p;
                bool foundIt = false;
                for( size_t i=0; i<VTt1.size(); i++ ) 
                {
                    p = VTt1[i];
// This list should be small. If not, consider more efficient methods.

                    // Pareto
                    if( p.first.d == trip.d && p.first.t == trip.t
                       && p.first.l == trip.l ) 
// non-pareto
//if( abs(p.first.cost-trip.cost) <= COST_DIFF  ) 
                    {
                        foundIt = true;
                        // mark as multiple
                        VTt1[i].second = -1;
                        break;
                    }
                }
                if( !foundIt ) {
                    // add it
                    pair<EventTriplet,int> p = make_pair(trip,idX);
                    VTt1.push_back( p );
                }
                if( !inVTt1all ) {
                    mVT1all[state.timeSlice][state.idU].push_back( idX );
                    inVTt1all = true;
                }
            } else if( br2 != -1  // second best is not optimal
                       && idX != br 
                       && !COST_GREATER( trip.cost, br2cost + mEpsilon ) )
                       //&& trip.cost <= br2cost + mEpsilon + COST_DIFF )
            {
                // add to second optimums list (mVT2) if within mEpsilon,
                mVTt2[state.timeSlice][state.idU].push_back( trip );
                if( !inVTt2all ) {
                    mVT2all[state.timeSlice][state.idU].push_back( idX );
                    inVTt2all = true;
                }
            }
        }
    }

    // stupid sort - fine for small lists sizes
    bool changes = true;
    while( changes ) {
        changes = false;
        for( size_t i = 1; i<VTt1.size(); i++ ) {
            if( VTt1[i] < VTt1[i-1] ) {
                // swap
                pair<EventTriplet,int> tmp = VTt1[i-1];
                VTt1[i-1] = VTt1[i];
                VTt1[i] = tmp;
                changes = true;
            }
        }
    }

    // remove duplicates from second optimums 
    cleanTriplets( mVTt2[state.timeSlice][state.idU], state );

    // remove duplicates
    removeDuplicates( mVT1all[state.timeSlice][state.idU] );
    removeDuplicates( mVT2all[state.timeSlice][state.idU] );
}
*/
/**
 * Triplet pair sort function. Sorts by cost.
 *
 * @return true if a < b
 */ 
bool sortByCost(
    pair< EventTriplet ,int> a, ///< first 
    pair< EventTriplet ,int> b ) ///< second
{
    return( a.first.cost < b.first.cost ); 
}

/**
 * Looks for trip in list, which is sorted by cost.
 *
 * @return true if triplet is in list
 */ 
bool DTLMatrixTriplets::tripInList(
        EventTriplet trip,
        vector< pair<EventTriplet,int> > list )
{
    // search backwards in list while cost is >=
    for( int i= (int) list.size() - 1; i>=0; i-- ) {
        if( trip.cost == list[i].first.cost ) {
            if( list[i].first.d == trip.d && list[i].first.t == trip.t
                && list[i].first.l == trip.l ) 
            {
                return true;
            }
        } else if( trip.cost > list[i].first.cost ) {
            return false;
        }
    }
    return false;
}


void DTLMatrixTriplets::computeVTs( 
    bool noSub, ///< no subdivision if true
    bool recompute, ///< replace existing values
    int timeSlice, ///< time slice (not used if noSub is true)
    int idU, ///< gene id
    vector<int> &speciesNodeIdsTS ) ///< ids of nodes to considier
{
    // create cost/id pairs from the optVector
    vector< pair< EventTriplet ,int> > tripletPairs;
    vector<int> idXs; // id lists ordered the same as speciesNodes
    for( size_t nodeIndex=0; nodeIndex<speciesNodeIdsTS.size(); nodeIndex++ ) {
        int idX = speciesNodeIdsTS[nodeIndex];
        idXs.push_back( idX );
        for( size_t t=0; t<mAllTripletsTS[nodeIndex].size(); t++ ) {
            EventTriplet trip = mAllTripletsTS[nodeIndex][t];
            tripletPairs.push_back( make_pair( trip, idX ) );
        }
    }

   if( !noSub && recompute ) {
        mVTt1[timeSlice][idU].clear();
        mVTt2[timeSlice][idU].clear();
        mVT1all[timeSlice][idU].clear();
        mVT2all[timeSlice][idU].clear();
    }

    // sort by cost 
    sort( tripletPairs.begin(), tripletPairs.end(), sortByCost );

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
            ts = mSpeciesTree->getTimeSlice( speciesNodeIdsTS[i] );
            brCost = mBestReceiverCost[ts][idU];
            if( recompute ) {
                mVTt1[ts][idU].clear();
                mVTt2[ts][idU].clear();
                mVT1all[ts][idU].clear();
                mVT2all[ts][idU].clear();
            }
        }

        bool firstTrip = true;
        EventTriplet curTrip;
        EventTriplet lastBestTrip;
        int bestId = -1;
        int costCount = 0;
        bool costGreater = false;
        for( size_t j=0; j<tripletPairs.size(); j++ ) {
            if( noSub && !mComparable[idXs[i]][tripletPairs[j].second] ) 
                continue;
     
            if( COST_GREATER( tripletPairs[j].first.cost, brCost + mEpsilon ) )
            {
                if( noSub ) 
                    break; // seen everything
               
                if( COST_GREATER( tripletPairs[j].first.cost, 
                                  br2cost + mEpsilon ))
                {
                    costGreater = true;
                    break; // second best finished, too
                } else if( br2 != -1  // second best is not optimal
                           && tripletPairs[j].second != br  )
                {
                    // second best receiver

                    // add cost to second optimums list (mVT2) 
                    // if within mEpsilon
                    if( !COST_EQUAL( tripletPairs[j].first.cost, 
                                     curTrip.cost ) ) 
                    {
                        mVTt2[ts][idU].push_back( tripletPairs[j].first );
                        curTrip = tripletPairs[j].first;
                    }
                    // all ids with a cost within epsilon for second best
                    // (except best receiver)
                    if( recompute )
                        mVT2all[ts][idU].push_back( tripletPairs[j].second );
                }
            } else { 
                // best receiver (keep one pair for each cost)
                if( firstTrip || !tripInList( tripletPairs[j].first,
                                              mVTt1[ts][idU] ) )
                {
                    if( firstTrip )
                        firstTrip = false;
                    else {
                        // cost is different, add the id for previous cost
                        pair<EventTriplet,int> p = make_pair(curTrip,bestId);
                        mVTt1[ts][idU].push_back( p );
                    }
                    bestId = tripletPairs[j].second;
                    curTrip = tripletPairs[j].first;
                    lastBestTrip = curTrip;
                    costCount = 1;
                } else {
                    // An id of -1 indicates that there are multiple species
                    // nodes with this cost.
                    costCount++;
                    bestId = -1;
                }
                // all ids with a cost within epsilon
                if( recompute )
                    mVT1all[ts][idU].push_back( tripletPairs[j].second );
            }
                // best receiver (keep one pair for each cost)
        }
        if( costCount > 0 ) 
            mVTt1[ts][idU].push_back( make_pair( lastBestTrip, bestId ) );

        // sort mVTt1[ts][idU] using triplet comparison
        // stupid sort - fine for small lists sizes
        vector< pair<EventTriplet,int> > &VTt1 = mVTt1[ts][idU]; 
        bool changes = true;
        while( changes ) {
            changes = false;
            for( size_t i = 1; i<VTt1.size(); i++ ) {
                if( VTt1[i] < VTt1[i-1] ) {
                    // swap
                    pair<EventTriplet,int> tmp = VTt1[i-1];
                    VTt1[i-1] = VTt1[i];
                    VTt1[i] = tmp;
                    changes = true;
                }
            }
        }

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

/*
void DTLMatrixTriplets::computeVTsNoSub( 
    vector<MySpeciesNode*> &speciesNodesTS, ///< nodes to considier
    DTLMatrixState &state ) ///< various current info
{
}
*/


/**
 * Add all possible transfer costs within suboptimal epsilon to
 * allTriplets.
 */
void DTLMatrixTriplets::computeTransferLossCostVT(
    int idx, ///< index into mAllTripletsTS
    DTLMatrixState &state ) ///< various current info
{
    if( mBestReceiver[state.timeSlice][state.idU] == -1 ) 
        return; // no best receiver

    int nonAlphaLoss = 0; // count a loww if this isn't alpha
    double otherCosts = state.hgtCost;
    if ( !mSpeciesTree->hasAlpha() || !mSpeciesTree->isAlpha( state.idX ) ) {
        nonAlphaLoss = 1;
        otherCosts += state.lossCost;
    }
   
    // for all VT, do a cost combination
    if( state.idX == mBestReceiver[state.timeSlice][state.idU] ) {
        // transfering to optimal species node
        BOOST_FOREACH( EventTriplet trip, mVTt2[state.timeSlice][state.idU] ) {
            double cost = trip.cost + otherCosts;
            EventTriplet newTrip( trip.d, trip.t+1, 
                                    trip.l + nonAlphaLoss, cost );
            if( !COST_GREATER( cost, state.curBestEpsCost) ) 
                mAllTripletsTS[idx].push_back( newTrip  );
        }

    } 

    // transfer all optimal to this node (except self)
    pair<EventTriplet,int> p;
    BOOST_FOREACH( p, mVTt1[state.timeSlice][state.idU]) {
        // don't add if it is this species node's cost
        if( p.second != state.idX ) {
            double cost = p.first.cost + otherCosts;
            if( !COST_GREATER( cost, state.curBestEpsCost) ) {
                EventTriplet newTrip( p.first.d, p.first.t+1,
                        p.first.l + nonAlphaLoss, cost );
                mAllTripletsTS[idx].push_back( newTrip );
            }
        }
    }

    // Transfer from alpha.
    if ( mSpeciesTree->hasAlpha() && !mSpeciesTree->isAlpha( state.idX ) ) {
        BOOST_FOREACH( EventTriplet alphaTrip, mAlphaTriplets ) {
            double cost = alphaTrip.cost + state.lossCost;
            if( !COST_GREATER( cost, state.curBestEpsCost) ) {
                EventTriplet newTrip( alphaTrip.d, alphaTrip.t, alphaTrip.l+1,
                               cost );
                mAllTripletsTS[idx].push_back( newTrip );
            }
        }
    }
}


void DTLMatrixTriplets::intersectionLines(
        double &intX,
        double &intY,
        double x1,
        double y1,
        double x2,
        double y2,
        EventTriplet line )
{    
    if( x1==x2 ) {
        intX = x1;
        intY = -(double)x1*line.d/line.l-(double)line.t/line.l;    
    } else {        
        double slope1 = (y2-y1)/(x2-x1);
        double in1 = (double)(y1*x2-y2*x1)/(x2-x1);
        if (line.l==0 ) 
            intX=-(double)line.t/line.d;
        else {          
            double slope2 = -(double)line.d/line.l;
            double in2 = -(double)line.t/line.l;      
            intX = (double)(in2-in1)/(slope1-slope2);
        }

        intY = intX*slope1+in1;
    }
} 


int DTLMatrixTriplets::value(
        EventTriplet line,
        double x,
        double y )
{
    if( line.d*x+line.l*y+line.t < 0) 
        return -1;

    if( line.d*x+line.l*y+line.t == 0) 
        return 0;

    return 1;
}

bool DTLMatrixTriplets::intersectionPolygonLine(
        list<double> &polygonX, 
        list<double> &polygonY,
        EventTriplet line )
{
    list<double> polyNewX;
    list<double> polyNewY;                
    list<double>::iterator iX = polygonX.begin();
    list<double>::iterator iY = polygonY.begin();
    double xFirst = *iX;
    double yFirst = *iY;
    double x = xFirst;
    double y = yFirst;
    if( value(line,x,y) <= 0 ) {
        polyNewX.push_back(x);
        polyNewY.push_back(y);
    }
    iX++;
    iY++;

    while( iX != polygonX.end() && iY != polygonY.end() ) {
        double xN = *iX;
        double yN = *iY;                           
        if( ( value(line,x,y)< 0 && value(line,xN,yN)> 0 ) 
            || ( value(line,x,y)>0 && value(line,xN,yN)<0 ) )
        {
            double xInt,yInt;                              
            intersectionLines( xInt, yInt, x, y, xN, yN, line );
            polyNewX.push_back(xInt);
            polyNewY.push_back(yInt);
        }              

        if( value(line,xN,yN) <=0 ) {
            polyNewX.push_back(xN);
            polyNewY.push_back(yN);
        }

        x = xN;
        y = yN;
        iX++;
        iY++;                      
    }

    if( (value(line,x,y) < 0 && value(line,xFirst,yFirst)>0 ) 
        || ( value(line,x,y)>0 && value(line,xFirst,yFirst)< 0 ))
    {
        if( polygonX.size() > 2 )  {
            double xInt, yInt;                      
            intersectionLines( xInt, yInt, x, y, xFirst, yFirst, line );
            polyNewX.push_back( xInt );
            polyNewY.push_back( yInt );
        }
    }               

    if( polyNewX.empty() ) 
        return false;
    else {
        polygonX = polyNewX;
        polygonY = polyNewY;  
        return true;
    }
}

/**
 * Remove non-minimal triplets, <d,t,l> such that another
 * triplet <d',t',l'> is in the list and d>=d' && t>=t' && l>=l'.
 */
void DTLMatrixTriplets::cleanTriplets( 
    vector<EventTriplet> &v, ///< the list
    DTLMatrixState &state ) ///< various current info
{
    if( v.size() <= 1 )
        return;

// non-pareto
/*
Non pareto sort needs to change in MyMatrixT.h also
    sort( v.begin(), v.end() ); 
    size_t pos=0;
    for( size_t i=1; i<v.size(); i++ ) {
        if( abs(v[i].cost-v[pos].cost) > COST_DIFF ) {
            pos++;
            if( pos != i )
                v[pos] = v[i];
        }
    }
    size_t p = v.size() - 1;
    while( p > pos ) {
        v.pop_back();
        p--;
    }
*/
    // pareto
    // sorts on the number of d's first, then t's, then l's
    sort( v.begin(), v.end() ); 


    size_t curPos = 0;
    for( size_t pos=1; pos<v.size(); pos++ ) {
        // scan backwards looking for duplicates and minimal
        bool minimal = true;
        size_t i = curPos+1;
        while( i > 0 ) {
            i--;
            // v[i].d <= v[pos].d  because of sorting
            if( v[i].t <= v[pos].t && v[i].l <= v[pos].l )
            {
                // not minimal
                minimal = false;
                break;
            }
        }
        if( minimal ) {
            curPos++; // update end of good list
            if( curPos != pos ) 
                v[curPos] = v[pos]; // move to the end of the good list
        }
    }
    // remove non-minimals not overwritten 
    size_t pos = v.size() - 1;
    while( pos > curPos ) {
        v.pop_back();
        pos--;
    }


    if( 0 ) { //mParetoMod == 2 || mParetoMod == 3 ) {
/* trim to just optimal for testing
        vector<EventTriplet> newV;
        double minCost;
        for( size_t i=0; i<v.size(); i++ ) {
            if( i==0 || v[i].cost < minCost )
                minCost = v[i].cost;
        }
        for( size_t i=0; i<v.size(); i++ ) {
            if( abs(v[i].cost-minCost) < COST_DIFF )
                newV.push_back( v[i] );
        }
        v = newV;
*/

        // keep dtl, s.t. there is no d't'l'
        // d'.(1+nD)(delta/tau)+l'.(1+nL)(lambda/tau)+t' 
        //   < d.(1-nD)(delta/tau)+l.(1-nL)(lambda/tau)+t  
        vector<double> upperRatioCosts( v.size() );
        vector<bool> keep ( v.size(), true );
        vector<EventTriplet> newV;
        for( size_t i=0; i<v.size(); i++ ) {
            double lowerRatioCost = v[i].d*mNDl + v[i].l*mNLl + v[i].t;
            for( size_t j=0; j<v.size(); j++ ) {
                if( i==0 ) {
                    upperRatioCosts[j] = v[j].d*mNDu + v[j].l*mNLu + v[j].t;
                }
                if( i!=j && 
                        COST_GREATER( lowerRatioCost, upperRatioCosts[j] ) ) 
                {
                    keep[i] = false;
                    if( i>0 )
                        break; // else continue upperRatioCosts calculation
                }
            }
            if( keep[i] ) 
                newV.push_back( v[i] );
        }
        v = newV;

    }

    
    if( 0 ) { //mParetoMod == 1 || mParetoMod == 3 || mParetoMod == 4 ) {
        vector<bool> keep ( v.size(), true );
        vector<EventTriplet> newV;
        //d\delta+l\lambda+t\tau < d'\delta+l'\lambda+(\delta+\lambda)t' 
        for( size_t i=0; i<v.size(); i++ ) {
            double iCost = v[i].d*state.dupliCost + v[i].l*state.lossCost
                         + v[i].t*(state.dupliCost+state.lossCost);
            for( size_t j=0; j<v.size(); j++ ) {
                if( i==j )
                    continue;
                // old critera
                if( v[i].d == v[j].d && v[i].t == v[j].t &&
                    v[i].l == v[j].l )
                {   
                    throw bpp::Exception("DTLMatrixTriplets::cleanTriplets:"
                            " Shouldn't have duplicates here.");
                    // only keep the first instance of this triplet
                    if( i > j ) {
                        keep[i] = false;
                        break;
                    }
                } else if( v[i].d >= v[j].d && v[i].t >= v[j].t &&
                           v[i].l >= v[j].l  ) 
                {
                    throw bpp::Exception("DTLMatrixTriplets::cleanTriplets:"
                            " Should be cleaned here.");
                    // not minimal
                    keep[i] = false;
                    break;
                }

                if( COST_GREATER( iCost, v[j].cost ) ) {
//cout << state.idU << "," << state.idX 
//    << " ============= " << v[j] << "=" << v[j].cost 
//    << " < " << v[i] << "=" << iCost << endl;
                    keep[i] = false;
                    break;
                }
            }
            if( keep[i] ) 
                newV.push_back( v[i] );
        }
        v = newV;
    }

    if( 0 ) { //mParetoMod == 4 ) {
        vector<bool> keep ( v.size(), true );
        vector<EventTriplet> newV;
        //d\delta+l\lambda+t\tau < d'\delta+l'\lambda+(\delta+\lambda)t' 
        for( size_t i=0; i<v.size(); i++ ) {
            double iCost = (v[i].d + v[i].l)*state.lossCost 
                         + v[i].t*state.hgtCost;
            for( size_t j=0; j<v.size(); j++ ) {
                if( i==j )
                    continue;

                if( COST_GREATER( iCost, v[j].cost ) ) {
//cout << state.idU << "," << state.idX 
//    << " ============= " << v[j] << "=" << v[j].cost 
//    << " < " << v[i] << "=" << iCost << endl;
                    keep[i] = false;
                    break;
                }
            }
            if( keep[i] ) 
                newV.push_back( v[i] );
        }
        v = newV;
    }


    //if ( mParetoMod == 5 ) {
    if ( mParetoMod == 1 ) {
        vector<bool> keep ( v.size(), true );  
        vector<EventTriplet> newV;       
        for( size_t i=0; i<v.size(); i++ ) {  
            list<double> polyX = mPolygonX;
            list<double> polyY = mPolygonY;
            for( size_t j=0; j<v.size(); j++ ) {
                if( j != i ) {                  
                    EventTriplet e( v[i].d-v[j].d,
                                    v[i].t-v[j].t,
                                    v[i].l-v[j].l,
                                    0);
                    if( !intersectionPolygonLine( polyX, polyY, e) ) {
                        keep[i]=false;
                        break;
                    }
                }              
            }
            if( keep[i] ) 
                newV.push_back( v[i] );            
        }
        v = newV;
    }    

    //if( mParetoMod == 6 ) {
    if( mParetoMod == 3 ) {
        vector<bool> keep ( v.size(), true );
        vector<EventTriplet> newV;
        for( size_t i=0; i<v.size(); i++ ) {             
            for( size_t j=0; j<v.size(); j++ ) {
                if( i==j )
                    continue;
                if( (v[j].t<v[i].t) 
                    && (v[j].d-v[i].d <= v[i].t-v[j].t) 
                    && (v[j].l<=v[i].l) ) 
                {
                    keep[i] = false;
                    break;
                }   
            }
            if( keep[i] ) 
                newV.push_back( v[i] );
        }

        v = newV;
    } 
}


/**
* Get event counts for each event.
*/
void DTLMatrixTriplets::eventCounts( 
    const char *event, ///< an event
    int &d,  ///< duplications for event
    int &t,  ///< transfers for event
    int &l ) ///< losses for event
{
    d = 0;
    t = 0;
    l = 0;
    if( !strcmp( event, "D" ) || !strcmp( event, "DD") )
        d = 1;
    else if( !strcmp( event, "T" ) || !strcmp( event, "TFD") )
        t = 1;
    else if( !strcmp( event, "TL" ) 
            || !strcmp( event, "TLFD" )
            || !strcmp( event, "TLTD" ) ) 
    {
        t = 1;
        l = 1;
    } else if( !strcmp( event, "SL" ) )
        l = 1;
    else if( !strcmp( event, "S" ) 
            || !strcmp( event, "TTD" ) 
            || !strcmp( event, "O" ) ) 
    {
        // no counts
    } else {
        cerr << "got event " << event << endl;
        throw bpp::Exception( "DTLMatrixTriplets::eventCounts: unknown event" );
    }
            
}


/**
 * Add graph vertices if computed costs match best costs.
 *
 * Create an event vertex off of the given vertex z and create
 * pair vertex from the given ids. The cost of this event is
 * computed. If the cost matches the current cost, add it to
 * the graph.
 */
void DTLMatrixTriplets::addVertices( 
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
	vector<EventTriplet> trips;
    int addD = 0;
    int addT = 0;
    int addL = 0;
    trips = mMatrixT->getValue(u,x);
    BOOST_FOREACH( EventTriplet t, trips )
        values.push_back( t.cost );
    eventCounts( event, addD, addT, addL );

    for( size_t i=0; i<values.size(); i++ ) {
        // check if this cost matches the matrix cost
	    double cost = values[i] + otherCost; 
	    if( !COST_EQUAL(cost, curCost) ) 
	        continue;

        // If triplets, make sure counts match
        if( trips[i].d+addD != zTrip.d 
                || trips[i].t+addT != zTrip.t
                || trips[i].l+addL != zTrip.l) 
        {
            continue;
        }

	    DTLGraph::MyGraph::Vertex eventVertex 
	                = graph.addEventVertex( z, event );
        DTLGraph::MyGraph::Vertex pairVertex;
        bool isNew = graph.addPairVertex( eventVertex, u, x, values[i], 
                        trips[i].d, trips[i].t, trips[i].l, pairVertex );
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
void DTLMatrixTriplets::addVertices( 
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
    vector<EventTriplet> trips1;
    vector<EventTriplet> trips2;
    int addD = 0;
    int addT = 0;
    int addL = 0;
    trips1 = mMatrixT->getValue(u1,x1);
    trips2 = mMatrixT->getValue(u2,x2);
    BOOST_FOREACH( EventTriplet t, trips1 )
        values1.push_back( t.cost );
    BOOST_FOREACH( EventTriplet t, trips2 )
        values2.push_back( t.cost );
    eventCounts( event, addD, addT, addL );
    
    // create an event vertex and 2 pair vertices for each combination of costs
    for( size_t i=0; i<values1.size(); i++ ) {
        for( size_t j=0; j<values2.size(); j++ ) {
            // check if this cost matches the matrix cost
            double cost = values1[i] + values2[j] + otherCost; 
            if( !COST_EQUAL(cost, curCost) ) 
                continue;

            // If triplets, make sure counts match
            if( trips1[i].d+trips2[j].d+addD != zTrip.d 
                    || trips1[i].t+trips2[j].t+addT != zTrip.t
                    || trips1[i].l+trips2[j].l+addL != zTrip.l )
            {
                continue;
            }

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
                                values2[j], trips2[j].d, trips2[j].t, 
                                trips2[j].l, pairVertex );
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
                isNew = graph.addPairVertex( eventVertex, u1, x1, 
                                values1[i], trips1[i].d, trips1[i].t, 
                                trips1[i].l, pairVertex ); 
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
			} else {
                isNew = graph.addPairVertex( eventVertex, u1, x1, 
                                values1[i], trips1[i].d, trips1[i].t, 
                                trips1[i].l, pairVertex ); 
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
                isNew = graph.addPairVertex( eventVertex, u2, x2, 
                                values2[j], trips2[j].d, trips2[j].t, 
                                trips2[j].l, pairVertex ); 
                if( isNew ) // add to queue if new
                    qList.push_back( pairVertex );
            }
		}
    }
}


/**
 * Get the root nodes of the reconciliation graph.
 *
 * The root nodes are all of the species nodes with minimal cost
 * at the clade root.
 *
 * @return root nodes
 */
vector<DTLGraph::MyGraph::Vertex> DTLMatrixTriplets::getRootNodes(
    DTLGraph &graph, ///< the graph
    bool keepRoot ) ///< use species root only
{
    /*
    vector<DTLGraph::MyGraph::Vertex> rootList;
	double bestCost = getBestCost() + mEpsilon;
    int rootClade = mCladesTrips->mClades.getRootClade();
    for ( int j=0; j<mSTnodes; j++ ) {	
		vector<double> values;
		vector<EventTriplet> trips = mMatrixT->getValueSure(rootClade,j);
        BOOST_FOREACH( EventTriplet t, trips ) 
            values.push_back( t.cost );
        for( size_t i=0; i<values.size(); i++ ) {
            if( values[i] <= bestCost + COST_DIFF ) {
                DTLGraph::MyGraph::Vertex vertex 
                    = graph.createPairVertex( rootClade, j, values[i],
                                    trips[i].d, trips[i].t, trips[i].l );
                rootList.push_back( vertex );
        	}
        }
    }
    */


    // make a list of all triplets
	double bestCost = getBestCost() + mEpsilon;
    int rootClade = mCladesTrips->mClades.getRootClade();
	vector<EventTriplet> allTrips;
    for ( int j=0; j<mSTnodes; j++ ) {	
		vector<EventTriplet> trips = mMatrixT->getValueSure(rootClade,j);
        BOOST_FOREACH( EventTriplet t, trips ) {
            if( !COST_GREATER( t.cost, bestCost ) ) 
                allTrips.push_back( t );
        } 
    }

    // clean entire list
    DTLMatrixState state;
    // assumes fixed costs
    state.dupliCost = mDupliCost;
    state.hgtCost = mHGTCost;
    state.lossCost = mLossCost;
    state.idU = rootClade;
    state.idX = -1;  // only used for debug printing
    cleanTriplets( allTrips, state );


    // get roots
    vector<DTLGraph::MyGraph::Vertex> rootList;
    int j=0;
    if( keepRoot )
        j = mSpeciesTree->getRootNode()->getId();
    for( ; j<mSTnodes; j++ ) {	
		vector<EventTriplet> trips = mMatrixT->getValueSure(rootClade,j);
        BOOST_FOREACH( EventTriplet t, trips ) {
            BOOST_FOREACH( EventTriplet tall, allTrips ) {
                if( t == tall ) {
//cout << "ROOT x=" << j << " trip=" << t << endl;
                    DTLGraph::MyGraph::Vertex vertex 
                        = graph.addRoot( rootClade, j, t.cost, t.d, t.t, t.l );
                    rootList.push_back( vertex );
                    break;
                }
        	}
        }
        if( keepRoot )
            break;
    }

    return rootList;
}

/**
 * Delete data structures for triplets calculations.
 */
void DTLMatrixTriplets::cleanUp() {
    if( mVTt1 != NULL ) {
        for( int ts=0; ts<=mMaxTS; ts++) {
            delete [] mVTt1[ts];
            delete [] mVTt2[ts];
            delete [] mVT1all[ts];
            delete [] mVT2all[ts];
        }
        delete [] mVTt1;
        delete [] mVTt2;
        delete [] mVT1all;
        delete [] mVT2all;
        mVTt1 = NULL;
    }


    if( mMatrixT != NULL ) {
        delete mMatrixT;
        mMatrixT = NULL;
    }
}

/**
 * Initialize data structures for suboptimal calculations.
 */
void DTLMatrixTriplets::initCalculation() 
{
#ifdef SPEED
    // triplets are not fast
    throw bpp::Exception( "DTLMatrix: triplets called with SPEED defined" );
#endif
    mTriplets = true;
    // remove the previous calculation if it exists
    cleanUp(); 

    mUseBestSplits = false;
    mMatrixT = new MyMatrixT();
    try{	
        // we can have at most 3GTnodes different clades all 
        // over all possible rootings of G
        mMatrixT->setDim(mCladeCount,mSTnodes); 
    }
    catch (exception& genericException){
        throw bpp::Exception("-----Low memory!!! The dataset is too big------");
    }

    mVTt1 = new vector< pair<EventTriplet,int> > *[mMaxTS+1];
    mVTt2 = new vector<EventTriplet> *[mMaxTS+1];
    mVT1all = new vector<int> *[mMaxTS+1];
    mVT2all = new vector<int> *[mMaxTS+1];
    for( int timeSlice=0; timeSlice<=mMaxTS; timeSlice++ ) {
        mVTt1[timeSlice] = new vector< pair<EventTriplet,int> > [mCladeCount];
        mVTt2[timeSlice] = new vector<EventTriplet> [mCladeCount];
        mVT1all[timeSlice] = new vector<int> [mCladeCount];
        mVT2all[timeSlice] = new vector<int> [mCladeCount];
    }
}


DTLMatrixTriplets::DTLMatrixTriplets( 
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
    double ilsCost,     ///< cost of an ils event 
    int paretoMod,      ///< modification to pareto optimal to use
    double nD,          ///< duplicate noise fraction
    double nL,         ///< loss noise fraction
    double nDL)         ///< DL noise fraction
    : DTLMatrix( speciesTree, cat, WGDCost, fixedCosts, computeT, computeTL,
               dupliCost, hgtCost, lossCost, maxTS, weight, useBestSplits,
               ilsCost )
{
    mMatrixT = NULL;
    mVTt1 = NULL;

    // 5-5-15 paretoMod 5->1, 0->2, 6->3 (commenting out the others)
    if( paretoMod < 1 || paretoMod > 3 )
        throw bpp::Exception("DTLMatrixTriplets: invalid mParetoMod option" );

    if( !fixedCosts ) {
        // not implmented, but can be added if needed
        throw bpp::Exception( "DLTMatrixTriplets called with non-fixed costs" );
    }

    mParetoMod = paretoMod;

    if( paretoMod == 1 ) {
    //if( paretoMod == 2 || paretoMod == 3 || paretoMod == 5) {
        if( nD < 0 || nD > 1 ) 
            throw bpp::Exception( "DTLMatrixTriplets: nD out of bounds [0,1]" );
        if( nL < 0 || nL > 1 ) 
            throw bpp::Exception( "DTLMatrixTriplets: nL out of bounds [0,1]" );

        mNDu = (1+nD)*dupliCost/hgtCost;
        mNDl = (1-nD)*dupliCost/hgtCost;
        mNLu = (1+nL)*lossCost/hgtCost;
        mNLl = (1-nL)*lossCost/hgtCost;
    //}
    //if( paretoMod == 5 ) {
        if( nDL < 0 || nDL > 1 ) 
            throw bpp::Exception( "DTLMatrixTriplets: nDL out of bounds [0,1]" );
        double nDLu  = (1+nDL)*lossCost/dupliCost; 
                        //(1+mNDl)*lossCost/dupliCost;
        double nDLl  = (1-nDL)*lossCost/dupliCost; 
                        //(1-mNDl)*lossCost/dupliCost;

        mPolygonX.push_back(mNDl);
        mPolygonY.push_back(mNLl);
        if( mNLl != mNLu ) {
            mPolygonX.push_back(mNDl);
            mPolygonY.push_back(mNLu);
        }
        if( mNLl != mNLu && mNDl != mNDu ) {
            mPolygonX.push_back(mNDu);
            mPolygonY.push_back(mNLu);
        }
        if( mNDl != mNDu ) {
            mPolygonX.push_back(mNDu);
            mPolygonY.push_back(mNLl);            
        }
        EventTriplet e1(nDLl*100,0,-100,0);
        EventTriplet e2(-nDLu*100,0,100,0);
        intersectionPolygonLine( mPolygonX, mPolygonY, e1 );
        intersectionPolygonLine( mPolygonX, mPolygonY, e2 );
        if( mPolygonX.empty() )
            throw bpp::Exception( "DTLMatrixTriplets:"
                    " The input cost range is empty." );
    }
}

DTLMatrixTriplets::~DTLMatrixTriplets() 
{
    cleanUp();
}


