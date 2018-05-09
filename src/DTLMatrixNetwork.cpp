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


#include <iostream>
#include <limits>

#include <Bpp/Text/TextTools.h>
#include <boost/foreach.hpp>

#include "DTLMatrixNetwork.h"

#ifdef SPEED
    #define UPDATE_OPT(optCost, bestSplit, cost, optimumSub, idX, c) if( cost < optCost) optCost = cost;
#else
    #define UPDATE_OPT(optCost, bestSplit, cost, optimumSub, idX, c); updateOpt(optCost, bestSplit, cost, optimumSub, idX, c);
#endif



/**
 * Calculate null cost for current matrix cell.
 */ 
void DTLMatrixNetwork::computeNullCost( 
        int idXl,     ///< species child
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    EventTriplet eventTrip;

    double otherCosts = 0;

    if(mSpeciesTree->getNodeById(idXl)->getInfos().secondaryFather!=NULL){
    	if( mSpeciesTree->getNodeById(idXl)->getInfos().secondaryFather->getId()== mSpeciesTree->getNodeById(idXl)->getFather()->getId() ){
        	otherCosts = state.hgtCost;  //needed of the switching cost
     	}
     }   	
        	
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



/**
 * Compute transfer cost with loss.
 *
 * @return true if computed cost is lower than opt
*/

bool DTLMatrixNetwork::computeTransferLossCost( 
        int idXt,     ///< species child
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{
    bool changed = false;
	EventTriplet eventTrip;

    double cost = computeCostList( state.hgtCost + state.lossCost, state.idU, idXt, 
                                    eventTrip, state ) ;	
                                    	
    UPDATE_OPT( optCost, bestSplit, cost, idXt, -1, 'x' );
        
    /*double bestRecCost = mBestReceiverCost[state.timeSlice][state.idU]; 
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
    */
    

    return changed;
    
    
}



/**
 * Compute transfer cost for the current split.
 */
void DTLMatrixNetwork::computeTransferCost( 
        int idUl,     ///< clade child
        int idUr,     ///< clade child
        int idXl,     ///< species child
        int idXr,     ///< species child
        double costThisSplit,  ///< current split cost
        DTLMatrixState &state, ///< various current info
        double &optCost, ///< current optimal cost
        BestSplit &bestSplit ) ///< current best split
{

	int idSecFatherSonXl =-2 ;
        if(mSpeciesTree->getNodeById(idXl)->getInfos().secondaryFather!=NULL)
        	idSecFatherSonXl =mSpeciesTree->getNodeById(idXl)->getInfos().secondaryFather->getId() ;
        	
//        int idSecFatherSonXr =-2;
//        if(mSpeciesTree->getNodeById(idXr)->getInfos().secondaryFather!=NULL)
//        	idSecFatherSonXr=mSpeciesTree->getNodeById(idXr)->getInfos().secondaryFather->getId() ;
          
		  
		int nodeWhereToTransfer;
		              
		if(idSecFatherSonXl==state.idX) //one of the two has to happen
			nodeWhereToTransfer= idXl;
		else
			nodeWhereToTransfer= idXr;	
			
			
		double cost = costThisSplit 
			+ mMatrix.getValueSure(idUl,state.idX) 
			+ mMatrix.getValueSure(idUr,nodeWhereToTransfer) 
			+ state.hgtCost; 
		
		UPDATE_OPT(optCost, bestSplit, cost, state.idX, nodeWhereToTransfer, 't');
 
		 cost = costThisSplit 
			+ mMatrix.getValueSure(idUr ,state.idX) 
			+ mMatrix.getValueSure(idUl,nodeWhereToTransfer) 
			+ state.hgtCost; 
		
		UPDATE_OPT(optCost, bestSplit, cost,  nodeWhereToTransfer, state.idX, 't');
		
//  BUG HERE? - transfer to self
		cost = costThisSplit 
			+ mMatrix.getValueSure(idUl,idXr) 
			+ mMatrix.getValueSure(idUr, nodeWhereToTransfer) 
			+ state.hgtCost; 
		
		UPDATE_OPT(optCost, bestSplit, cost, idXr, nodeWhereToTransfer, 't');
 
		 cost = costThisSplit 
			+ mMatrix.getValueSure(idUr ,idXr) 
			+ mMatrix.getValueSure(idUl,nodeWhereToTransfer) 
			+ state.hgtCost; 
		
		UPDATE_OPT(optCost, bestSplit, cost,  nodeWhereToTransfer, idXr, 't');
             
   /* int idUsub = idUl; 
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
                BOOST_FOREACH( double c, mVT2[state.timeSlice][idUsub] ) {
                    computeCostList( costThisSplit + state.hgtCost + c,
                                     idUother, state.idX, eventTrip, state ); 
                }
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
    */
}



//#define PRINT_COSTS 1;
/**
 * Compute costs for the current split.
 */
void DTLMatrixNetwork::computeOptimaForCladeSplit(
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

// T can be compute as S, see below
    // maxTS has one node (therefore no transfer) normally, but with
    // ILS, there can be multiple idX at maxTS
//    if( mComputeT && ( state.timeSlice != mMaxTS || mUseILS )) 
//        computeTransferCost(idUl, idUr, costThisSplit, state, 
 //               optCost, bestSplit );
//#ifdef PRINT_COSTS
//cout << "  T: " << optCost << endl;
//#endif

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
            
        int idSecFatherSon1 =-2 ;
        if(mSpeciesTree->getNodeById(speciesSplits[i].first)->getInfos().secondaryFather!=NULL)
        	idSecFatherSon1 =mSpeciesTree->getNodeById(speciesSplits[i].first)->getInfos().secondaryFather->getId() ;
        	
        int idSecFatherSon2 =-2;
        if(mSpeciesTree->getNodeById(speciesSplits[i].second)->getInfos().secondaryFather!=NULL)
        	idSecFatherSon2=mSpeciesTree->getNodeById(speciesSplits[i].second)->getInfos().secondaryFather->getId() ;
                
		if(idSecFatherSon1==state.idX || idSecFatherSon2==state.idX){
			#ifdef PRINT_COSTS
			cout << "Speciation not possible, T edge " << speciesSplits[i].first << "," << speciesSplits[i].second <<  " vs "  <<state.idX  <<" \n";
        	#endif  
        	computeTransferCost(idUl, idUr, speciesSplits[i].first, speciesSplits[i].second, costThisSplit, state,  optCost, bestSplit );				
			#ifdef PRINT_COSTS
			cout << "  T : " << optCost << endl;
			#endif    
		} 
		else{
        	computeSpeciationCost( idUl, idUr, speciesSplits[i].first,
            speciesSplits[i].second, costThisSplit, 
            state, optCost, bestSplit );
			#ifdef PRINT_COSTS
			cout << "  S: " << optCost << endl;
			#endif
		}
    }
    

}


/**
 * Compute costs for the current split.
 */
void DTLMatrixNetwork::computeOptimaForSpeciesSplit(
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
        
          int idSecFatherSon1 =-2 ;
          if(mSpeciesTree->getNodeById(speciesSplits[i].first)->getInfos().secondaryFather!=NULL)
        	idSecFatherSon1 =mSpeciesTree->getNodeById(speciesSplits[i].first)->getInfos().secondaryFather->getId() ;
        	
          int idSecFatherSon2 =-2;
          if(mSpeciesTree->getNodeById(speciesSplits[i].second)->getInfos().secondaryFather!=NULL)
        	idSecFatherSon2=mSpeciesTree->getNodeById(speciesSplits[i].second)->getInfos().secondaryFather->getId() ;
          
		                
		  if(idSecFatherSon1==state.idX || idSecFatherSon2==state.idX){
		  
		    if(idSecFatherSon1==state.idX){
            	computeNullCost( speciesSplits[i].second, state, optCost, bestSplit);
            	
            	#ifdef PRINT_COSTS
				cout << "  N primary edge: " << optCost << endl;
				
				#endif  
				
            	computeTransferLossCost( speciesSplits[i].first, state, optCost, bestSplit);

				#ifdef PRINT_COSTS
				cout << "  TL secondary edge: " << optCost << endl;
				
				#endif  
  
			}
			else{
            	computeNullCost( speciesSplits[i].first, state, optCost, bestSplit);
				#ifdef PRINT_COSTS
				cout << "  N primary edge: " << optCost << endl;
				#endif   
            	computeTransferLossCost( speciesSplits[i].second, state, optCost, bestSplit);

				#ifdef PRINT_COSTS
				cout << "  TL secondary edge: " << optCost << endl;
				
				#endif  
				 
			}  
			 #ifdef PRINT_COSTS
			 cout << "Speciation loss not possible, T edge " << speciesSplits[i].first << "," << speciesSplits[i].second <<  " vs "  <<state.idX  <<" \n";
			 #endif
			 continue;
		  }
		
            computeSpeciationPlusLossCost( speciesSplits[i].first, 
                        speciesSplits[i].second, state, optCost, bestSplit );
			#ifdef PRINT_COSTS
			cout << "  SL: " << optCost << endl;
			#endif
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
 * Calculate matrix for species trees with no subdivisions.
 *
 * Species Tree has no subdivision and time slice = idX 
 * (set this here or before?).
 */
void DTLMatrixNetwork::calculateMatrixNoSub( 
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
//    double alphaCost = 0;
    vector<double> alphaCosts;
	initTS(); // clear triplets
    
    for( int ts=0; ts<=mMaxTS; ts++ ) {
        state.timeSlice = ts;
		//cout << ts << endl;
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

		//not used
        /*if( mSubOpt ) 
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
                suboptSetFinalCosts( optVector.back(), allCosts.back(), state);
            else if( mTriplets ) 
                setFinalCosts( ts, state );

            // Save value for setting cost list.
            // - pre TL value need for speciation calculations
            curBestCosts.push_back( mMatrix.getValueSure( idU, state.idX ) );
        }
        */
		
        
    	mMatrix.setValue( idU, state.idX, optVector.back() ); 
        speciesNodeIds.push_back( singleNode[0] );
        
        
    }
    
	/*
	
	
	//the code is not used because the bestReceiver is computed in a different way 
	
	
    // set best receiever and cost for each node
    if( mSubOpt ) 
        computeVTs( true, -1, idU, speciesNodeIds, allCosts ); 
    else if( mTriplets ) 
        computeVTs( true, -1, idU, speciesNodeIds ); 
    else
    	computeBestReceiversNoSub( idU, speciesNodeIds, optVector );
    */
        




    /*if( !mComputeTL && !mSubOpt && !mTriplets )
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
                computeTransferLossCostVT( alphaCosts,
                                        allCosts[ts], state );
                suboptSetFinalCosts( curBestCosts[ts], allCosts[ts], state );
            } else if( mTriplets ) {
                computeTransferLossCostVT( ts, state );
                setFinalCosts( ts, state );
            } else {
                computeTransferLossCost( alphaCost, opt, state );
                mMatrix.setValue( idU, state.idX, opt ); 
            }
        }

        // make sure values are reset, event if TL is not computed
        if( mSubOpt || mTriplets ) 
            mMatrix.setValue( idU, state.idX, curBestCosts[ts] ); 
    }*/
    
    
}

/**
 * Return 1 if idX equals the secondary father of idXl.
 * Return 2 if idX equals the secondary father of idXr.
 * Else return 0.
 */
int DTLMatrixNetwork::secondaryFather(
    int idX,
    int idXl,
    int idXr )
{
    MySpeciesNode *secFatherL = 
        	    mSpeciesTree->getNodeById(idXl)->getInfos().secondaryFather;
    MySpeciesNode *secFatherR = 
        	    mSpeciesTree->getNodeById(idXr)->getInfos().secondaryFather;
        
    if( secFatherL != NULL && idX == secFatherL->getId() )
        return 1;
    if( secFatherR != NULL && idX == secFatherR->getId() )
        return 2;
    return 0;
}


/**
 * Add all possible vertices to the given vertex.
 *
 * Check all possible events and add vertices for events
 * that have the same cost as the vertex.
 */
void DTLMatrixNetwork::createVertices(
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
/* To allow polytomy graphs
  for( int splitIdx=0; splitIdx<splitCount; splitIdx++ ) {
    pair<int,int> cladeSplit = mCladesTrips->getCladeSplit( idU, splitIdx );
*/

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
            if( split.second == -1 )
                continue;

            int sf = secondaryFather( idX, split.first, split.second );
            if( sf == 0 ) {
                // speciation 
                addVertices( graph, pairVertex, qList, idUl, split.first, 
                             idUr, split.second, "S", cost, 0, trip );
                addVertices( graph, pairVertex, qList, idUl, split.second, 
                             idUr, split.first, "S", cost, 0, trip );
            } else {
                int toX = split.first; // sf == 1
                if( sf == 2 )
                    toX = split.second;
                addVertices( graph, pairVertex, qList, idUl, idX,
                    idUr, toX, "T", cost, mHGTCost, trip ); 
                addVertices( graph, pairVertex, qList, idUr, idX,
                    idUl, toX, "T", cost, mHGTCost, trip ); 
                addVertices( graph, pairVertex, qList, idUl, split.second,
                    idUr, toX, "T", cost, mHGTCost, trip ); 
                addVertices( graph, pairVertex, qList, idUr, split.second,
                    idUl, toX, "T", cost, mHGTCost, trip ); 
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
/*       
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
*/
    }
//}
   
    // Speciation plus loss costs
    BOOST_FOREACH( split, speciesSplits ) {
        if( split.second != -1 ) {
            int sf = secondaryFather( idX, split.first, split.second );

            if( sf == 1 ) {
                // NULL cost for second child
                addVertices( graph, pairVertex, qList, idU, split.second, "O", 
                             cost, 0, trip );
                addVertices( graph, pairVertex, qList, idU, split.first, "TL",
                             cost, mLossCost+mHGTCost, trip ); 
            } else if( sf == 2 ) {
                // NULL cost for first child
                addVertices( graph, pairVertex, qList, idU, split.first, "O", 
                             cost, 0, trip );
                addVertices( graph, pairVertex, qList, idU, split.second, "TL",
                             cost, mLossCost+mHGTCost, trip ); 
            }
            if( sf != 2 ) {
                addVertices( graph, pairVertex, qList, idU, split.first, "SL", 
                             cost, mLossCost, trip );
                addVertices( graph, pairVertex, qList, idU, split.second, "SL", 
                             cost, mLossCost, trip );
            }
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
/*
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
*/
}
