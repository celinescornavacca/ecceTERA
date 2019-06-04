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


#include <float.h>
#include <sstream>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/algorithm/string.hpp>
#include <Bpp/Text/TextTools.h>
#include "DTLGraph.h"

// Rounding error tolerance for comparing costs.
const double DTLGraph::SCORE_DIFF = 0.00000000001; 
const int DTLGraph::NAME_PRECISION = 12; // precision for cost in vertex names

/**
 * Score comparison accounting for rounding.
 *
 * @return true if a is greater than b
 */
bool DTLGraph::SCORE_GREATER( double a, double b ) {
    if( a > (b+SCORE_DIFF) ) 
        return true;
    return false;
}
/**
 * Score comparison accounting for rounding.
 *
 * @return true if a is equal to b
 */
bool DTLGraph::SCORE_EQUAL( double a, double b ) {
    if( abs(a-b) <= SCORE_DIFF )
        return true;
    return false;
}



////////////////////////////////////////////////////////////////////////
/////////////////  Graph Traversal Functions //////////////////////////
////////////////////////////////////////////////////////////////////////


/**
 * Recursive part of depthFirstTraversal.
 */
void DTLGraph::depthFirstTraversalAux( 
    MyGraph::Vertex z,  ///<currently visited vertex
    ArgBase &args, ///< structure containing traversal specific variables
    void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&), 
        ///< function called before visiting children
    void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&) ) 
        ///< function called after visiting children
{
    if( mGraph.properties(z).visits > 0 ) 
        return;
    
    mGraph.properties(z).visits++;
    
    if( discoverVertexPtr != NULL ) 
        (this->*discoverVertexPtr)( z, args );

    MyGraph::adjacency_vertex_range_t sons =  mGraph.getAdjacentVertices( z );
    BOOST_FOREACH( MyGraph::Vertex son, sons ) 
        depthFirstTraversalAux( son, args, discoverVertexPtr, finishVertexPtr );
   
    if( finishVertexPtr != NULL ) 
        (this->*finishVertexPtr)( z, args );
}

/**
 * Depth first traversal from each root. 
 *
 * The given functions are called at each vertex (if not null). 
 */
void DTLGraph::depthFirstTraversal( 
    ArgBase &args, ///< structure containing traversal specific variables
    void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&),
        ///< function called before visiting children
    void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&) ) 
        ///< function called after visiting children
{
    resetVisits();
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices )
        depthFirstTraversalAux(root, args, discoverVertexPtr, finishVertexPtr);
}



/**
 * Recursive part of breadthFirstTraversal.
 *
 * A vertex is visited only after all of it's parents have been visited.
 */
void DTLGraph::breadthFirstTraversalAux( 
    MyGraph::Vertex z,  ///<currently visited vertex
    bool isRoot,        ///<vertex is visited as root
    ArgBase &args, ///< structure containing traversal specific variables
    void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&),
        ///< function called before visiting children
    void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&) ) 
        ///< function called after visiting children
{
    if( !isRoot ) // no in degree if isRoot is true
        mGraph.properties(z).visits++;
   
    if( mGraph.properties(z).visits < mGraph.getVertexInDegree(z) ) 
        return;

    if( mGraph.properties(z).visits > mGraph.getVertexInDegree(z) ) {
        cout << "vertex " << mGraph.properties(z).name << " visits "
            << mGraph.properties(z).visits << " > degree "
            << mGraph.getVertexInDegree(z) << endl;
        throw bpp::Exception( "DTLGraph::breadthFirstTraversalAux:"
                " too many visits" );
    }
    
    if( discoverVertexPtr != NULL ) 
        (this->*discoverVertexPtr)( z, args );

    MyGraph::adjacency_vertex_range_t sons =  mGraph.getAdjacentVertices( z );
    BOOST_FOREACH( MyGraph::Vertex son, sons ) 
        breadthFirstTraversalAux( son, false, args, discoverVertexPtr, 
                                  finishVertexPtr );
    
    if( finishVertexPtr != NULL ) 
        (this->*finishVertexPtr)( z, args );
}


/**
 * Breadth first traversal from each root. 
 *
 * The given functions are called at each vertex (if not null). 
 */
void DTLGraph::breadthFirstTraversal( 
    ArgBase &args, ///< structure containing traversal specific variables
    void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&),
        ///< function called before visiting children
    void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&) ) 
        ///< function called after visiting children
{
    resetVisits();
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        if(mGraph.getVertexInDegree(root) == 0 ) 
            breadthFirstTraversalAux( root, true, args, discoverVertexPtr, 
                                  finishVertexPtr );
    }
}



/** 
 * Find number of reconciliations in sub-graph rooted at each vertex
 * (recNumber property). 
 *
 * Traverse the graph and label each vertex with the number of 
 * reconciliations possible for the sub-graph at that vertex by
 * multiplying the reconciliations of children for event vertices
 * and summing the reconciliations of children for mapping vertices.
 * Put the value in recNumber.
 *
 * @return number of reconciliations for the entire graph
*/
double DTLGraph::countSubReconciliations() {

    ArgBase args; // No extra variables needed.
    depthFirstTraversal( args, NULL,
            &DTLGraph::countSubReconciliationsFinishVertex );

    // sum root sons (except counts from null sons),
    // and remove null sons counts from roots with no parents if
    // only canonical
    double rootSum = 0;
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {

        MyGraph::adjacency_vertex_range_t sons 
            = mGraph.getAdjacentVertices(root);

        BOOST_FOREACH( MyGraph::Vertex son, sons) {
            if( !mOnlyCanonical || mGraph.properties(son).name[0] != 'O' ) 
                rootSum += mGraph.properties(son).recNumber;
            else if( mGraph.getVertexInDegree(root) == 0 ) 
                // for canonical, remove 'O' event counts to
                // get the right root recNumber
                mGraph.properties(root).recNumber 
                        -= mGraph.properties(son).recNumber;
        }
    }

    return rootSum;
}


/**
 * Traversal function for countSubReconciliations
 *
 * Leaves have a single reconciliation. Mapping nodes sum sons.
 * Event nodes with one children sum valid grandchilren. Event
 * nodes with two children sum all combinations of valid cousins.
 *
 */
void DTLGraph::countSubReconciliationsFinishVertex( 
        MyGraph::Vertex z, ///< currently visited vertex
        ArgBase &baseArgs ) ///< not used 
{

    int sonCount = mGraph.getVertexOutDegree(z);
    if( sonCount == 0 ) {  // leaf 
        mGraph.properties(z).recNumber = 1;
        return;
    }

    MyGraph::adjacency_vertex_range_t sons = mGraph.getAdjacentVertices(z);

    // mapping nodes sum sons
    if( mGraph.properties(z).isMapping ) {
        BOOST_FOREACH( MyGraph::Vertex son, sons)
            mGraph.properties(z).recNumber += mGraph.properties(son).recNumber;
        return;
    }

    // EVENT NODE - visit children and propogate recNumbers
    MyGraph::adjacency_vertex_range_t eventSons1 
                    = mGraph.getAdjacentVertices( *(sons.first) );
    MyGraph::adjacency_vertex_range_t eventSons2; 
    if( sonCount == 2 ) 
        eventSons2 = mGraph.getAdjacentVertices( *(sons.first+1) );

    double recNumber = 0;
    BOOST_FOREACH( MyGraph::Vertex eventSon1, eventSons1 ) {
        if( sonCount == 1 ) {
            if( !mOnlyCanonical || validEvent( z, eventSon1 ) ) 
                recNumber += mGraph.properties(eventSon1).recNumber;
        } else {
            BOOST_FOREACH( MyGraph::Vertex eventSon2, eventSons2 ) 
                recNumber += validRecNumber( z, eventSon1, eventSon2 );
        }
    }
    mGraph.properties(z).recNumber = recNumber;

    if( recNumber < 0 ) // sanity check 
        throw bpp::Exception( "DTLGraph::CountSubReconciliationsObject:"
                          " negative recNumber" );
}





/**
 * Label each vertex with the number of reconciliations it is in for
 * each root (support or supportVector property).
 *
 * The support for each vertex for the subgraph
 * of each root is stored in supportVector, indexed by the order of
 * the root vertices.
 *
 * The computeSupportDiscoverVertex function is used in a traversal
 * to do the calculation once the root sons are intialized.
 */
void DTLGraph::computeSupport( 
        double epsilon ) ///< Epsilon of costs (used if suboptimal)
{

    MyGraph::vertex_iter iter, v_end;

    if( mWeighted ) {
        // initialize supportVector for all events
        int costCount = mRootVertices.size();
        for(boost::tie(iter,v_end) = mGraph.getVertices(); 
             iter != v_end; iter++) 
	{
            if( !mGraph.properties(*iter).isMapping ) 
                mGraph.properties(*iter).supportVector.resize( costCount, 0 );
        }
    }

    // initialize root event sons
    int idx = 0;
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        BOOST_FOREACH( MyGraph::Vertex son, mGraph.getAdjacentVertices(root) ) {
            if( mOnlyCanonical && mGraph.properties(son).name[0] == 'O' ) 
                continue; // non-canonical, skip
            if( mWeighted ) 
                mGraph.properties(son).supportVector[idx]
                    = mGraph.properties(son).recNumber;
            else
                mGraph.properties(son).support
                    = mGraph.properties(son).recNumber;
        }
        idx++;
    }
    
    // call the traversal function to calculate support.
    ArgBase args;
    breadthFirstTraversal( args, 
            &DTLGraph::computeSupportDiscoverVertex, NULL );

    if( mSuboptimal ) 
        generateMergedSupport( epsilon );
}


/**
 * Traversal function for computeSupport.
 *
 * Calculates support for event nodes. Each event vertice's reconciliations
 * are assigned proportionally to the grand children. For weighted support
 * (mWeighted), the support for each root's subgraph is calculated
 * separately.
 */
void DTLGraph::computeSupportDiscoverVertex( 
        MyGraph::Vertex z, ///< currently visited vertex
        ArgBase &baseArgs ) ///< not used 
{
    if( mGraph.getVertexOutDegree(z) == 0 )
        return; // leaf - processed by parent

    if( mGraph.properties(z).isMapping ) 
        return;
   
    if( mGraph.properties(z).recNumber == 0 ) 
        return; // no point in continuing

    computeR( z ); // compute r for grandchildren

    int rootCount = mRootVertices.size();
    MyGraph::adjacency_vertex_range_t sons = mGraph.getAdjacentVertices(z);
    BOOST_FOREACH( MyGraph::Vertex son, sons ) {
        MyGraph::adjacency_vertex_range_t eventSons 
                = mGraph.getAdjacentVertices(son);
        BOOST_FOREACH( MyGraph::Vertex eventSon, eventSons ) {
            if( mWeighted ) {
                // loop over each separate subgraph (different roots)
                for( int idx=0; idx<rootCount; idx++ ) {
                    mGraph.properties(eventSon).supportVector[idx] 
                        += mGraph.properties(z).supportVector[idx]
                            * mGraph.properties(eventSon).r 
                            / mGraph.properties(z).recNumber;
                }
            } else {
                mGraph.properties(eventSon).support 
                    += mGraph.properties(z).support
                       * mGraph.properties(eventSon).r 
                       / mGraph.properties(z).recNumber;
            }
        }
    }
}


/**
 * Compute r values for grand children.
 *
 * Used by computeSupportDiscoverVertex to calculate valid
 * reconciliations for the grandchilren of z. The r value
 * is the recNumber for valid grandchildren if the z has
 * one child. If z has two children, r is the product
 * of valid cousins (calculated in validRecNumber).
 */
void DTLGraph::computeR( 
        MyGraph::Vertex z ) ///< compute r for event grandsons of this vertex 
{
    int sonCount = mGraph.getVertexOutDegree(z);
    if( sonCount == 0 ) 
        throw bpp::Exception( "DTLGraph::computeR given a leaf");

    if( mGraph.properties(z).isMapping ) 
        throw bpp::Exception( "DTLGraph::computeR: called on pair vertex" );

    MyGraph::adjacency_vertex_range_t sons = mGraph.getAdjacentVertices(z);

    MyGraph::adjacency_vertex_range_t eventSons2;
    if( sonCount == 2 ) {
        eventSons2 = mGraph.getAdjacentVertices( *(sons.first+1) );
        BOOST_FOREACH( MyGraph::Vertex eventSon2, eventSons2) 
            mGraph.properties(eventSon2).r = 0; // initialize
    }

    MyGraph::adjacency_vertex_range_t eventSons1 
                    = mGraph.getAdjacentVertices( *(sons.first) );
    BOOST_FOREACH( MyGraph::Vertex eventSon1, eventSons1) {
        if( sonCount == 1 ) {
            if( !mOnlyCanonical || validEvent( z, eventSon1 ) ) {
                mGraph.properties(eventSon1).r = 
                    mGraph.properties(eventSon1).recNumber;
            } else 
                mGraph.properties(eventSon1).r = 0;
        } else {
            mGraph.properties(eventSon1).r = 0;

            BOOST_FOREACH( MyGraph::Vertex eventSon2, eventSons2) {
                double recNumber = validRecNumber( z, eventSon1, eventSon2 );
                mGraph.properties(eventSon1).r += recNumber; 
                mGraph.properties(eventSon2).r += recNumber; 
            }
        }
    }
}


/**
 * Merge support for similar events for suboptimal costs.
 *
 * For the suboptimal case, event supports (support)
 * must be merged for events with identical parent and
 * children clades and species (u,x), where the costs
 * can be different.
 *
 */
void DTLGraph::generateMergedSupport( 
        double epsilon ) ///< Epsilon of cost range for weighted calculation.
{

    // create map to identify similar events
    map<string,vector<MyGraph::Vertex> > descriptors
        = createEventDescriptorMap(); 

    // loop over all unique events and merge support for the simliar events
    map<string,vector<MyGraph::Vertex> >::const_iterator mapIter;
    for( mapIter = descriptors.begin(); mapIter != descriptors.end(); mapIter++)
    {
        vector<MyGraph::Vertex> events = mapIter->second;
        if( mWeighted ) {
            for( size_t idx=0; idx<mRootVertices.size(); idx++ ) {

                // get sum for this event
                double supportSum = 0;
                BOOST_FOREACH( MyGraph::Vertex event, events ) 
                    supportSum += mGraph.properties(event).supportVector[idx];

                // now set the value for each event
                BOOST_FOREACH( MyGraph::Vertex event, events )
                    mGraph.properties(event).supportVector[idx] = supportSum;
            }
        } else {
            // get sum for this event
            double supportSum = 0;
            BOOST_FOREACH( MyGraph::Vertex event, events ) 
                supportSum += mGraph.properties(event).support;

            // now set the value
            BOOST_FOREACH( MyGraph::Vertex event, events ) 
                mGraph.properties(event).support = supportSum;
        }
    }

    if( mWeighted ) 
        weightSupport( epsilon );   
}

/**
 * Create string identifiers for events and map them to the events.
 * 
 * Generate a string identifer for each event 
 * of the form u,x|u1,x1;u2,x2 and insert event into a map
 */
map<string,vector<DTLGraph::MyGraph::Vertex> > 
DTLGraph::createEventDescriptorMap() 
{ 
    map<string,vector<MyGraph::Vertex> > descriptorMap; 

    MyGraph::vertex_iter iter, v_end;
    for(boost::tie(iter,v_end) = mGraph.getVertices(); iter != v_end; iter++) 
    {
        if( mGraph.properties(*iter).isMapping ) 
            continue; // event nodes only

        if( mGraph.properties(*iter).name[0] == 'C' ) 
            continue;  //leaf

        MyGraph::in_edge_range_t fathers = mGraph.getInEdges(*iter);
        int fatherU = mGraph.getEdgeSource(*fathers.first).id_u; 
        int fatherX = mGraph.getEdgeSource(*fathers.first).id_x; 

        ostringstream eventString;
        eventString << fatherU << "," << fatherX << "|";

        MyGraph::adjacency_vertex_range_t sons 
                = mGraph.getAdjacentVertices(*iter);

        int son1u = mGraph.properties(*(sons.first)).id_u;
        int son1x = mGraph.properties(*(sons.first)).id_x;

        int sonCount = mGraph.getVertexOutDegree(*iter);
        if( sonCount == 1 ) {
            eventString << son1u << "," << son1x;
        } else {
            int son2u = mGraph.properties(*(sons.first+1)).id_u;
            int son2x = mGraph.properties(*(sons.first+1)).id_x;
            if( son2u < son1u || (son1u == son2u && son2x < son1x) ) //swap
                eventString << son2u<<","<<son2x <<";"<< son1u<<","<<son1x;
            else
                eventString << son1u<<","<<son1x <<";"<< son2u<<","<<son2x;
        }
        descriptorMap[eventString.str()].push_back( *iter );
    }

    return descriptorMap;
}

/**
 * Weight the support.
 *
 * Finish the weighted support calculation and save it in
 * the support property for each event vertex.
 */
void DTLGraph::weightSupport( 
        double epsilon ) ///< Epsilon of cost range for weighted calculation.
{
    // get optimal cost
    double optCost = 0;
    bool first = true;
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        if( first || mGraph.properties(root).cost < optCost ) {
            optCost = mGraph.properties(root).cost;
            first = false;
        }
    }

    // calculate p(c)
    vector<double> costAdjs; // p(c) in algorithm
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        double ratio = 0; 
        if( epsilon > 0 )
           ratio =  (mGraph.properties(root).cost - optCost) / epsilon;
        costAdjs.push_back( exp( -ratio ) );
    }

    double costAdjsTotal = 0;
    for( size_t idx=0; idx<costAdjs.size(); idx++ ) 
        costAdjsTotal += costAdjs[idx];

    // find weighted number of solutions
    int idx = 0;
    mWeightedNumberSolutions = 0;
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) 
        mWeightedNumberSolutions +=
            mGraph.properties(root).recNumber*costAdjs[idx++] / costAdjsTotal;


    MyGraph::vertex_iter iter, v_end;
    for(boost::tie(iter,v_end) = mGraph.getVertices(); iter != v_end; iter++) {
        if( mGraph.properties(*iter).isMapping )
            continue; // event nodes only

        mGraph.properties(*iter).support = 0;
        for( size_t idx=0; idx<costAdjs.size(); idx++ ) 
            mGraph.properties(*iter).support 
                += mGraph.properties(*iter).supportVector[idx]
                 *costAdjs[idx] / costAdjsTotal;
    }
}



/**
 * Remove non-canonical vertices. 
 *
 * Find and remove non-canoical vertices. These will be
 * vertices with recNumber=0 if reconciliations were counted
 * with mOnlyCanonical=true.
 *
 * This removes vertices with no reconciliations because there
 * are no canonical reconciliations. However, non-conical reconciliations
 * will still exist in the graph, especially if some of the more
 * complicated non-canonical conditions are use, e.g. condition7.
 */
void DTLGraph::removeNoncanonicalVertices() {

    // With the boost graph library, removing vertices invalidates
    // vertice references. This makes removing them costly. It is easier
    // just to create a new graph.

    mVertexMap.clear();
    MyGraph newGraph;

    RemoveArgs args( &newGraph) ;
    depthFirstTraversal( args, NULL,
            &DTLGraph::removeNoncanonicalFinishVertex );

    // redo root vertices
    vector<MyGraph::Vertex> newRoots;
    for( size_t i=0; i<mRootVertices.size(); i++ ) {
        MyGraph::Vertex root = mRootVertices[i];
        if( mGraph.properties(root).recNumber == 0 ) 
            continue; // skip

        map<string,MyGraph::Vertex>::iterator iter 
                = mVertexMap.find( mGraph.properties(root).name );
        if( iter == mVertexMap.end() ) {
            cout << "COULD NOT FIND " << mGraph.properties(root).name << ", "
                 << mGraph.properties(root).recNumber << endl;
            throw bpp::Exception( "DTLGraph::removeNoncanonicalVertices: "
                        " could not find created root vertex" );
        }
        newRoots.push_back( iter->second );
    }

    mRootVertices = newRoots;
    mGraph = newGraph;
}


/**
 * Traversal function for removeNoncanonicalVertices
 *
 * Copy current vertex and it's out edges to the new graph.
 * The baseArgs variable is cast to RemoveArgs, which contains
 * the new graph.
 *
 */
void DTLGraph::removeNoncanonicalFinishVertex( 
        MyGraph::Vertex z, ///< currently visited vertex
        ArgBase &baseArgs ) ///< cast to RemoveArgs
{

    RemoveArgs *args = (RemoveArgs*) &baseArgs;

    // non-canonical node, don't add
    if( mGraph.properties(z).recNumber == 0 )
        return;

    // create vertex
    MyGraph::Vertex newVertex = args->mNewGraph->addVertex(
                                    mGraph.properties(z) );
    string vertexId = mGraph.properties(z).name;
    mVertexMap.insert( make_pair(vertexId, newVertex) );

    // create edges
    MyGraph::adjacency_vertex_range_t sons = mGraph.getAdjacentVertices(z);
    BOOST_FOREACH( MyGraph::Vertex son, sons ) {
    
        if( mGraph.properties(son).recNumber == 0 ) 
            continue;

        // get the son
        string sonVertexId = mGraph.properties(son).name;
        map<string,MyGraph::Vertex>::iterator iter 
                = mVertexMap.find( sonVertexId );	
        if( iter == mVertexMap.end() ) {
            cout << "Couldn't find " << sonVertexId << endl;
            throw bpp::Exception( "DTLGraph::removeNoncanonicalObject: "
                        " could not find created vertex" );
        }

        // create edge to new son
        EdgeProperties e_properties_temp; //not used right now
        e_properties_temp.propertyNotUsed = "";
        args->mNewGraph->addEdge( newVertex, iter->second, e_properties_temp );
    } 
}


/**
 * Traversal function for calculating the median score (score property).
 *
 * For event nodes, sum the support with the children scores
 * (summing sons differently depending on the problem type). 
 * 
 *
 */
void DTLGraph::bestScoreFinishVertex( 
        MyGraph::Vertex z, ///< currently visited vertex
        ArgBase &baseArgs ) ///< cast to BestScoreArgs
{
    BestScoreArgs *args = (BestScoreArgs*) &baseArgs;

    if( mGraph.properties(z).isMapping ) 
        return;
    
    // support for event vertex
    mGraph.properties(z).score = 0;
    if( mScoredProblem != 5 || mGraph.properties(z).name[0] != 'S'  
             || mGraph.properties(z).name[1] != '_' )
    {
        mGraph.properties(z).score = mGraph.properties(z).support
                                   + args->mScoreMod;
    }

    if( mScoredProblem == 4 && mGraph.properties(z).score < 0 ) 
        mGraph.properties(z).score = 0;

    int sonCount = mGraph.getVertexOutDegree( z );
    if( sonCount == 0 ) // leaf
        return;

    MyGraph::adjacency_vertex_range_t sons = mGraph.getAdjacentVertices(z);
    MyGraph::adjacency_vertex_range_t eventSons1 
                    = mGraph.getAdjacentVertices( *(sons.first) );

    MyGraph::adjacency_vertex_range_t eventSons2;
    if( sonCount == 2 ) 
        eventSons2 = mGraph.getAdjacentVertices( *(sons.first+1) );

    double sMax = 0;
    bool first = true;
    BOOST_FOREACH( MyGraph::Vertex eventSon1, eventSons1) {
    
        double sonScore1 = mGraph.properties(eventSon1).score;
        if( sonCount == 1 ) {
            // take max valid grand son
            if( (first || sonScore1 > sMax)
                && (!mOnlyCanonical || validEvent( z, eventSon1 ))) 
            {
                sMax = sonScore1;
                first = false;
            }
        } else {
            BOOST_FOREACH( MyGraph::Vertex eventSon2, eventSons2 ) {
                double sum = sonScore1 + mGraph.properties(eventSon2).score;
                if( ( first || sum > sMax )
                    && validRecNumber( z, eventSon1, eventSon2 ) != 0 )
                {
                    sMax = sum;
                    first = false;
                }
            }
        }
    }

    mGraph.properties(z).score += sMax;
}

 
/**
 * Traversal function for printGraph.
 */
void DTLGraph::printDiscoverVertex( 
        MyGraph::Vertex z, ///< currently visited vertex
        ArgBase &baseArgs ) ///< cast to PrintArgs
{
    PrintArgs *args = (PrintArgs*) &baseArgs;

    MyGraph::adjacency_vertex_range_t sons =  mGraph.getAdjacentVertices( z );
    BOOST_FOREACH( MyGraph::Vertex son, sons ) {

        // print edge and visit son
        if( !args->mCleanGraph || mGraph.properties(son).recNumber !=0 ) 
        { 

            string v1 = mGraph.properties(z).name;
            if( !mGraph.properties(z).isMapping ) {
                v1 += ","  
                   + bpp::TextTools::toString(mGraph.properties(z).support); 
            } else if( !args->mUseInternalIds ) {
                v1 = createNameWithExternalIds( z, args->mCladeToPOrd );
            }

            string v2 = mGraph.properties(son).name;
            if( !mGraph.properties(son).isMapping ) {
                v2 += ","  
                 + bpp::TextTools::toString(mGraph.properties(son).support); 
            } else if( !args->mUseInternalIds ) {
                v2 = createNameWithExternalIds( son, args->mCladeToPOrd );
            }

            (*args->mOut) << v1 << " -> " << v2 << endl;
        }
    }
}


/**
 * Traversal function for countTLs.
 */
void DTLGraph::countTLsDiscoverVertex( 
        MyGraph::Vertex z, ///< currently visited vertex
        ArgBase &baseArgs ) ///< cast to CountTLsArgs
{
    CountTLsArgs *args = (CountTLsArgs*) &baseArgs;

    if( mGraph.properties(z).name[0] == 'T'
        && mGraph.properties(z).name[1] == 'L' )
    {
        args->tlSupport += mGraph.properties(z).support;
    }

    if( !mGraph.properties(z).isMapping 
        && mGraph.properties(z).name[0] != 'C' ) 
    {
        args->eventSupport += mGraph.properties(z).support;
    }
}

/**
 * Average number of TLs per reconciliation.
 *
 * @return The support for TLs.
 */
double DTLGraph::countTLs( 
        double &frequency )  ///< frequency of TL (among non-C events)
{
    // Use graph traversal.
    CountTLsArgs args;
    args.tlSupport = 0;
    args.eventSupport = 0;
    depthFirstTraversal( args, &DTLGraph::countTLsDiscoverVertex, NULL );
    double tlCount = args.tlSupport/getNumberSolutions();
    frequency = args.tlSupport/ args.eventSupport;
    return tlCount;
}





/**
 * Get events for this vertex.
 */
void DTLGraph::getEvents( 
    MyGraph::Vertex z, ///< vertex 
    string &d,  ///< duplication count
    string &t,  ///< transfer count
    string &l ) ///< loss count
{
    vector<string> tokens;
    boost::split( tokens, mGraph.properties(z).name, 
                  boost::is_any_of("(,)") );
    if( tokens.size() != 8 ) {
        cout << "token count = " << tokens.size() << endl;
        BOOST_FOREACH( string t, tokens )
            cout << " <" << t << ">";
        cout << endl;
        cout << "name=" << mGraph.properties(z).name << endl;
        throw bpp::Exception( "DTLGraph::getEvents:"
                " triplets not found in vertex name" );
    }

    d = tokens[4];
    t = tokens[5];
    l = tokens[6];
}

/**
 * Get vertex properties.
 */
void DTLGraph::getVertexIdentfiers( 
        MyGraph::Vertex vertex,  ///< the vertex
        int &id_u, ///< clade id 
        int &id_x, ///< species id
        double &cost, ///< cost
        int &d,  ///< duplication count
        int &t,  ///< transfer count
        int &l ) ///< loss count
{
    id_u = mGraph.properties(vertex).id_u;
    id_x = mGraph.properties(vertex).id_x;
    cost = mGraph.properties(vertex).cost;

    string dupl, trans, loss;
    getEvents( vertex, dupl, trans, loss );

    d = bpp::TextTools::toInt( dupl );
    t = bpp::TextTools::toInt( trans );
    l = bpp::TextTools::toInt( loss );
}

/**
 * Create a name using external ids (clade post order and 
 * species real post order).
 *
 * The name is of the format "(idU,idX,cost)".
 *
 * @return name
 */
string DTLGraph::createNameWithExternalIds(
    MyGraph::Vertex z, ///< vertex 
    vector<int> &cladeToPOrd ) ///< mapping for clade ids
{
    int pOrd = cladeToPOrd[mGraph.properties(z).id_u];

    int idX = mGraph.properties(z).id_x;
    int rpo = mSpeciesTree->getRPO(idX);

    double cost = mGraph.properties(z).cost;

    if( mPareto ) {
        // get event triplets from name
        string duplCount, transCount, lossCount;
        getEvents( z, duplCount, transCount, lossCount );
        return  "(" + bpp::TextTools::toString(pOrd) + "," 
            + bpp::TextTools::toString(rpo) + ","
            + bpp::TextTools::toString(cost,NAME_PRECISION) 
            + "," + duplCount + "," + transCount + "," + lossCount + ")";
    } else
        return  "(" + bpp::TextTools::toString(pOrd) + "," 
            + bpp::TextTools::toString(rpo) + ","
            + bpp::TextTools::toString(cost,NAME_PRECISION) + ")";
}



/**
 * Print the edges of the graph (pairs of vertices) with the number
 * of reconciliations each vertex is in (support).
 */
void DTLGraph::printGraph( 
    string path, ///< File name
    bool useInternalIds,
    bool cleanGraph ) ///< If true, only print canonical vertices.
{
    vector<int> cladeToPOrd;
    if( !useInternalIds )
        cladeToPOrd = mCladesTrips->getPostOrderMapping();

    filebuf file;
    file.open(path.c_str(), ios::out);
    ostream out(&file);

    // Use graph traversal.
    PrintArgs args( cleanGraph, useInternalIds, cladeToPOrd, &out );
    depthFirstTraversal( args, &DTLGraph::printDiscoverVertex, NULL );
    
    file.close();
}




/**
 * Remove all vertices with non-optimal costs.
 */
void DTLGraph::pruneNonoptimal()
{
    // Use removeNoncanonicalVertices by setting recNumber= 0
    // for non-optimal nodes

    // Get all best costs for u,x
    map<string,double> bestCosts; 
    MyGraph::vertex_iter iter, v_end;
    for(boost::tie(iter,v_end) = mGraph.getVertices(); iter != v_end; iter++) {
        if( mGraph.properties(*iter).isMapping ) {
            int idU = mGraph.properties(*iter).id_u;
            int idX = mGraph.properties(*iter).id_x;
            double cost = mGraph.properties(*iter).cost;
            string pairStr = bpp::TextTools::toString(idU) + "," 
                             + bpp::TextTools::toString(idX);
            map<string,double>::iterator iter = bestCosts.find( pairStr );	
            if ( iter == bestCosts.end() ) {
                // insert
                bestCosts.insert( make_pair(pairStr, cost) );
            } else if( cost < iter->second ) {
                iter->second = cost;
            }
        }
    }

    // mark for removal vertices without best costs
    for(boost::tie(iter,v_end) = mGraph.getVertices(); iter != v_end; iter++) {
        if( mGraph.properties(*iter).isMapping ) {
            int idU = mGraph.properties(*iter).id_u;
            int idX = mGraph.properties(*iter).id_x;
            double cost = mGraph.properties(*iter).cost;
            string pairStr = bpp::TextTools::toString(idU) + "," 
                             + bpp::TextTools::toString(idX);
            map<string,double>::iterator costIter = bestCosts.find( pairStr );	
            if ( costIter == bestCosts.end() ) {
                // shouldn't happen
                throw bpp::Exception( "DTLGraph::pruneNonoptimal:"
                        " vertex pair not found" );
            } else if( !SCORE_GREATER( cost, costIter->second )) {
            //} else if( cost <= costIter->second + SCORE_DIFF ) {
                mGraph.properties(*iter).recNumber = 1;
            }
        } else {
            // include all event vertices, if they are disconnected because
            // the mapping parent is set to 0, then they should not be
            // put into the new graph with removeNoncanonical
            mGraph.properties(*iter).recNumber = 1;
        }
    }

    // remove non-optimal roots
    double bestCost = std::numeric_limits<double>::max();
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        double cost = mGraph.properties(root).cost;
        if( cost < bestCost )
            bestCost = cost;
    }
    vector<MyGraph::Vertex> newRoots;
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        double cost = mGraph.properties(root).cost;
        if( !SCORE_GREATER( cost, bestCost ) ) 
        //if( cost <= bestCost + SCORE_DIFF ) 
            newRoots.push_back( root );
    }
    mRootVertices = newRoots;
        
    // remove vertices
    removeNoncanonicalVertices();

    // reset recNumber
    for(boost::tie(iter,v_end) = mGraph.getVertices(); iter != v_end; iter++)
        mGraph.properties(*iter).recNumber = 0;
}

/**
 * Get event supports as a fracation of the number of solutions.
 */
void DTLGraph::getEventSupports( 
        map<string,double> &eventSupports )  ///< event supports to return
{
    MyGraph::vertex_iter iter, v_end;
    for(boost::tie(iter,v_end) = mGraph.getVertices(); iter != v_end; iter++) {
        if( !mGraph.properties(*iter).isMapping ) {
            string name = mGraph.properties(*iter).name;
            double support 
                = mGraph.properties(*iter).support/getNumberSolutions();
            eventSupports.insert( make_pair(name, support) );
        }
    }
}







////////////////////////////////////////////////////////////////////////
/////////////////  Graph Building Functions ////////////////////////////
////////////////////////////////////////////////////////////////////////


/**
 * Creates vertex names.
 *
 * Used in the map from names to vertices to retreive
 * vertices by name
 */
string DTLGraph::mappingFromIds(
    int id_u, ///< clade identifier
    int id_x, ///< species identifier
    double cost ) ///< cost for this id_u,id_x pair
{
    return  "(" + bpp::TextTools::toString(id_u) + "," 
        + bpp::TextTools::toString(id_x) + ","
        + bpp::TextTools::toString(cost,NAME_PRECISION) + ")";
}

/**
 * Creates vertex names with event counts.
 *
 * Used in the map from names to vertices to retreive
 * vertices by name
 */
string DTLGraph::mappingFromIds(
    int id_u, ///< clade identifier
    int id_x, ///< species identifier
    double cost, ///< cost for this id_u,id_x pair
    int d, ///< duplication counts
    int t, ///< transfer counts
    int l) ///< loss counts
{
    return  "(" + bpp::TextTools::toString(id_u) + "," 
        + bpp::TextTools::toString(id_x) + ","
        + bpp::TextTools::toString(cost,NAME_PRECISION)  + ","
        + bpp::TextTools::toString(d)  + ","
        + bpp::TextTools::toString(t)  + ","
        + bpp::TextTools::toString(l)  + ")";
}

/**
 * Add an edge between eventVertex and the vertex identified by
 * by id_u, id_x, and cost. 
 *
 * A vertex is created if it doesn't 
 * exist, using id_u, id_x, and cost as identifiers.
 *
*  @return true if it is novel (added)
*/
bool DTLGraph::addPairVertex( 
    MyGraph::Vertex eventVertex, ///< the parent vertex of the new vertex
    int id_u,   ///< clade identifier
    int id_x,   ///< species node identifier 
    double cost, ///< node cost 
    MyGraph::Vertex &pairVertex ) ///< returns the vertex identifed by id_u,
    ///< id_x, and cost - either newly created or found in the current set 
{

    string pairVertexName = mappingFromIds( id_u, id_x, cost );
    
    EdgeProperties e_properties_temp; //not used right now
    e_properties_temp.propertyNotUsed ="";

    map<string, MyGraph::Vertex>::iterator iter =
                            mVertexMap.find( pairVertexName );	
    bool added = false;
    if ( iter == mVertexMap.end() ) {
        pairVertex = createPairVertex( id_u, id_x, cost );
        mGraph.addEdge( eventVertex, pairVertex, e_properties_temp );
        added = true;
    } else {
        mGraph.addEdge( eventVertex, iter->second, e_properties_temp );
    }

    return added;
}

/**
 * Add an edge between eventVertex and the vertex identified by
 * by id_u, id_x, cost, and d, t, and l counts. 
 *
 * A vertex is created if it doesn't 
 * exist, using id_u, id_x, cost, d, t, an l as identifiers.
 *
*  @return true if it is novel (added)
*/
bool DTLGraph::addPairVertex( 
    MyGraph::Vertex eventVertex, ///< the parent vertex of the new vertex
    int id_u,   ///< clade identifier
    int id_x,   ///< species node identifier 
    double cost, ///< node cost 
    int d,      ///< duplication counts
    int t,      ///< transfer counts
    int l,      ///< loss counts
    MyGraph::Vertex &pairVertex ) ///< returns the vertex 
        ///< - either newly created or found in the current set 
{

    string pairVertexName = mappingFromIds( id_u, id_x, cost, d, t, l );
    
    EdgeProperties e_properties_temp; //not used right now
    e_properties_temp.propertyNotUsed ="";

    map<string, MyGraph::Vertex>::iterator iter =
                            mVertexMap.find( pairVertexName );	
    bool added = false;
    if ( iter == mVertexMap.end() ) {
        pairVertex = createPairVertex( id_u, id_x, cost, d, t, l );
        mGraph.addEdge( eventVertex, pairVertex, e_properties_temp );
        added = true;
    } else {
        mGraph.addEdge( eventVertex, iter->second, e_properties_temp );
    }

    return added;
}

   		

/** 
 * Create an event vertex and add an edge from z to it.
 *
 * @param event 
 */
DTLGraph::MyGraph::Vertex DTLGraph::addEventVertex( 
        MyGraph::Vertex z, ///< connecting vertex
        string event )  ///< type of event for the new vertex
{

    EdgeProperties e_properties_temp;//not used right now
    e_properties_temp.propertyNotUsed ="";

    string tempEvent = event + "_" + bpp::TextTools::toString(mEventNumber++);
    VertexProperties v_properties;
    v_properties.name = tempEvent;
    v_properties.isMapping = false; 
    v_properties.visits = 0;
    v_properties.recNumber = 0;
    v_properties.support = 0;

    MyGraph::Vertex eventVertex = mGraph.addVertex( v_properties );

    mGraph.addEdge( z, eventVertex, e_properties_temp );

    mVertexMap.insert( make_pair(tempEvent, z) );

    return eventVertex;
}

/** 
 * Create a root vertex and return it.
 *
 * @param new vertex
 */
DTLGraph::MyGraph::Vertex DTLGraph::addRoot( 
        int id_u, ///< clade identifier
        int id_x, ///< species node identifier
        double cost, ///< node cost
        int d, ///< duplication counts
        int t, ///< transfer counts
        int l) ///< loss counts
{
    DTLGraph::MyGraph::Vertex v = createPairVertex( id_u, id_x, cost, d, t, l );

    mRootVertices.push_back( v );

    return v;
}



/**
 * Create vertex using id_u, id_x, and cost as the vertex
 * identifiers and add it to vertex map.
 *
 * @return new vertex
 */
DTLGraph::MyGraph::Vertex DTLGraph::createPairVertex( 
        int id_u, ///< clade identifier
        int id_x, ///< species node identifier
        double cost, ///< node cost
        int d, ///< duplication counts
        int t, ///< transfer counts
        int l) ///< loss counts
{
    if( id_u > mMaxIdU )
        mMaxIdU = id_u;

    string pairVertexName;
    if( d == -1 )
        pairVertexName = mappingFromIds( id_u, id_x, cost );
    else
        pairVertexName = mappingFromIds( id_u, id_x, cost, d, t, l );
       
    VertexProperties v_properties;
    v_properties.isMapping = true; 
    v_properties.name = pairVertexName;
    v_properties.id_u = id_u;
    v_properties.id_x = id_x;
    v_properties.cost = cost;
    v_properties.visits = 0;
    v_properties.recNumber = 0;
    v_properties.support = 0;

    MyGraph::Vertex vertex = mGraph.addVertex( v_properties );

    mVertexMap.insert( make_pair(pairVertexName, vertex) );

    return vertex;
}


		
		











////////////////////////////////////////////////////////////////////////
/////////////////  Graph Property Functions ////////////////////////////
////////////////////////////////////////////////////////////////////////




     
/**
 * Label the graph with the number of reconciliations each vertex
 * is part of. 
 *
 * If mOnlyCanonical is true, remove any vertices not
 * in any canonical reconciliation (recNumber=0).
 *
 * @return True if the number of solutions is finite.
 */
bool DTLGraph::countReconciliationNumberAndCheck( 
    double epsilon, ///< Epsilon, for suboptimal graphs.
    bool onlyCanonical, ///< Remove non-canonical vertices
    bool verbose, ///< print debugging information
    bool weighted ) ///< weight the support for each vertex (suboptimal)
{
    clock_t start = clock();
    mOnlyCanonical = onlyCanonical;

    // count the number of reconciliations from each subtree
    mNumberSolutions = countSubReconciliations();
    if( !boost::math::isfinite( mNumberSolutions ) ) 
        return false;


    if( mSuboptimal ) //&& weighted ) {
        mWeighted = weighted;

    // label graph with number reconciliations in which each vertex
    // is in (support)
    computeSupport( epsilon );
    if( mWeighted ) {
        if( !boost::math::isfinite( mWeightedNumberSolutions ) ) 
            return false;
        if( verbose )
            cout << "weighted solutions " << mNumberSolutions 
                 << " -> " << mWeightedNumberSolutions << endl;
    }


    // Remove non-canonical vertices
    if( mOnlyCanonical ) 
        removeNoncanonicalVertices();

    
    double time = (double) (clock()-start) / CLOCKS_PER_SEC * 1000.0;
    if( verbose ) 
        cout << "count reconciliation time = " << time << endl;
    
    return true;

}




/**
 * Calculates reconciliation score (4 types).
 *
 * @return the score
 */
double DTLGraph::getBestScore( 
    int problem ) ///< The type of score.
{

    if( problem < 0 || problem > 5 )
        throw bpp::Exception("DTLGraph::getBestScoringReconciliation:"
                "Unknown problem");

    if( problem == 2 || problem == 4 )
        throw bpp::Exception("DTLGraph::getBestScoringReconciliation:"
                "Problem not implemented");


    double scoreMod = 0;
    if( problem > 2 && problem < 5 )
        scoreMod = -getNumberSolutions()/2;

    if( mScoredProblem != problem ) { // else already done

        mScoredProblem = problem;

        // score the graph
        BestScoreArgs args( scoreMod );
        depthFirstTraversal( args, NULL, &DTLGraph::bestScoreFinishVertex );

        // now all sons, ignoring grandson null events for canonical
        double maxScore = 0;
        bool first = true;
        BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
            MyGraph::adjacency_vertex_range_t eventSons
                    = mGraph.getAdjacentVertices(root);
            BOOST_FOREACH( MyGraph::Vertex eventSon, eventSons) {
                double score = mGraph.properties(eventSon).score;
                if( problem == 5 ) {
                    int u, x, d, t, l;
                    double cost;
                    getVertexIdentfiers( root, u, x, cost, d, t, l );
                    int sum = d + t + l;
                    if( sum > 0 )
                        score /= sum;
                }
                if( (first || score > maxScore)
                    && ( !mOnlyCanonical 
                        || mGraph.properties(eventSon).name[0] != 'O')) 
                {
                    maxScore = score;
                    first = false;
                }
            }
        }
        mProblemScore = maxScore;
    }
    
    return mProblemScore;
}


/**
 * Get the reconciliation score.
 *
 * @return reconciliation score
 */
double DTLGraph::scoreReconciliation(
    vector< vector<MyGraph::Vertex> > &reconciliation )
        ///< Reconciliation being built
{
    double score = 0;
    BOOST_FOREACH( vector<MyGraph::Vertex> vv, reconciliation )
        BOOST_FOREACH( MyGraph::Vertex z, vv ) 
            if( !mGraph.properties(z).isMapping )
                score += mGraph.properties(z).support;
    return score;
}


/**
 * Recursive function to implment getScoreReconciliation.
 *
 * Fill and return a vector with nodes from a single reconciliation,
 * corresponding to the best score calculated.
 * The tree is created by choosing just one child for each
 * mapping node. The choice is either the child corresponding
 * to the currently score problem or a random child.
 *
 * This function handles the mapping nodes and sets up the call
 * to the backtrack through events.
 *
 * Note: this function recursive, but does not use the graph traversal
 *   functionality since it is selectively visiting vertices,
 *   and not all vertices.
 *
 */
void DTLGraph::backtrack( 
    MyGraph::Vertex vertex, ///< current event vertex in the recursion
    vector< vector<MyGraph::Vertex> > &reconciliation, 
        ///< Reconciliation being built
    bool random, ///< Build a random reconciltaion if true.
    double scoreMod ) ///< Value to add to score for certain problems.
{
    // make sure vertices are only visited once
    if( mGraph.properties(vertex).visits > 0 )
        return;
    mGraph.properties(vertex).visits++;
    if( mGraph.properties(vertex).isMapping ) 
        throw bpp::Exception( "DTLGraph::backtrack visited mapping node" );
    
    // leaf
    if( mGraph.getVertexOutDegree( vertex ) == 0 ) 
        return;

    MyGraph::adjacency_vertex_range_t sons = mGraph.getAdjacentVertices(vertex);
    MyGraph::Vertex son1 = *(sons.first);
    MyGraph::adjacency_vertex_range_t eventSons1 
                    = mGraph.getAdjacentVertices(son1);

    // call backtrack through events
    if( mGraph.getVertexOutDegree( vertex ) == 1 ) {
        backtrackEvents( vertex, son1, eventSons1, reconciliation, random, 
                         scoreMod );
    } else { 
        // son count == 2
        MyGraph::Vertex son2 = *(sons.first+1);
        MyGraph::adjacency_vertex_range_t eventSons2 
                    = mGraph.getAdjacentVertices(son2);
        backtrackEvents( vertex, son1, son2, eventSons1, eventSons2, 
                         reconciliation, random, scoreMod );
    }
}

/**
 * Choose a support randomly, weighting the supports.
 *
 * @return weighted random number
 */
size_t DTLGraph::chooseRandomSupport( 
        vector<double> &supports ) ///< list of supports
{
    double supportSum = 0;
    BOOST_FOREACH( double d, supports )
        supportSum += d;
    vector<double> fracSupports;
    BOOST_FOREACH( double d, supports )
        fracSupports.push_back( d/supportSum );

    double randSupport = rand() / RAND_MAX; 
                // slightly biased towards first bin
    double fracSum = 0;
    for( size_t i=0; i< fracSupports.size(); i++ ) {
        fracSum += fracSupports[i];
        if( !SCORE_GREATER(randSupport, fracSum) ) 
        //if( randSupport <= fracSum+SCORE_DIFF ) 
            return i;
    }


    // failed because supports are infinite, choose random bin
    return rand() % supports.size();
}


/**
 * Backtrack through events for vertices with one son.
 *
 * Follow sons whose score matches the best score.
 */
void DTLGraph::backtrackEvents( 
    MyGraph::Vertex parentVertex, ///< parent event vertex
    MyGraph::Vertex mappingVertex, ///< parent mapping vertex
    MyGraph::adjacency_vertex_range_t eventSons, 
        ///< current event verticies in the recursion
    vector< vector<MyGraph::Vertex> > &reconciliation, 
        ///< Reconciliation being built
    bool random, ///< build a random reconciltaion if true
    double scoreMod ) ///<Value to add to score for certain problems.
{
    double fz = mGraph.properties(parentVertex).support + scoreMod;

    // find chosen event
    bool foundIt = false;
    vector<MyGraph::Vertex> validSons;
    vector<double> validSupports;
    BOOST_FOREACH( MyGraph::Vertex eventSon1, eventSons) {
        if( !mOnlyCanonical || validEvent( parentVertex, eventSon1 ) ) {
            if( random ) {
                validSons.push_back( eventSon1 );
                validSupports.push_back( mGraph.properties(eventSon1).support );
            } else {
                double diff = abs( mGraph.properties(parentVertex).score
                              -(fz+mGraph.properties(eventSon1).score) );
                if( ( SCORE_EQUAL( mGraph.properties(parentVertex).score,
                              fz+mGraph.properties(eventSon1).score) ||
                //if( ( diff < SCORE_DIFF ||
                    diff/mGraph.properties(parentVertex).score < SCORE_DIFF )
                 )   //&& (!mPareto || sameEvents( parentVertex, eventSon1 ))
                {
                    int idU = mGraph.properties(mappingVertex).id_u;
                    reconciliation[idU].push_back( mappingVertex );
                    reconciliation[idU].push_back( eventSon1 );
                    backtrack( eventSon1, reconciliation, random, scoreMod );
                    foundIt = true;
                    break;
                }
            } 
        }
    }

    if( random ) {
        // pick a random valid son
        //int randIdx = rand() % validSons.size();
        int randIdx = chooseRandomSupport( validSupports );

        int idU = mGraph.properties(mappingVertex).id_u;
        reconciliation[idU].push_back( mappingVertex );
        reconciliation[idU].push_back( validSons[randIdx] );
        backtrack( validSons[randIdx], reconciliation, random, scoreMod );
    } else if( !foundIt ) {
        cout << mGraph.properties(parentVertex).name 
             << " score = " << mGraph.properties(parentVertex).score << endl;
        throw bpp::Exception( "DTLGraph::backtrackEvents1 - son not found" );
    }
}



/**
 * Backtrack through events for vertices with one son.
 *
 * Follow sons whose score matches the best score.
 */
void DTLGraph::backtrackEvents( 
    MyGraph::Vertex parentVertex, ///< parent mapping vertex
    MyGraph::Vertex mappingVertex1, ///< parent mapping vertex 1
    MyGraph::Vertex mappingVertex2, ///< parent mapping vertex 2
    MyGraph::adjacency_vertex_range_t eventSons1, 
        ///< current event verticies in the recursion from son1
    MyGraph::adjacency_vertex_range_t eventSons2, 
        ///< current event verticies in the recursion from son2
    vector< vector<MyGraph::Vertex> > &reconciliation, 
        ///< Reconciliation being built
    bool random, ///< built a random reconciltaion
    double scoreMod ) ///<Value to add to score for certain problems.
{
    double fz = 0;
    if( mScoredProblem != 5 || mGraph.properties(parentVertex).name[0] != 'S'  
             || mGraph.properties(parentVertex).name[1] != '_' )
    {
        fz = mGraph.properties(parentVertex).support + scoreMod;
    }

    // double loop over all pairs of cousins
    bool foundIt = false;
    vector<pair<MyGraph::Vertex,MyGraph::Vertex> > validSonPairs;
    vector<double> validSupports;
    BOOST_FOREACH( MyGraph::Vertex eventSon1, eventSons1 ) {
        BOOST_FOREACH( MyGraph::Vertex eventSon2, eventSons2 ) {
            if( validRecNumber( parentVertex, eventSon1, eventSon2 ) != 0 )
            {
                if( random ) {
                    validSonPairs.push_back( 
                            make_pair( eventSon1, eventSon2 ));
// Is this right? What is the support of a pair of events?
                    validSupports.push_back(
                            mGraph.properties(eventSon1).support 
                            *mGraph.properties(eventSon2).support );
                } else {
                    // this calculation must be done in the same
                    // order as in bestScoreFinishVertex due
                    // to rounding errors
                    double sonSum = mGraph.properties(eventSon1).score
                                    + mGraph.properties(eventSon2).score;
                    double score = fz + sonSum;
                    double diff = abs(mGraph.properties(parentVertex).score
                                      - score );
                    if( ( SCORE_EQUAL( mGraph.properties(parentVertex).score,
                                       score ) ||
                    //if( ( diff < SCORE_DIFF ||
                       diff/mGraph.properties(parentVertex).score < SCORE_DIFF )
                      ) //&& (!mPareto 
                        // || sameEvents( parentVertex, eventSon1, eventSon2 )))
                    {
                        int idU = mGraph.properties(mappingVertex1).id_u;
                        reconciliation[idU].push_back( mappingVertex1 );
                        reconciliation[idU].push_back( eventSon1 );
                        backtrack(eventSon1, reconciliation, random, scoreMod );
                        idU = mGraph.properties(mappingVertex2).id_u;
                        reconciliation[idU].push_back( mappingVertex2 );
                        reconciliation[idU].push_back( eventSon2 );
                        backtrack(eventSon2, reconciliation, random, scoreMod );
                        foundIt = true;
                        break;
                    }
                } 
            }
        }
        if( foundIt )
            break;
    }

    if( random ) {
        // pick a random valid pair of sons
        //int randIdx = rand() % validSonPairs.size();
        int randIdx = chooseRandomSupport( validSupports );

        int idU = mGraph.properties(mappingVertex1).id_u;
        reconciliation[idU].push_back( mappingVertex1 );
        reconciliation[idU].push_back( validSonPairs[randIdx].first );
        backtrack( validSonPairs[randIdx].first, reconciliation, random,
                    scoreMod );
        idU = mGraph.properties(mappingVertex2).id_u;
        reconciliation[idU].push_back( mappingVertex2 );
        reconciliation[idU].push_back( validSonPairs[randIdx].second );
        backtrack( validSonPairs[randIdx].second, reconciliation, random,
                    scoreMod );
    } else if( !foundIt ) {
        cout << mGraph.properties(parentVertex).name 
             << " score = " << mGraph.properties(parentVertex).score << endl;
        throw bpp::Exception( "DTLGraph::backtrackEvents2 - son not found" );
    }
}




/**
 * Get the reconciliation corresponding to the best score.
 *
 * The graph is traversed, following only vertices corresponding
 * to the best score.
 *
 * @return False if the score is infinite.
 */
bool DTLGraph::getScoredReconciliation( 
        int problem, ///< Which score to calculate.
        size_t cladeCount, ///< Number of clades.
        vector< vector<DTLGraph::MyGraph::Vertex> > &reconciliation,
            ///< The reconciliation as a list of clades (id_u) with associated 
            ///< specie nodes (id_x).
        bool random )   ///< Create a randome reconciliation if true.
{
    // score the graph if necessary
    if( !random && mScoredProblem != problem ) 
        getBestScore( problem );

    if( !random && !boost::math::isfinite( mProblemScore ) ) 
        return false;

    double scoreMod = 0;
    if( problem > 2 && problem < 5 )
        scoreMod = -getNumberSolutions()/2;

    resetVisits(); 
    reconciliation.resize( cladeCount, vector<MyGraph::Vertex>() );

    double optimalCost = 0;
    if( random ) {
        // find the optimal cost
        bool first = true;
        BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
            if( first || mGraph.properties(root).cost < optimalCost ) {
                optimalCost = mGraph.properties(root).cost;
                first = false;
            }
        }
    }

    // Call backtracking on each root, avoiding non-canoical if mOnlyCanonical.
    bool foundIt = false;
    vector<MyGraph::Vertex> validEventSons;
    vector<double> validSupports;
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        MyGraph::adjacency_vertex_range_t eventSons
                = mGraph.getAdjacentVertices(root);
        BOOST_FOREACH( MyGraph::Vertex eventSon, eventSons) {
            if( !mOnlyCanonical || mGraph.properties(eventSon).name[0] != 'O')
            {
                if( random ) {
                    if( mGraph.properties(root).cost == optimalCost ) {
                        validEventSons.push_back( eventSon );
                        validSupports.push_back( 
                                mGraph.properties(eventSon).support );
                    }
                } else {
                    double score = mGraph.properties(eventSon).score;
                    if( problem == 5 ) {
                        int u, x, d, t, l;
                        double cost;
                        getVertexIdentfiers( root, u, x, cost, d, t, l );
                        int sum = d + t + l;
                        if( sum > 0 )
                            score /= sum;
                    }
                    double diff = abs(mProblemScore-score);
                    if( ( SCORE_EQUAL( mProblemScore, score ) 
                    //if( (diff < SCORE_DIFF 
                        || diff/mProblemScore < SCORE_DIFF )
                      ) // && (!mPareto || sameEvents( root, eventSon )))
                    {
                        int idU = mGraph.properties(root).id_u;
                        reconciliation[idU].push_back( root );
                        reconciliation[idU].push_back( eventSon );
                        backtrack( eventSon, reconciliation, random, scoreMod );
                        foundIt = true;
                        break;
                    }
                }
            }
        }
        if( foundIt ) 
            break;
    }

    if( random ) {
        // pick a random valid son
        //int randIdx = rand() % validEventSons.size();
        int randIdx = chooseRandomSupport( validSupports );

        MyGraph::Vertex randVertex = validEventSons[randIdx];
        // get it's father 
        MyGraph::in_edge_range_t fathers = mGraph.getInEdges(randVertex);
        MyGraph::Vertex root = mGraph.returnSource(*fathers.first);

        int idU = mGraph.properties(root).id_u;
        reconciliation[idU].push_back( root );
        reconciliation[idU].push_back( randVertex );
        backtrack( randVertex, reconciliation, random, scoreMod );
    } else if( !foundIt ) {
        cout << "score=" << mProblemScore << endl;
        throw bpp::Exception( "DTLGraph::getScoredReconciliation - son not found" );
    }

    return true;
}













////////////////////////////////////////////////////////////////////////
/////////////////  Utility Functions ///////////////////////////////////
////////////////////////////////////////////////////////////////////////

/**
 * For the given node, find the species node associated with the
 * vertex (mapping vertex) and return the post order id of that node's parent
 * in the species tree.
 *
 * @return parent id or -1 if no father
 */
int DTLGraph::getSpeciesParentId( 
        MyGraph::Vertex z ) ///< A mapping vertex.
{
    MySpeciesNode *x = mSpeciesTree->getNodeById( mGraph.properties(z).id_x ); 
    if( !x->hasFather() )
        return -1;
    return x->getFather()->getId();
}


/**
 * For a given event vertex, get the id of its parent mapping vertex.

 * @return id of parent
 */
int DTLGraph::getFatherXId( 
        MyGraph::Vertex event ) ///< an event vertex
{
    MyGraph::in_edge_range_t fathers = mGraph.getInEdges(event);
    if( mGraph.getVertexInDegree( event ) == 0 ) 
        throw bpp::Exception( "DTLGraph::getFatherXId: called with no in edges" );
    return mGraph.getEdgeSource(*fathers.first).id_x; 
}




/**
 * Validity check for event vertices with one son (e.g. SL events).
 *
 * If onlyCanonical is true, do not count multiple copies of
 * events added because of artificial nodes. Specifically, do not count the
 * following:
 *  1. Duplication events with two null grandchildren.
 *  2. Transfer events with two null grandchildren.
 *  3. Transfer loss events with a null grandchild, and the 
 *     mapping vertex preceding the TL event is associated
 *     with an artifical node (single child in species tree).
 *  4. Null child events of the roots. 
 */
bool DTLGraph::validEvent( 
        MyGraph::Vertex event, ///< an event vertex
        MyGraph::Vertex eventSon ) ///< the mapping vertex son of event
{

    bool valid = true;

    string eventName = mGraph.properties(event).name;

    // null events are always valid (except root sons, which
    // are handled outside of this function)
    if( eventName[0] == 'O' ) 
        return true;

    string eventSonName = mGraph.properties(eventSon).name;

    if( eventName[0] == 'T' && eventName[1] == 'L'
        && eventSonName[0] == 'O' ) 
    {
        int id_x = getFatherXId( event );
        if( !mSpeciesTree->isReal( id_x ) ) {
            valid = false;
        }
    }
/*
    if( !mSuboptimal )
        return valid;


IF THIS IS TURNED ON, THE PRINT AND ALL RECONCILIATIONS FUNCTION NEED TO BE REDONE SINCE removenNonCanonical MIGHT NOT CLEAN CORRECTLY.
    //////////////////// SUBOPTIMAL CHECKS ///////////////////////
   
    ////// SL ///////////////
    if( eventName[0] == 'S' && eventName[1] == 'L' && eventSonName[0] == 'T' ) 
    {
        if( eventSonName[1] == 'L' ) {
            // SL followed by a TL, the son of TL cannot 
            // be the sibling of the son of SL (transfer to sibling)
            int eventFatherId = getFatherXId( event );

            MyGraph::Vertex TLson
                = *((mGraph.getAdjacentVertices(eventSon)).first);
            int fatherTLsonId = getSpeciesParentId(TLson);
            
cout << "SL+TL " << eventFatherId << " vs " << fatherTLsonId << endl;
            //  check if TLson is also a child of z (eventFather)
            if( eventFatherId == fatherTLsonId ) 
                valid = false;

        } else if( eventSonName[1] != 'L' ) {

            // SL followed by T, children of T (T1 and T2) cannot be siblings
            MyGraph::adjacency_vertex_range_t sons 
                    = mGraph.getAdjacentVertices(eventSon);

            int fatherT1_x = getSpeciesParentId(*(sons.first) );
            int fatherT2_x = getSpeciesParentId(*(sons.first+1) );
cout << "SL+T_ " << fatherT1_x << " vs " << fatherT2_x << endl;
            // check if T sons are siblings
            if( fatherT1_x == fatherT2_x ) 
                valid = false;
        }
    }


    ////// TL ///////////////
    // TL followed by T, a son of T (T1 and T2) cannot be father of TL 
    // (transfer back to original)
    if( eventName[0] == 'T' && eventName[1] == 'L' 
        && eventSonName[0] == 'T' && eventSonName[1] == '_' ) 
    {
        int eventFather_x = getFatherXId( event );

        MyGraph::adjacency_vertex_range_t sons 
                = mGraph.getAdjacentVertices(eventSon);
        int T1_x = mGraph.properties(*(sons.first)).id_x;
        int T2_x = mGraph.properties(*(sons.first+1)).id_x;
cout << "TL+T_ " << eventFather_x << " == " << T1_x << " or " << T2_x << endl;
        // check if T1_x or T2_x equals event father id_x
        if( T1_x == eventFather_x || T2_x == eventFather_x ) {
            valid = false;
        }
    }
    
*/
    return valid;
}


/**
 * Find the number of valid reconciliations from the two sons of
 * the given event.
 *
 * If not only canoncial, the sons reconciliation are multiplied.
 * For canonical, D or T events with null children are not valid,
 * resulting in zero reconciliations.
 *
 * Furthermore, if suboptimal, certain conditions (e.g. 7),
 * restrict the number of reconciliations (not currently implemented).
 *
 * @return number of reconciliations or zero for invalid
 */
double DTLGraph::validRecNumber( 
        MyGraph::Vertex event, ///< an event vertex
        MyGraph::Vertex eventSon1, ///< the first mapping vertex son of event
        MyGraph::Vertex eventSon2 ) ///< the second mapping vertex son of event
{

    double recNumber = mGraph.properties(eventSon1).recNumber
                       * mGraph.properties(eventSon2).recNumber;

    if( !mOnlyCanonical )
        return recNumber;

    string eventName = mGraph.properties(event).name;
    if( eventName[0] == 'S' ) // two son speciation events are always valid
        return recNumber;

    string eventSon1Name = mGraph.properties(eventSon1).name;
    string eventSon2Name = mGraph.properties(eventSon2).name;

    if( eventName[0] == 'D' || eventName[0] == 'T' ) 
        if( eventSon1Name[0] == 'O' && eventSon2Name[0] == 'O' ) {
//if( mGraph.properties(event).recNumber > 0 )
//cout << "TOO - " << eventName << endl;
            return 0;
        }
/*
    if( !mSuboptimal )
        return recNumber;


    //////////////////// SUBOPTIMAL CHECKS ///////////////////////

    ////// D ///////////////
    if( eventName[0] == 'D' ) {
        // No duplication followed by a TL.
        if( (eventSon1Name[0] == 'T' && eventSon1Name[1] == 'L') 
           || (eventSon2Name[0] == 'T' && eventSon2Name[1] == 'L') )
        {
cout << "D+TL " << eventName << endl;
            return 0;   
        }

        // No duplications followed by two SLs,
        // with children SL as siblings (ok if children are identical).
        if( (eventSon1Name[0] == 'S' && eventSon1Name[1] == 'L') 
            && (eventSon2Name[0] == 'S' && eventSon2Name[1] == 'L') )
        {
            MyGraph::Vertex firstSLson = 
                *((mGraph.getAdjacentVertices(eventSon1)).first);
            MyGraph::Vertex secondSLson = 
                *((mGraph.getAdjacentVertices(eventSon2)).first);

cout << "D+2SL: " << mGraph.properties(firstSLson).id_x << " and " 
<< mGraph.properties(secondSLson).id_x << endl;
            // check if SL sons don't have the same id_x
            if( mGraph.properties(firstSLson).id_x 
                    != mGraph.properties(secondSLson).id_x )
                return 0;
        }
    }

    ////// T ///////////////
    if( eventName[0] == 'T' ) {

        // The father(x) of TL must be transfered to,
        // i.e., not the same as the original x.
        int eventFather_x = getFatherXId( event );

        if( eventSon1Name[0] == 'T' && eventSon1Name[1] == 'L' )  {
            // check if eventSon1 father id_x != event father id_x 
            // (or equivalently eventSon2 father id_x = event father id_x)
            int eventSon1Father_x = getFatherXId( eventSon1 );
            if( eventSon1Father_x != eventFather_x ) {
cout << "T+TL 1 " << eventName << endl;
                return 0;
            }
        }
        if( eventSon2Name[0] == 'T' && eventSon2Name[1] == 'L' ) { 
            // check if eventSon2 father id_x != event father id_x 
            // (or equivalently eventSon1 father id_x = event father id_x)
            int eventSon2Father_x = getFatherXId( eventSon2 );
            if( eventSon2Father_x != eventFather_x ) {
cout << "T+TL 2 " << eventName << endl;
                return 0;
            }
        }


        // condition7 check
        double e21RecNum = condition7( eventSon1, eventSon2,);
        recNumber -= e21RecNum*mGraph.properties(eventSon1).recNumber;
        double e11RecNum = condition7( eventSon2, eventSon1 );
        recNumber -= e11RecNum*mGraph.properties(eventSon2).recNumber;


    }
*/
    return recNumber;
}


/*
// SL and null children of T
// for each child of null, if type is TL
// return recNumber if son of TL is a sibling of the son of SL
double DTLGraph::condition7( MyGraph::Vertex eventSon1, 
            MyGraph::Vertex eventSon2 )
{

    string eventSon1Name = mGraph.properties(eventSon1).name;
    string eventSon2Name = mGraph.properties(eventSon2).name;

    double recNumber = 0;

    if( eventSon1Name[0] == 'S' && eventSon1Name[1] == 'L'
         && eventSon2Name[0] == 'O' ) 
    {

        MyGraph::adjacency_vertex_range_t SLsons 
                        = mGraph.getAdjacentVertices( eventSon1 );
        int SLsonParentX = getSpeciesParentId( *(SLsons.first) );

        MyGraph::adjacency_vertex_range_t sons 
                = mGraph.getAdjacentVertices(eventSon2);
        MyGraph::Vertex nullSon = *(sons.first); // son of null event
        // check z21 sons for a TL
        MyGraph::adjacency_vertex_range_t nullGrandsons 
                = mGraph.getAdjacentVertices(nullSon);
        BOOST_FOREACH( MyGraph::Vertex nullGrandson, nullGrandsons) {
            string nullGrandsonName = mGraph.properties(nullGrandson).name;
            if( nullGrandsonName[0] == 'T' && nullGrandsonName[1] == 'L' )
            {
                MyGraph::adjacency_vertex_range_t TLsons 
                        = mGraph.getAdjacentVertices(nullGrandson);
                int TLsonParentX = getSpeciesParentId( *(TLsons.first) );

                // check if TLson id_x is a sibling of eventSon1 father
                if( SLsonParentX == TLsonParentX ) {
cout << "COND7 " << SLsonParentX << " vs " << TLsonParentX << endl;
                    recNumber = mGraph.properties(nullGrandson).recNumber;
                    break;
                }
            }
        }
        
    }


    return recNumber;
}

*/

/**
 * Recursively find orthologous gene pairs.
 *
 * @return List of descedent genes by id.
 */
vector<int> DTLGraph::orthologyOutputAux(
    int idU, ///< current gene id
    vector< vector<MyGraph::Vertex> > &reconciliation, ///< reconciliation
    vector< pair<int,int> > &orthologs ) ///< resulting pairs
{
    int lastIdx = reconciliation[idU].size() - 1;
    string eventStr = mGraph.properties( reconciliation[idU][lastIdx] ).name;
    pair<int,int> p = mCladesTrips->getCladeSplit( idU, 0 );
    vector<int> children;
    if( p.first != -1 ) {
        children = orthologyOutputAux( p.first, reconciliation, orthologs );
        vector<int> children2 = 
                orthologyOutputAux( p.second, reconciliation, orthologs );

        if( eventStr[0] == 'S' || eventStr[0] == 'I' ) {
            BOOST_FOREACH( int child1, children )
                BOOST_FOREACH( int child2, children2 ) {
                    orthologs.push_back( make_pair( child1, child2 ) );
                }
        }

        BOOST_FOREACH( int child, children2 )
            children.push_back( child );
    }
    children.push_back( idU );

    return children;
}


/**
 * Print pairs of orthologous gene ids to a file.
 *
 * Orthologous gene pairs are those whose least common ancestor
 * map to a speciation event.
 */
void DTLGraph::orthologyOutput(
    string fileName ) ///< base output file name
{
    bool random = false;
    int problem = 3; // asymmetric
    vector< vector<MyGraph::Vertex> > reconciliation;
    if( !getScoredReconciliation( problem, 
            mCladesTrips->mClades.getCladeCount(), reconciliation, random ))
    {
        // score infinite, no reconciliation
        return;
    }

    vector< pair<int,int> > orthologs;
    orthologyOutputAux( mCladesTrips->mClades.getRootClade(), reconciliation,
                        orthologs );

    ofstream outFile;
    outFile.open(fileName.c_str() );
    vector<int> cladeToPOrd = mCladesTrips->getPostOrderMapping();
    pair<int,int> p;
    // only leaf pairs!
    BOOST_FOREACH( p, orthologs ) {
        // by id
        //outFile << cladeToPOrd[p.first] 
        //     << "," << cladeToPOrd[p.second] << endl;
        // leaf pairs, by name
        if( mCladesTrips->mClades.isLeaf( p.first ) 
            && mCladesTrips->mClades.isLeaf( p.second ) )
        {
            outFile << mCladesTrips->mClades.getLeafName(p.first) 
                    << "," << mCladesTrips->mClades.getLeafName(p.second)  
                    << endl;
        }
    }
    outFile.close();
}


/**
 * Print one or all problem reconciliations to a file.
 */
void DTLGraph::printReconciliation( 
    string problemStr, 
        ///< which reconciliation: symmetric, asymmetric, random, or all
    string fileName, ///< base output file name
    bool sylvxFormat, ///< use the Sylvx format for the reconciliation
    bool recPhyloXMLFormat, ///< use the recPhyloXMLFormat for the reconciliation
    bool checkConsistent, ///< check consistency of reconciliation
    bool &isConsistent,  ///< returns true if at least one is consistent
    map<string,double> &eventSupports ) ///< event supports to use
{
    vector<string> problems;
    if( problemStr == "all" || problemStr == "allTriplets" ) {
        problems.push_back( "symmetric" );
        problems.push_back( "asymmetric" );
        problems.push_back( "random" );
    } else 
        problems.push_back( problemStr );

    if( problemStr == "allTriplets" ) 
        problems.push_back( "maxAvg" );

    isConsistent = false;
    BOOST_FOREACH( string curProblemStr, problems ) {
        bool random = false;
        int problem = 3;
        if( curProblemStr == "random" ) 
            random = true;
        else if( curProblemStr == "asymmetric" ) 
            problem = 1;
        else if( curProblemStr == "maxAvg" ) 
            problem = 5;
        else if( curProblemStr != "symmetric" ) {
            cout << "Given problem " << curProblemStr << endl;
            throw bpp::Exception("DTLGraph::printReconciliation:"
                    " Given unknown problem" );
        }

        // create reconcilation
        vector< vector<MyGraph::Vertex> > reconciliation;
        string otherStr = "";
        if( getScoredReconciliation( problem, 
                mCladesTrips->mClades.getCladeCount(), reconciliation, random ))
        {
            vector<string> reconStrings; 
            string problemFileName = fileName + "_" + curProblemStr;
// move consistency check here and pass ordering to sylvx if not dated?
            if( sylvxFormat ) {
                reconStrings = getSylvxReconciliation( reconciliation );
                reverse( reconStrings.begin(), reconStrings.end() );
                problemFileName += ".mr";
            } 
            else if (recPhyloXMLFormat) {
                reconStrings = getRecPhyloXMLReconciliation( reconciliation );
                problemFileName += ".recPhyloXML";
            }
            else {
                reconStrings = printReconciliationAux( 
                                reconciliation, eventSupports );
                problemFileName += ".txt";
            }

            if( checkConsistent ) {
                if( isTimeConsistent( reconciliation ) ) {
                    isConsistent = true;
                    otherStr += ", consistent";
                } else {
                    otherStr += ", not consistent";
                }
            }

            // merge clade strings using pOrd ordering and print them to file
            if( fileName != "" ) {
                ofstream outFile;
                outFile.open( problemFileName.c_str() );
                for( int i = reconStrings.size()-1; i>=0; i-- ) {
                    if( reconStrings[i] != "" ) 
                        outFile << reconStrings[i] << endl;
                }
                outFile.close();
            }
        } else {
            // score infinite, no reconciliation
            if( checkConsistent ) 
                otherStr += ", consistency unknown";
        }

        double score = -1;
        if( random )
            score = scoreReconciliation( reconciliation );
        else
            score = getBestScore( problem );

        if( curProblemStr != "random" )
            curProblemStr += " median";
        cout << curProblemStr << " score: " << score 
             << otherStr << endl;
    }  
}
/*
void printGTree(
    vector<int> cladeToPOrd,
    int idU = 0,
    int depth = 0 )
{
    for( int i=0; i<depth; i++ )
        cout << " ";
    int pOrd = cladeToPOrd[idU];

    pair<int,int> cladeSplit = mCladesTrips->getCladeSplit( idU, 0 );
    int idUl = cladeSplit.first;
    int idUr = cladeSplit.second;
    cout << idU << " (" << pOrd << ") --> " << idUl << ", " << idUr << endl;

    if( idUl != -1 ) {
        printGTree( cladeToPOrd, idUl, depth+1 );
        printGTree( cladeToPOrd, idUr, depth+1 );
    }
}*/

bool debugSylvx = false;

double DTLGraph::getSpeciesDates(
        MySpeciesNode *node,    ///< species node
        vector<double> &startDates, ///< start date (earliest/highest)
        vector<double> &endDates ) ///< end dates for species
{
    int idX = node->getId();
    double len = 0;
    for( size_t i=0; i<node->getNumberOfSons(); i++ )
        // all child lengths should be the same for an ultrametric tree
        len = getSpeciesDates( node->getSon( i ), startDates, endDates );

    endDates[idX] = len;
    if( node->hasDistanceToFather() )
        len += node->getDistanceToFather();
    startDates[idX] = len;

    if( mSpeciesTree->isAlpha( idX ) ) 
        endDates[idX] = 0;
if( debugSylvx && node->getNumberOfSons() != 1 )
cout << idX << " " << startDates[idX] << "-" << endDates[idX]  << endl;

    return len;
}



int DTLGraph::getChildIdx(
        int idU,
        vector< vector<MyGraph::Vertex> > &reconciliation, 
        size_t reconciliationIdx ) ///< current event in reconciliation
{
    MyGraph::Vertex mappingV = 0;
    string name;
    while( reconciliationIdx < reconciliation[idU].size() ) {
        MyGraph::Vertex v = reconciliation[idU][reconciliationIdx++];
        if( mGraph.properties(v).isMapping ) {
            mappingV = v; // save mapping vertex 
        } else {
            name = mGraph.properties(v).name;
            if( name[0] != 'O' )
                break; // skip artificial nodes
        }
    }

    // get real idX (could be an artificial node id)
    int idX = mGraph.properties(mappingV).id_x;
    MySpeciesNode *node = mSpeciesTree->getNodeById( idX );
    return node->getId();
}

string DTLGraph::makeIntervals(
        int idU,
        vector< vector<MyGraph::Vertex> > &reconciliation,
        vector<double> &speciesStartDates, ///< input species start times
        vector<double> &speciesEndDates, ///< input species end times
        map<string,double> &eventDates, ///< output event end times
        size_t reconciliationIdx, ///< current event in reconciliation
        double parentStartDate, ///< defaults to -1 for root
        string parentEvent,  ///< defaults to R for root
        int prevRpo, ///< previous species real post order
        int seqNum ) ///< number of event in current species
{
    // event name and mapping vertex
    MyGraph::Vertex mappingV = 0;
    string name;
    while( reconciliationIdx < reconciliation[idU].size() ) {
        MyGraph::Vertex v = reconciliation[idU][reconciliationIdx++];
        if( mGraph.properties(v).isMapping ) {
            mappingV = v; // save mapping vertex 
        } else {
            name = mGraph.properties(v).name;
            if( name[0] != 'O' )
                break; // skip artificial nodes
        }
    }
   
    // get real idX (could be an artificial node id)
    int idX = mGraph.properties(mappingV).id_x;
    MySpeciesNode *node = mSpeciesTree->getNodeById( idX );
    idX = node->getId();

    // leaf, return name
    if( name[0] == 'C' ) {
        eventDates[name] = 0;
if( debugSylvx )
cout << name << " (" << parentEvent << ")   (" 
     << speciesStartDates[idX] << "-" << speciesEndDates[idX] << ") " 
     << node->getName() << endl;
        return name;
    }

    // update number events in this species (seqNum)
    double endDate =  speciesEndDates[idX];
    int rpo = node->getInfos().realPostOrder;
    if( rpo == prevRpo )
        seqNum++;
    else
        seqNum = 1; 

    // initialize startDate to species start or parent start
    double startDate = speciesStartDates[idX];
    if( parentStartDate == -1 )  // for root
        eventDates[parentEvent] = startDate;
    else if( startDate > parentStartDate ) 
        startDate = parentStartDate;
    

    // find gene id of next event
    int nextU = idU;
    pair<int,int> p;
    if( reconciliationIdx == reconciliation[idU].size() ) {
        // split event (S,D,T)
        p = mCladesTrips->getCladeSplit( idU, 0 );
        nextU = p.first;
        reconciliationIdx = 0;
    }

    double maxToChild = speciesStartDates[idX];
    if( maxToChild > startDate ) 
        maxToChild = startDate;

    if( name[0] == 'T' ) {
        // get start dates of children - must go to at least there
        int child1idX = getChildIdx( nextU, reconciliation, 
                                    reconciliationIdx );
        if( maxToChild > speciesStartDates[child1idX] ) 
            maxToChild = speciesStartDates[child1idX];
        if( reconciliationIdx == 0 ) {
            int child2idX = getChildIdx( p.second, reconciliation,
                                         reconciliationIdx);
            if( maxToChild > speciesStartDates[child2idX] ) 
                maxToChild = speciesStartDates[child2idX];
        }
        startDate =  maxToChild;
    }

    string child1Name = makeIntervals( nextU, reconciliation, 
        speciesStartDates, speciesEndDates,
        eventDates, 
        reconciliationIdx, maxToChild, name,
        rpo, seqNum );

    string child2Name;
    if( reconciliationIdx == 0 ) {
        nextU = p.second;
        child2Name = makeIntervals( nextU, reconciliation, 
            speciesStartDates, speciesEndDates,
            eventDates, 
            0, maxToChild, name,
            rpo, seqNum );
    }

/*
    if( seqNum > 1 && name[0] == 'T' ) {
// This should be a percentage.
// Make sure it doesn't move out of species.
// OR - put as part of a collision detection step.
        endDate -= seqNum*0.2;    
cout << name << " SAME TUBE endDate-> " << endDate << endl; } */
    
    if( name[0] == 'T' || name[0] == 'D' || name[0] == 'I' ) {
        double childEnd = eventDates[child1Name];
        if( childEnd < eventDates[child2Name] )
            childEnd = eventDates[child2Name];
        if( childEnd < speciesEndDates[idX] )
            childEnd = speciesEndDates[idX];
        double dist = startDate-childEnd;
if( dist < 0 ) {
    cout << "====Impossible interval for T or D event =====" << endl;
    cout << seqNum << " " << name << endl;
    cout << "species interval (" << idX << ") = " << speciesStartDates[idX]
        << " - " << speciesEndDates[idX] << endl;
    cout << "startDate = " << startDate << endl;
    cout << "child1 event = " << eventDates[child1Name] << endl;
    cout << "child2 event = " << eventDates[child2Name] << endl;
    throw bpp::Exception( "DTLGraph::makeIntervals: impossible interval" );
}
        endDate = dist/(seqNum+1) + childEnd;
if( debugSylvx )
cout << " " << seqNum << " " << endDate << " from "
    << parentStartDate << " - " << childEnd << endl;
    }

if( debugSylvx )
cout << name << " (" << parentEvent << ") " << endDate 
    << "  (" << speciesStartDates[idX] << "-" << speciesEndDates[idX] 
        << ")" << endl;

    eventDates[name] = endDate;

    return name;
}

void DTLGraph::printEvents(
        int idU, ///< current gene id
        vector<int> &cladeToPOrd, ///< clade mapping
        vector< vector<MyGraph::Vertex> > &reconciliation, 
            ///< the reconciliation
        vector<double> &speciesEndDates, ///< input species end times
        map<string,double> &eventDates, ///< output event end times
        vector< vector<string> > &eventsByIdU, ///< output strings
        vector<string> &losses, ///< loss strings
        size_t reconciliationIdx, ///< current event in reconciliation
        int fatherX,
        string fatherEvent,  ///< defaults to R for root
        string prevRealX, ///< previous species id
        string prevRealFatherX ) ///< previous father species id
{
    // get the next non-artificial mapping vertex and event name
    MyGraph::Vertex mappingV = 0;
    string name;
    while( reconciliationIdx < reconciliation[idU].size() ) {
        MyGraph::Vertex v = reconciliation[idU][reconciliationIdx++];
        if( mGraph.properties(v).isMapping ) {
            mappingV = v; // save mapping vertex 
        } else {
            name = mGraph.properties(v).name;
            if( name[0] != 'O' )
                break; // skip artificial nodes
        }
    }

    // get the species info
    int idX = mGraph.properties(mappingV).id_x;
    MySpeciesNode *node = mSpeciesTree->getNodeById( idX);
    idX = node->getId(); // real post order

    double startDate = eventDates[fatherEvent];
    double endDate =  eventDates[name];

if( fatherEvent != "ROOT" && startDate < endDate ) {
    cout << name << " invalid dates [" << startDate << "-"
         << endDate << "]" << endl;
    throw bpp::Exception( "DTLGraph::printEvents: invalid dates" );
}
    string realX;
    if( node->isLeaf() ) 
        realX = node->getName();
    else 
        realX = bpp::TextTools::toString( node->getInfos().realPostOrder );

    string realFatherX;
    if( realX == "OUTGROUP" ) {
        realFatherX = "OUTGROUP";
    } else {
        MySpeciesNode *node = mSpeciesTree->getNodeById( 
                                       mGraph.properties(mappingV).id_x ); 
        if( node->hasFather() )
            node = node->getFather();

        realFatherX = bpp::TextTools::toString( 
                                node->getInfos().realPostOrder );
    }

    // real gene id
    string pOrd; 
    if( mCladesTrips->mClades.isLeaf( idU ) )
        pOrd = mCladesTrips->mClades.getLeafName( idU );
    else
        pOrd = bpp::TextTools::toString( cladeToPOrd[idU] );


    // print gene line and recursion
    if( name[0] != 'C' ) {
        int nextU = idU;
        pair<int,int> p;
        if( reconciliationIdx == reconciliation[idU].size() ) {
            // split event (I,S,D,T)
            p = mCladesTrips->getCladeSplit( idU, 0 );
            nextU = p.first;
            reconciliationIdx = 0;
        }
        string childPord;
        if( mCladesTrips->mClades.isLeaf( nextU ) )
            childPord = mCladesTrips->mClades.getLeafName( nextU );
        else
            childPord = bpp::TextTools::toString( cladeToPOrd[nextU] );

        // recursion 
        printEvents( nextU, cladeToPOrd, reconciliation, 
            speciesEndDates, eventDates, eventsByIdU, losses,
            reconciliationIdx, idX, name, realX, realFatherX );

        // print the gene id lines if this is the first event
        if( reconciliationIdx == 0 ) {
            eventsByIdU[nextU].push_back( "{" + pOrd + "," + childPord + "}" );
if( debugSylvx )
cout << "{" << pOrd << "," << childPord << "}  idU=" << nextU << endl;
        }


        if( reconciliationIdx == 0 ) {
            // recursion on second gene child
            nextU = p.second;
            if( mCladesTrips->mClades.isLeaf( nextU ) )
                childPord = mCladesTrips->mClades.getLeafName( nextU );
            else
                childPord = bpp::TextTools::toString( cladeToPOrd[nextU] );
            printEvents( nextU, cladeToPOrd, reconciliation, 
                speciesEndDates, eventDates, eventsByIdU, losses,
                0, idX, name, realX, realFatherX );
            eventsByIdU[nextU].push_back( "{" + pOrd + "," + childPord + "}" );
if( debugSylvx )
cout << "{" << pOrd << "," << childPord << "}  idU=" << nextU << endl;
        }
    }

    // Add events to make up for ILS species skipping.
    string realNextX = realX;
    string realNextFatherX = realFatherX;
    string nextEnd =  bpp::TextTools::toString( endDate );
    string nextStart =  bpp::TextTools::toString( startDate );
    int nextX = idX;
    if( node->hasFather() )
        nextX = node->getFather()->getId();
    if( fatherEvent != "ROOT" && fatherEvent[0] != 'T' && idX != fatherX && nextX != fatherX ) {
        MySpeciesNode *nextNode = node;
        while( nextNode->getId() != fatherX && nextNode->hasFather() ) {
            if( !nextNode->hasFather() ) 
                throw bpp::Exception( 
                        "DTLGraph::printEvents: father not found for ILS" );
            // get father info 
            nextNode = nextNode->getFather();
            string fatherStart = bpp::TextTools::toString( 
                                    speciesEndDates[nextNode->getId()] ); 
            string realNextFatherX = bpp::TextTools::toString( 
                                nextNode->getInfos().realPostOrder );

            string event = "\t(" + realNextFatherX + "," + realNextX + ")"
                        + " [" + fatherStart
                        + "-" + nextEnd + "]";
if( debugSylvx )
cout << "ILS ADD: " << event << endl;
            eventsByIdU[idU].push_back( event );
            realNextX = realNextFatherX;
            nextEnd = fatherStart;
        }
//realNextFatherX = prevRealFatherX;
    } else {
        string event = "\t(" + realNextFatherX + "," + realNextX + ")"
                    + " [" + nextStart + "-" + nextEnd + "]";
        eventsByIdU[idU].push_back( event );
if( debugSylvx )
cout << "  " << name << " "  << event << endl;
    }
    


    if( fatherEvent[0] == 'T' && realX != prevRealX ) {
        string event = "\t(" + prevRealFatherX + "," + prevRealX 
                    + ") [" + bpp::TextTools::toString( startDate ) 
                    + "-" + bpp::TextTools::toString( startDate ) 
                    + ":" + bpp::TextTools::toString( startDate ) + "]";
        eventsByIdU[idU].push_back( event );
if( debugSylvx )
cout << "  " << name << "-Ta " << event << endl;
    }

    if( fatherEvent[1] == 'L' ) {
        string lostXstr;
        string fatherXstr;
        int lostX;
        if( fatherEvent[0] == 'S' || fatherEvent[0] == 'I' ) {
            fatherXstr = realFatherX;
            lostX = getSibling( fatherX, idX, true );
            MySpeciesNode *node = mSpeciesTree->getNodeById( lostX );
            if( node->isLeaf() ) 
                lostXstr = node->getName();
            else 
                lostXstr = bpp::TextTools::toString( 
                                        node->getInfos().realPostOrder );
        } else {
            lostXstr = prevRealX;
            fatherXstr = prevRealFatherX;
            lostX = fatherX;
        }
if( debugSylvx )
cout << name << " Loss from father " << fatherEvent 
     << " lostX=" << lostX  << " lostXstr=" << lostXstr
     << " startDate = " << startDate
     << " endDate=" << speciesEndDates[lostX] << endl;
        double lossEndDate = 4*(startDate - speciesEndDates[lostX])/5
                           + speciesEndDates[lostX];

        // add losses
        string event =  "{" + pOrd + "," + pOrd + "-" + lostXstr + "}";
if( debugSylvx )
cout << "loss: " << event << endl;
        losses.push_back( event );
        //losses.push_back( "{" + pOrd + "," + pOrd + "-}" );
        event = "\t(" + fatherXstr + "," + lostXstr + ")"
                    + " [" + bpp::TextTools::toString( startDate )
                    + "-" + bpp::TextTools::toString( lossEndDate ) + "]";
if( debugSylvx )
cout << "loss: " << event << endl;
        losses.push_back( event );
    }
}



/**
 * Create a string representation of a reconciliation.
 * in the RecPhyloXML format (see http://phylariane.univ-lyon1.fr/recphyloxml/
 * for more details.
 */


vector<string> DTLGraph::getRecPhyloXMLReconciliation(
    vector< vector<MyGraph::Vertex> > &reconciliation
        ///< list of vertices, by clade in the reconcilation
        )
{
	vector<string> temp;
	vector<int> cladeToPOrd = mCladesTrips->getPostOrderMapping();
	string geneTree_recPhyloXML_format= mSpeciesTree->MySpeciesTree::toRecPhyloXML() ;
	geneTree_recPhyloXML_format= geneTree_recPhyloXML_format + "<recGeneTree>\n<phylogeny rooted=\"true\">\n";
	geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + DTLGraph::getRecPhyloXMLReconciliation(reconciliation,-1, mCladesTrips->mClades.getRootClade(), cladeToPOrd);
	geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "</phylogeny>\n</recGeneTree>\n</recPhylo>";
	temp.push_back(geneTree_recPhyloXML_format);
	//cout << temp[0] << endl;
	//cout << geneTree_recPhyloXML_format << endl;
	return temp; 
};

string DTLGraph::getRecPhyloXMLReconciliation(
    vector< vector<MyGraph::Vertex> > &reconciliation,
        ///< list of vertices, by clade in the reconcilation
        int idUP,
        int idU,
     vector< int>  &cladeToPOrd  ///< map of clade ids to DF post-order numbering
    )
{    
	
	//    int idX = mGraph.properties(mappingV).id_x;
   // MySpeciesNode *node = mSpeciesTree->getNodeById( idX );
    //idX = node->getId();
	// string eventName = mGraph.properties(reconciliation[idU][z+1]).name;
	//            string eventStr=getEventString(idU, z, reconciliation, idUl, idUr,
    //                                        eventSupports );
	
	int reconciliationIdx =0;
	MyGraph::Vertex mappingV = 0; 
	int reconciliationIdxMappingV=-1;
	int eventNumber=0;

	string name;
	
	bool isLeaf = mCladesTrips->mClades.isLeaf( idU );
	string geneTree_recPhyloXML_format;
	
	int lossId=0;
			
	bool doneWithIdU=false;
	//cout << "idU " << idU << endl;
	while (!doneWithIdU){ 
		while( reconciliationIdx < reconciliation[idU].size() ) {
			MyGraph::Vertex  v = reconciliation[idU][reconciliationIdx++];
			if( mGraph.properties(v).isMapping ) {
				mappingV = v; // save mapping vertex 
				reconciliationIdxMappingV = reconciliationIdx-1;
			} else {
				name = mGraph.properties(v).name;
				if( name[0] != 'O' )
					break; // skip artificial nodes	
			}
		}
		eventNumber++;
   
		string idToPrint;
	
		string event;
		
		if (name.size()>2 && name[2]=='D')
			event = name.substr(0,3); //TFD TTD
		else if (name.size()>3 && name[3]=='D')
			event = name.substr(0,4); //TLFD TLTD
		else
			event = name.substr(0,2);
	
				
		geneTree_recPhyloXML_format= geneTree_recPhyloXML_format + "<clade>\n<name>";
	    
		if(event=="C_" || event=="S_" || event=="DD" || event=="D_" || event=="T_" || event=="TTD" || event=="TFD" ){
			
        	if( isLeaf )
            	geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + mCladesTrips->mClades.getLeafName(idU) ;
       		else
            	geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + bpp::TextTools::toString( cladeToPOrd[idU] ) ;   
                           
		}
		//else if (event=="SL"){
		//	lossId ++;  
		//	geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + event + bpp::TextTools::toString(lossId) + "_" + bpp::TextTools::toString(cladeToPOrd[idU]) + "_" + bpp::TextTools::toString(mSpeciesTree->getNodeById(mGraph.properties(mappingV).id_x)->getId()) ;
		//	
		//}
		else if (event=="SL" || event=="TL" || event=="TLFD" || event=="TLTD"){
			lossId ++;  
			geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + event + "_" + bpp::TextTools::toString(lossId) + "_" + bpp::TextTools::toString(cladeToPOrd[idU]) + "_" + bpp::TextTools::toString(mSpeciesTree->getNodeById(mGraph.properties(mappingV).id_x)->getId()) ;	
		}
		
		geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "</name>\n";
	
		geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "<eventsRec>\n";
				
		if(mSpeciesTree->getNodeById(mGraph.properties(mappingV).id_x)->isLeaf())
			idToPrint=mSpeciesTree->getNodeById(mGraph.properties(mappingV).id_x)->getName() ;
		else
			idToPrint=bpp::TextTools::toString(mSpeciesTree->getNodeById(mGraph.properties(mappingV).id_x)->getId()) ;

		
		bool previousEventIsTransfert=false;
		string previousEvent;
		int keptSonT =-1;
		
		if(idUP!=-1 && eventNumber==1){ // the previous event is associated to the parent of idU, ie idUP
			string namePreviousEvent = mGraph.properties(reconciliation[idUP][ reconciliation[idUP].size()-1]).name;
			if (namePreviousEvent.size()>2 && namePreviousEvent[2]=='D')
				previousEvent = namePreviousEvent.substr(0,3); //TFD TTD
			else if (namePreviousEvent.size()>2 && namePreviousEvent[3]=='D')
				previousEvent = namePreviousEvent.substr(0,4); //TLFD TLTD
			else
				previousEvent = namePreviousEvent.substr(0,2);
			
			keptSonT = mGraph.properties(reconciliation[idUP][reconciliation[idUP].size()-2]).id_x; 

		}	
		else if (reconciliationIdxMappingV-2 >=0){ // the previous event is associated to idU
			
			string namePreviousEvent = mGraph.properties(reconciliation[idU][reconciliationIdxMappingV-1]).name;
			if (namePreviousEvent[2]=='D')
				previousEvent = namePreviousEvent.substr(0,3); //TFD TTD
			else if (namePreviousEvent[3]=='D')
				previousEvent = namePreviousEvent.substr(0,4); //TLFD TLTD
			else
				previousEvent = namePreviousEvent.substr(0,2);	
		}
		
		if(previousEvent=="TL"  || previousEvent=="TLFD"  || event=="TLTD" || ((previousEvent=="T_" || previousEvent=="TTD" || previousEvent=="TFD") && mGraph.properties(mappingV).id_x != keptSonT) ){// in TL, the keptSon is always id_x while in T is the sibling of the keptOne
			geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "<transferBack destinationSpecies=\""  + idToPrint + "\"></transferBack>\n";
		}	
			
		
		//cout << "previousEvent " << previousEvent << endl;
	
				 		
		 if(event=="C_" ){
			geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "<leaf speciesLocation=\"" + idToPrint + "\" geneName=\"" +mCladesTrips->mClades.getLeafName(idU)  + "\"></leaf>";
		}
		 else if(event=="S_" || event=="SL"){
			geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "<speciation speciesLocation=\""  + idToPrint + "\"></speciation>";
		}
		 else if(event=="D_" || event=="DD"){
			geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "<duplication speciesLocation=\""  + idToPrint + "\"></duplication>";
		}
		else if(event=="T_" || event=="TL" || event=="TTD" || event=="TFD" || event=="TLTD" || event=="TLFD"){
			geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "<branchingOut speciesLocation=\""  + idToPrint + "\"></branchingOut>";
		}
		
		geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "\n</eventsRec>\n";
		
		 if(event=="SL" || event=="TL" || event=="TLFD" || event=="TLTD" ){ // add the loss
			int lostSonX =-1;
		 	if(event=="SL" ){
				MyGraph::Vertex vNext = reconciliation[idU][reconciliationIdxMappingV+2];
        		int keptSonX = mGraph.properties(vNext).id_x; 
        		lostSonX = getSibling( mGraph.properties(reconciliation[idU][reconciliationIdxMappingV]).id_x, keptSonX, false );
            }
            else
            	lostSonX = mGraph.properties(mappingV).id_x;
                        
            geneTree_recPhyloXML_format= geneTree_recPhyloXML_format + "<clade>\n<name>LOST</name>\n<eventsRec>\n<loss speciesLocation=\"";
			string idToPrintLoss;
			
			if(mSpeciesTree->getNodeById(lostSonX)->isLeaf())
				idToPrintLoss=mSpeciesTree->getNodeById(lostSonX)->getName() ;
			else
				idToPrintLoss=bpp::TextTools::toString(mSpeciesTree->getNodeById(lostSonX)->getId()) ;
			
			geneTree_recPhyloXML_format = geneTree_recPhyloXML_format + idToPrintLoss;

            geneTree_recPhyloXML_format= geneTree_recPhyloXML_format + "\"></loss>\n</eventsRec>\n</clade>\n";
		}
		
		if(reconciliationIdx== reconciliation[idU].size())
			doneWithIdU=true;
			
	}
	

	if( ! isLeaf ){
		pair<int,int> cladeSplit = mCladesTrips->getCladeSplit( idU, 0 );
    	int idUl = cladeSplit.first;
    	int idUr = cladeSplit.second;

		geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + DTLGraph::getRecPhyloXMLReconciliation(reconciliation, idU, idUl, cladeToPOrd);
		geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + DTLGraph::getRecPhyloXMLReconciliation(reconciliation, idU, idUr, cladeToPOrd);
	}
		
	geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "</clade>\n";

	for (int l=1; l <=lossId ; l++){
		geneTree_recPhyloXML_format=geneTree_recPhyloXML_format + "</clade>\n"; //closing the clade for the SL and TL events
	}
	return geneTree_recPhyloXML_format;
};




vector<string> DTLGraph::getSylvxReconciliation(
    vector< vector<MyGraph::Vertex> > &reconciliation )
{
    int nodeCount = mSpeciesTree->getNumberOfNodes();
    vector<double> speciesStartDates( nodeCount );
    vector<double> speciesEndDates( nodeCount );
    if( mSpeciesTree->isSubdivided() ) {
        getSpeciesDates( mSpeciesTree->getRootNode(), speciesStartDates, 
                         speciesEndDates );
    } else {
        // 1. check consistency
// OR DO THIS BEFORE CALLING THIS FUNCTION (need ordering)
        // 2. create a node ordering
        // 3. set species start end dates
        // 4. print a species tree (what file name? or return a string)
// What about partial orderings?
// For all reconciliations, need to print all species trees.
    }

    // makeIntervals: connected transfers and event intervals
    vector<int> cladeToPOrd = mCladesTrips->getPostOrderMapping();
    int cladeCount = cladeToPOrd.size();
    map<string,double> eventDates;
    int rootClade = mCladesTrips->mClades.getRootClade();
    makeIntervals( rootClade, reconciliation, speciesStartDates, 
        speciesEndDates, eventDates );


    vector< vector<string> > eventsByIdU( cladeCount );
    vector<string> losses;
    printEvents( rootClade, cladeToPOrd, reconciliation,
        speciesEndDates, eventDates, eventsByIdU, losses );



    vector<int> pOrdToClade( cladeCount );
    for( int i=0; i<cladeCount; i++ ) 
        pOrdToClade[cladeToPOrd[i]] = i;

    // create a string for each clade
    vector<string> lines;
    lines.push_back( "------------------------" );
// add more info?
    lines.push_back( "Rec teraRec" );
	for( int pOrd = cladeCount-1; pOrd >=0; pOrd-- ) {
        int idU = pOrdToClade[pOrd];
        for( size_t i=eventsByIdU[idU].size(); i>0; i-- ) {
            string str = eventsByIdU[idU][i-1];
//cout << str << endl;
if( eventsByIdU[idU].size() > 1 )
            lines.push_back( str );
else if( idU != rootClade ) {
    cerr << "*** ERROR NO EVENTS" << endl;
    cout << "*** NO EVENTS FOR idU " << idU << endl;
    exit(1);
}
        }

    }

// use a vector merge
    BOOST_FOREACH( string str, losses ) 
        lines.push_back( str );

    return lines;
}




/**
 * Create a string representation of a reconciliation.
 *
 * Each line of the reconciliation is a clade. Ids are postorder and 
 * not internal ids. For clade leaves, the leaf name is given
 * rather than the postorder id.
 *
 * The input reconciliation are lists of vertices, by clade.
 *
 * Format:
 * cladeId: [event];...;[event]: cladeSonId1, cladeSonId2
 * [event] = "speciesId,eventName,[aux]@eventSupport"
 * [aux] = speciesSonId1,speciesSonId2 for S and T events
 *       = lost species id for SL events
 *
 * @return A reconciliation as a list of strings.
 */
vector<string> DTLGraph::printReconciliationAux(
    vector< vector<MyGraph::Vertex> > &reconciliation, 
        ///< list of vertices, by clade in the reconcilation,
    map<string,double> &eventSupports )
        ///< supports to use if this is not empty
{
    vector<int> cladeToPOrd = mCladesTrips->getPostOrderMapping();
    int cladeCount = cladeToPOrd.size();
    vector<int> pOrdToClade( cladeCount );
    for( int i=0; i<cladeCount; i++ ) 
        pOrdToClade[cladeToPOrd[i]] = i;
    int rootClade = mCladesTrips->mClades.getRootClade();

    // create a string for each clade
	vector<string> outStrings( cladeCount, "" ); 
    vector<string> firstEvents( cladeCount, "" );
	for( int pOrd = cladeCount-1; pOrd >=0; pOrd-- ) {
        int idU = pOrdToClade[pOrd];

        // start with clade name
        string cladeString;
        bool isLeaf = mCladesTrips->mClades.isLeaf( idU );
        if( isLeaf )
            cladeString = mCladesTrips->mClades.getLeafName(idU) + ":";
        else
            cladeString = bpp::TextTools::toString( pOrd ) + ":";

        pair<int,int> cladeSplit = mCladesTrips->getCladeSplit( idU, 0 );
        int idUl = cladeSplit.first;
        int idUr = cladeSplit.second;

        // add last split event if not root
        string semicolon = "";
        if( idU != rootClade ) {
            if( firstEvents[idU] == "" ) {
                cout << idU << " not set" << endl;
                throw bpp::Exception("DTLGraph::printReconciliationAux: "
                        "firstEvents not set" );
            }
            cladeString += firstEvents[idU];
            semicolon = ";";
        }

        // create event strings
        int eventsAdded = 0;
	    for( size_t z=0; z<reconciliation[idU].size(); z+=2 ) {
            string eventName = mGraph.properties(reconciliation[idU][z+1]).name;
            // sanity checks (last event must be C,D,T, or S )
            if( z+2 == reconciliation[idU].size() ) {
                if( eventName[1] == 'L' || eventName[0] == 'O' )
                    throw bpp::Exception( "DTLGraph::printReconciliationAux: "
                            " last event is not a leaf or split" );
            } else {
                if( eventName[1] != 'L' && eventName[0] != 'O' ) {
                    cout << "Bad mid event " << eventName 
                         << " for idU=" << idU << endl;
                    throw bpp::Exception( "DTLGraph::printReconciliationAux: "
                            "mid event is a leaf or split" );
                }
            }

            if( eventName[0] == 'O' || eventName[0] == 'C'  
              ) //|| (eventName[0] == 'D' && eventName[1] == 'D') )
                continue;

            // event description
            string eventStr=getEventString(idU, z, reconciliation, idUl, idUr,
                                            eventSupports );
            cladeString += semicolon + eventStr;
            eventsAdded++;
            semicolon = ";"; // separate following events

            // last split is first split for children
            if( idUl != -1 ) {
                firstEvents[idUl] = eventStr;
                firstEvents[idUr] = eventStr;
            }
        }
        if( eventsAdded == 0 )
            continue; // no events, skip

        // add child clades
        if( !isLeaf ) {
            if( mCladesTrips->mClades.isLeaf( idUl ) )
                cladeString += ":" + mCladesTrips->mClades.getLeafName( idUl );
            else
                cladeString += ":" 
                            + bpp::TextTools::toString( cladeToPOrd[idUl] );

            if( mCladesTrips->mClades.isLeaf( idUr ) )
                cladeString += "," + mCladesTrips->mClades.getLeafName( idUr );
            else
                cladeString += "," 
                            + bpp::TextTools::toString( cladeToPOrd[idUr] );
        }
        outStrings[pOrd] = cladeString;
    }

    return outStrings;
}

/**
 * Get leaf name for a node by id or it's post order id.
 */
string DTLGraph::getStringId( 
    int idX ) ///< species postorder id
{
    MySpeciesNode *node = mSpeciesTree->getNodeById( idX );
    if( node->isLeaf() ) 
        return "'" + node->getName() + "'";
    else  
        return bpp::TextTools::toString( node->getInfos().realPostOrder );
}

/**
 * For the given id, search the sons to find the son id and return
 * the other one.
 *
 * @return The sibling of sonX.
 */
int DTLGraph::getSibling(
    int idX, ///< species id
    int sonX, ///< kept son species id
    bool realX ) ///< use real ids
{
    vector< pair<int,int> > splits = mSpeciesTree->getSplits( idX );
    pair<int,int> split;
    for( int i=0; i<2; i++ ) {
        BOOST_FOREACH( split, splits ) {
            if( realX ) {
                int realId1 = mSpeciesTree->getNodeById(split.first)->getId();
                int realId2 = mSpeciesTree->getNodeById(split.second)->getId();
                if( sonX == realId1 )
                    return realId2;
                if( sonX == realId2 )
                    return realId1;
            } else {
                if( sonX == split.first )
                    return split.second;
                if( sonX == split.second )
                    return split.first;
            }
        }

        if( mSpeciesTree->hasILS() ) 
            splits = mSpeciesTree->getIlsSplits( idX ); // try ils splits
        else
            break;
    }

    throw bpp::Exception("DTLGraph::getSibling: sibling not found");
}

/**
 * Create an event string.
 */
string DTLGraph::getEventString(
    int idU, ///< current clade
    int z, ///< index in reconciliation of current mapping vertex
    vector< vector<MyGraph::Vertex> > &reconciliation, 
        ///< list of vertices, by clade in the reconcilation
    int idUl, ///< current left clade child
    int idUr, ///< current right clade child
    map<string,double> &eventSupports )
        ///< supports to use if this is not empty
{
    // event indices follow their mapping vertex
    MyGraph::Vertex event = reconciliation[idU][z+1];

    // remove event number for output
    string eventOutName = "";
    string eventName = mGraph.properties(event).name;
    for( size_t i=0; i<eventName.length(); i++ ) {
        if( eventName[i] == '_' )
            break;
        eventOutName += eventName[i];
    }

    // output post order species id, rather than internal species id
    int idX = mGraph.properties(reconciliation[idU][z]).id_x; 
    string idXstr = getStringId( idX );

    double support;
    if( eventSupports.size() == 0 ) {
        support = mGraph.properties(event).support/getNumberSolutions();
    } else {
        map<string,double>::iterator iter = eventSupports.find( eventName );	
        if ( iter == eventSupports.end() ) 
            throw bpp::Exception( "DTLGraph::getEventString:"
                    " Could not find event in supports map" );
        support = iter->second;
    }

    if( support < 0.0001 )
        support = 0;
    string supportStr = bpp::TextTools::toString( support, 2 );

    string auxStr;
    if( eventName[1] == 'L' ) { 
        // put lost id in aux, which is the son not in the next event
        MyGraph::Vertex vNext = reconciliation[idU][z+2];
        int keptSonX = mGraph.properties(vNext).id_x; 
        if( eventOutName == "SL" || eventOutName == "IL" ) {
            int lostSonX = getSibling( idX, keptSonX, false );

            string lostId;
            if( mSpeciesTree->isAlpha( lostSonX ) )
                lostId = "-1";
            else
                lostId = getStringId( lostSonX );

            string keptId;
            if( mSpeciesTree->isAlpha( keptSonX ) )
                keptId = "-1";
            else
                keptId = getStringId( keptSonX );

            auxStr = lostId + "," + keptId;
        } else if( eventOutName == "TL" ) {
            auxStr = idXstr + "," + getStringId( keptSonX );
        } else if( eventOutName == "TLFD" ) {
            idXstr = "-1";
            auxStr = "-1," + getStringId( keptSonX );
        } else if( eventOutName == "TLTD" ) {
            auxStr = idXstr + ",-1";
        } else {
            cout << "Unknown event: " << eventOutName << endl;
            throw bpp::Exception("DTLGraph::getEventString: unknown event" );
        }
    } else if( eventOutName == "S" || eventOutName == "T" || eventOutName == "I"
               || eventOutName == "TFD"  || eventOutName == "TTD" ) 
    {
        // put species children in aux, ordered as ul, ur
//if( reconciliation[idUl].size() > 0 && reconciliation[idUr].size() > 0 ) {
        MyGraph::Vertex vl = reconciliation[idUl][0];
        int idXl = mGraph.properties(vl).id_x; 
        MyGraph::Vertex vr = reconciliation[idUr][0];
        int idXr = mGraph.properties(vr).id_x; 
        if( eventOutName == "S" || eventOutName == "T" || eventOutName == "I" ) 
        {
            auxStr = getStringId( idXl ) + "," + getStringId( idXr );
        } else {
            if( eventOutName == "TFD" )
                idXstr = "-1";
            // transfer to/from dead, replace alpha (dead) by -1
            if( mSpeciesTree->isAlpha( idXl ) )
                auxStr = "-1," + getStringId( idXr );
            else if( mSpeciesTree->isAlpha( idXr ) )
                auxStr = getStringId( idXl ) + ",-1";
            else 
                throw bpp::Exception("DTLGraph::getEventString: alpha not in "
                                "TTD or TFD" );
        }
//}
    } else if( eventOutName == "D" ) {
        auxStr = idXstr + "," + idXstr;
    } else if( eventOutName == "DD" ) {
        idXstr = "-1";
        auxStr = "-1,-1";
    } else {
        cout << "Unknown event: " << eventOutName << endl;
        throw bpp::Exception("DTLGraph::getEventString: unknown event" );
    }

    string eventString = idXstr + "," + eventOutName + "," + auxStr
                       + "@" + supportStr;
    return eventString;
}





/**
 * Call the given function on each reconciliation.
 */
void DTLGraph::traverseAllReconciliations( 
    ArgBase &args, ///< structure containing traversal specific variables
    bool (DTLGraph::*handleReconPtr)
        (int, vector< vector<MyGraph::Vertex> > &, ArgBase&) )
        ///< function to call on each reconciliation
{
    resetVisits();
    long reconCount=0;
    BOOST_FOREACH( MyGraph::Vertex root, mRootVertices ) {
        MyGraph::adjacency_vertex_range_t rootSons 
                    = mGraph.getAdjacentVertices( root );

        int sonCount = -1;
        BOOST_FOREACH( MyGraph::Vertex rootSon, rootSons ) {
            sonCount++;
            string sonName = mGraph.properties(rootSon).name;

            if( mOnlyCanonical && sonName[0] == 'O' ) 
                continue;

            DecisionNode rootDecision( root );
            // set visits to current event son
            mGraph.properties(root).visits = sonCount;

            // get reconciliations for this root son
            while( true ) {
                fillRecon( rootDecision );
                // convert to vertex list format and check validity
                vector< vector<MyGraph::Vertex> > recon( mMaxIdU+1 );
                bool valid = convertRecon( rootDecision, recon );
                if( valid ) {
                    reconCount++;
                    bool continueTraverse = 
                        (this->*handleReconPtr)( reconCount, recon, args );
                    if( !continueTraverse )
                        return;
                }

                // increment decision tree
                depthFirstOrderIncrement( rootDecision );
                if( rootDecision.mChanged )
                    break; // root incremented (handled by containing loop)
            }
        }
    }

    if( reconCount != mNumberSolutions ) 
        throw( "DTLGraph::traverseAllReconciliations:"
                "wrong number of reconciliations" );
}


/**
 * Convert a decision tree into a reconcilation structure.
 *
 * A reconciliation structure is list of vertices for each clade.
 * Each vertex list is  ordered by the traversal order of the graph,
 * which keeps species parent nodes before the children.
 *
 * @return True if reconcialition is valid.
 */ 
bool DTLGraph::convertRecon(
    DecisionNode &decision, ///< current decision node in recursion
    vector< vector<MyGraph::Vertex> > &recon ) ///< reconcilation structure
{
    int idU = -1;
    BOOST_FOREACH( MyGraph::Vertex v, decision.mReconVertices )
    {
        if( mGraph.properties(v).isMapping ) 
            idU = mGraph.properties(v).id_u;
        else if( mOnlyCanonical && !checkEventValidity( v ) ) 
            return false;
        recon[idU].push_back( v );
    }

    for( size_t i=0; i<decision.mChildren.size(); i++ ) 
        if( !convertRecon( decision.mChildren[i], recon ) )
            return false;

    return true;
}

/**
 * Check if the given event vertex is canonical for the current
 * reconciliation.
 *
 * If non-canoical nodes were removed, this might not do anything
 * unless the more complex criteria are used.
 *
 * The visits vertex property is used to choose mapping sons
 * if the event has two sons.
 *
 * @return True if the event is valid.
 */
bool DTLGraph::checkEventValidity(
        MyGraph::Vertex z ) ///< vertex to check
{
    if( mGraph.properties(z).isMapping ) 
        throw bpp::Exception( "DTLGraph::checkValidity called on mapping vertex" );

    int sonCount = mGraph.getVertexOutDegree(z);
    if( sonCount == 0 )   // leaf 
        return true;

    // mapping node sons (1 or 2 of them)
    MyGraph::adjacency_vertex_range_t sons = mGraph.getAdjacentVertices(z);
    MyGraph::Vertex son1 = *(sons.first);

    // get chosen mapping sons
    MyGraph::adjacency_vertex_range_t eventSons1 
                    = mGraph.getAdjacentVertices( son1 );
    int chosenSon1 = mGraph.properties(son1).visits;
    MyGraph::Vertex eventSon1 = *(eventSons1.first + chosenSon1);

    bool valid = true;
    if( sonCount == 1 ) {
        valid = validEvent( z, eventSon1 );
    } else {
        // get second event son
        MyGraph::Vertex son2 = *(sons.first+1);
        MyGraph::adjacency_vertex_range_t eventSons2 
                    = mGraph.getAdjacentVertices( son2 );
        int chosenSon2 = mGraph.properties(son2).visits;
        MyGraph::Vertex eventSon2 = *(eventSons2.first + chosenSon2);
        if( validRecNumber( z, eventSon1, eventSon2 ) == 0 )
            valid = false;
    }

    return valid;
}

/**
 * Traverse the decision tree and fill it with vertices
 * in the current reconciliation.
 */
void DTLGraph::fillRecon( 
        DecisionNode &decision ) ///< current decision node
{
    // only fill if changed
    if( decision.mChanged ) {
        decision.mChildren.clear();
        decision.mReconVertices.clear();

        int chosenSon = mGraph.properties(decision.mVertex).visits;

        // visit mapping nodes
        vector<MyGraph::Vertex> vertexStack;
        vertexStack.push_back( decision.mVertex );
        bool first = true;
        while( vertexStack.size() > 0 ) {
            MyGraph::Vertex v = vertexStack.back();
            vertexStack.pop_back();
            int sonCount = mGraph.getVertexOutDegree(v);
            if( first || sonCount == 1 ) {
                // add vertex to recon
                MyGraph::adjacency_vertex_range_t sons 
                        = mGraph.getAdjacentVertices( v );
                // only one son unless this is the decision
                MyGraph::Vertex event = *(sons.first + chosenSon);
                decision.mReconVertices.push_back( v );
                decision.mReconVertices.push_back( event );
                MyGraph::adjacency_vertex_range_t eventSons
                        = mGraph.getAdjacentVertices( event );
                BOOST_FOREACH( MyGraph::Vertex eventSon, eventSons ) 
                    vertexStack.push_back( eventSon );
                chosenSon = 0;
                first = false;
            } else {
                // new decision (mapping node with multiple sons)
                // Pushes a copy of childDecision unto the mChildren vector.
                // The original is deleted after this block.
                DecisionNode childDecision( v );
                decision.mChildren.push_back( childDecision );
            }
        }
    }

    // recursion
    for( size_t i=0; i<decision.mChildren.size(); i++ ) 
        fillRecon( decision.mChildren[i] );
}

/**
 * Print decision tree.
 *
 * For debugging.
 */
void DTLGraph::printDecisionTree( 
        DecisionNode &decision, ///< current node
        int depth ) ///< depth of recursion
{
    for( int i=0; i<depth; i++ ) 
        cout << "  ";
    cout << mGraph.properties(decision.mVertex).name 
         << " visits=" << mGraph.properties(decision.mVertex).visits
         << "/" << mGraph.getVertexOutDegree(decision.mVertex)
         << " changed=" << decision.mChanged
         << " children=" << decision.mChildren.size() 
         << endl;
    for( size_t i=0; i<decision.mChildren.size(); i++ ) 
        printDecisionTree( decision.mChildren[i], depth+1 );
}

/**
 * Visits nodes of a decision tree in depth first order, 
 * incrementing the visits property.
 *
 * The first node is visited. If it can be incremented, it is 
 * marked as changed and the recursion ends.
 *
 * If not, the node is reset, marked as changed, and the next node in DFS is
 * visted, repeating the procedure.
 *
 * All combinations have been visited when the first call returns
 * true (the root reset).
 *
 * @return True if this vertex was reset.
 */
bool DTLGraph::depthFirstOrderIncrement( 
        DecisionNode &decision )
{
    decision.mChanged = false;

    // pre order visit
    for( size_t i=0; i<decision.mChildren.size(); i++ ) {
        if( !depthFirstOrderIncrement( decision.mChildren[i] ) )
            // done because a descendant was incremented
            return false;
    }

    // all children changed, increment or reset this node
    decision.mChanged = true;

    if( mGraph.properties(decision.mVertex).visits+1 
            < mGraph.getVertexOutDegree(decision.mVertex) ) {
        // increment 
        mGraph.properties(decision.mVertex).visits++;
        return false;
    }

    // reset to first
    mGraph.properties(decision.mVertex).visits = 0;
    return true;
}




/**
 * Print all reconciliations to the given file.
 *
 * @return The number of reconciliations printed.
 */
long DTLGraph::printAllReconciliations(
    string path, ///< File name
    bool sylvxFormat, ///< output in the Sylvx format
    bool recPhyloXMLFormat, ///< use the recPhyloXMLFormat for the reconciliation
    bool checkConsistent, ///< only print consistent reconciliations
    int limit ) ///< maximum number of reconciliations to print
{
    ofstream outFile;
    outFile.open( path.c_str() );

    PrintReconsArgs args( &outFile, sylvxFormat, checkConsistent, limit );
    traverseAllReconciliations( args, &DTLGraph::printHandleRecon );

    outFile.close();
   
    return args.mPrintCount;
}

/**
 * Print all reconcilations function, used with reconciliation
 * traversal.
 *
 * @return True to continue printing.
 */ 
bool DTLGraph::printHandleRecon(
        int reconNumber, ///< ordinal of reconciliation in traversal
        vector< vector<MyGraph::Vertex> > &recon,  ///< a reconciliation
        ArgBase &baseArgs )   ///< class containing other arguments
{
    PrintReconsArgs *args = (PrintReconsArgs*) &baseArgs;

    vector<string> reconStrings;
    if( args->mSylvxFormat ) {
        reconStrings = getSylvxReconciliation( recon );
        reverse( reconStrings.begin(), reconStrings.end() );
    } else {
        map<string,double> eventSupports; // not used
        reconStrings = printReconciliationAux( recon, eventSupports );
    }
    bool printThisOne = true;
    if( args->mCheckConsistent && !isTimeConsistent( recon ) ) 
        printThisOne = false;

    if( printThisOne ) {
        if( !args->mSylvxFormat ) 
            *(args->mOutFile) << reconNumber
                << " =============================================" << endl;
        for( int i = reconStrings.size()-1; i>=0; i-- ) {
            if( reconStrings[i] != "" ) 
                *(args->mOutFile) << reconStrings[i] << endl;
        }
        args->mPrintCount++;
    }

    if( args->mLimit == 0 || reconNumber < args->mLimit )
        return true;
    return false;
}


/**
 * Check time consistencies of all reconciliations.
 *
 * @return 0 for inconsistent, 1 for consistent, 2 for unknown
 */
int DTLGraph::checkTimeConsistencies(
        int limit, ///< maximum number of reconciliations to check
        bool checkAll ) ///< check all, else until found consistent or limit
{
    TimeConsistentReconsArgs args( limit, checkAll );
    traverseAllReconciliations( args, &DTLGraph::timeConsistentHandleRecon );

    if( args.mLimitReached )
        return 2;

    long cycleCount = 0; 
    int recNum = 1;
    bool foundConsistent = 0;
    BOOST_FOREACH( bool result, args.results ) {
        if( !result ) {
            cycleCount++;
//            cout << recNum << " cycle" << endl;
        } else {
            foundConsistent = 1;
            if( !checkAll )
                return true; // found one consistent, return true
        }
        recNum++;
    }
   
    if( checkAll && cycleCount > 0 ) {
        long consistentCount = args.results.size() - cycleCount;
        cout << consistentCount << "/" << args.results.size()
             << " of the reconciliations are consistent" << endl;
    }

    return foundConsistent;
}

/**
 * All reconciliation traversal function that checks if the current
 * reconciliation is consistent (no cycles).
 *
 * A reconciliation is consistent if there no cycles in the species
 * tree when transfers are considered.
 *
 * @return True if not consistent or not searching all (true=continue).
 */
bool DTLGraph::timeConsistentHandleRecon(
        int reconNumber, ///< ordinal of reconciliation in traversal
        vector< vector<MyGraph::Vertex> > &recon,  ///< a reconciliation
        ArgBase &baseArgs )   ///< class containing other arguments
{
    TimeConsistentReconsArgs *args = (TimeConsistentReconsArgs*) &baseArgs;

    bool isConsistent = isTimeConsistent( recon );
    args->results.push_back( isConsistent ); // results are consistent

    if( isConsistent && !args->mCheckAll )
        return false; // stop traversal

    // check if the maximum number to check has been reached
    if( reconNumber == args->mLimit ) {
        args->mLimitReached = true;
        return false;
    }

    return true; // continue traversal
}


/**
 * Get all transfers and return parallel vectors with their
 * source and target vertices. Also return a vector with 
 * all sources (by transferFrom index) from the same gene.
 */
void DTLGraph::getAllTransfers(
    vector< vector<MyGraph::Vertex> > &recon,  ///< a reconciliation
    vector<MyGraph::Vertex> &transferFrom, ///< transfer sources
    vector<MyGraph::Vertex> &transferTo,  ///< transfer targets
    vector< vector<int> > &transfersForGene ) ///< index for each gene
{
    MyGraph::Vertex mappingV;
    BOOST_FOREACH( vector<MyGraph::Vertex> vect, recon ) {
        BOOST_FOREACH( MyGraph::Vertex v, vect ) {
            string name = mGraph.properties(v).name;
            if( mGraph.properties(v).isMapping ) 
                // save mapping vertex for pair (parent of event)
                mappingV = v; 
            else if( name[0] == 'T' ) {
                // ignore alpha transfers (TTD, TLFD, TLTD)
                if( name[1] == 'T' 
                    || (name[1] == 'L' 
                        && (name[2] == 'F' || name[2] == 'T' )))
                {
                    continue;
                }

                transferFrom.push_back( mappingV );
                MyGraph::adjacency_vertex_range_t sons 
                        = mGraph.getAdjacentVertices( v );
                MyGraph::Vertex son1 = *(sons.first);
                int toX = mGraph.properties(son1).id_x;
                if( toX == mGraph.properties(mappingV).id_x ) {
                    if( name[1] == 'L' ) 
                        throw bpp::Exception( "DTLGraph::getAllTransfers:"
                                " TL to same species" );
                    MyGraph::Vertex son2 = *(sons.first+1);
                    transferTo.push_back( son2 );

                } else {
                    transferTo.push_back( son1 );
                }
                // all transfers (by index) for idU
                size_t index = transferFrom.size() - 1;
                int idU = mGraph.properties(mappingV).id_u;
                transfersForGene[idU].push_back( index );
            }
        }
    }
}

/**
 * Recursive function to check for cycles in an id tree.
 *
 * The given tree is a vector of vector ids, one for each node in
 * the tree, where each node's vector contains a list of children
 * ids. The seenIt vector contains the list of parents of this node.
 * There is a cycle if a parent is visited a second time.
 *
 * @return True if there is a cycle.
 */
bool DTLGraph::hasCycles(
    int id, ///< current id in the recursion
    vector< vector<int> > &tree, ///< a tree of ids
    vector<bool> &seenIt ) ///< list of ids seen above this node
{
    if( seenIt[id] ) 
        return true;

    seenIt[id] = true; // mark this node as seen and descend
   
    BOOST_FOREACH( int childId, tree[id] ) {
        if( hasCycles( childId, tree, seenIt ) ) 
            return true;
    }

    // undo since it is not a cycle if another branch revisits this node
    seenIt[id] = false; 
   
    return false;
}


/**
 * Checks if the reconciliation is consistent (no cycles).
 *
 * A reconciliation is consistent if there no cycles in the species
 * tree when transfers are considered.
 *
 * @return True if consistent.
 */
bool DTLGraph::isTimeConsistent(
        vector< vector<MyGraph::Vertex> > &recon)  ///< a reconciliation
{
    // get all tranfers
    vector<MyGraph::Vertex> transferFrom;
    vector<MyGraph::Vertex> transferTo;
    int cladeCount = mCladesTrips->mClades.getCladeCount(); 
    vector< vector<int> > transfersForGene( cladeCount );
    getAllTransfers( recon, transferFrom, transferTo, transfersForGene );

    // get tree
    vector< vector<int> > cladeParents = mCladesTrips->getCladeParents();

    // create graph of species tree vertices with transfer or transfer
    // parents, using postorder ids (id_x) as indices
    vector< vector<int> > graph( mSpeciesTree->getNumberOfNodes() ); 

    // add edges
    vector<bool> inTree( graph.size(), false );
    for( size_t i=0; i<transferFrom.size(); i++ ) {
        int fromX = mGraph.properties(transferFrom[i]).id_x;
        int toX = mGraph.properties(transferTo[i]).id_x;
        MySpeciesNode *fromNode = mSpeciesTree->getNodeById( fromX ); 
        MySpeciesNode *toNode = mSpeciesTree->getNodeById( toX ); 

        // 1. add all edges between source/target and root
        MySpeciesNode *node = fromNode; // sources first
        int nodeId = fromX;
        for( int j=0; j<2; j++ ) {
            while( !inTree[nodeId] && node->hasFather() ) {
                inTree[nodeId] = true;
                MySpeciesNode *father = node->getFather();
                int fatherId = father->getId();
                graph[fatherId].push_back( nodeId );
                node = father;
                nodeId = fatherId;
            }
            // repeat with targets 
            node = toNode;
            nodeId = toX;
        }

        // 2. Cross link transfer with its parent nodes
        int parentFromX = getSpeciesParentId( transferFrom[i] );
        if( parentFromX != -1 )
            graph[parentFromX].push_back( toX );
        int parentToX = getSpeciesParentId( transferTo[i] );
        if( parentToX != -1 )
            graph[parentToX].push_back( fromX );

        // 3. condition 2 from Stolzer
        int idU = mGraph.properties(transferFrom[i]).id_u;
        while( cladeParents[idU].size() > 0 ) {
            idU = cladeParents[idU][0]; 
            // find each ancestor transfer
            BOOST_FOREACH( int parentIdx, transfersForGene[idU] ) {
                // cross link ancestor transfer parents to this transfer
                parentFromX = getSpeciesParentId( transferFrom[parentIdx] );
                if( parentFromX != -1 ) {
                    graph[parentFromX].push_back( toX );
                    graph[parentFromX].push_back( fromX );
                }
                parentToX = getSpeciesParentId( transferTo[parentIdx] );
                if( parentToX != -1 ) {
                    graph[parentToX].push_back( fromX );
                    graph[parentToX].push_back( toX );
                }
            }
        }
    }
        

    // check for cycles in species tree
    vector<bool> seenIt( graph.size(), false );
    int rootId = mSpeciesTree->getRootNode()->getId();
    bool foundCycle = hasCycles( rootId, graph, seenIt );
/*
// -put this in another function
// -pass this in (use another function if it isn't required since can't default)
vector<int> ordering( mSpeciesTree->getNumberOfNodes(), -1 ); 
    if( ordering.size() > 0 ) {
        vector<MySpeciesNode> leaves = mSpeciesTree->getSortedLeaves();
        vector<MySpeciesNode> curNodes;
        BOOST_FOREACH( MySpeciesNode *node, leaves ) {
            if( node->hasFather() )
                curNodes.push_back( node->getFather() );
        }



    }
*/


    return !foundCycle;
}






 
/////////////////////////////////////////////
///////////// Hali's Code ////////////////////
/////////////////////////////////////////////

/**
 * Return a reconciliation in Hali's format. DEPRECATED
 */
void DTLGraph::printReconciliationHali( 
        MyGeneTree &geneTree, vector<int> &pOrdToClade, 
        string outDir, string extensionId, bool triplets ) 
{



    vector<string> fileNames; 
    fileNames.push_back( "Rec_AllEventMapping_randomMPR.txt" );
    fileNames.push_back( "Rec_AllEventMapping_problem1.txt" );
    fileNames.push_back( "Rec_AllEventMapping_problem3.txt" );

    if( triplets )
        fileNames.push_back( "Rec_AllEventMapping_problem5.txt" );

    string header[4] = {
        "RandomRec", "Problem1Rec", "Problem3Rec", "Problem5Rec" };


    // first reconciliation (just any reconciliation will do)
    for( size_t i=0; i<fileNames.size(); i++ ) {
        ofstream outFileAll;
        outFileAll.open( string(outDir+fileNames[i]+extensionId).c_str() );
        string recContent= "******"
            + header[i]
            + "******\n";

cout << recContent;

        bool random = false;
        int problem = 1;
        if( i==0 ) {
            random = true;
        } else {
            if( i==2 )
                problem = 3;
            else if( i==3 )
                problem = 5;
            double score = getBestScore( problem );
            cout << "Problem " << problem << " score: " << score << endl;
        }
        vector< vector<MyGraph::Vertex> > reconciliation;
        if( getScoredReconciliation( problem, geneTree.getNumberOfNodes(), 
                                     reconciliation, random ) ) 
        {
            recContent += getReconciliationFormatHali( 
                                    geneTree, pOrdToClade, reconciliation );
            outFileAll << recContent;
            outFileAll.close();
        } // else score infinite, no reconciliation
    }
}

/**
 * Sort a list of strings and returned them as one comma separated string.
 *
 * @return sorted string
 */
string getSortedStrings( vector<string> & myStrings ) ///< a list of strings
{
    string strTemp = "";
    if( myStrings.size() < 1 )
        return strTemp;

    if (myStrings.size() >=2){
        // sort the elements of vector of string
        sort( myStrings.begin(), myStrings.end() );
   
        for( size_t i = 0; i <= myStrings.size() - 2; i++ ) 
            strTemp += myStrings[i] + "," ;
    }

    // the last element
    strTemp += myStrings[ myStrings.size() - 1 ];
    
    return strTemp;
}

/**
 * Format a support string.
 * @return formatted string.
 */
string supportString( double support ) ///< a support number
{
    if( support < 0.0001 )
        support = 0;
    return bpp::TextTools::toString( support, 2 );
}
    
/**
 * Return a reconciliation in Hali's format. DEPRECATED
 */
string DTLGraph::getReconciliationFormatHali(
        MyGeneTree &geneTree, 
        vector<int> &pOrdToClade, 
        vector< vector<MyGraph::Vertex> > &reconciliation ) 
{
    double numberSolutions = getNumberSolutions();
    int cladeCount = pOrdToClade.size();

	string recContent = "";
	for (int pOrd = cladeCount -1 ; pOrd >=0 ; pOrd-- ) {
        MyGeneNode *u = geneTree.getNodeById(pOrd);
	    vector <string> stru = bpp::TreeTemplateTools::getLeavesNames(*u);
        string struName = getSortedStrings( stru );
        int id_u = pOrdToClade[pOrd];
//cout << "pOrd=" << pOrd << " id_u=" << id_u << " name=" << struName << endl;


	    int size = reconciliation[id_u].size();
	    for (int z=0; z < size; z+=2 ) {
	     
            MyGraph::Vertex v = reconciliation[id_u][z];
            int id_x = mGraph.properties(v).id_x;
            vector< pair<int,int> > splits = mSpeciesTree->getSplits( id_x );
            if( splits.size() > 1 ) 
                throw bpp::Exception( "DTLGraph::getReconciliationFormatHali:"
                            " multiple species splits not support" );
//            vector<int> xSons = mSpeciesTree->getSplit( id_x );
            int timeSlice = mSpeciesTree->getTimeSlice( id_x );
            int realDonorId= mSpeciesTree->getRPO( id_x );

//cout << " " << pOrd << ", " << realDonorId << " z=" << z << "/" << size<< endl;
           
            //obtaining event name
            MyGraph::Vertex eventVertex = reconciliation[id_u][z+1];
            string supportStr = supportString( 
                    mGraph.properties(eventVertex).support /numberSolutions );

            string event="";
            string outString="";

            if( z==(size-2) ) {
                if(( splits.size() ==0) && (u->getNumberOfSons() ==0)){ 
                    MySpeciesNode* x = mSpeciesTree->getNodeById(id_x);
		            if(u->getName().substr(0,x->getName().size())==x->getName())
		                event="Extant";
		        } else if(u->getNumberOfSons() ==2) {
                    int pOrd_l = u->getSon(0)->getId();
                    int pOrd_r = u->getSon(1)->getId();
                    int id_u_l = pOrdToClade[pOrd_l];
                    int id_u_r = pOrdToClade[pOrd_r];
		            //int id_alpha1_ul= reconciliation[id_u_l][0];
                    MyGraph::Vertex vl = reconciliation[id_u_l][0];
                    int id_alpha1_ul = mGraph.properties(vl).id_x; 
                    //int id_alpha1_ur= reconciliation[id_u_r][0];
                    MyGraph::Vertex vr = reconciliation[id_u_r][0];
                    int id_alpha1_ur = mGraph.properties(vr).id_x; 
	
                    MyGeneNode* ul= u->getSon(0);
                    MyGeneNode* ur= u->getSon(1);
	
                    vector <string>strul =
                        bpp::TreeTemplateTools::getLeavesNames(*ul);
                    vector <string>strur = 
                        bpp::TreeTemplateTools::getLeavesNames(*ur);
                    string struNameL = getSortedStrings( strul );
                    string struNameR = getSortedStrings( strur );
                    //if(xSons.size() ==2){
                    if( splits.size() == 1 && splits[0].second != -1 ) {
                        // stored in the graph mapping
                        if((id_alpha1_ul==splits[0].first 
                                    && id_alpha1_ur == splits[0].second)
                                || (id_alpha1_ul==splits[0].second
                                    && id_alpha1_ur==splits[0].first))
                        {
                            event = "Spec";
                            outString = supportStr
                                //supportString( specEvents/numberSolutions ) 
                                + "@" + bpp::TextTools::toString(pOrd) 
                                + "@Spec@" 
                                + bpp::TextTools::toString(realDonorId) 
                                +"@NULL@" + struNameL + "@" + struNameR + "@\n";
		                }
		            }
		  
                    if(id_alpha1_ul==id_x  && id_alpha1_ur==id_x){
                        event = "Dup";
                        outString = 
                            //supportString( dupEvents/numberSolutions ) 
                            supportStr
                            + "@" + bpp::TextTools::toString(pOrd) + "@Dup@" 
                            + bpp::TextTools::toString(realDonorId) +"@NULL@" 
                            + struNameL + "@" + struNameR + "@\n";
		            }
                    if( (id_alpha1_ul==id_x  
                         && id_alpha1_ur!=id_x 
                         && mSpeciesTree->getTimeSlice(id_alpha1_ur)
                             == timeSlice ) 
                         || (id_alpha1_ur==id_x 
                             && id_alpha1_ul!=id_x 
                             && mSpeciesTree->getTimeSlice(id_alpha1_ul)
                                    == timeSlice))
                    {
                        int recipientId;
                        if (id_alpha1_ul == id_x)
                            recipientId = id_alpha1_ur;
                        else
                            recipientId = id_alpha1_ul;

                        int realRecipientId = mSpeciesTree->getRPO(recipientId);
                        event ="Tran";
                        outString =  supportStr
                                //supportString( transEvents/numberSolutions ) 
                                + "@" + bpp::TextTools::toString(pOrd) +"@Tran@"
                                + bpp::TextTools::toString(realDonorId) +"@" 
                                + bpp::TextTools::toString(realRecipientId) 
                                + "@" 
                                + struNameL + "@" + struNameR + "@\n";
                    }								
                }	    
            } else {
                //int id_alphaiplus1_u= reconciliation[id_u][z+1];	
                MyGraph::Vertex vNext = reconciliation[id_u][z+2];
                int id_alphaiplus1_u = mGraph.properties(vNext).id_x; 
                if(splits.size() == 1 && splits[0].second == -1) {
                    // stored in the graph mapping
                    if(id_alphaiplus1_u == splits[0].first)
                        event="NoEvent";								
		        } else if(splits.size() ==1 && splits[0].second != -1) {
                    // stored in the graph mapping
                    if (id_alphaiplus1_u == splits[0].first 
                            || id_alphaiplus1_u == splits[0].second) 
                    {
                        event="SpecLoss"; // SpecLoss
		                outString = supportStr
                            //supportString( specEvents/numberSolutions ) 
                            + "@" + bpp::TextTools::toString(pOrd) + "@Spec@" 
                                  + bpp::TextTools::toString(realDonorId) 
                                  +"@NULL@"
                                  + struName + "@@\n" ;
		    
                        int donorLossId= (id_alphaiplus1_u == splits[0].first) 
                                           ? splits[0].second :splits[0].first ;
                        int realDonorLossId = mSpeciesTree->getRPO(donorLossId);
                        outString += supportStr
                            //supportString( lossEvents/numberSolutions ) 
                            + "@" + bpp::TextTools::toString(pOrd) + "@Loss@" 
                            + bpp::TextTools::toString(realDonorLossId) 
                            + "@NULL@@@\n";
                    }
		        }
		 
                if( mSpeciesTree->getTimeSlice(id_alphaiplus1_u) == timeSlice) {
		            event="TranLoss";
                    int realRecipientId = mSpeciesTree->getRPO(
                                            id_alphaiplus1_u );
                    outString = supportStr
                        //supportString( transEvents/numberSolutions ) 
                        + "@" + bpp::TextTools::toString(pOrd) +"@Tran@" 
                                + bpp::TextTools::toString(realDonorId) +"@"  
                                + bpp::TextTools::toString(realRecipientId) 
                                + "@" 
                                +  struName +"@@\n";
                    // Loss at id_x
                    outString += supportStr
                        //supportString( lossEvents/numberSolutions ) 
                            + "@" + bpp::TextTools::toString(pOrd) + "@Loss@" 
                              + bpp::TextTools::toString(realDonorId) 
                              +"@NULL@@@\n";
		        }
            }
            recContent += outString;
        }
    }
    return recContent;
}
 


