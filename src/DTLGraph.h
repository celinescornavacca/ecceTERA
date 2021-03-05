#ifndef DTLGRAPH_H_
#define DTLGRAPH_H_

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
 

#include <string>
#include <vector>


#include "Graph.h"
#include "MySpeciesTree.h"
#include "CladesAndTripartitions.h"

using namespace std;

/**
A DTLGraph represents all possible solutions to a most parsimonious
reconciliation. There are two types of vertices: pairs and events.
A pair vertex contains a clade (subtree of genes) and 
a species node which are reconciled.
Additionally, a pair (mapping) vertex contains a cost for suboptimal 
reconciliations, which can contain the same clade/species node 
reconciliations are different costs. Event vertices represent the 
possible DTL (duplications, transfer, loss) 
events that can occur from a vertex pair. A traversal following all
vertex pairs and a single event vertex from each node will result in
a single solution. The countReconciliationNumberAndCheck function calculates
the number of solutions and finds the number of solutions each event
vertex is in. The getBestScore calculates a graph score.

A DTLGraph is a forest (mRootVertices). Each root represents a differnt 
mapping from the species tree to the root of the gene tree.


The DTLGraph has three types of public functions:

1. Graph creation: addPairVertex, addEventVertex, addRoot, getPairVertexByName.

2. Property calculation: countReconciilationNumberAndCheck and 
   getBestScore.

3. Property retreival: getScoredReconciliation,
   getVertexIdentfiers, and getEventCounts.

For debugging, there is also the printGraph function.

Finally, the getAlphas returns a represenation of all solutions,
but is no longer supported and should be redone if it is to be used.


DTLGraph is modeled on the BOOST graph library, and therefore uses
generic vertices with a vertex property structure attached to each
node. The BOOST library functions are restricted to the Graph class,
which should help if the Graph class needs to be rewritten without BOOST.


As many functions as possible have been written as graph traversal
functions to simplify the code. The functions depthFirstTraversal
and breadthFirstTraversal take two functions pointers as input, which
are run pre-order and post-order at each vertex.
*/
class DTLGraph
{
public:

    // Define MyGraph and vertex properties
	struct VertexProperties {
        string name;    ///< vertex name 
        bool isMapping; ///< an event node if false
       
        // mapping vertex info
	    int id_u;       ///< clade id for node
        int id_x;       ///< species id for node
        double cost;    ///< node cost used to diferentiate suboptimal vertices

		int visits;     ///< used to track graph traversals

        // computed values
        double recNumber;       ///< number of reconciliations in sub-graph,
            ///< also used by vertex removal function (removeNoncanonical)
		double support;  ///< number of reconciltions vertex is in

        //! weighted support for vertex
        /*! number of reconciliations vertex is in for each root 
         (for weighted support) */
        vector<double> supportVector;
        double score;           ///< various median scores
        double r;               ///< used by computeR

        int count; ///< whatever count for any function 
	};

	struct EdgeProperties {
		string propertyNotUsed; ///< not used right now
	};

	typedef Graph<VertexProperties,EdgeProperties> MyGraph; 
        ///< graph type
       
private:
    static const double SCORE_DIFF;     ///<Minimal difference between scores 
            ///< at which they are considered different.
    static const int NAME_PRECISION; ///< Precision of costs in vertex names.
    bool SCORE_EQUAL( double a, double b );
    bool SCORE_GREATER( double a, double b );

    bool mOnlyCanonical;  ///< remove non-canonical vertices if true 
     int mEventNumber;     ///< counter for numbering events 
    int mScoredProblem;   ///< which problem was scored in vertices
    int mMaxIdU;          ///< maximum idU value
    double mProblemScore; ///< median score for current problem (mScoredProblem)
    double mNumberSolutions; ///< number of solutions (traversals of graph)
    MySpeciesTree *mSpeciesTree;  ///< species tree used to calculate matrix
    CladesAndTripartitions *mCladesTrips; 
        ///< clades and tripartitions used in matrix calculation
    MyGraph mGraph;       ///< the graph
    vector<MyGraph::Vertex> mRootVertices; ///< graph is forest
    map<string,MyGraph::Vertex> mVertexMap; ///< vertex name -> vertex



    // The base class containing no variables. If a traversal function
    // needs variables, a child class is created with the variables,
    // and the ArgBase class is cast to the child in the function.
    class ArgBase {};  


    // graph creation
    MyGraph::Vertex createPairVertex( int id_u, int id_x, double cost,
                                      int d = -1, int t = -1, int l = -1 );

    ///////////////////////////////////////////
    //////  Graph Traversal Functions /////////
    ///////////////////////////////////////////


    // breadth and depth first traversal
    void depthFirstTraversalAux( MyGraph::Vertex z, ArgBase &args,
        void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&),
        void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&));
    void depthFirstTraversal( ArgBase &args,
        void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&),
        void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&));
    void breadthFirstTraversalAux( MyGraph::Vertex z, bool isRoot, 
            ArgBase &args,
        void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&),
        void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&));
    void breadthFirstTraversal( ArgBase &args,
        void (DTLGraph::*discoverVertexPtr)(MyGraph::Vertex, ArgBase&),
        void (DTLGraph::*finishVertexPtr)(MyGraph::Vertex, ArgBase&));



    // reconciliation counting
    double countSubReconciliations();
    void countSubReconciliationsFinishVertex( 
            MyGraph::Vertex z, ArgBase &baseArgs );

    void computeSupport();
    void computeSupportDiscoverVertex( MyGraph::Vertex z, ArgBase &baseArgs );
    void computeR( MyGraph::Vertex z );
    map<string,vector<MyGraph::Vertex> > createEventDescriptorMap(); 
    void weightSupport( double epsilon );

    /** Arguments for remove non-canonical vertices traversal function. */
    class RemoveArgs : public ArgBase {
        public:
        MyGraph *mNewGraph; ///< the new canonical graph
        /** Constructor */
        RemoveArgs( MyGraph *newGraph ) ///< the new graph
            : mNewGraph(newGraph) {};
    };
    void removeNoncanonicalFinishVertex( 
            MyGraph::Vertex z, ArgBase &baseArgs );
    void removeNoncanonicalVertices();


    /** Arguments for median scoring traversal function. */
    class BestScoreArgs : public ArgBase {
        public:
        double mScoreMod;    ///< amount to add to score
        /** Constructor */
        BestScoreArgs( double scoreMod ) : ///< amount to add to score
            mScoreMod(scoreMod) {};
    };
    void bestScoreFinishVertex( MyGraph::Vertex z, ArgBase &baseArgs );
    
    
    /** Arguments for printing traversal function. */
    class PrintArgs : public ArgBase {
        public:
        bool mCleanGraph; ///< the graph
        bool mUseInternalIds; ///< use internal ids rather than matching
                              ///< Those in the reconciliation
        vector<int> mCladeToPOrd;    ///< mapping for external clade ids
        std::ostream *mOut; ///< output stream
        /** Constructor */
        PrintArgs( bool cleanGraph,  ///< the graph
                bool useInternalIds, ///< use internal ids
                vector<int> cladeToPOrd, ///< mapping for external clade ids
                std::ostream *out) ///< output stream
            : mCleanGraph(cleanGraph), mUseInternalIds(useInternalIds),
                mCladeToPOrd(cladeToPOrd), mOut(out) {}
    };
    void printDiscoverVertex( MyGraph::Vertex z, ArgBase &args );

    /** Arguments for countTLs traversal function. */
    class CountTLsArgs : public ArgBase {
        public:
        double tlSupport;
        double eventSupport;
    };
    void countTLsDiscoverVertex( MyGraph::Vertex z, ArgBase &args );



    ////////////////////////////////////////////////////
    //////  Reconciliation Traversal Functions /////////
    ////////////////////////////////////////////////////

    void traverseAllReconciliations( ArgBase &args, 
            bool (DTLGraph::*handleReconPtr) 
                (int, vector< vector<MyGraph::Vertex> >&, ArgBase&) );

    // Print all reconciliations
    /** Arguments for printing traversal function. */
    class PrintReconsArgs : public ArgBase {
        public:
        std::ofstream *mOutFile; ///< file output stream
        bool mSylvxFormat;       ///< print in the Sylvx format
        bool mCheckConsistent; ///< only print time consistent recons
        long mLimit; ///< maximum number of reconciliations to print
        long mPrintCount; ///< number printed

        /** Constructor */
        PrintReconsArgs( 
            std::ofstream *outFile, ///< file output stream
            bool sylvxFormat,       ///< print in the Sylvx format
            bool checkConsistent, ///< only print consistent recons
            int limit ) ///< maximum number of reconciliation to print
            : mOutFile(outFile), mSylvxFormat(sylvxFormat),
              mCheckConsistent(checkConsistent), mLimit(limit) 
        {
            mPrintCount = 0;
        }
    };
    bool printHandleRecon( int reconNumber,
            vector< vector<MyGraph::Vertex> > &recon, ArgBase &baseArgs );

    /** Arguments for time consistency traversal function. */
    class TimeConsistentReconsArgs : public ArgBase {
        public:
        int mLimit;
        bool mCheckAll;
        bool mLimitReached;
        vector<bool> results; ///< result for each reconcisliation
        /** Constructor */
        TimeConsistentReconsArgs( int limit, bool checkAll )
            : mLimit(limit), mCheckAll(checkAll) {
                mLimitReached = false;
            }
    };
    bool timeConsistentHandleRecon( int reconNumber,
            vector< vector<MyGraph::Vertex> > &recon, ArgBase &baseArgs );
    void getAllTransfers( vector< vector<MyGraph::Vertex> > &recon, 
                          vector<MyGraph::Vertex> &transferFrom,
                          vector<MyGraph::Vertex> &transferTo,
                          vector< vector<int> > &transfersForGene );
    bool isTimeConsistent( vector< vector<MyGraph::Vertex> > &recon);


    ///////////////////////////////////////////
    //////  Other Functions /////////
    ///////////////////////////////////////////

    // non-canonical validation function
    int getSpeciesParentId( MyGraph::Vertex z );
    int getFatherXId( MyGraph::Vertex event );
    bool validEvent( MyGraph::Vertex event, MyGraph::Vertex eventSon );
    double validRecNumber( MyGraph::Vertex event, MyGraph::Vertex eventSon1,
                     MyGraph::Vertex eventSon2 );
//    double condition7( MyGraph::Vertex eventSon1, MyGraph::Vertex eventSon2 );


    /** various utility functions */
    void resetVisits() {
        MyGraph::vertex_iter iter, v_end;
        for(boost::tie(iter,v_end) = mGraph.getVertices(); 
            iter != v_end; iter++) 
	{
            mGraph.properties(*iter).visits = 0;
	}
    }
    string mappingFromIds(int id_u, int id_x, double cost );
    string mappingFromIds(int id_u, int id_x, double cost, int d, int t, int l);
    double scoreReconciliation(
            vector< vector<MyGraph::Vertex> > &reconciliation );

    

    // find reconciliation associated with a scored problem
    size_t chooseRandomSupport( vector<double> &supports );
    void backtrack( MyGraph::Vertex vertex,
                    vector< vector<MyGraph::Vertex> > &reconciliation, 
                    bool random, double scoreMod );
    void backtrackEvents( MyGraph::Vertex vertex,
                    MyGraph::Vertex mappingVertex,
                    MyGraph::adjacency_vertex_range_t eventSons, 
                    vector< vector<MyGraph::Vertex> > &reconciliation, 
                    bool random, double scoreMod );
    void backtrackEvents( MyGraph::Vertex vertex,
                    MyGraph::Vertex mappingVertex1,
                    MyGraph::Vertex mappingVertex2,
                    MyGraph::adjacency_vertex_range_t eventSons1, 
                    MyGraph::adjacency_vertex_range_t eventSons2, 
                    vector< vector<MyGraph::Vertex> > &reconciliation, 
                    bool random, double scoreMod );



    // output
    vector<string> printReconciliationAux( 
            vector< vector<MyGraph::Vertex> > &reconciliation, 
            map<string,double> &eventSupports );
   vector<int> orthologyOutputAux( int idU,
                vector< vector<MyGraph::Vertex> > &reconciliation,
                vector< pair<int,int> > &orthologs );
    double getSpeciesDates( MySpeciesNode *node, vector<double> &startDates,
                            vector<double> &endDates );
    int getChildIdx( int idU,
        vector< vector<MyGraph::Vertex> > &reconciliation,
        size_t reconciliationIdx );

    string makeIntervals( int idU,
        vector< vector<MyGraph::Vertex> > &reconciliation,
        vector<double> &speciesStartDates, vector<double> &speciesEndDates,
        map<string,double> &eventEndDates,
        size_t reconciliationIdx = 0,
        double parentStartDate = -1, string fatherEvent = "ROOT",
        int prevRpo = -1,
        int seqNum = 0 );

    void printEvents( int idU, vector<int> &cladeToPOrd,
        vector< vector<MyGraph::Vertex> > &reconciliation,
        vector<double> &speciesEndDates,
        map<string,double> &eventDates,
        vector< vector<string> > &eventsByIdU,
        vector<string> &losses,
        size_t reconciliationIdx = 0,
        int fatherX = -1,
        string fatherEvent = "ROOT",
        string prevRealX = "ROOT", 
        string prevRealFatherX = "ROOT" );

        vector<string> getRecPhyloXMLReconciliation(
            vector< vector<MyGraph::Vertex> > &reconciliation );

        string getRecPhyloXMLReconciliation(
            vector< vector<MyGraph::Vertex> > &reconciliation, int idUP, int idU,  vector< int>  &cladeToPOrd );   

    string getStringId( int idX );

    int getSibling( int idX, int sonX, bool realX );

    string getEventString( int idU, int z, 
            vector< vector<MyGraph::Vertex> > &reconciliation, 
            int idUl, int idUr, map<string,double> &eventSupports );
            
    // used by traverseAll
    class DecisionNode {
        public:
        bool mChanged;
        MyGraph::Vertex mVertex;
        vector<DecisionNode> mChildren;
        vector<MyGraph::Vertex> mReconVertices;
        /** Constructor */
        DecisionNode( MyGraph::Vertex v ) 
            ///< mapping vertex with multiple events
            : mChanged(true), mVertex(v) {}
    };
    bool depthFirstOrderIncrement( DecisionNode &decision);
    void combineEventLists( MyGraph::Vertex z, MyGraph::Vertex son1,
                            MyGraph::Vertex son2 );
    void allReconsFinishVertex( MyGraph::Vertex z, ArgBase &args );

    void printDecisionTree( DecisionNode &decision, int depth=0 );
    void fillRecon( DecisionNode &decision );
    bool convertRecon( DecisionNode &decision,
        vector< vector<MyGraph::Vertex> > &recon );
    bool checkEventValidity( MyGraph::Vertex z );
    string createNameWithExternalIds( 
            MyGraph::Vertex z, vector<int> &cladeToPOrd );
    void getEvents( MyGraph::Vertex z, string &d, string &t, string &l );


    bool hasCycles( int id, vector< vector<int> > &tree, vector<bool> &seenIt );

public: 
    /** DTL Graph Constructor */
    DTLGraph( 
            MySpeciesTree *speciesTree, 
                    ///< Species tree used to construct matrix.
            CladesAndTripartitions *cladesTrips ) 
                ///< Clades/tripartions for the matrix.
        : mOnlyCanonical(false),
        mEventNumber(0), mScoredProblem(-1), mMaxIdU(-1),
        mSpeciesTree(speciesTree), mCladesTrips(cladesTrips) {}


    // graph building functions
    bool addPairVertex( MyGraph::Vertex eventVertex, int id_u, 
                int id_x, double cost, MyGraph::Vertex &pairVertex );
    bool addPairVertex( MyGraph::Vertex eventVertex, int id_u, 
                int id_x, double cost, int d, int t, int l, 
                MyGraph::Vertex &pairVertex );
    MyGraph::Vertex addEventVertex( MyGraph::Vertex z, string event );
    MyGraph::Vertex addRoot( int id_u, int id_x, double cost,
                                      int d = -1, int t = -1, int l = -1 );



    // compute graph properties
    bool countReconciliationNumberAndCheck(
                    bool onlyCanonical, bool verbose);
    void pruneNonoptimal();


    // return graph properties
    double getBestScore( int problem );
    bool getScoredReconciliation( int problem, size_t cladeCount, 
            vector< vector<MyGraph::Vertex> > &reconciliation,
            bool random=false );
    double countTLs( double &frequency );
    int checkTimeConsistencies( int limit, bool checkAll=false );



    // graph output
    void printGraph( string path, bool useInternalIds, bool cleanGraph=false);
    void orthologyOutput( string fileName );
    void printReconciliation( string problemStr, string fileName,
                bool sylvxFormat,
                bool recPhyloXMLFormat, 
                bool checkConsistent, bool &isConsistent,
                map<string,double> &eventSupports );
    long printAllReconciliations( string path, bool sylvxFormat, bool recPhyloXMLFormat, 
                                  bool checkConsistent, int limit=0 );



    // other functions
    void getEventSupports( map<string,double> &eventSupports );
    void getVertexIdentfiers( MyGraph::Vertex vertex, 
                              int &id_u, int &id_x, double &cost,
                              int &d, int &t, int &l ); 


    /**
     * Return the number of solutions for graph.
     */
    double getNumberSolutions() { 
        return mNumberSolutions; 
    } 

    /**
     * Return ids and cost for vertex
     */
    void getVertexIdentfiers( MyGraph::Vertex vertex, 
                              int &id_u, int &id_x, double &cost ) 
    {
        id_u = mGraph.properties(vertex).id_u;
        id_x = mGraph.properties(vertex).id_x;
        cost = mGraph.properties(vertex).cost;
    }

    /**
     * Return vertex name.
     *
     * @vertex name
     */
    string getVertexName( MyGraph::Vertex vertex ) 
    {
        return mGraph.properties(vertex).name;
    }

};


#endif
