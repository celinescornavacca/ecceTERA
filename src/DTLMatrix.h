#ifndef DTLMATRIX_H_ 
#define DTLMATRIX_H_ 

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

/**
Functions to compute amalgamation matrix wrapped in a class in order
to have state variables accesible from sub-functions.  
The constructors are used to set the input and parameters before
calling the function calculateMatrix (or recalculateMatrix).

To calculate the matrix, clades are visited from leaves to root.
For each time slice, costs for each clade species node combination
are calculated, initially without transferLoss, which is calculated
after the bestReciever (lowest cost species node) is found.

The printCSVmatrix outputs the matrix for subsequent recalculations.
The matrix and bestReceiver data can be written to a file and
read for a fast recalculation for changes in the species tree.

The backtrack function returns the consensus gene tree. To calculate
it, all of the bestSplits must be known, which isn't the case
for recalculated matrices.
**/

#include "MyMatrix.h"
#include "CladesAndTripartitions.h"
#include "DTLGraph.h"


class DTLMatrix {
protected:
    struct DTLMatrixState {
        // state variables for matrix computation
        int timeSlice;            ///< current time slice
        int idU;                  ///< current clade id
        int idX;                  ///< current species id
        double curBestEpsCost;    ///< subopt best cost
        vector<double> costList;  ///< subopt costs
        double dupliCost;       ///< duplication cost
        double hgtCost;         ///< transfer cost
        double lossCost;        ///< loss cost
    };

    // split information
    struct BestSplit {
        char event;                 ///< event for best split
        pair<int,int> cladeSplit;   ///< clade split for event
        pair<int,int> speciesSplit; ///< species split for event
    };

    // compute costs except for TL
    virtual void computeOptimaForCladeSplit( int toCompute, int splitIdx,
            pair<int,int> cladeSplit, DTLMatrixState &state, 
            double &optCost, BestSplit &bestSplit );
    virtual void computeOptimaForSpeciesSplit( int toCompute,
            pair<int,int> cladeSplit, DTLMatrixState &state, 
            double &optCost, BestSplit &bestSplit );
    

    virtual void pushList( DTLMatrixState &state ) {}
    
            
    void init(); // helper for constructors
	virtual void initTS() { }

    virtual void setAlpha( int alphaIdx ) {}
  
                      
    void computeBestReceiversNoSub( int idU, bool recompute,
            vector<int> &speciesNodeIds, vector<double> &optVector );
        


    static const double COST_DIFF;     ///<Minimal Difference between costs
            ///< at which they are considered different.
    bool COST_EQUAL( double a, double b );
    bool COST_GREATER( double a, double b );

    
    bool mUseILS;

    bool mUseBestSplits; ///< save split information
    BestSplit **mBestSplits; ///< backtrack information
      
      
    // input parameters (or derived from input)
    int mCladeCount; ///< number of clades
	int mSTnodes;    ///< number of species nodes for a time slice
    int mMaxTS;      ///< maximum time slice
    MySpeciesTree *mSpeciesTree; ///< species tree
    CladesAndTripartitions *mCladesTrips; ///< clades and tripartitions

    bool mRecal; ///< true if this is a recalculation

    bool mMatrixCalculated; ///< true after calculateMatrix called



    // input parameters (or derived from input)
    bool mFixedCosts;    ///< used fixed costs instead of node specific
    bool mComputeT;      ///< compute transfer cost
    bool mComputeTL;     ///< compute transfer loss cost
    double mDupliCost;   ///< duplication cost
    double mHGTCost;     ///< transfer cost
    double mLossCost;    ///< loss cost
    double mSplitWeight; ///< split weighting
    double mWGDCost;     ///< whole genome duplication costs


    MyMatrix mMatrix; ///< the matrix


    // computed variables  
	int **mBestReceiver;  ///< all best receivers
    double **mBestReceiverCost;  ///< all best receiver costs
	int **mSecondBestReceiver;   ///< all second best receivers
    double **mSecondBestReceiverCost;  ///< all second best receiver costs
    vector <int> mAllOtherBestReceivers; 

    vector< vector<bool> > mComparable; ///< comparable species nodes


    vector< vector<int> > mSubTreeIdLists;
        ///< lists of ids in each subtree


    void removeDuplicates( vector<int> &v );
    double computeCostList( double otherCost, int idUa, int idXa, 
            int idUb, int idXb, DTLMatrixState &state );
    double computeCostList( double otherCost, int idUa, int idXa,
            DTLMatrixState &state );

    void computeTransferCost( int idUl, int idUr, double costThisSplit, 
            DTLMatrixState &state, double &optCost, BestSplit &bestSplit );
    virtual void computeTransferCostSub( int optimumSub, int idUsub,    
        int idUother, double costThisSplit, DTLMatrixState &state ) {}
    void computeDuplicationCost( int idUl, int idUr, double costThisSplit, 
            DTLMatrixState &state, double &optCost, BestSplit &bestSplit );
    void computeSpeciationCost( int idUl, int idUr, int idXl, int idXr,
            double costThisSplit, 
            DTLMatrixState &state, double &optCost, BestSplit &bestSplit );
    void computeIlsCost( int idUl, int idUr, int idXl, int idXr,
            double costThisSplit, 
            DTLMatrixState &state, double &optCost, BestSplit &bestSplit );
    void computeSpeciationPlusLossCost( int idXl, int idXr, 
            DTLMatrixState &state, double &optCost, BestSplit &bestSplit );
    void computeIlsPlusLossCost( int idXl, int idXr, 
            DTLMatrixState &state, double &optCost, BestSplit &bestSplit );
    void computeNullCost( int idXl, DTLMatrixState &state, double &optCost, 
            BestSplit &bestSplit );
    double computeOptimaForCell( int speciesNodeId, 
            int events, DTLMatrixState &state );
    bool computeTransferLossCost( double alphaCost, double &opt, 
            DTLMatrixState &state );

    virtual void addZeroCost() {}

    int computeBestReceivers( int timeSlice, int idU, 
            vector<int> &speciesNodeIdsTS, vector<double> &optVectorTS );

    void updateOpt( double &optCost, BestSplit &bestSplit, double cost, 
            int firstSplit, int secondSplit, char event );

    virtual void addVertices( DTLGraph &graph, DTLGraph::MyGraph::Vertex z, 
            vector<DTLGraph::MyGraph::Vertex> &qList, int u, int x,
            const char*event, double curCost, double otherCost );
    virtual void addVertices( DTLGraph &graph, DTLGraph::MyGraph::Vertex z, 
            vector<DTLGraph::MyGraph::Vertex> &qList, 
            int u1, int x1, int u2, int x2, const char*event, 
            double curCost, double otherCost);

private:
    // difference at which costs are considered different;
    // used for estimating new costs
    static const int PRINT_MATRIX_PRECISION = 15; ///< precision for output

                          


    // cost functions
    virtual void computeList( int idUa, int idXa, int idUb, int idXb,
            double lowestCost, double curBestEpsCost, 
            double otherCost ) {}
    virtual void computeList( int idUa, int idXa, 
            double lowestCost, double curBestEpsCost, 
            double otherCost) {}





    // compute best receivers and transfer loss
    void initComparable( bool partialDates );

    int saveBestReceivers( int timeSlice, int idU, int bestReceiver, 
            int secondBestReceiver, double bestCost, double secondBestCost, 
            vector<int> &otherBests );

    double findIdsWithCurrentCost( vector<int> &speciesNodeIdsTS, 
        vector< vector<double> > &allCostsTS,
        double lastValue, vector<size_t> &indices, vector<int> &idList );

    // create graph functions
    vector<int> ***getAllBestReceivers();
    virtual vector<DTLGraph::MyGraph::Vertex> getRootNodes( DTLGraph &graph,
                                bool keepRoot );
    virtual void createVertices( DTLGraph &graph, 
                DTLGraph::MyGraph::Vertex pairVertex, 
                vector<DTLGraph::MyGraph::Vertex> &qList, 
                vector<int> *** allBestReceivers );

    void clearBestReceivers(); // sets all best receiver info to -1


    void calculateMatrixTS( int idU, int timeSlice );
    virtual void calculateMatrixNoSub( int idU );

    // various functions
    bool updateCosts( bool weightOnly, bool verbose );
    void checkSolution( int splitNum, int id_u, int z, 
            vector< vector<int> > &recAlphaVector,
            double &costRec, bool &canonical, bool &isRec );
    double countEventsAux( int idU, int idX, int &duplications, int &transfers, 
            int &losses, int &ils, double &cost );
    void printAllCostsError( double optVal, vector<double> &allCosts,
                             DTLMatrixState &state ); 

    // ILS
    double mIlsCost;    ///< cost of an ILS event
public:
    DTLMatrix( MySpeciesTree *speciesTree, CladesAndTripartitions *cat, 
               double WGDCost, bool fixedCosts, bool computeT, bool computeTL,
               double dupliCost, double hgtCost, double lossCost, 
               int maxTS, double weight, bool mUseBestSplits,
               double ilsCost );

    DTLMatrix() {
        mUseBestSplits = false;
    }
    virtual ~DTLMatrix();

    // the main functions
    void calculateMatrix( bool verbose=false, int maxIterations=1, 
            bool updateWeightOnly=false, int dated=2,
            int tsStart=0, int tsEnd=-1 );

    // lowest cost associated with gene root (all clade)
    double getBestCost( bool keepRoot=false );
    double getBestCost( int &u, int &x );

    void printMatrixCSV( const char *fileName,
            const char *speciesTreeFile, const char *genesTreeFile );

    /**
     * Change the species tree.
     *
     * for strale
     */
    void updateSpeciesTree( 
            MySpeciesTree *speciesTree )  ///< a species tree
    {
        mSpeciesTree = speciesTree;
    }

    // find best splits and events (consensus gene tree)
    MyGeneNode *backtrack( bool printTree=false, bool keepRoot=false,
                           int idU=-1, int idX=-1, int level=0 );
    double countEvents( int &duplications, int &transfers, int &losses,
                        int &ils );


    // graph
    DTLGraph constructGraph( bool verbose=false, bool keepRoot=false );

};


#endif 
