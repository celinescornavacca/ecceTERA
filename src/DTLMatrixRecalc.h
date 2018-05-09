#ifndef DTLMATRIX_RECALC_H_ 
#define DTLMATRIX_RECALC_H_ 

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

**/

#include "DTLMatrix.h"


class DTLMatrixRecalc : public DTLMatrix {
private:

    // matrix reading state
    static const int READING_MATRIX = 0;    ///< reading the matrix
    static const int READING_BEST_REC = 1;  ///< reading best receivers
    static const int READING_BEST_REC_COST = 3; ///< reading br costs

    // recalculation lists - maps pair<idU,idX> to an event list
    // Stores what needs to be recalculated.
    boost::unordered_map<pair<int,int>,int> mRecalQ; 
        ///< current time slice events to be recalculated
    boost::unordered_map<pair<int,int>,int> mRecalNextQ; 
        ///< next time slice events to be recalculated

    virtual void computeOptimaForCladeSplit( int toCompute, int splitIdx,
            pair<int,int> cladeSplit, DTLMatrixState &state, 
            double &optCost, BestSplit &bestSplit );
    virtual void computeOptimaForSpeciesSplit( int toCompute,
            pair<int,int> cladeSplit, DTLMatrixState &state, 
            double &optCost, BestSplit &bestSplit );

    bool recomputeOptima( bool &changed, bool recomputeAll,
            double &cost, DTLMatrixState &state );
    void redoOptVector( 
            boost::unordered_map<int,double> &calculatedSpeciesIds,
            vector<double> &optVectorTS, DTLMatrixState &state );
    void updateEventQueue( bool changed, vector<int> &parentsU, 
                           DTLMatrixState &state );

    void readMatrixVariable( string line );
    void readMatrixLine( string line, int spNodeCnt, int &brTs, int &state );

    void recomputeBestReceivers( int brChanged, 
        vector<int> &speciesNodeIdsTS, vector<double> &optVectorTS, 
        boost::unordered_map<int,double> &calculatedSpeciesIds,
        DTLMatrixState &state );
    bool recomputeOptimaForTS( vector<int> &speciesNodeIdsTS,
        vector<double> &optVectorTS,
        boost::unordered_map<int,double> &calculatedSpeciesIds,
        double &alphaCost, bool recomputeAll, DTLMatrixState &state );
    void recalculateMatrixTS( int idU, int timeSlice, vector<int> &idUparents,
                              bool recomputeAll );

public:
    // reads matrix file for recalculation
    DTLMatrixRecalc( MySpeciesTree *speciesTree, CladesAndTripartitions *cat,  
               const char *matrixFileName, int maxTS );

    // reads matrix file parameters
    static void getCSVparams( const char *matrixFileName, bool &transferDead );


    void recalculateMatrix( vector<int> tsChangeList, bool verbose = false );
};

#endif
