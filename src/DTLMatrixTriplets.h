#ifndef DTLMATRIX_TRIPLETS_H_ 
#define DTLMATRIX_TRIPLETS_H_ 

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


class DTLMatrixTriplets : public DTLMatrix {
private:
    int mParetoMod;
    double mNDu;
    double mNDl;
    double mNLu;
    double mNLl;
    list<double> mPolygonX; 
    list<double> mPolygonY;

    MyMatrixT *mMatrixT;    ///< suboptimal matrix (dtl triplets) 

    vector<EventTriplet> mTripletList;  ///< subopt triplets 
    vector< vector<EventTriplet> > mAllTripletsTS;
    vector<EventTriplet> mAlphaTriplets; ///< alpha triplets

    vector< pair<EventTriplet,int> > **mVTt1; 
            ///< triplet suboptimal best receivers
    vector<EventTriplet> **mVTt2; 
            ///< triplet suboptimal second best receivers


    bool tripInList( EventTriplet trip, vector< pair<EventTriplet,int> > list );

    virtual void computeList( int idUa, int idXa, int idUb, int idXb,
            double lowestCost, double curBestEpsCost, 
            double otherCost, EventTriplet eventTrip );
    virtual void computeList( int idUa, int idXa, 
            double lowestCost, double curBestEpsCost, 
            double otherCost, EventTriplet eventTrip );
    virtual void addZeroCost() {
                EventTriplet newTrip;
                mTripletList.push_back( newTrip );
    }
    virtual void computeTransferCostSub( int optimumSub, int idUsub,    
        int idUother, double costThisSplit, DTLMatrixState &state );
    virtual void pushList( DTLMatrixState &state ) {
        // remove duplicates and non-minimal
//if( state.idU==1 && (state.idX==46 || state.idX==40 )) 
//    cout << state.idU << "," << state.idX << " pushList: ";
        cleanTriplets( mTripletList, state );
// replace mTripletList by index into mAllTripletsTS?

//cout << "PUSHING " << state.idU << "," << state.idX << " size=" 
//    << mAllTripletsTS.size() << endl;
//for( size_t i=0; i<mTripletList.size(); i++ ) {
//    cout << "pushing " << mTripletList[i] << endl;
//}
        mAllTripletsTS.push_back( mTripletList );
        mTripletList.clear();
    }
    virtual void setAlpha( int alphaIdx ) {
        mAlphaTriplets = mAllTripletsTS[alphaIdx];
    }
    virtual void computeTransferLossCostVT( int idx, DTLMatrixState &state );
    virtual void setFinalCosts( int index, DTLMatrixState &state ) {
//if( state.idU==1 && (state.idX==46 || state.idX==40 )) 
//    cout << state.idU << "," << state.idX << " FINAL: ";
        cleanTriplets( mAllTripletsTS[index], state );
//cout << "SETTING index=" << index << ": " << state.idU << "," << state.idX << " size=" 
//    << mAllTripletsTS[index].size() << endl;
//for( size_t i=0; i<mAllTripletsTS[index].size(); i++ ) 
//   cout << " " << mAllTripletsTS[index][i];
//cout << endl;
        mMatrixT->setValues( state.idU, state.idX, mAllTripletsTS[index] );
    }
    virtual void initTS() {
        mAllTripletsTS.clear();
    }
    void cleanUp();
    virtual void initCalculation();

//    virtual void computeVTs( vector<MyNode*> &speciesNodesTS, 
//            DTLMatrixState &state );
//    virtual void computeVTsNoSub( vector<MyNode*> &speciesNodesTS, 
//            DTLMatrixState &state );
    virtual void computeVTs( bool noSub, bool recompute, int timeSlice, 
                             int idU, vector<int> &speciesNodeIdsTS );

    void intersectionLines( double &intX, double &intY,
        double x1, double y1, double x2, double y2, EventTriplet line );
    int value( EventTriplet line, double x, double y );
    bool intersectionPolygonLine( list<double> &polygonX, 
                list<double> &polygonY, EventTriplet line );
    void cleanTriplets( vector<EventTriplet> &v, DTLMatrixState &state );

    void eventCounts( const char *event, int &d, int &t, int &l );
    virtual void addVertices( DTLGraph &graph, DTLGraph::MyGraph::Vertex z, 
            vector<DTLGraph::MyGraph::Vertex> &qList, int u, int x,
            const char*event, double curCost, double otherCost, 
            EventTriplet zTrip );
    virtual void addVertices( DTLGraph &graph, DTLGraph::MyGraph::Vertex z, 
            vector<DTLGraph::MyGraph::Vertex> &qList, 
            int u1, int x1, int u2, int x2, const char*event, 
            double curCost, double otherCost, EventTriplet zTrip );
    virtual vector<DTLGraph::MyGraph::Vertex> getRootNodes( DTLGraph &graph,
                                                            bool keepRoot );

public:
    DTLMatrixTriplets( MySpeciesTree *speciesTree, CladesAndTripartitions *cat, 
               double WGDCost, bool fixedCosts, bool computeT, bool computeTL,
               double dupliCost, double hgtCost, double lossCost, 
               int maxTS, double weight, bool mUseBestSplits,
                double ilsCost, 
                int paretoMod, double nD, double nL, double nDL );
    ~DTLMatrixTriplets();
};
#endif
