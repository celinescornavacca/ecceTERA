#ifndef DTLMATRIXNETWORK_H_ 
#define DTLMATRIXNETWORK_H_ 

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
#include "CladesAndTripartitions.h"
#include "DTLGraph.h"
#include "MyNetwork.h"


class DTLMatrixNetwork : public DTLMatrix {
	
    int secondaryFather( int idX, int idXl, int idXr );
	
	public: 
	
    DTLMatrixNetwork( MyNetwork *speciesNetwork, CladesAndTripartitions *cat, 
               double WGDCost, bool fixedCosts, bool computeT, bool computeTL,
               double dupliCost, double hgtCost, double lossCost, 
               int maxTS, double weight, bool mUseBestSplits,
               double ilsCost ): DTLMatrix( (MySpeciesTree*) speciesNetwork,  cat, 
                WGDCost,  fixedCosts,  computeT,  computeTL,
                dupliCost,  hgtCost,  lossCost, 
                maxTS,  weight,  mUseBestSplits,
                ilsCost ){}
  
    void calculateMatrixNoSub( int idU );

    void computeNullCost( int idXl, DTLMatrixState &state, double &optCost, 
            BestSplit &bestSplit );
            
    bool computeTransferLossCost( int idXt, DTLMatrixState &state, double &optCost, 
            BestSplit &bestSplit );
               
    virtual void computeOptimaForCladeSplit( int toCompute, int splitIdx,
            pair<int,int> cladeSplit, DTLMatrixState &state, 
            double &optCost, BestSplit &bestSplit );
    virtual void computeOptimaForSpeciesSplit( int toCompute,
            pair<int,int> cladeSplit, DTLMatrixState &state, 
            double &optCost, BestSplit &bestSplit );
    
    void computeTransferCost( int idUl, int idUr, int idXl, int idXr, double costThisSplit, 
            DTLMatrixState &state, double &optCost, BestSplit &bestSplit );

    virtual void createVertices( DTLGraph &graph, 
                DTLGraph::MyGraph::Vertex pairVertex, 
                vector<DTLGraph::MyGraph::Vertex> &qList, 
                vector<int> *** allBestReceivers );
            

};

#endif 
