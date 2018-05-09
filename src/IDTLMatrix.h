#ifndef IDTLMATRIX_H_ 
#define IDTLMATRIX_H_ 

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

class IDTLMatrix : public DTLMatrix {
private:
    int mU1;
    int mX1;
    int mU2;
    int mX2;

    double mILScost;    ///< cost of an ILS event


    void setBacktrack( int u1, int x1, int u2, int x2 );

    double computeTransferCost( int idUsource, int idUtarget,
                                int timeSlice, int idX, int &targetSpecies );
    double computeTL( int idU, int timeSlice, int idX );
    void computeBestReceivers( int timeSlice, int idU,
                vector<int> &speciesInTS );
    string computeOptimaForCell( int idU, int timeSlice, int idX,
                                 bool backtrack=false );
    void calculateMatrixTS( int idU, int timeSlice );

public:
    IDTLMatrix( MySpeciesTree *speciesTree, CladesAndTripartitions *cat, 
               bool computeT, bool computeTL,
               double dupliCost, double hgtCost, double lossCost, 
               double ilsCost );
    void calculateMatrix( bool verbose=false, bool unDated=false, 
            bool partialDates=false);
    void getReconciliation( int u=-1, int x=-1 );
};

#endif
