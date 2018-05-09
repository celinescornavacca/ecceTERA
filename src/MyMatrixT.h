
#ifndef MYMATRIXT_H_
#define MYMATRIXT_H_

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

*/

#include <vector>
using namespace std;

struct EventTriplet
{
    int d, t, l;
    double cost;

    // for sorting (d first, then t, then l)
    bool operator < (const EventTriplet& trip) const {
        // non-pareto sort
        //return cost < trip.cost;
       
        // pareto sort
        if( d == trip.d ) {
            if( t == trip.t ) 
                return (l < trip.l);
            else
                return (t < trip.t);
        }
        return (d < trip.d);
    }

    bool operator == (const EventTriplet& trip) const {
        if( d == trip.d && l == trip.l && t == trip.t )
            return true;
        return false;
    }

    EventTriplet() : d(0), t(0), l(0), cost(0) {}
    EventTriplet( int dIn, int tIn, int lIn, double costIn ) 
        : d(dIn), t(tIn), l(lIn), cost(costIn) {}

    friend ostream& operator<<(ostream& os, EventTriplet const & t) {
        return os << "<" << t.d << "," << t.t << "," << t.l 
                  << ":" << t.cost << ">";
    }
};


class MyMatrixT
{	
	private:
	int dim1;   ///< first dimension of matrix
	int dim2;   ///< second dimension of matrix
	vector<EventTriplet> *triplets; ///< matrix with vector cells
	
	public:
    MyMatrixT() {
        triplets = NULL;
    }

    ~MyMatrixT() {
        if( triplets != NULL ) 
            delete [] triplets;			
    }

	void setDim(int d1, int d2);
	vector<EventTriplet> getValue(int j,int z);
	vector<EventTriplet> getValueSure(int j,int z);
	void setValues(int j, int z, vector<EventTriplet> &values);

};	

#endif /*MYMATRIXV_H_*/
