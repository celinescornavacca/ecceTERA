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


@section DESCRIPTION
Wrapper for a matrix with each cell containing a list of doubles.

*/
#include "MyMatrixV.h"
#include <Bpp/Exceptions.h>

		
/**
 * This function allocates memory for a matrix of size d1 x d2.
 */
void MyMatrixV::setDim(int d1, int d2)
{
    dim1 = d1;
    dim2 = d2; 

    costs = new vector<double>[dim1*dim2];
}


/**
 * Return the matrix values at the given coordinantes, checking
 * that the coordinantes are valid and the cell has been assigned. 
 *
 * @return list of matrix values
 */
vector<double> MyMatrixV::getValue( int j, int z ) 
{
    if( j>dim1 || j<0 || z>dim2 || z<0 ) {
        cout << "coordinates=" << j << "," << z << endl;
        throw bpp::Exception("MyMatrixV::getValue: invalid coordinate");
    }

    int i = j*dim2 + z;
/* This is possible with suboptimal.
   if(costs[i].size() == 0) {
        cout << "error: matrix cell " << j << ", " << z << " not initialized"
            << endl;
        throw Exception( "MyMatrixV::getValue: tryng to access a cell" 
                " that hasn't initialized yet!" );
    }
*/            
    return costs[i];
}


/**
 * Return the matrix values at the given coordinantes, without
 * checking validity.
 *
 * @return list of matrix values
 */
vector<double> MyMatrixV::getValueSure( int j, int z ) 
{	
    //without checking that it is not -1, 
    //because in this case it may be -1			
    int i = j*dim2 + z;
    return costs[i];
}


/**
 * Set the matrix values.
 */
void MyMatrixV::setValues( int j, int z, vector<double> &values ) 
{
    int i = j*dim2 + z;
    costs[i] = values;
}
