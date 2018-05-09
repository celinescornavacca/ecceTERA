
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
Wrapper for a matrix of doubles.

*/


#include "MyMatrix.h"
#include <iostream>
#include <Bpp/Exceptions.h>

		
/**
 * Create vectors for the matrix.
 *
 * All values are initialized to -1.
 */
void MyMatrix::setDim( int d1, int d2 )
{
    mDim1 = d1;
    mDim2 = d2; 
    
    costs.resize( mDim1 );
    for (int j=0; j<mDim1; j++)
        costs[j].resize( mDim2, -1 );
}


/**
 * Return the matrix value at the given coordinantes, checking
 * that the coordinantes are valid and the cell has been assigned 
 * (not -1).
 *
 * @return matrix value
 */
double MyMatrix::getValue( int j, int z ) 
{
    return costs[j][z];
}

/**
 * Return the matrix value at the given coordinantes, without
 * checking validity.
 *
 * @return matrix value
 */
double MyMatrix::getValueSure( int j, int z ) 
{	
    if( j>mDim1 || j<0 || z>mDim2 || z<0 ) {
        cout << "coordinates=" << j << "," << z << endl;
        throw bpp::Exception("MyMatrix::getValue: invalid coordinate");
    }

    if( costs[j][z] == -1 ) {
        cout << "error: tryng to access cell (" << j 
            << "," << z << ") which hasn't initialized yet!" << endl;
        throw bpp::Exception("MyMatrix::getValue: invalid coordinate");
    }

    return costs[j][z]; 
}

/**
 * Set the matrix value.
 */
void MyMatrix::setValue( int j, int z, double value ) 
{
     costs[j][z] = value;	
}
