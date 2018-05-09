#ifndef ADJ_SOLUTION_H_
#define ADJ_SOLUTION_H_


/*

This file contains several classes intended to represent adjacency scores and solutions

Created the: 09-12-2015
by: Wandrille Duchemin

Last modified the: 29-04-2016
by: Wandrille Duchemin

*/

#include <limits>
#include <cstddef>
#include <vector>
#include <iostream>
#include <Bpp/Exceptions.h>


using namespace std;
using namespace bpp;

//This class represents a c1 or c0 score in the AdjMatrix
class AdjScore
{
public:
	bool C1;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	int id1;//id of a node in the first reconciled gene tree
	int id2;//id of a node in the second reconciled gene tree

	AdjScore(bool c1, int Id1, int Id2) : C1(c1), id1(Id1), id2(Id2)
	{}

	AdjScore(){}
	~AdjScore(){}
};

//This class represent a potential solution to an adjacency history problem as a set of AdjScore and a number of gains / losses
class AdjSolution
{
public:
	vector<AdjScore> components; //the different scores composing the solution
	int NbGain; //number of gain
	int NbBreak; //number of loss

	double score; // the score of that solution

	bool coevent; // wether this specific solution is implying a co-event or not


	AdjSolution(int Nbgain,int Nbbreak,bool coev) : NbGain(Nbgain), NbBreak(Nbbreak), coevent(coev)
	{}

	AdjSolution()
	{}

	~AdjSolution()
	{
		components.clear();
	}
};

class MyMatrixAdjSolution
{	
	private:
	int dim1;   ///< first dimension of matrix
	int dim2;   ///< second dimension of matrix
	vector<AdjSolution> *solutions; ///< matrix with vector cells
	
	public:
	MyMatrixAdjSolution() {
		solutions = NULL;
	}

	~MyMatrixAdjSolution() {
		if( solutions != NULL )
		{
			//solutions->clear();
			//cout << "plop mymatrixadjsolution"<< endl;

			for(int i = 0 ; i < dim1; i++ )
			{
				for(int j = 0 ; j < dim2; j++ )
				{
					solutions[i*dim2 + j].clear();
				}
			}


			delete [] solutions;
		}
						
	}

	
	int getDim1() const{return dim1;}
	int getDim2() const{return dim2;}
	vector<AdjSolution> * getSolutions() const{return solutions;}

	void setDim(int d1, int d2);
	vector<AdjSolution> getValue(int j,int z) const;
	vector<AdjSolution> getValueSure(int j,int z) const;
	void setValues(int j, int z, vector<AdjSolution> &values);

	vector<int> getIndexNonNull(vector<AdjSolution> Vsolution,vector<int> tochoosefrom, double nullScore);
	vector<int> getIndexOnlyBest(vector<AdjSolution> Vsolution,vector<int> tochoosefrom, bool minimum);
	vector<int> getIndexCoevStatus(vector<AdjSolution> Vsolution,vector<int> tochoosefrom, bool wishedCoev);

	void printMe();

	MyMatrixAdjSolution& operator=(const MyMatrixAdjSolution& rhs)
	{
		if( solutions != NULL )
			delete [] solutions;

		if( rhs.getSolutions() == NULL ) // array not set -> nothing to copy
			return *this;

		setDim( rhs.getDim1(), rhs.getDim2()); // initialize solutions

		for(int i = 0 ; i < dim1; i++ )
		{
			for(int j = 0 ; j < dim2; j++ )
			{
				vector<AdjSolution> sol = rhs.getValue(i,j);
				setValues(i,j, sol);
			}
		}
		//copy(rhs.getSolutions(), rhs.getSolutions() + dim1*dim2, solutions); // copy solutions

		return *this;
	}


};	

#endif