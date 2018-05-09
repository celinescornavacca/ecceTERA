#ifndef MYCLADESANDTRIPARTITIONS_H_
#define MYCLADESANDTRIPARTITIONS_H_


/**

@file
@author Wandrille Duchemin
@created 06-Nov-2015
@modified 09-Jun-2016


**/

#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "MyGeneTree.h"

#include "CladesAndTripartitions.h"



class MyCladesAndTripartitions : public CladesAndTripartitions
{
	protected:
		//pair <int,int> getNextOrRandomClade(CladeForest * CF , int idU , double probaFixed= 1.0);

		pair<int,int> getRandomClade(int idU);
		pair<int,int> getMaxClade(int idU);
		MyGeneNode *getRandomTreeAux(int &pOrd,int idU );
		MyGeneNode *getMaxTreeAux(int &pOrd,int idU );
		pair<double,vector<string> > getTreeLikelihoodAux(MyGeneNode *node, double default_proba);

		double getMaxLikelihoodAux(int idU);

	public:
		MyCladesAndTripartitions( char charSep, vector<MyGeneTree*> &geneTrees, bool verbose, bool &overflow, string &errStr, bool polytomy = false );

	    // rooted version
//	    CladesAndTripartitions( char charSep, MyTree &geneTree );

	    // Ale reader
		MyCladesAndTripartitions( char charSep, string aleFileName, string &errStr, bool verbose );


		MyGeneTree *getRandomTree();
		MyGeneTree *getMaxTree();


		double getTreeLikelihood(MyGeneNode *node, double default_proba);

		bool isCompatible(MyGeneNode * node);
		
		double getMaxLikelihood();


};

#endif
