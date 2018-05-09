#ifndef RECONCILIATION_EVENT_H_
#define RECONCILIATION_EVENT_H_



/*

This file contains a class for reconciliation events

Created the: 28-10-2015
by: Wandrille Duchemin

Last modified the: 20-01-2016
by: Wandrille Duchemin

*/

#include "CladesAndTripartitions.h"
#include "DTLGraph.h"
#include "MySpeciesTree.h"

#include "XMLUtils.h"

#include <string>
#include <map>
#include <Bpp/Exceptions.h>

using namespace std;
using namespace bpp;


class ReconciliationEvent
{
	protected:

		int idX; //postOrder id in the species tree
		int idXl; //postOrder id in the species tree of the left child of idX
		int idXr; //postOrder id in the species tree of the right child of idX
		//for postOrder id, -1 means that it is in a dead/unsampled lineage

		string evtName; //one of : S, SL, D, T, TL, TFD, TTD, TLTD, TLFD

		float support;

		int timeSlice;

		int getIdXForNullEvent(MySpeciesTree * t, int idx);

	public:

		int getidX();
		int getidXl();
		int getidXr();
		string getevtName();
		float getsupport();
		int gettimeSlice();


		void setidX( int idx);
		void setidXl( int idxl);
		void setidXr( int idxr);
		void setevtname( string evtname);
		void setsupport( float supp);
		void settimeSlice( int ts);

		ReconciliationEvent()
		{evtName=""; timeSlice = -1;}

		ReconciliationEvent(CladesAndTripartitions * CandT, DTLGraph *graph, MySpeciesTree * speciesTree, int recIndex, int evtIndex, vector< vector<DTLGraph::MyGraph::Vertex> > *reconciliation, bool VERBOSE);

		ReconciliationEvent(string recstr, map<string, int> speciesLeafNameToId);

		ReconciliationEvent(string XMLline);

		~ReconciliationEvent()
		{}


		bool hasTimeSlice();
		bool is_empty();
		void printMe();

};




#endif