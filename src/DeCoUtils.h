#ifndef DECO_UTIL_H_
#define DECO_UTIL_H_

/*

This file contains various functions used by DeCo

Created the: 02-03-2016
by: Wandrille Duchemin

Last modified the: 13-06-2016
by: Wandrille Duchemin

*/

#include "MyGeneTree.h"
#include "MySpeciesTree.h"
#include "MyCladesAndTripartitions.h"
#include "ReconciliationEvent.h"
#include "CladeReconciliation.h"
#include "ReconciledTree.h"
#include "GeneFamily.h"
#include "EquivalenceClass.h"
#include "EquivalenceClassFamily.h"
#include "AdjMatrix.h"
#include "DeCoOutputManager.h"
#include "CoEvent.h"

#include "CladesAndTripartitions.h"
#include "DTLMatrix.h"
#include "DTLGraph.h"


#include <map>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>




#define AFOREST vector <AdjTree *>
#define ECFsample vector < AFOREST *>
#define NECFsample vector < ECFsample * >





void printNewick( MyGeneNode *node, string fileName );
void printNewick( MySpeciesNode *node, string fileName );

double computeTopoScore( vector <GeneFamily *> GeneFamilyList, double TopoWeight = 1);
double computeReconciliationScore(vector <GeneFamily *> GeneFamilyList, double ReconWeight = 1);
double computeAdjacenciesScore(vector <EquivalenceClass> * RefinedEquivalenceClasses, double AGainCost, double ABreakCost, double AdjWeight = 1);
double computeAdjacenciesScore(vector <EquivalenceClassFamily> * ECFams, double AGainCost, double ABreakCost, double AdjWeight = 1);
double computeCoEventScore(vector <CoEvent> Lcoevent, double DupCost, double LossCost, double HGTCost, double ReconWeight = 1 );
double computeSystemScore( vector <GeneFamily *> GeneFamilyList, vector <EquivalenceClass> *RefinedEquivalenceClasses, vector <CoEvent> Lcoevent, double AGainCost, double ABreakCost, double DupCost, double LossCost, double HGTCost, double TopoWeight = 1, double ReconWeight = 1, double AdjWeight = 1 );
double computeSystemScore( vector <GeneFamily *> GeneFamilyList, vector <EquivalenceClassFamily> * ECFams, vector <CoEvent> Lcoevent, double AGainCost, double ABreakCost, double DupCost, double LossCost, double HGTCost, double TopoWeight = 1, double ReconWeight = 1, double AdjWeight = 1 );


//not used anymore - could still be useful
vector < EquivalenceClass > CreateEquivalenceClasses(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList, bool DeCoLTrefine,bool Verbose, bool SuperVerbose);
vector < EquivalenceClass > CreateAllPairEquivalenceClasses(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList,bool Verbose, bool SuperVerbose);

//new version
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList, MySpeciesTree *speciesTree, map < string, int > speciesChrNb, map < int,vector <float> > &speciesC0C1, map<int,string> &species_id_name, map<int, map<string,int> > &speGeneAdjNb, float Break, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose);
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose);

void ComputeEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList, map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,double WDupCost, double WLossCost, double WHgtCost, bool Verbose, bool SuperVerbose, double boltzmannTemp , double absencePenalty );
void backtrackInPlaceEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);

void ComputeOneEquivalenceClassFamily( EquivalenceClassFamily * ECF, vector <GeneFamily *> * GeneFamilyList, map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,double WDupCost, double WLossCost, double WHgtCost, bool Verbose, bool SuperVerbose, double boltzmannTemp , double absencePenalty );
void backtrackInPlaceOneEquivalenceClassFamily(EquivalenceClassFamily * ECF, vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);

void backtrackNtimesOneEquivalenceClassFamily(NECFsample * Samples, EquivalenceClassFamily * ECF, int N , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);
void backtrackNtimesNequivalenceClassFamilies(vector <NECFsample * > * AllSamples, vector <EquivalenceClassFamily > * ECFams, int N , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage);

void ComputeAndBacktrackEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList , map<int,vector<float> > speciesC0C1, double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,double WDupCost, double WLossCost, double WHgtCost, bool Verbose, bool SuperVerbose, bool galwaysGain, double gC1Advantage, double boltzmannTemp = 1, double absencePenalty = -1);
void PopulatesCoeventsFromAdjForest( vector <AdjTree * > * AForest, vector <CoEvent> * CoEventSet, int EclassId , int gfam1, int gfam2, ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool ignoreTime = false);


MySpeciesTree * getSpeciesTree(string speciesFile, bool dateAsBootstrap = false);
int processSpeciesTree( MySpeciesTree *speciesTree , bool dateAsBootstrap, bool dated, bool computeT, bool computeTD , double hgtCost, double lossCost ,bool verbose );

void checkGeneTree( int counter, MyGeneTree *geneTree);

vector<string> readGeneDistributionsFile(string ListGeneFile);
void readGeneDistributions(vector <GeneFamily *> * GeneFamilyList, MySpeciesTree *speciesTree ,vector<string> geneFiles, bool ale, bool reconciled, char charSep ,bool verbose, bool superverbose);

vector< pair <string,string > > readAdjacencies(string AdjacencyFile);



void WriteSpeciestree(MySpeciesTree * speciesTree, bool newick, string prefix);
void WriteGeneFamilyReconciledTrees( vector < GeneFamily * > * GeneFamilyList, bool newick, bool hideLosses,string prefix);

void WriteECFamTrees(EquivalenceClassFamily * ECF, bool newick,bool hideLosses, string prefix);
void WriteECFsample(ECFsample * sample,int g1, int g2, int sampleIndex, bool newick,bool hideLosses, string prefix);

void WriteAdjacencies(vector < EquivalenceClassFamily > * ECFams, string prefix);

map <string,int> storeSpeciesChrNumber(string chrNumberFile);

#endif