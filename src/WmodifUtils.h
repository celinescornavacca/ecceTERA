#ifndef WMODIF_UTIL_H_
#define WMODIF_UTIL_H_


/*

This file contains utils specific to Wmodif stuff like redrawing reconciliations , ...

Created the: 24-03-2016
by: Wandrille Duchemin

Last modified the: 06-06-2016
by: Wandrille Duchemin

*/


#include "GeneFamily.h"
#include "EquivalenceClassFamily.h"
#include "CoEvent.h"
#include "DeCoUtils.h"
#include "XMLUtils.h"



class GfamSave
{
public:

    ReconciledTree * recTree;
    double treeLkh;
    double RecScore;
    map <int , EquivalenceClassFamily > ECFamsMap;
    double AdjScore;
    
    vector <CoEvent> * CoEventSet;
    double CoEventScore;

    GfamSave(){}
    ~GfamSave(){}

    void clearP()
    {
        //delete recTree;
        CoEventSet->clear();
        delete CoEventSet;
    }

};

void prepRecRedrawEqClass(int GeneFamilyId, vector < EquivalenceClassFamily > * ECFams);
void prepRecRedrawCoEvent(int GeneFamilyId, vector <CoEvent> * CoEventSet, bool verbose);

void prepECFRedrawCoEvent(int ECFidToReset, vector <CoEvent> * CoEventSet, bool verbose);

void refiningAndComputingECFamily(EquivalenceClassFamily * ECF, 
									ReconciledTree * Rtree1, ReconciledTree * Rtree2,
									double Again, double Bgain, 
									bool absencePenalty, bool doAllPair, bool withTransfer, bool substractRecoToAdj, 
									double weightedDupCost, double weightedLossCost, double weightedHGTCost, 
									bool alwaysGainAtTop, double C1Advantage,bool superverbose);


void refiningAndComputingECFamilyList(vector <EquivalenceClassFamily> * ECFams, vector < GeneFamily * > * GeneFamilyList,
										int GfamIdToReset,
										double Again, double Bgain, 
										bool absencePenalty, bool doAllPair, bool withTransfer, bool substractRecoToAdj, 
										double weightedDupCost, double weightedLossCost, double weightedHGTCost, 
										bool alwaysGainAtTop, double C1Advantage, bool verbose, bool superverbose);


void GettingCoEvents(vector <CoEvent> * CoEventSet, vector <EquivalenceClassFamily> * ECFams, vector < GeneFamily * > * GeneFamilyList);


void ReDrawGeneFamily(int GfamIdToReset, bool resetTopo, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                        MySpeciesTree * speciesTree, bool datedSpeciesTree, bool boundedTS, bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        double Again, double Abreak, bool absencePenalty, bool doAllPair, bool substractRecoToAdj,
                        double weightedDupCost, double weightedLossCost, double weightedHGTCost, bool alwaysGainAtTop, double C1Advantage,
                        bool verbose, bool superverbose );


void ReDrawGeneFamily(int GfamIdToReset, bool resetTopo, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                        MySpeciesTree * speciesTree, bool datedSpeciesTree, bool boundedTS, bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        double Again, double Abreak, bool absencePenalty, bool doAllPair, bool substractRecoToAdj,
                        double weightedDupCost, double weightedLossCost, double weightedHGTCost, bool alwaysGainAtTop, double C1Advantage,
                        bool verbose, bool superverbose ,
                        bool trySeveralECFSolutions, double Temperature, double RecWeight, bool deleteCoev);


void ReplaceGeneFamilyRecTree(int GfamIdToReset, ReconciledTree * newRTree, double treeLkh, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                        bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        double Again, double Abreak, bool absencePenalty, bool doAllPair, bool substractRecoToAdj,
                        double weightedDupCost, double weightedLossCost, double weightedHGTCost, bool alwaysGainAtTop, double C1Advantage,
                        bool verbose, bool superverbose );


/*void ReDrawECFamily(int ECFidToReset , 
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool alwaysGainAtTop, double C1Advantage,bool superverbose);*/

void ReplaceECFamily(int ECFidToReset , vector < vector< AdjTree *> * > * NewAdjacencyTreesVector, vector < CoEvent > * NewCoEventSet,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool superverbose);


bool acceptProposition(double scoreInit, double scoreProposed, double Temperature=1);


bool tryECFreDraw(int ECFidToReset , double Temperature, double oldCoEventScore,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool alwaysGainAtTop, double C1Advantage,bool superverbose,
                    double DupliCost, double LossCost, double HGTCost, double RecWeight);

bool tryGfamRedraw(int GfamIdToReset, bool resetTopo, bool trySeveralECFSolutions, double Temperature,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * * CoEventSetPointer, MySpeciesTree * speciesTree,
                    int maxTS, bool datedSpeciesTree, bool boundedTS, bool withTransfer, 
                    double DupliCost, double HGTCost, double LossCost,
                    double TopoWeight, double RecWeight, double AdjWeight,
                    double Again, double Abreak,
                    bool absencePenalty, bool doAllPair, bool substractRecoToAdj, 
                    bool alwaysGainAtTop, double C1Advantage,
                    bool verbose, bool superverbose);


bool updateGfam(int GfamIdToReset, int maxNbTry, bool resetTopo, bool trySeveralECFSolutions, double Temperature,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * * CoEventSetPointer, MySpeciesTree * speciesTree,
                    int maxTS, bool datedSpeciesTree, bool boundedTS, bool withTransfer, 
                    double DupliCost, double HGTCost, double LossCost,
                    double TopoWeight, double RecWeight, double AdjWeight,
                    double Again, double Abreak,
                    bool absencePenalty, bool doAllPair, bool substractRecoToAdj, 
                    bool alwaysGainAtTop, double C1Advantage,
                    bool verbose, bool superverbose);




vector <AdjTree * > * readAdjForest(string filename, int gfam1, int gfam2, bool gainAtRoot, bool VERBOSE);

void LoadGeneFamilyReconciledTrees( vector < GeneFamily * > * GeneFamilyList, string prefix, MySpeciesTree * Stree, bool VERBOSE );
void LoadECFamTrees(EquivalenceClassFamily * ECF, string prefix,bool gainAtRoot, bool VERBOSE);


#endif