// Created by: Celine Scornavacca


/**

@file ecceTERA.cpp
@author Celine Scornavacca
@author Edwin Jacox
@version 1.2.4 
@date 14/09/2016

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

**/

#include <sstream>
#include <sys/stat.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>


#include "CladesAndTripartitions.h"
#include "DTLMatrix.h"
#include "DTLMatrixTriplets.h"
#include "DTLGraph.h"
#include <Bpp/Phyl/Io/Newick.h>
#include "parameters.h"

			
#include <Bpp/App/ApplicationTools.h>




/**
 * Return the memory used in MB.
 */
void printMemory( string msg ){

#ifdef __linux__
    // This gets the value in KB.
    // VmPeak = peak virtual memory (all memory requested)
    // VmRSS = resident memory (currently in RAM)
    FILE* file = fopen("/proc/self/status", "r");
    char line[128];
    while( fgets(line, 128, file) != NULL ) {
        if (strncmp(line, "VmPeak:", 6) == 0)
            break;
    }
    fclose(file);

    // parse line
    int len = strlen(line);
    int result = 0;
    for( int i=0; i<len; i++ ) {
        if( line[i] >= '0' && line[i] <= '9' ) {
            result = result*10 + line[i]-'0';
        }
    }
    result = round( result/1000 );
   
    if( gBoolParams.find("verbose")->second )
        cout << msg << " memory usage (MB): " << result << endl;

#elif __APPLE__
/*
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    if( KERN_SUCCESS != task_info(mach_task_self(),
        TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count) )
    {
            return -1;
    }
    // resident size is in t_info.resident_size;
    // virtual size is in t_info.virtual_size;
*/
#endif
}


/**
 * Read a file that maps gene names to species and add them
 * to the clades and tripartitions.
 */
void readGeneMappingFile(
    CladesAndTripartitions *cladesAndTripartitions ) ///< clades structure
{
    const char *mappingFileName 
        = gStringParams.find("gene.mapping.file")->second.c_str();

    ifstream fileStream( mappingFileName, ios::in );  
    if( !fileStream) { 
        cerr << "Failed to open file: "
            << gStringParams.find("gene.mapping.file")->second << endl;
        exit(1);
    }

    string description;
    vector<string> speciesNames;
    vector<string> geneNames;
	while( !fileStream.eof() ) {
        string temp;
        // Copy current line to temporary string
        getline(fileStream, temp, '\n');  
        if( temp == "" ) // blank line
            continue;

        boost::char_separator<char> sep(" ");
        Tokenizer tok( temp, sep );
        Tokenizer::iterator iter=tok.begin();

        geneNames.push_back( *iter );
        iter++;
        if( iter==tok.end() ) {
            cerr << "Invalid line, missing species name in "
                << gStringParams.find("gene.mapping.file")->second << endl;
            exit( 1 );
        }
        speciesNames.push_back( *iter );
    }

    string errStr;
    cladesAndTripartitions->mClades.mapSpeciesNames( 
                                speciesNames, geneNames, errStr );
    if( errStr != "" ) {
        delete cladesAndTripartitions;
        cerr << "Error in file " 
             << gStringParams.find("gene.mapping.file")->second 
             << ": " << errStr << endl;
        exit( 1 );
    }
}


/**
 * Create the clades and tripartitions strucutre.
 */
CladesAndTripartitions *getCladesAndTripartitions(
    bool amalgamation, ///< do amalgamation if true
    vector<MyGeneTree*> geneTrees, ///< gene trees
    MyGeneTree *bootstrapGeneTree, ///< gene tree with bootstrap values
    bool polytomic ) ///< 
{
    if( gBoolParams.find("clades.only")->second ) // output clade count
        gBoolParams.find("verbose")->second = true;

    CladesAndTripartitions *cladesAndTripartitions;
    if( amalgamation && gBoolParams.find("ale")->second ) {
        // read a saved version from a file
        string errStr = "";
        cladesAndTripartitions = 
            new CladesAndTripartitions( 
                    gStringParams.find("char.sep")->second[0], 
                    gStringParams.find("gene.file")->second, 
                    errStr, gBoolParams.find("verbose")->second );
        if( errStr != "" ) {
            delete cladesAndTripartitions;
            cout << errStr << endl;
            exit(1);
        }
    } else if( polytomic || amalgamation ) 
    {
        // polytomic or amalgamated
        bool overflow = false;
        string errStr = "";
        bool polytomicUnrooted = false;
        if( polytomic && gIntParams.find("resolve.trees")->second == 1 )
            polytomicUnrooted = true;
        cladesAndTripartitions = 
            new CladesAndTripartitions( 
                    gStringParams.find("char.sep")->second[0], 
                    geneTrees, gBoolParams.find("verbose")->second, 
                    overflow, errStr, 
                    polytomic, bootstrapGeneTree, polytomicUnrooted, 
                    gBoolParams.find("print.clade.info")->second );
        if( overflow ) 
            cerr << "Too many possible gene trees due to polytomies" << endl;
        if( errStr != "" ) 
            cerr << errStr << endl;
        if( overflow || errStr != "" ) {
            delete cladesAndTripartitions;
            exit( 1 );
        }
        if( gBoolParams.find("verbose")->second ) {
            bpp::ApplicationTools::displayTime("Computing ALE done:");
            printMemory( "Computing ALE" );
        }
    } else {
        // read a simple tree
        cladesAndTripartitions = new CladesAndTripartitions( 
                   gStringParams.find("char.sep")->second[0], *(geneTrees[0]) );
    }

    if( gStringParams.find("gene.mapping.file")->second != "none" ) 
        readGeneMappingFile( cladesAndTripartitions );

    if( gStringParams.find("print.ale")->second != "none" ) {
        cladesAndTripartitions->saveALE( 
                gStringParams.find("print.ale")->second );           
    }

    // just want clade and tripartition count
    if( gBoolParams.find("clades.only")->second ) {
        delete cladesAndTripartitions;
        exit(1);
    }

    // experimental option to skip some clades
    if( gIntParams.find("skip.clades")->second > 0 ) {
        int skipCount = cladesAndTripartitions->excludeClades( 
                                gIntParams.find("skip.clades")->second );
        if( skipCount == -1 ) {
            cerr << "Invalid skip.clades - too many clades removed." << endl;
            exit( 1 );
        }
        if( gBoolParams.find("verbose")->second )
            cout << "Skipping " << skipCount << " clades." << endl;
    }

    return cladesAndTripartitions;
}

/**
 * Create matrix and run it.
 */
DTLMatrix *createAndRunMatrix(
        bool returnTree,  ///< return backtrack gene tree
        MySpeciesTree *speciesTree, ///< species tree
        CladesAndTripartitions *cladesAndTripartitions, ///< gene clades
        int maxTS,  ///< maximum time slice
        vector<int> &changedTimeSlices, ///< list of changed time slices
        double &inEps ) ///< suboptimal epsilon value used
{
    inEps = 0;
    DTLMatrix *dtlMatrix = NULL;

    // record best splits if needed later
    bool useBestSplits = false;
    if( gBoolParams.find("verbose")->second 
        || gBoolParams.find("print.newick")->second
        || gStringParams.find("print.reroot.file")->second != "none" 
        || returnTree )
    {
        useBestSplits = true; // required to print event counts or newick
    }

    // Create matrix
    if( gIntParams.find("pareto.mod")->second > 0 ) {
        // create triplets matrix
        if( gBoolParams.find("verbose")->second )
            cout << "============ pareto mod: " 
                    << gIntParams.find("pareto.mod")->second << endl;
        dtlMatrix = new DTLMatrixTriplets( speciesTree, 
                        cladesAndTripartitions, 
                        gDoubleParams.find("WGD.cost")->second, 
                        gFixedCosts, 
                        gBoolParams.find("compute.T")->second,
                        gBoolParams.find("compute.TL")->second,
                        gDoubleParams.find("dupli.cost")->second, 
                        gDoubleParams.find("HGT.cost")->second, 
                        gDoubleParams.find("loss.cost")->second, 
                        maxTS, 
                        gDoubleParams.find("weight.amalgamation")->second, 
                        useBestSplits, 
                        gDoubleParams.find("ils.cost")->second, 
                        gIntParams.find("pareto.mod")->second, 
                        gDoubleParams.find("nD")->second, 
                        gDoubleParams.find("nL")->second, 
                        gDoubleParams.find("nDL")->second );
    } else
        dtlMatrix = new DTLMatrix( speciesTree, 
                cladesAndTripartitions, 
                gDoubleParams.find("WGD.cost")->second, 
                gFixedCosts,
                gBoolParams.find("compute.T")->second, 
                gBoolParams.find("compute.TL")->second, 
                gDoubleParams.find("dupli.cost")->second, 
                gDoubleParams.find("HGT.cost")->second, 
                gDoubleParams.find("loss.cost")->second, 
                maxTS, 
                gDoubleParams.find("weight.amalgamation")->second, 
                useBestSplits,
                gDoubleParams.find("ils.cost")->second );

    // run calculation
    dtlMatrix->calculateMatrix( 
            gBoolParams.find("verbose")->second, 
            gIntParams.find("max.iterations")->second, 
            gBoolParams.find("fix.dtl.costs")->second,
            gIntParams.find("dated")->second );

    if( gFixedCosts && gBoolParams.find("verbose")->second ) {
    // count events needs to be updated to handle variable costs
        int duplications = 0;
        int transfers = 0; 
        int losses = 0;
        int ils = 0;
        double N_am = dtlMatrix->countEvents(duplications, 
                                                transfers, losses, ils );
        cout << "nAm=" << N_am << ", duplications=" << duplications 
                << ", transfers=" << transfers << ", losses=" << losses;
        if( speciesTree->hasILS() )
            cout << ", ils=" << ils;
        cout <<endl;
    }

    if( gDoubleParams.find("suboptimal.epsilon")->second > 0 
        || gIntParams.find("pareto.mod")->second > 0 
        || ( gBoolParams.find("print.info")->second 
                && gDoubleParams.find("min.recs")->second > 0 ) )
    {
        // rerun for suboptimal
        inEps = gDoubleParams.find("suboptimal.epsilon")->second;
        if( gNoEpsilon ) {
            inEps = std::numeric_limits<double>::max();
        } else if( !gBoolParams.find("real.epsilon")->second )
            inEps *= dtlMatrix->getBestCost(
                        gBoolParams.find("gene.origination.species.root")
                                    ->second ) /100;

        dtlMatrix->calculateMatrix( gBoolParams.find("verbose")->second, 
                    1, gBoolParams.find("fix.dtl.costs")->second, 
                    gIntParams.find("dated")->second, inEps );
    }



    printMemory( "Matrix" );

    return dtlMatrix;
}

/**
 * Print newick gene trees if requested.
 */
MyGeneTree *printNewickTrees(
    bool returnTree,  ///< return backtrack gene tree
    bool append,     ///< append to existing if true
    bool polytomic, ///< gene tree is polytomic
    bool amalgamation, ///< gene trees are amalgamated
    vector<MyGeneTree*> geneTrees, ///< gene trees
    CladesAndTripartitions *cladesAndTripartitions, ///< gene clades
    DTLMatrix *dtlMatrix ) ///< the matrix
{
    // newick gene tree
    MyGeneTree *treeToReturn = NULL;
    if( returnTree || gBoolParams.find("print.newick")->second ) {
        MyGeneTree *tree = NULL;
        MyGeneNode *node = NULL;
        if( polytomic || amalgamation ) {
            node = dtlMatrix->backtrack( false,
                gBoolParams.find("gene.origination.species.root")->second );
            // create a tree to delete node structre easily
            tree = new MyGeneTree(*node); 
        } else {
            tree = cladesAndTripartitions->getTree();
            node = tree->getRootNode();
        }

        if( gBoolParams.find("print.newick")->second ) {
            string pathName = gPathPrefix
                + gStringParams.find("print.newick.gene.tree.file")->second;
            tree->printNewick( pathName, append );
        }

        if( !returnTree )
            delete tree;
        else
            treeToReturn = tree;
    }

    if( gStringParams.find("print.reroot.file")->second != "none" 
            && amalgamation && geneTrees.size() == 1 ) 
    {
        MyGeneNode *newRoot = dtlMatrix->backtrack( false,
                gBoolParams.find("gene.origination.species.root")->second );
        vector<string> sonLeaves = bpp::TreeTemplateTools::getLeavesNames( 
                                                *(newRoot->getSon(0)) );
        MyGeneTree *geneTreeCopy = new MyGeneTree( *(geneTrees[0]) );
        bool success = geneTreeCopy->reroot( sonLeaves, 
                            gIntParams.find("reroot.proportion")->second );
        if( success ) {
            string pathName = gPathPrefix
                + gStringParams.find("print.reroot.file")->second;
            geneTreeCopy->printNewick( pathName, append );
        } else 
            cout << "Failed to reroot tree. =======================" << endl;
    }

    return treeToReturn;
}

/**
 * Redo matrix and graph while number of solutions is less than min.recs.
 *
 * @return True if continue to loop.
 */
bool minRecsLoop(
        bool &notTooBig,    ///< too many solutions if true
        double &numberSolutions, ///< number of reconciliations
        double inEps,   ///< epsilon value, if used
        DTLMatrix *&dtlMatrix, ///< current matrix
        DTLGraph &graph) ///< current graph
{
    if( numberSolutions >= 0 ) { // Isn't this always true?
        gDoubleParams["suboptimal.epsilon"] = 
            gDoubleParams.find("suboptimal.epsilon")->second 
            + gDoubleParams.find("min.recs.increment")->second;
        if( gDoubleParams.find("suboptimal.epsilon")->second 
                 > gDoubleParams.find("max.epsilon")->second ) 
        {
            cout << "Surpassed maximum epsilon (" 
                 << gDoubleParams.find("max.epsilon")->second
                 << ")." << endl;
            return false;
        }

        inEps = gDoubleParams.find("suboptimal.epsilon")->second;
        if( !gBoolParams.find("real.epsilon")->second )
            inEps *= dtlMatrix->getBestCost(
                            gBoolParams.find("gene.origination.species.root")
                                        ->second ) /100;
        cout << numberSolutions
             << " reconciliations less than min, trying epsilon="
             << inEps << endl;
        dtlMatrix->calculateMatrix( gBoolParams.find("verbose")->second,
                1, gBoolParams.find("fix.dtl.costs")->second, 
                gIntParams.find("dated")->second, 
                inEps );
    }
    graph = dtlMatrix->constructGraph( 
                gBoolParams.find("verbose")->second );
    notTooBig = graph.countReconciliationNumberAndCheck( inEps,
      gBoolParams.find("keep.only.canonical.reconciliations")->second, 
                gBoolParams.find("verbose")->second, 
                gBoolParams.find("weighted.support")->second );
    numberSolutions = graph.getNumberSolutions( false ); // unweighted

    return true;
}

/**
 */
void printReconciliations(
    double numberSolutions, ///< number of reconciliations
    string ext,     ///< file name extension
    vector<MyGeneTree*> geneTrees, ///< gene trees
    map<string,double> &eventSupports, ///< event supports
    CladesAndTripartitions *cladesAndTripartitions, ///< gene clades
    DTLMatrix *dtlMatrix, ///< the matrix
    MySpeciesTree *speciesTree, ///< species tree
    DTLGraph &graph ) ///< the current graph
{

    // which reconciliations to print
    string whichRecs = "all";
    if( gIntParams.find("pareto.mod")->second > 0 )
        whichRecs = "allTriplets";

    // print median reconciliations
    bool isConsistent = false;
    if( gBoolParams.find("subopt.support")->second ) {
        // make new graph and remove non-optimal nodes
        DTLGraph optGraph = dtlMatrix->constructGraph( 
                                gBoolParams.find("verbose")->second );
        optGraph.pruneNonoptimal();

        bool notTooBig = optGraph.countReconciliationNumberAndCheck( 0,
             gBoolParams.find("keep.only.canonical.reconciliations")
                    ->second,
                    gBoolParams.find("verbose")->second, 
                    gBoolParams.find("weighted.support")->second );
        if( !notTooBig ) {
            cout << "Optimal umber of reconciliations bigger than "
                << "std::numeric_limits<double>::max()\n";
            cout << "How does this happen?" << endl;
            exit(1);
        }

        // create map of eventName-> support/numSolutions with graph
        string pathName = "";
        if( gIntParams.find("print.reconciliations")->second == 1 ) {
            pathName = gPathPrefix + 
                gStringParams.find("print.reconciliations.file")->second
                + ext;
            graph.getEventSupports( eventSupports );
        }
        optGraph.printReconciliation( whichRecs, pathName,
                gBoolParams.find("sylvx.reconciliation")->second, 
                gBoolParams.find("recPhyloXML.reconciliation")->second, 
                gBoolParams.find("check.time.consistency")->second, 
                isConsistent, eventSupports );
    } else {
        string pathName = "";
        if( gIntParams.find("print.reconciliations")->second == 1 ) 
            pathName = gPathPrefix + 
                gStringParams.find("print.reconciliations.file")->second
                + ext;
        graph.printReconciliation( whichRecs, pathName,
                gBoolParams.find("sylvx.reconciliation")->second, 
                gBoolParams.find("recPhyloXML.reconciliation")->second, 
                gBoolParams.find("check.time.consistency")->second, 
                isConsistent, eventSupports );
    }

    if( gIntParams.find("print.reconciliations")->second > 0
        && gBoolParams.find("sylvx.reconciliation")->second ) 
    {
        // print a sylvx species tree
        string fileName = gPathPrefix  
                    + gStringParams.find("print.reconciliations.file")->second
                    + "_sylvx_species.txt";
        MySpeciesTree *tree = speciesTree->getPostorderTree( true );
        tree->printNewick( fileName );
        delete tree;
    }

    // check time consistency
    if( gBoolParams.find("check.time.consistency.all")->second )  {
        if( !graph.checkTimeConsistencies( 
                 gIntParams.find("consistency.limit")->second, true ) ) 
            cout << "FOUND TIME INCONSISTENCIES" << endl;
    } else if( gBoolParams.find("check.time.consistency")->second 
               && !isConsistent )  
    {
        int result = graph.checkTimeConsistencies( 
                gIntParams.find("consistency.limit")->second, false );
        if( result == 1 )
            cout << "time consistent" << endl;
        else if( result == 2 )
            cout << "time consistency unknown, limit reached (" 
                 << gIntParams.find("consistency.limit")->second 
                 << ")" << endl;
        else
            cout << "NOT TIME CONSISTENT" << endl;
    }


    // print all reconciliations
    if( gIntParams.find("print.reconciliations")->second == 2 ) {
        string pathName = gPathPrefix + 
          gStringParams.find("print.reconciliations.file")->second;
        if( gBoolParams.find("check.time.consistency")->second )
            pathName += "_consistent";
        long printCount = graph.printAllReconciliations( pathName, 
                        gBoolParams.find("sylvx.reconciliation")->second, 
                    	gBoolParams.find("recPhyloXML.reconciliation")->second, 
                        gBoolParams.find("check.time.consistency")->second,
                        gIntParams.find("solution.limit")->second );
        if( printCount != numberSolutions )
            cout << "Printed " << printCount << "/" <<
                numberSolutions << " reconciliations" << endl;
    }

}

/**
 * Create graph and count reconciliations. Print them if requested.
 */
void makeGraph(
    double inEps,   ///< epsilon value, if used
    string ext,     ///< file name extension
    MyGeneTree *polytomyTree, ///< polytomy tree to expand
    vector<MyGeneTree*> geneTrees, ///< gene trees
    CladesAndTripartitions *cladesAndTripartitions, ///< gene clades
    DTLMatrix *&dtlMatrix, ///< the matrix
    bool rerunParetoMod, ///< set to return if needed
    MySpeciesTree *speciesTree ) ///< species tree
{ 
    DTLGraph graph = dtlMatrix->constructGraph( 
                        gBoolParams.find("verbose")->second );
//graph.checkScore( gDoubleParams.find("dupli.cost")->second, gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second );

    // print a non-canonical graph
    if( !gBoolParams.find("keep.only.canonical.reconciliations")->second 
        && gBoolParams.find("print.graph")->second ) 
    {
        string pathName = gPathPrefix + 
                gStringParams.find("print.graph.file")->second + ".sif";
        graph.printGraph( pathName,
                gBoolParams.find("internal.graph.ids")->second, false );
    }

    // get number of solutions
    bool notTooBig = graph.countReconciliationNumberAndCheck(
            inEps, 
            gBoolParams.find("keep.only.canonical.reconciliations") ->second, 
            gBoolParams.find("verbose")->second, 
            gBoolParams.find("weighted.support")->second );
    double numberSolutions = graph.getNumberSolutions();

    // repeat calculations if not enough solutions
    if( gDoubleParams.find("min.recs")->second > 0 ) 
    while( notTooBig && !gNoEpsilon  
           && numberSolutions < gDoubleParams.find("min.recs")->second )
    {
        if( !minRecsLoop( notTooBig, numberSolutions, inEps, dtlMatrix, graph ))
            break;
    }

    // print a canonical graph
    if( gBoolParams.find("keep.only.canonical.reconciliations")->second ) {
        ext += "_canonical";
        if( gBoolParams.find("print.graph")->second ) {
            string pathName = gPathPrefix + 
                gStringParams.find("print.graph.file")->second + ext+".sif";
            graph.printGraph( pathName,
                    gBoolParams.find("internal.graph.ids")->second, true );
        }

    }

    // print reconciliations and do other operations on reconciliations
    map<string,double> eventSupports;
    if( !notTooBig ) {
        cout << "Number of reconciliations bigger than "
            << "std::numeric_limits<double>::max()\n";
        string pathName = "";
        if( gIntParams.find("print.reconciliations")->second > 0 ) {
            pathName = gPathPrefix + 
                gStringParams.find("print.reconciliations.file")->second
                + ext;
        }
        bool isConsistent = false;
        graph.printReconciliation( "random", pathName,
            gBoolParams.find("sylvx.reconciliation")->second, 
            gBoolParams.find("recPhyloXML.reconciliation")->second, 
            gBoolParams.find("check.time.consistency")->second, 
            isConsistent, eventSupports);
    } else if( numberSolutions == 0 ) {
        cout << "Number of reconciliations: 0" << endl;
        if( gIntParams.find("pareto.mod")->second > 0 )
            rerunParetoMod = true;
        return;
    } else if( notTooBig ) {
        cout << "Number of reconciliations: " << numberSolutions <<endl;
//cout << "Count consistent: " << graph.countConsistent() << endl;
        printReconciliations( numberSolutions, ext,
                              geneTrees, eventSupports, cladesAndTripartitions,
                              dtlMatrix, speciesTree, graph );

        // count transfer losses (experimental option)
        if( gBoolParams.find("tls")->second ) {
            double tlFrequency;
            double tlCount = graph.countTLs( tlFrequency );
            tlFrequency *= 100;
            cout << "Average TLs: " << tlCount 
                  << " (" << tlFrequency << "%)" << endl;
        }
    }

    if( gStringParams.find( "orthology.output" )->second != "none" ) {
        string pathName = gPathPrefix + 
                gStringParams.find( "orthology.output" )->second;
        graph.orthologyOutput( pathName );
    }
}

/**
 * Run MPR and output results
 */
void run( 
    bool amalgamation, ///< do amalgamation if true
    MySpeciesTree *speciesTree, ///< species tree
    vector<MyGeneTree*> geneTrees, ///< gene trees
    int maxTS,  ///< maximum time slice
    vector<int> &changedTimeSlices, ///< list of changed time slices
    int counter, ///< which gene tree this is in the gene tree file
    bool constructGraph,  ///< constuct graph if true
    bool backtrackTree,  ///< create and run backtrack gene tree
    bool printCost=true, ///< print cost 
    MyGeneTree *bootstrapGeneTree=NULL, ///< gene tree with bootstrap values
    MyGeneTree *polytomyTree=NULL ) ///< polytomy tree to expand
{
    //compute the clades and tripartitions
    bool polytomic = false; 
    // don't create polytomic tree if graph is being constructed
    if( !constructGraph && gIntParams.find("resolve.trees")->second != -1 )
        polytomic = true;
    CladesAndTripartitions *cladesAndTripartitions = 
        getCladesAndTripartitions( amalgamation, geneTrees, bootstrapGeneTree,
                                   polytomic );
    // Create the matrix and print the best cost
    double inEps;
    DTLMatrix *dtlMatrix = createAndRunMatrix( backtrackTree, speciesTree,
                                    cladesAndTripartitions, maxTS, 
                                    changedTimeSlices, inEps );
    double bestCost = dtlMatrix->getBestCost(
                gBoolParams.find("gene.origination.species.root")->second );
    if( printCost ) {
        if( bestCost == std::numeric_limits<double>::max() ) 
            cout << "Cost of a most parsimonious reconciliation: OVERFLOW" 
                 << endl;
        else
            cout << "Cost of a most parsimonious reconciliation: " 
                 << bestCost << endl;
    }
    
    string ext = "";
    if( !gBoolParams.find("amalgamate")->second && counter > 0 )
        ext = "_" + bpp::TextTools::toString( counter );

    // matrix - only print once
    string matrixFileStr = gStringParams.find("print.matrix.file")->second;
    if( matrixFileStr != "none" ) {
        dtlMatrix->printMatrixCSV( 
                (matrixFileStr+ext+".csv").c_str(),
                gStringParams.find("species.file")->second.c_str(),
                gStringParams.find("gene.file")->second.c_str() );
    }

    // print resulting gene tree 
    bool append = false; // append newick or rerooted if not the first
    if( counter > 1 )
        append = true;
    MyGeneTree *treeToBacktrack = printNewickTrees( backtrackTree, append, 
                                        polytomic, amalgamation, geneTrees, 
                                        cladesAndTripartitions, dtlMatrix );
    // Graph 
    bool rerunParetoMod = false;
    if( constructGraph ) 
        makeGraph( inEps, ext, polytomyTree, geneTrees, 
                    cladesAndTripartitions, dtlMatrix, rerunParetoMod,
                    speciesTree );

    printMemory( "Final" );

    delete cladesAndTripartitions;
    delete dtlMatrix;

    if( rerunParetoMod ) {
        // no renconciliations with pareto mod, run again without
        gIntParams.find("pareto.mod")->second = 0;
        cout << "Rerunning without pareto mod." << endl;
        run( amalgamation, speciesTree, geneTrees, maxTS, 
                    changedTimeSlices, counter, constructGraph, backtrackTree,
                    printCost, bootstrapGeneTree, polytomyTree );
    }

    if( backtrackTree ) {
        // construct graph from backtracked tree
        vector<MyGeneTree*> backtrackGeneTree;
        backtrackGeneTree.push_back( treeToBacktrack );
        MyGeneTree *geneTree = NULL;
        //if( gIntParams.find("resolve.trees")->second != -1 ) 
        //    geneTree = geneTrees[0]; 
        //    // for experimental polytomic gene tree in graph
        run( false, speciesTree, backtrackGeneTree, maxTS, 
                changedTimeSlices, 0, true, false,
                false, NULL, geneTree );
    }

    if( treeToBacktrack != NULL )
        delete treeToBacktrack;
}

/**
 * Print occurrences of clades.
 */
void printCladeSupport( 
    vector<MyGeneTree*> geneTrees ) ///< gene trees
{

    bool overflow = false;
    string errStr = "";
    CladesAndTripartitions *cladesAndTripartitions = 
          new CladesAndTripartitions( gStringParams.find("char.sep")->second[0],
                            geneTrees, gBoolParams.find("verbose")->second, 
                            overflow, errStr );
    if( errStr != "" ) {
        cerr << errStr << endl;
        exit( 1 );
    }
    if( overflow ) {
        cerr << "Too many possible gene trees due to polytomies" << endl;
        exit( 1 );
    }

    int cladeCount = cladesAndTripartitions->mClades.getCladeCount();
    vector<int> bins( 10, 0 );
    vector<int> lowBins( 10, 0 );
    int treeCount = (int) geneTrees.size();
    for( int i=0; i<cladeCount; i++ ) { 
        if( i == cladesAndTripartitions->mClades.getRootClade() )
            continue;   // skip root
        double occ = cladesAndTripartitions->getCladeOccurrences( i );
        int bin = int(occ/treeCount*10);
        if( bin == 10 )
            bin--;
        bins[bin]++;

        if( bin == 0 ) {
            int lowBin = int(occ/treeCount*100);
            if( lowBin == 10 )
                lowBin--;
            lowBins[lowBin]++;
        }
    }
   
    cout << "support: 10%, 20% ... " << endl;
    for( size_t i=0; i<bins.size(); i++ ) {
        cout << bins[i] << " ";
    }
    cout << endl;
    cout << "support: 1%, 2% ... " << endl;
    for( size_t i=0; i<lowBins.size(); i++ ) {
        cout << lowBins[i] << " ";
    }
    cout << endl;
         
    delete cladesAndTripartitions;
}

/**
 * Make a polytomic tree.
 */
void collapseTree(
        MyGeneTree *&geneTree, ///< gene tree to collapse 
        MyGeneTree *&bootstrapGeneTree ) ///< gene tree with bootstrap values
{
    if( gIntParams.find("resolve.trees")->second == 0  
        || !gBoolParams.find("use.bootstrap.weighting")->second ) 
    {
        geneTree->rootTree(); // randomly root if unrooted
    }

    // copy original tree in case 
    // gDoubleParams.find("tree.limit")->second exceeded
    MyGeneTree *geneTreeCopy = NULL;
    if( //gDoubleParams.find("tree.limit")->second != 0 
        //|| 
        gIntParams.find("degree.limit")->second != 0 ) 
    {
         geneTreeCopy = new MyGeneTree( *geneTree );
    }

    double threshold = gDoubleParams.find("collapse.threshold")->second;
    double usedThreshold = threshold;
    int maxDegree = 0;
    bool overflow;
    bool firstLoop = true;
    bool continueLoop = true;
    unsigned long long numberOfTrees = 0;
    
    
    vector<double> valuesForCollapsing;
    vector<MyGeneNode*> nodes= geneTree->getNodes();
	for (unsigned int n=0; n <nodes.size() ; n++){
		if(nodes[n]->getNumberOfSons()!=0){
			if(gIntParams.find("collapse.mode")->second==0 && nodes[n]->hasDistanceToFather())
				valuesForCollapsing.push_back(nodes[n]->getDistanceToFather()+0.00001);
			else 
				if(gIntParams.find("collapse.mode")->second==1 && nodes[n]->hasBranchProperty(bpp::TreeTools::BOOTSTRAP)){
				valuesForCollapsing.push_back(dynamic_cast<const bpp::Number<double> *> (nodes[n]->getBranchProperty(bpp::TreeTools::BOOTSTRAP))->getValue()+0.00001);          
				}   
		}
            
    }
       
    std::sort( valuesForCollapsing.begin(), valuesForCollapsing.end() );
     
    int nextTresholdToTry=valuesForCollapsing.size()-1;

    //valuesForCollapsing.push_back(valuesForCollapsing[valuesForCollapsing.size()-1]+0.00001);
        
    //std::sort( valuesForCollapsing.begin(), valuesForCollapsing.end() );

    bool foundNextTresholdToTry=false;
    for (int n=valuesForCollapsing.size()-1; n>=0 ; n--){
       if(!foundNextTresholdToTry && (valuesForCollapsing[n])<threshold){
          nextTresholdToTry=n;
          foundNextTresholdToTry=true;
      }
    }

    while( firstLoop || continueLoop ) {
    
        if( !firstLoop ) {
            //if( gDoubleParams.find("tree.limit")->second!=0 ) 
            //    cout << "number of trees (" << numberOfTrees
            //        << ") greater than limit, trying threshold "
            //        << threshold << endl;
            if( gIntParams.find("degree.limit")->second!=0 ) 
                cout << "max out degree (" << maxDegree 
                    << ") greater than limit, trying threshold "
                    << threshold << endl;
            // revert to copy
            delete geneTree;
            geneTree = new MyGeneTree( *geneTreeCopy );
        }

        // save the gene tree for bootstrap weighting
        if( gBoolParams.find("use.bootstrap.weighting")->second ) {
            // root the tree if necessary
            if( gIntParams.find("resolve.trees")->second == 1 ) {
                if( !geneTree->rootWithThreshold( threshold,
                                gIntParams.find("collapse.mode")->second ) )
                {
                    cerr << "Unrooted tree has unequal values at root" << endl;
                    exit(1);
                }
                // update bootstrap with new rooting
                if( bootstrapGeneTree != NULL ) {
                    delete bootstrapGeneTree;
                    bootstrapGeneTree = NULL;
                }
            }
            if( bootstrapGeneTree == NULL )
                bootstrapGeneTree = new MyGeneTree( *geneTree );
        }  
        
        // collapse tree to make it polytomic
        maxDegree = geneTree->collapseOnTree( threshold, 
                                gIntParams.find("collapse.mode")->second );
                       
        usedThreshold = threshold;
        numberOfTrees = CladesAndTripartitions::findPolytomies( 
                                geneTree->getRootNode(), overflow );
        continueLoop = false;
        if( gIntParams.find("degree.limit")->second != 0 )
           // || gDoubleParams.find("tree.limit")->second != 0 ) 
        {
            if( 
                //(gDoubleParams.find("tree.limit")->second != 0 
                //   && (overflow
                //    || numberOfTrees>gDoubleParams.find("tree.limit")->second)) ||
               gIntParams.find("degree.limit")->second != 0 
                   && maxDegree>gIntParams.find("degree.limit")->second) 
            {
                continueLoop = true;
            }
        //} else if( overflow ) {
        //    cerr << "Overflowed number of possible polytomic"
        //         << " trees - exiting." << endl;
        //    exit(1);
        }

		if(threshold==0){
			continueLoop = false;
			
		}
		else{		

        	threshold = valuesForCollapsing[nextTresholdToTry] ;

			if (!firstLoop)
				nextTresholdToTry--;
		}		
		
        if(nextTresholdToTry ==-2 && continueLoop){
            cout << "Impossible to get a max out degree of " << maxDegree << " for this tree. Collapsing the branches with the smallest threshold (";
            cout << valuesForCollapsing[0] << ") makes the tree reaches the chosen max out degree!"<< endl;
            cout << "The gene tree is thus unchanged\n";
            threshold=0;
            //exit(-1);
        }
        firstLoop = false;
    }

    if( gBoolParams.find("verbose")->second ) {
        cout << "collapsed gene tree using a threshold of " << usedThreshold 
             << endl;
        if( numberOfTrees == 0 )
            cout << "   number of trees overflowed, degree.limit is too high" << endl;
        else
            cout << "   number of trees = " << numberOfTrees << endl;
        cout << "max out degree = " << maxDegree << endl;
    }

    if( geneTreeCopy != NULL )
        delete geneTreeCopy;
}

/**
 * Read the species tree and do basic checks.
 */
MySpeciesTree *getSpeciesTree()
{
    // read species tree
    string errString = "";
    MySpeciesTree* speciesTree = MySpeciesTree::readMySpeciesTree( 
            gStringParams.find("species.file")->second.c_str(),
            errString, 
            gBoolParams.find("dates.as.bootstraps")->second );
                        
    if( errString != "" || speciesTree == NULL ) {
        cerr << "Error reading species tree: " << errString << endl;
        exit(1);
    }

    //to have random dates
    //TreeTools::computeBranchLengthsGrafen( *speciesTree, 1, true);
    //TreeTools::convertToClockTree( *speciesTree, 
    //                          speciesTree->getRootNode()->getId());
    //Newick * print = new Newick(false,false);
    //print->write( *speciesTree, cout);		
    
    // species tree must be binary
    int nonBinaryCount = speciesTree->getBinaryCount();
    if( nonBinaryCount != 0 ) {
        cerr << "ERROR: Species tree is not binary." << endl;
        cerr << nonBinaryCount << " nodes are not binary." << endl;
        exit(1);
    }

    // species tree must have unique leaves
    string dupName;
    if( !speciesTree->uniqueLeaves( dupName ) ) {
        cerr << "ERROR: Species tree has duplicate leaf names: <" 
             << dupName << ">" << endl;
        exit(1);
    }

    return speciesTree;
}


/** 
 * Process species tree: trimming, costs, subdivision, date changes
 *
 * @return Maximum time slice.
 */
int processSpeciesTree(
        MySpeciesTree *speciesTree, ///< the species tree
        vector<MyGeneTree*> geneTrees, ///< the gene tree
        vector<int> &changedTimeSlices ) ///< return changed time slices
{
    // parse data changing parameters
    vector< pair<int, int> > dateMap;
    string errStr = "empty error";

    // read variable costs
    if( gStringParams.find("costs.file")->second != "none" ) {
        if( !speciesTree->assignCosts( 
                    gStringParams.find("costs.file")->second ) )
            exit(1);
        gFixedCosts = false;
    }

    // to add the outGroup to transfer from the dead 
    // create a sibling of the root and a new root
    if( gBoolParams.find("compute.TD")->second )  
        speciesTree->addAlphaForDeadTransfer( 
                gBoolParams.find("dates.as.bootstraps")->second,
                gDoubleParams.find("HGT.cost")->second, 
                gDoubleParams.find("loss.cost")->second );
   
    // must be done before subdivision
    speciesTree->compute_RealPostOrder();  

    // trim tree data
    boost::unordered_map<string, int> taxaNamesGenes;
    if( gBoolParams.find("trim.species.tree")->second ) {
        // create map of all taxa name in the gene trees
        BOOST_FOREACH( MyGeneTree *geneTree, geneTrees ) {
            vector<string> leaves = geneTree->getLeavesNames();
            BOOST_FOREACH( string leafName, leaves) {
                size_t pos = leafName.find( 
                            gStringParams.find("char.sep")->second[0] ); 
                string taxaName = leafName.substr(0,pos);
                taxaNamesGenes.insert( make_pair(taxaName,1) );
            }
        }

        // moved undated/partially dated tree trimming from here
    }

    // subdivide tree if dated
    if( gIntParams.find("dated")->second == 2 ) {
        bool success = speciesTree->computeSubdivision( 
                        gBoolParams.find("dates.as.bootstraps")->second, 
                        gBoolParams.find("ultrametric.only")->second,
                        errStr ); 
        if( !success ) {
            cerr << "ERROR: " << errStr << endl;
            exit(1);
        }
        if( gBoolParams.find("verbose")->second ) {

            bpp::ApplicationTools::displayTime(
                    "Computing the dated subdivision done:");
            cout << "leaves and nodes: " << speciesTree->getNumberOfLeaves()
                 << " " << speciesTree->getNumberOfNodes() << endl;
        }

        // dated tree trimming
        if( gBoolParams.find("trim.species.tree")->second ) {
            if( !gBoolParams.find("compute.TD")->second ) {
                cerr << "Trimming with dated species tree requires "
                    << "transfer from the dead." << endl;
                // The dated tree can have transfers at each timeslice,
                // which is a node in the original species tree. Trimming
                // the species tree removes these timeslices and therefore
                // possiblities to transfer at these times. To compensate,
                // transfer from the dead can be used.
                exit(1);
            }

            // trim species tree
            if( !speciesTree->trimTree( taxaNamesGenes, 
                gBoolParams.find("verbose")->second ) ) 
            {
                cerr << "Species tree has none of the gene taxa." << endl;
                exit(1);
            }
        }
    } else {
        if( gBoolParams.find("trim.species.tree")->second ) {
            // do undated or partially dated trimming here
            if( gIntParams.find("dated")->second != 2 ) {
                // trim species tree
                if( !speciesTree->trimTree( taxaNamesGenes, 
                            gBoolParams.find("verbose")->second ) ) 
                {
                    cerr << "Species tree has none of the gene taxa." << endl;
                    exit(1);
                }
            }
        }
        speciesTree->assignNoSubdivisionTimeSlices();
    }

    // assign ids (correspondance)
    speciesTree->assignPostOrderIds();

    int maxTS = 0;
    if( gDoubleParams.find("ils.cutoff")->second > 0 ) {
        speciesTree->computeSpeciesCladesAndSplits( 
            gDoubleParams.find("ils.cutoff")->second,
            gIntParams.find("ils.max.cluster.size")->second,
            gBoolParams.find("verbose")->second );
        maxTS = speciesTree->getRootNode()->getInfos().timeSlice;
    } else {
        // create vector of nodes by time slices
        maxTS = speciesTree->setVectorTimeSlices();
    }

//speciesTree->printIds();
//speciesTree->printTreeInfo(); // debugging
    if( !gFixedCosts ) {
        // check input date costs
        vector<MySpeciesNode*> nodes = speciesTree->getNodes();
        BOOST_FOREACH( MySpeciesNode* node, nodes ) {
            if( node->getInfos().isAlpha ) {
                if( node->getInfos().hgtCost == 0
                    || node->getInfos().lossCost == 0 ) 
                {
                    cout << "BAD ALPHA COST=========" << endl;
                    exit(1);
                }
            } else if( !node->hasFather() ) {
                // root doesn't matter
            } else if( node->getInfos().duplicationCost == 0 
                    || node->getInfos().hgtCost == 0
                    || node->getInfos().lossCost == 0 ) 
            {
                cout << "BAD COST=========" 
                     << node->getId() << endl;
                exit(1);
            }
        }
    }

    return maxTS;
}

/**
 * Check if gene tree has unique leaves that correspond to species names.
 * Collapse tree if requested or make the gene binary and rooted if it not.
 */
bool checkGeneTree(
        int counter, ///< counter for gene tree
        MyGeneTree *&geneTree, ///< gene tree to check
        MyGeneTree *&bootstrapGeneTree, ///< gene tree with bootstrap values
        boost::unordered_map<string, int> &taxaNamesSpecies,  ///< specie names
        boost::unordered_map<string, string> &mapNames //map for gene and species name 
        ) 
{
    if( !geneTree->restrictTreeToASetOfTaxa( 
                taxaNamesSpecies, mapNames,
                gStringParams.find("char.sep")->second[0], 
                gBoolParams.find("verbose")->second) ) 
    {
        cerr << "ERROR: Gene tree " << counter
             << " has no valid leaves." << endl;
        exit(1);
    }
    
    // needed by the algoritm, otherwise it is seg fault!!!	
    if( geneTree->getNumberOfLeaves()<3 &&  (gBoolParams.find("amalgamate")->second || gIntParams.find("resolve.trees")->second >-1)) {
        cerr << "ERROR: Gene tree " << counter
             << " has only " << geneTree->getNumberOfLeaves();
        //exit(1);
        if(gBoolParams.find("amalgamate")->second){      
            cerr << " leaves. At least three are required to amalgamate trees." << endl;
			exit(1);
        }
     	else{
        	cerr << " leaves. At least three are required to resolve trees." << endl;
        	return false;

        }      
             
    }

    // gene trees must have unique leaves
    string dupName;
    if( !geneTree->uniqueLeaves( dupName ) ) {
        cerr << "ERROR: A gene tree " << counter
            << " has duplicate leaf names: <" 
             << dupName << ">" << endl;
        exit(1);
    }

    if( gIntParams.find("resolve.trees")->second != -1 ) {
        if( gIntParams.find("collapse.mode")->second != -1 ) 
            collapseTree( geneTree, bootstrapGeneTree );
        else if( gBoolParams.find("use.bootstrap.weighting")->second )
            bootstrapGeneTree = new MyGeneTree( *geneTree );
    } else {
        // Randomly root tree if it is unrooted (3 root sons)
        geneTree->rootTree();

        // randomly make tree binary if necessary.
        if( !geneTree->isBinary() ) {
            geneTree->makeBinary();
            //if( gBoolParams.find("verbose")->second ) 
            //    cout << "Randomly made gene tree " << counter
            //         << " binary." << endl;
        }
    }

    return true;
}

/////////////////////////////////////////////////
// Main 
/////////////////////////////////////////////////

int main(int args, char ** argv)
{
	if(args == 1)
	{
		help();
		exit(0);
	}
	
	try {
        // fill global parameter variables
        map<string, string> params = 
            bpp::AttributesTools::parseOptions(args, argv);
        readParameters( params );

        if( gBoolParams.find("verbose")->second ) {
/*
	        cout << "**********************************************"
                    "********************" << endl;
	        cout << "*                 ecceTERA, version " 
                 << version << "          *" << endl;
	        cout << "* Authors: C. Scornavacca, E. Jacox         "
                    "Created     11/12/13 *" << endl;
	        cout << "*                                           "
                    "Last Modif. 01/10/15 *" << endl;
	        cout << "**********************************************"
                    "********************" << endl;
	        cout << endl;
*/
	        cout << "ecceTERA, version " << version << endl;
        }

        if( gBoolParams.find("verbose")->second ) 
            bpp::ApplicationTools::startTimer();
            
		boost::unordered_map<string, string> mapNames; //map for gene and species name 

		
		if( gStringParams.find("gene.mapping.file")->second != "none" ) {
			const char *mappingFileName  = gStringParams.find("gene.mapping.file")->second.c_str();

			ifstream fileStream( mappingFileName, ios::in );  
			if( !fileStream) { 
				cerr << "Failed to open file: "
					<< gStringParams.find("gene.mapping.file")->second << endl;
				exit(1);
			}

			string description;
	
			while( !fileStream.eof() ) {
				string temp;
				// Copy current line to temporary string
				getline(fileStream, temp, '\n');  
				if( temp == "" ) // blank line
					continue;

				boost::char_separator<char> sep(" ");
				Tokenizer tok( temp, sep );
				Tokenizer::iterator iter=tok.begin();

				string geneName =  *iter;
				iter++;
				if( iter==tok.end() ) {
					cerr << "Invalid line, missing species name in "
						<< gStringParams.find("gene.mapping.file")->second << endl;
					exit( 1 );
				}
				string speciesName =  *iter;
	
				mapNames.insert( make_pair(geneName,speciesName) );
					
			}

		}

        // read the species tree
        MySpeciesTree *speciesTree = getSpeciesTree();
        
        if( speciesTree->getNumberOfLeaves()<=1) {
            cerr << "Error: species tree has less than 2 leaves\n";
            exit(1);
        }		
        


        // map of taxaNames required for the function restrictTreeToASetOfTaxa
        boost::unordered_map<string, int> taxaNamesSpecies;
        vector<string> leaves = speciesTree->getLeavesNames();
        BOOST_FOREACH( string leaf, leaves) 
            taxaNamesSpecies.insert(make_pair(leaf,1));

        //read the gene trees
        string errString = "";
        vector<MyGeneTree*> geneTrees;
        if( !gBoolParams.find("ale")->second ) {
            bool readBootstrap = false;
            if( gIntParams.find("collapse.mode")->second == 1 
                || gBoolParams.find("use.bootstrap.weighting")->second ) 
            {
                readBootstrap = true;
            }
            geneTrees = MyGeneTree::readMyGeneTrees( 
                    gStringParams.find("gene.file")->second.c_str(), errString,
                    readBootstrap );
                                
            if( errString != "" ) {
                cerr << "Error reading gene trees: " << errString << endl;
                exit(1);
            }
            if( gBoolParams.find("verbose")->second ) 
                cout << geneTrees.size() << " gene trees" << endl;
        }

        // process the species tree
        vector<int> changedTimeSlices;
        int maxTS = processSpeciesTree( speciesTree, geneTrees, 
                                        changedTimeSlices );
        // Gene trees loop
        int counter = 0;
        for( size_t i=0; i<geneTrees.size(); i++ ) {
            counter++;

            MyGeneTree *bootstrapGeneTree = NULL;
            if( !checkGeneTree( counter, geneTrees[i], bootstrapGeneTree,
                                taxaNamesSpecies, mapNames ) )
                continue; //ship these gene trees 
            if(geneTrees[i]->getNumberOfLeaves()<=1){
                    cout << "Gene tree number " << i+1 << ", once restricted to the species present in the species "
                    << "tree, has only one leaf, skipped\n";
                continue;
            }
            
            // do each tree individually if this isn't an amalgamation
            if( !gBoolParams.find("amalgamate")->second ) {
                vector<MyGeneTree*> singleGeneTree;
                singleGeneTree.push_back( geneTrees[i] );
                if( gBoolParams.find("verbose")->second )
                    cout << "Gene tree " << counter << endl;

                bool backtrackTree = false; // run backtrack tree
                bool constructGraph = gBoolParams.find("print.info")->second;
                if( gIntParams.find("resolve.trees")->second  != -1
                    && gBoolParams.find("print.info")->second ) 
                {
                    backtrackTree = true;
                    constructGraph = false;  // not this time, after
                }

                // run the calculations
                if( geneTrees.size() == 1 )
                    counter = 0; // don't add extension
                run( false, speciesTree, singleGeneTree, maxTS, 
                        changedTimeSlices, counter, constructGraph, 
                        backtrackTree, true, bootstrapGeneTree );
                if( gBoolParams.find("verbose")->second )
                    cout << endl; // add a gap in output
            }

            if( bootstrapGeneTree != NULL )
                delete bootstrapGeneTree;
        }

        if( !gBoolParams.find("ale")->second && geneTrees.size() == 0 ) {
            cout << "No gene trees found." << endl;
            exit(1);
        }

        if( gBoolParams.find("verbose")->second ) 
            bpp::ApplicationTools::displayTime("Reading the gene trees done:");

        if( gBoolParams.find("print.support")->second ) {
            // print clade support and exit 
            printCladeSupport( geneTrees );
            return 0;
        } 

    
        // run the amalgamation
        if( gBoolParams.find("amalgamate")->second ) {
            run( true, speciesTree, geneTrees, maxTS, 
                    changedTimeSlices, 0, false, 
                    gBoolParams.find("print.info")->second );
        }

        // print a species file with post order ids 
        if( gBoolParams.find("print.newick")->second ) {
            string pathName = gPathPrefix + 
                gStringParams.find("print.newick.species.tree.file")->second;
            MySpeciesTree *tree = speciesTree->getPostorderTree();
            tree->printNewick( pathName );
            delete tree;
        }
    
        // clean up
        BOOST_FOREACH( MyGeneTree *tree, geneTrees ) 
            delete tree;
        delete speciesTree;
    
		
        if( gBoolParams.find("verbose")->second )
		    bpp::ApplicationTools::displayTime("Done:");

	}  catch(exception & e) {
		cerr << e.what() << endl;
		exit(-1);
	}
	
	return 0;
}








