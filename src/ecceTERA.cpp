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

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Utils/AttributesTools.h>

#include "CladesAndTripartitions.h"
#include "DTLMatrix.h"
#include "DTLMatrixTriplets.h"
#include "DTLMatrixRecalc.h"
#include "DTLGraph.h"
#include <Bpp/Phyl/Io/Newick.h>
#include "MyNetwork.h"
#include "NetAlg.h"
#include "DTLMatrixNetwork.h"

			
string version = "1.2.4";

//name, type, default, description
std::map<string,bool> gBoolParams;
std::map<string,int> gIntParams;
std::map<string,double> gDoubleParams;
std::map<string,string> gStringParams; // string, char, and path 


// Ordered for the help message.
// Parameters with an empty description are not displayed in the help.
const int gParameterCount = 82; // This must matche the size of gParameters
string const gParameters[gParameterCount][4] = {
// basic options
 {"species.file", "path", "required", "species tree file (newick)" },
 {"gene.file", "path", "required", "gene trees file (newick)"},
 {"char.sep", "char", "_", "char separating gene names in gene tree files"},
 {"parameter.file", "path", "none", "a file with input parameters" },
 {"dupli.cost", "double", "2", "cost of a duplication"},
 {"HGT.cost", "double", "3", "cost of an HGT"}, 
 {"loss.cost", "double", "1", "cost of a loss"},

// advanced options
 {"check.time.consistency.all", "bool", "false", 
    "check if all reconciliations are time consistent (no subdivision only)"},
 {"check.time.consistency", "bool", "false", 
 "check if there exists a time consistent reconciliation (no subdivison only)"},
 {"collapse.mode", "int", "-1", 
     "collapse nodes below collapse.threshold to create polytomic tree:"
         " -1=no collapse 0=distances, 1=bootstrap values"},
 {"collapse.threshold", "double", "0.5", 
     "collapse trees if nodes are below threshold"},
 {"consistency.limit", "int", "10000", 
     "maximum number of reconciliations to check for consistency"},
 {"compute.T", "bool", "true", "compute transfers"},
 {"compute.TL", "bool", "true", "use transfer losses"},
 {"dates.as.bootstraps", "bool", "false", 
     "read the species node ordering directly from the bootstrap values"},
 {"dated", "int", "2", "species tree dating: 0=undated, 1=partially, 2=fully"}, 
 {"gene.mapping.file", "path", "none", 
     "file with mapping between gene names and species names, space separated"},
 {"output.dir", "string", "none", "directory for printed files"}, 
 {"output.prefix", "string", "none", 
     "A prefix to prepend to all output files."}, 
 {"orthology.output", "string", "none", "print pairs of orthologous genes" 
    " to the given file"},
 {"print.newick", "bool", "false", 
    "print resulting newick gene tree and a species tree with reconciliation"
        " ids"}, 
 {"print.newick.gene.tree.file", "string", "geneTree", 
    "file name for the printed gene tree"}, 
 {"print.newick.species.tree.file", "string", "speciesTree", 
    "file name for the printed species tree"}, 
 {"resolve.trees", "int", "-1", 
     "-1=don't resolve, 0,1=create combination clades from polytomies and amalgamate them, consider trees to be 0=rooted, 1=unrooted"},
 {"compute.TD", "bool", "true", "allow transfers from the dead"}, 
 {"trim.species.tree", "bool", "false", 
        "Remove elements of the species tree above LCA of shared taxa"
        " with the gene trees."}, 
 {"ultrametric.only", "bool", "true", 
    "return an error if the dated species tree is not ultrametric"},
 {"gene.origination.species.root", "bool", "false", 
    "force the reconciliation to use the species root"},
 {"verbose", "bool", "false", "show progress and timing"},
{"use.bootstrap.weighting", "bool", "false", "weight splits using bootstrap values. Default is true if collapse.mode is set to 1"},
{"degree.limit", "int", "12", "maximum out degree of collapsed polytomic trees (0=no limit)"},


 // graph options
 {"keep.only.canonical.reconciliations", "bool", "true", 
     "simplify the reconciliation graph"},
 {"print.graph", "bool", "false", "print the graph"}, 
 {"print.graph.file", "string", "DTLGraph", "file for graph output"}, 
 {"print.info", "bool", "false", "print graph information."}, 
 {"print.reconciliations", "int", "0", 
     "0=no reconciliations, 1=print median reconciliations,"
     " 2=print all reconciliations"}, 
 {"print.reconciliations.file", "string", "reconciliationsFile", 
     "base file name for printed reconciliations"}, 
 {"solution.limit", "int", "100000", 
      "maximum number of solutions for print.reconciliations=2"},
 {"sylvx.reconciliation", "bool", "false", "print the reconciliations using the Sylvx format"},

 // amalgamation options
 {"ale", "bool", "false", "gene.file is an ALE file"},
 {"amalgamate", "bool", "false", "try all amalgamations"}, 
 {"fix.dtl.costs", "bool", "false", 
     "when iterating, only update weight, not costs"}, 
 {"max.iterations", "int", "1", 
     "max iterations to estimate costs (amalgamation)"}, 
 {"weight.amalgamation", "double", "0", 
     "weight multipler for clade costs (amalgamation)"}, 


 // suboptimal parameters
 {"min.recs", "double", "0", 
     "epsilon increased until this number of reconciliations is reached"},
 {"min.recs.increment", "double", "1", "increments of epsilon with min.recs"},
 {"max.epsilon", "double", "30", "maximum epsilon with min.recs"},
 {"nD", "double", "0", "optional parameter for some triplets algorithms"},
 {"nL", "double", "0", "optional parameter for some triplets algorithms"},
 {"nDL", "double", "0", "optional parameter to some triplets algorithms"},
 {"pareto.mod", "int", "0", "triplets algorithm to use (1-3)"},
 {"real.epsilon", "bool", "false", 
     "epsilon is a real value, rather than a percentage if true)"},
 {"suboptimal.epsilon", "double", "0", "compute suboptimal costs for epsilon)"},

 // ils options
 {"ils.cost", "double", "1", ""}, //"cost of a incomplete lineage sorting"},
 {"ils.cutoff", "double", "0", ""},
     //"branch length cutoff for incomplete lineage sorting (0 is disabled)"},
 {"ils.max.cluster.size", "int", "15", ""},
     //"maximum size of an ils cluster (abort if exceeded)"},

//netTERA
{"input.network", "bool", "false", "The input species phylogeny is a network (in extended newick)"},
{"best.switch", "bool", "true", "Run minimum switching algorithm"},
{"min.recon", "bool", "false", "Run minimum network reconciliation algorithm"},
 {"print.newick.best.switching.file", "string", "bestSwitching", 
    "file name for the printed species switching"}, 
   
   
 // hidden
 {"weighted.support", "bool", "false", ""}, //"use weighted event support)"},
 {"print.reroot.file", "string", "none", ""},
  //"output a rerooted input gene tree based on the amalgamation to this file"},
 {"reroot.proportion", "int", "2", ""},
     //"ratio of distances of root sons in rerooted tree"},
 //{"collapse.decrement", "double", "0.1", ""},
 //    //"amount to decrement collapse.threshold of tree.limit is surpassed"},
 //{"tree.limit", "double", "0", ""},
 //    //"maximum number of possible polytomic trees (0=no limit)"},
 {"print.clade.info", "bool", "false", ""},
     //"output clade and tripartition information"},
 {"run.brute", "bool", "false", ""}, 
 //Run brute minimum switching algorithm


 // Utility options
 {"print.matrix.file", "string", "none", ""}, //"file for output matrix"}, 
 {"print.support", "bool", "false", ""}, 
    //"exit after printing clade support stats"},
 {"clades.only", "bool", "false", ""}, //"print the number of clades and exit"},
 {"internal.graph.ids", "bool", "false", ""}, 
        // "otherwise they match the reconciliation"},
 {"print.ale", "string", "none", ""}, // "output ALE file name"}


 // Experimental options
 {"skip.clades", "int", "0", ""}, //"skip low occuring clades"},
 {"subopt.support", "bool", "false", ""}, 
   //"optimal reconciliation (epsilon=0), but with supports from non-optimal"},
 {"tls", "bool", "false", ""}, 
    //"print average number of transfer losses per reconciliation"},


 // strale parameters
 {"costs.file", "path", "none", "" },
 {"other.species.file", "path", "none", "" },
 {"other.costs.file", "path", "none", "" },
 {"date.change.matrix", "path", "none", ""},
    //"matrix csv file from a previous run"},
 {"date.change.swap", "string", "", ""}, 
    //"comma-separated list of 'dates' to swap"},
 {"date.changes", "string", "", ""},
     //"comma-separated list of time slices that changed"},


 // TO DEPRECATE
 {"eNewick", "bool", "false", ""}, //"allow eNewick trees)"},
 {"hali", "bool", "false", ""}, // "output reconciliations in Hali's format"},
 {"extension.id", "string", "", ""}, // "ID added to the end of file"}
};



void help() {
    cout << "The program outputs the cost of the most parsimonious"
            " reconcilation between the species tree and the ALE"
            " corresponding to the given gene trees." << endl;
    cout << "version: " << version << endl;
    cout << "__________________________________________________________"
            "_____________________" << endl;
    cout << "___________________________________INPUT__________________"
            "_____________________" << endl;
    cout << "format: parameter.name | [type,default] <description>" << endl;

    for( size_t i=0; i<gParameterCount; i++ ) {
        if( gParameters[i][3] != "" )
            cout << gParameters[i][0] << " | [" << gParameters[i][1] << ","
                 << gParameters[i][2] << "] " << gParameters[i][3] << endl;
    }

}


// parameter variables
bool gNoEpsilon;
bool gFixedCosts = true;
string gPathPrefix = ""; // path and prefix for all output files



// change to boost::filesystem?
static bool do_mkdir( const char *path ) {
    struct stat st;
    if( stat(path, &st) != 0 ) {
        if( mkdir( path, 0777 ) != 0 && errno != EEXIST)
            return false;
    } else if( !S_ISDIR(st.st_mode) ) { // check if it is a directory
        return false;
    }

    return true;
}
/**
 ** mkpath - ensure all directories in path exist
 ** Algorithm takes the pessimistic view and works top-down to ensure
 ** each directory in path exists, rather than optimistically creating
 ** the last element and working backwards.
 **/
int mkpath( string path ) // mode_t mode )
{
    char *copypath = strdup(path.c_str());
    char *sp;
    int status = 0;
    char *pp = copypath;
    while( status == 0 && (sp = strchr(pp, '/')) != 0 ) {
        if (sp != pp) { 
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath ); //, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    //if (status == 0)
        status = do_mkdir(path.c_str()); //, mode);
    free(copypath);
    return status;
}

struct MatchPathSeparator {
        bool operator()( char ch ) const { return ch == '/'; }
};



/**
 * read input parameters
 */
void readParameters( 
        map<string,string> &params ) ///< input paramsters
{
    map<string,string>::iterator iterStr;
    bool problem = false;

    // check for unknown parameters
    map<string,bool> allParams;
    for( size_t i=0; i<gParameterCount; i++ ) 
        allParams[gParameters[i][0]] = true;
    for( iterStr = params.begin(); iterStr != params.end(); ++iterStr ) {
        map<string,bool>::iterator iterBool = allParams.find( iterStr->first );
        if( iterBool == allParams.end() ) {
            problem = true;
            cerr << "Unknown parameter: " << iterStr->first << endl;
        }
    }

    // check for a parameter file
    iterStr = params.find( "parameter.file" );	
    if( iterStr != params.end() ) {
        map<string,string> fileParams;
        bpp::AttributesTools::getAttributesMapFromFile( iterStr->second,
                                                       fileParams, "=" );
        // add file parameters if not already there
        map<string,string>::iterator fIter;
        for( fIter=fileParams.begin(); fIter!=fileParams.end(); ++fIter) 
        {
            iterStr = params.find( fIter->first );
            if( iterStr == params.end() ) 
                params[fIter->first] = fIter->second;
        }

    }
        
    bool bootstrapWeightingSet=false;
    
    for( size_t i=0; i<gParameterCount; i++ ) {
        string name = gParameters[i][0];
        string type = gParameters[i][1];
        string value;

             
        // get value 
        iterStr = params.find( name );	
        if( iterStr == params.end() ) {
            if( gParameters[i][1] == "required" ) {
                problem = true;
                cerr << name << " must be specified" << endl;
            }
            value = gParameters[i][2]; // default
        } else {
            value = iterStr->second;
            if(name=="use.bootstrap.weighting") // "use.bootstrap.weighting" has been set
                bootstrapWeightingSet=true;
        }

        if( type == "bool" ) {
            if ((value == "true") 
                || (value == "TRUE")
                || (value == "t")
                || (value == "T")
                || (value == "yes")
                || (value == "YES")
                || (value == "y")
                || (value == "Y")
                || (value == "1") ) 
            {
                gBoolParams[name] = true;
            } 
            else if ((value == "false") 
                || (value == "FALSE")
                || (value == "f")
                || (value == "F")
                || (value == "no")
                || (value == "NO")
                || (value == "n")
                || (value == "N")
                || (value == "0") ) 
            {
                gBoolParams[name] = false;
            }
            else {
                cerr << "Invalid boolean value (" << value << ") for " 
                     << name << endl;
                problem = true;
            }
        } else if( type == "string" ) {
            gStringParams[name] = value;
        } else if( type == "char" ) {
            if( value.size() != 1 ) {
                cerr << "Invalid char value (" << value << ") for " 
                     << name << endl;
                problem = true;
            } else {
                gStringParams[name] = value;
            }

        } else if( type == "int" ) {
            gIntParams[name] = bpp::TextTools::toInt( value );
        } else if( type == "double" ) {
            gDoubleParams[name] = bpp::TextTools::toDouble( value );
        } else if( type == "path" ) {
            if( value == "required" ) {
                cerr << name << " is a required parameter." << endl;
                problem = true;
            } else if( value != "none" ) {
                if( !bpp::FileTools::fileExists( value ) ) {
                    cerr << value << " does not exist." << endl;
                    problem = true;
                } else {
                    gStringParams[name] = value;
                }
            } else 
                gStringParams[name] = "none";
        } else {
            cerr << "Found type <" << type << "> for " << name << endl;
            cerr << "Is the gParameterCount variable up to date?" << endl;
            throw bpp::Exception( "readParameters: unknown type" );
        }
    }

    if( problem )
        exit(1);

    // PARAMETER CHECKS
       
    // print the newick tree if a file name is given
    if( gStringParams.find("print.newick.species.tree.file")->second 
                != "speciesTree" 
       || gStringParams.find("print.newick.gene.tree.file")->second 
            != "geneTree" ) 
    {
        if( !gBoolParams.find("print.newick")->second ) {
            cout << "Forcing print.newick=true" << endl;
            gBoolParams["print.newick"] = true;
        }
    }

    // print graph if a file name is given
    if( gStringParams.find("print.graph.file")->second != "DTLGraph" ) { 
        if( gBoolParams.find("verbose")->second
            && !gBoolParams.find("print.graph")->second ) 
        {
            cout << "Forcing print.graph=true" << endl;
            gBoolParams["print.graph"] = true;
        }
    }

    // force graph construction if orthology output
    if( gStringParams.find("orthology.output")->second != "none" ) { 
        gBoolParams["print.info"] = true;
    }

    if( gIntParams.find("print.reconciliations")->second < 0 
        || gIntParams.find("print.reconciliations")->second > 2 ) 
    {
        cout << "ERROR: Valid values for print.reconciliations are 0, 1, and 2"
             << endl;
        exit(1);
    } else if( gIntParams.find("print.reconciliations")->second > 0 ) 
        gBoolParams["print.info"] = true;

    // make graph printing it
    if( gBoolParams.find("print.graph")->second )
        gBoolParams["print.info"] = true;

    // print reconciliation if a file name is given
    if( gIntParams.find("print.reconciliations")->second == 0 
        && gStringParams.find("print.reconciliations.file")->second 
            != "reconciliationsFile" )
    {
        if( gBoolParams.find("verbose")->second
            && gIntParams.find("print.reconciliations")->second != 1 ) 
        {
            cout << "Forcing print.reconciliations=1" << endl;
            gIntParams["print.reconciliations"] = 1;
        }
        gBoolParams["print.info"] = true;
    }

    // create path prefix for output files
    if( gStringParams.find("output.prefix")->second != "none" ) 
        gPathPrefix = gStringParams.find("output.prefix")->second;
    string outputDir = gStringParams.find("output.dir")->second;
    if( outputDir != "none" ) {
        if( !mkpath( outputDir ) ) {
            cout << "ERROR: " << outputDir << " is not a directory." 
                 << endl;
            exit(1);
        }
        gPathPrefix = outputDir + "/" + gPathPrefix;
    }

    // amalgamate if reading ale file
    if( gBoolParams.find("verbose")->second
        && gBoolParams.find("ale")->second ) 
    {
        cout << "Forcing amalgamate=true" << endl;
        gBoolParams["amalgamate"] = true;
    }

    if( gIntParams.find("skip.clades")->second >= 100 ) {
        cerr << "Can't skip all clades." << endl;
        exit(1);
    }

    // date change specific options
    if( gStringParams.find("date.change.matrix")->second != "none" ) {
        // Date changing parameteres
        if( gBoolParams.find("verbose")->second ) 
            bpp::ApplicationTools::displayResult(
                    "Changing dates using matrix csv file", 
                    gStringParams.find("date.change.matrix")->second);

        // Get parameters necessary for preparing species and gene trees
        // from the matrix file.
        // Read other parameters from the matrix file later.
        bool td;
        DTLMatrixRecalc::getCSVparams( 
                gStringParams.find("date.change.matrix")->second.c_str(), td );
        gBoolParams["compute.TD"] = td;
       
        // read date function 
        if( gStringParams.find("other.species.file")->second == "none" 
            && gStringParams.find("date.change.swap")->second == "" 
            && gStringParams.find("date.changes")->second == "" ) 
        {
            cerr << "other species tree, date changes or swap "
                    "required with date.change.matrix"
                << endl;
            exit(1);
        }
    }

    if( gIntParams.find("max.iterations")->second < 1 ) {
        cerr << "max.iterations must be greater than 0" << endl;
        exit(1);
    }

    // use resolve trees if colapsing nodes
    if( gIntParams.find("collapse.mode")->second < -1 
        || gIntParams.find("collapse.mode")->second > 2 ) 
    {
        cerr << "Valid values for collapse.mode are -1, 0, and 1." << endl;
        exit(1);
    } else if( gIntParams.find("collapse.mode")->second != -1
              && gBoolParams.find("verbose")->second
               && gIntParams.find("resolve.trees")->second == -1 )
    {
        cout << "Forcing resolve.trees=0" << endl;
        gIntParams["resolve.trees"] = 0;
    }

    // dates must be in bootstrap values for partially dated trees 
    if( gIntParams.find("dated")->second == 1 
        && gBoolParams.find("verbose")->second
        && !gBoolParams.find("dates.as.bootstraps" )->second )
    {
        cout << "Forcing dates.as.bootstraps=true" << endl;
        gBoolParams["dates.as.bootstraps"] = true;
    }

    if( gBoolParams.find("check.time.consistency")->second 
        || gBoolParams.find("check.time.consistency.all")->second ) 
    {
        gBoolParams["print.info"] = true;
    }

    // no consistency check if dated
    if( gIntParams.find("dated")->second == 2 ) {
        if( gBoolParams.find("check.time.consistency")->second ) {
            gBoolParams["check.time.consistency"] = false;
            cout << "WARNING: check.time.consistency set to false"
                 << " since the tree is fully dated." << endl;
        }
        if( gBoolParams.find("check.time.consistency.all")->second ) {
            gBoolParams["check.time.consistency.all"] = false;
            cout << "WARNING: check.time.consistency.all set to false"
                 << " since the tree is fully dated." << endl;
        }
    }

    if( gDoubleParams.find("suboptimal.epsilon")->second == 0 )
        gNoEpsilon = true;
    else 
        gNoEpsilon = false;

    if( gIntParams.find("pareto.mod")->second < 0 
        || gIntParams.find("pareto.mod")->second > 3 )
    {
        cout << "ERROR: pareto.mod option valid values are 0 to 3" << endl;
        exit(1);
    }

    if( gIntParams.find("pareto.mod")->second != 1 ) {
        if( gDoubleParams.find("nD")->second != 0 
            || gDoubleParams.find("nDL")->second != 0 
            || gDoubleParams.find("nL")->second != 0 ) 
        {
            cout << "WARNING: setting pareto.mod=1 because " 
                << "it's paramaters (nD,nDL,nL) are non-zero." << endl;
            gIntParams["pareto.mod"] = 1;
        }
    }

    if( gIntParams.find("pareto.mod")->second == 1 ) {
        if( gDoubleParams.find("nD")->second < 0 
            || gDoubleParams.find("nD")->second > 1 )  
        {
            cout << "ERROR: given nD=" << gDoubleParams.find("nD")->second 
                 << ", valid values are [0,1]" << endl;
            exit(1);
        }
        if( gDoubleParams.find("nL")->second < 0 
            || gDoubleParams.find("nL")->second > 1 ) 
        {
            cout << "ERROR: given nL=" << gDoubleParams.find("nL")->second 
                 << ", valid values are [0,1]" << endl;
            exit(1);
        }

        if( gDoubleParams.find("nDL")->second < 0 
            || gDoubleParams.find("nDL")->second > 1 ) {
            cout << "ERROR: given nDL=" << gDoubleParams.find("nDL")->second 
                 << ", valid values are [0,1]" << endl;
            exit(1);
        }

        gNoEpsilon = true;
    }

    if( gIntParams.find("pareto.mod")->second > 0 
        && gDoubleParams.find("suboptimal.epsilon")->second < 0 ) 
    {
        cerr << "Epsilon cannot be negative." << endl;
        exit(1);
    }

    if( gIntParams.find("pareto.mod")->second > 0 
        && gDoubleParams.find("ils.cutoff")->second > 0 ) 
    {
        cerr << "ILS not supported with pareto.mod" << endl;
        exit(1);
    }

    if( gDoubleParams.find("ils.cutoff")->second > 0  
        && gIntParams.find("dated")->second != 2 ) 
    {
        cerr << "ILS not supported with undated or partially"
             << " dated species trees" << endl;
        exit(1);
    }

    if( gBoolParams.find("hali")->second ) 
        gBoolParams["print.info"] = true;
        
 if( gIntParams.find("collapse.mode")->second ==1 && (!gBoolParams.find("use.bootstrap.weighting")->second) &&  !bootstrapWeightingSet){
        cout << "Setting use.bootstrap.weighting=true" << endl;
        gBoolParams["use.bootstrap.weighting"] = true;   //if "use.bootstrap.weighting" is false and was not set, if collapse.mode we force it to true 
    }    


//        string WGDPath;
//        if(gConsiderWGD){
//            throw bpp::Exception( "WGD not implmented" );
//        NOT IMPLEMENTED 
//            WGDPath = bpp::ApplicationTools::getAFilePath("WGD.file", params,
//                                                         false, true );
//			bpp::ApplicationTools::displayResult("WGD.file", WGDPath);
//       	if(WGDPath == "none") throw bpp::Exception(
//            "You must provide a file containing the information about the WGD.");
//        }

}


/** 
 * read date swap, which is a comma separated list of
 * pairs of integers, separated by a -, e.g. 7-8,15-24
 *
 * @return true if no errors
 */ 
bool parseDateSwap( 
        vector<pair<int, int> > &dateMap, ///< the swap pairs in a map
        string &errStr ) ///< description of any errors
{

    // loop over changes
    string::size_type index = -1;
    while( true ) {
        string val;
        string::size_type prevIndex = index;
        string swapStr = gStringParams.find("date.change.swap")->second;
        index = swapStr.find(",", index+1 );
        if( index == string::npos ) 
            val = swapStr.substr(prevIndex+1); // last token 
        else 
            val = swapStr.substr(prevIndex+1, index-prevIndex-1);
        size_t index2 = val.find("-");
        if( index2 == string::npos ) {
            errStr = "Date function missing separator (-)";
            return false;
        }
        try {
            int ts1 = bpp::TextTools::toInt( (val.substr(0,index2)).c_str() );
            int ts2 = bpp::TextTools::toInt( (val.substr(index2+1)).c_str() );
            dateMap.push_back(make_pair(ts1,ts2));

//            dateMap.insert(make_pair(ts1,ts2));
//            dateMap.insert(make_pair(ts2,ts1));
// Not catching non-numerics (change in ComputeClades too)
// Use boost cast
//cout << "reorder: <" << ts1 << "/" << ts2 << ">" << endl;
        } catch( bpp::Exception e ) {
              errStr = "Date swap contains non-numeric time slices.";
              return false;
        }
        if( index == string::npos ) 
            break;
    }

    return true;
}

/** 
 * read date changes, which is a comma separated list of time slices
 *
 * @return true if no errors
 */ 
bool parseDateChanges( 
        vector<int> &changedTimeSlices,
        string &errStr ) ///< description of any errors
{
    // loop over changes
    string::size_type index = -1;
    while( true ) {
        string val;
        string::size_type prevIndex = index;
        string changeStr = gStringParams.find("date.changes")->second;
        index = changeStr.find(",", index+1 );
        if( index == string::npos ) 
            val = changeStr.substr(prevIndex+1); // last token 
        else 
            val = changeStr.substr(prevIndex+1, index-prevIndex-1);
        try {
            int ts = bpp::TextTools::toInt( val );
            changedTimeSlices.push_back( ts );
        } catch( bpp::Exception e ) {
              errStr = "Date changes contains non-numeric time slices.";
              return false;
        }
        if( index == string::npos ) 
            break;
    }

    return true;
}

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

        boost::char_separator<char> sep(" ");
        Tokenizer tok( temp, sep );
        Tokenizer::iterator iter=tok.begin();
        if( *iter == "" ) // blank line
            continue;
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
DTLMatrixNetwork *createAndRunMatrix(
        bool returnTree,  ///< return backtrack gene tree
        MyNetwork *speciesTree, ///< species tree
        CladesAndTripartitions *cladesAndTripartitions, ///< gene clades
        int maxTS,  ///< maximum time slice
        vector<int> &changedTimeSlices, ///< list of changed time slices
        double &inEps ) ///< suboptimal epsilon value used
{
    inEps = 0;
    DTLMatrixNetwork *dtlMatrix = NULL;

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
		dtlMatrix = new DTLMatrixNetwork( speciesTree, 
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
			false );


    printMemory( "Matrix" );

    return dtlMatrix;
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
    if( gStringParams.find("date.change.matrix")->second == "none" ) {

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

    } else {
        // Date changing algorithm, which reads existing algorithm
        DTLMatrixRecalc *dtlMatrixR = new DTLMatrixRecalc( speciesTree,
                cladesAndTripartitions, 
                gStringParams.find("date.change.matrix")->second.c_str(),
                maxTS );
        if( changedTimeSlices.size() > 0 ) {
        //cout << "=== changed time slices:";
        //BOOST_FOREACH( int ts, changedTimeSlices ) 
        //    cout << " " << ts;
        //cout << endl;
    
            dtlMatrixR->recalculateMatrix( 
                    changedTimeSlices, gBoolParams.find("verbose")->second );
           // For testing:
//            dtlMatrix->calculateMatrix( gBoolParams.find("verbose")->second, 
//                   1, false, 2, 0, changedTimeSlices[0], maxTS);
//            dtlMatrix->calculateMatrix( gBoolParams.find("verbose")->second );
        } else {
            cout << "No time slice changes were given." << endl;
        }

        dtlMatrix = dtlMatrixR;
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
    ofstream &statFile, ///< hali's stat file
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
                        gBoolParams.find("check.time.consistency")->second,
                        gIntParams.find("solution.limit")->second );
        if( printCount != numberSolutions )
            cout << "Printed " << printCount << "/" <<
                numberSolutions << " reconciliations" << endl;
    }

    // print Hali's reconciliations
    if( gBoolParams.find("hali")->second ) {
        statFile << numberSolutions << endl;
        statFile << "totalWeight:" << numberSolutions;
        statFile.close();

        MyGeneTree *geneTree = geneTrees[0];
        geneTree->assignPostOrderIds();
        vector<int> cladeToPOrd =
            cladesAndTripartitions->getPostOrderMapping();
        vector<int> pOrdToClade( cladeToPOrd.size() );
        for( size_t id_u=0; id_u<cladeToPOrd.size(); id_u++ ) 
            pOrdToClade[cladeToPOrd[id_u]] = id_u;
        bool triplets = false;
        if( gIntParams.find("pareto.mod")->second > 0 )
            triplets = true;
        graph.printReconciliationHali( *geneTree, pOrdToClade, 
                gPathPrefix, 
                gStringParams.find("extension.id")->second, 
                triplets );
    }
}

/**
 * Create graph and count reconciliations. Print them if requested.
 */
void makeGraph(
    double inEps,   ///< epsilon value, if used
    string ext,     ///< file name extension
    ofstream &statFile, ///< hali's stat file
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
        printReconciliations( numberSolutions, ext, statFile,
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

    // Hali's stat output
    ofstream statFile;
    if( gBoolParams.find("hali")->second ) {
        statFile.open( string(gPathPrefix + "statistics.txt" 
                       + gStringParams.find("extension.id")->second).c_str(), 
                        ios::out );
        statFile << "\"runID\", \"optCost\", \"dupCost\", \"tranCost\","
                    " \"lossCost\", \"nbRecs\", \"nbCanonicals\"" << endl; 
        statFile << 1 << "," << bestCost<< "," 
                 << gDoubleParams.find("dupli.cost")->second << "," 
                 << gDoubleParams.find("HGT.cost")->second<< "," 
                 << gDoubleParams.find("loss.cost")->second << ",,";
    }

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
        makeGraph( inEps, ext, statFile, polytomyTree, geneTrees, 
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
        
        vector<double> valuesForCollapsing;
        vector<MyGeneNode*> nodes= geneTree->getNodes();
        for (unsigned int n=0; n <nodes.size() ; n++){
            if(nodes[n]->getNumberOfSons()!=0){
                if(gIntParams.find("collapse.mode")->second==0 && nodes[n]->hasDistanceToFather())
                    valuesForCollapsing.push_back(nodes[n]->getDistanceToFather()+0.00001);
                else if(gIntParams.find("collapse.mode")->second==1 && nodes[n]->hasBranchProperty(bpp::TreeTools::BOOTSTRAP))
                    valuesForCollapsing.push_back(dynamic_cast<const bpp::Number<double> *> (nodes[n]->getBranchProperty(bpp::TreeTools::BOOTSTRAP))->getValue()+0.00001);             
            }
            
        }
        
        std::sort( valuesForCollapsing.begin(), valuesForCollapsing.end() );

        //valuesForCollapsing.push_back(valuesForCollapsing[valuesForCollapsing.size()-1]+0.00001);
        
        std::sort( valuesForCollapsing.begin(), valuesForCollapsing.end() );

        int nextTresholdToTry=valuesForCollapsing.size()-1;
        bool foundNextTresholdToTry=false;
        for (int n=valuesForCollapsing.size()-1; n>=0 ; n--){
            if(!foundNextTresholdToTry && (valuesForCollapsing[n])<threshold){
                nextTresholdToTry=n;
                foundNextTresholdToTry=true;
            }
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

        threshold = valuesForCollapsing[nextTresholdToTry] ;

        nextTresholdToTry--; // -= gDoubleParams.find("collapse.decrement")->second; 
        if(nextTresholdToTry ==-1 && continueLoop){
            cout << "Impossible to get a max out degree of " << maxDegree << " for this tree. Collapsing the branches with the smallest threshold (";
            cout << valuesForCollapsing[0] << ") makes the tree reaches the chosen max out degree!"<< endl;
            exit(-1);
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
 * Use another species tree to calculate changed time slices.
 */
void processOtherSpeciesTree(
        MySpeciesTree *speciesTree, ///< the species tree
        vector< pair<int, int> > &dateMap, ///< input date changes
        vector<int> &changedTimeSlices ) ///< return changed time slices
{
    //cout << "OTHER SPECIES TREE: " 
    //     << gStringParams.find("other.species.file")->second << endl;
    if( changedTimeSlices.size() > 0 ) {
        cout << "ERROR: Cannot do other species tree, already have "
                " changed time slices." << endl;
        exit(1);
    }
    if( gIntParams.find("dated")->second != 2 ) {
        cout << "ERROR: other.species.file only applies to"
                " subdivided trees" << endl;
        exit(1);
    }

    // read tree
    string errString = "";
    MySpeciesTree* otherSpeciesTree = MySpeciesTree::readMySpeciesTree(
               gStringParams.find("other.species.file")->second.c_str(),
               errString,
               gBoolParams.find("dates.as.bootstraps")->second );
    if( errString != "" || otherSpeciesTree == NULL ) {
        cerr << "Error reading other species tree: " << errString 
              << endl;
        exit(1);
    }

    if( gStringParams.find("other.costs.file")->second != "none" ) 
        otherSpeciesTree->assignCosts( 
                gStringParams.find("other.costs.file")->second );

    // transfer dead
    if( gBoolParams.find("compute.TD")->second )  
        otherSpeciesTree->addAlphaForDeadTransfer( 
                gBoolParams.find("dates.as.bootstraps")->second, 
                gDoubleParams.find("HGT.cost")->second, 
                gDoubleParams.find("loss.cost")->second );

    otherSpeciesTree->compute_RealPostOrder();  

    // subdivision
    otherSpeciesTree->computeSubdivision( dateMap, 
                     gBoolParams.find("dates.as.bootstraps")->second, 
                     gBoolParams.find("ultrametric.only")->second,
                     changedTimeSlices, errString ); 

    otherSpeciesTree->assignPostOrderIds();

    // find the changed time slices
    changedTimeSlices = speciesTree->findChangedTimeSlices( 
                            otherSpeciesTree );
    delete otherSpeciesTree;

    cout << "CHANGED TIME SLICES: ";
    BOOST_FOREACH( int i, changedTimeSlices ) 
        cout << i << " ";
    cout << endl;
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
    if( gStringParams.find("date.change.swap")->second != "" ) {
        if( !parseDateSwap( dateMap, errStr ) ) {
            cerr << "ERROR: " << errStr << endl;
            exit(1);
        }
    }
    if( gStringParams.find("date.changes")->second != "" ) {
        if( !parseDateChanges( changedTimeSlices, errStr ) ) {
            cerr << "ERROR: " << errStr << endl;
            exit(1);
        }
    }

    // for strale ordering to identify nodes to swap
    if( dateMap.size() != 0 ) 
        speciesTree->breadthFirstreNumber();

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
        bool success = speciesTree->computeSubdivision( dateMap, 
                        gBoolParams.find("dates.as.bootstraps")->second, 
                        gBoolParams.find("ultrametric.only")->second,
                        changedTimeSlices, errStr ); 
        if( !success ) {
            cerr << "ERROR: " << errStr << endl;
            exit(1);
        }
        if( gBoolParams.find("verbose")->second ) {
            BOOST_FOREACH(int ts, changedTimeSlices ) 
                cout << "dates changed: " << ts << endl;

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

    // process another species file
    if( gStringParams.find("other.species.file")->second != "none" ) 
        processOtherSpeciesTree( speciesTree, dateMap, changedTimeSlices );

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
        boost::unordered_map<string, int> &taxaNamesSpecies ) ///< specie names
{
    if( !geneTree->restrictTreeToASetOfTaxa( 
                taxaNamesSpecies, 
                gStringParams.find("char.sep")->second[0], 
                gBoolParams.find("verbose")->second ) ) 
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

		if(! gBoolParams.find("input.network")->second ){
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
									taxaNamesSpecies ) )
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
		}
		else{
			// read species tree
			string errString = "";
			MyNetwork* speciesNetwork = MyNetwork::readMyNetwork( 
					gStringParams.find("species.file")->second.c_str(),
										 errString );
			if( errString != "" || speciesNetwork == NULL ) {
				cerr << "Error reading species network: " << errString << endl;
				exit(1);
			}

			speciesNetwork->assignNetworkPostOrderIds();
					

			vector<MySpeciesNode*> allNodes = speciesNetwork->getNodes();
		
			string errStr;
		

				BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
				int sonCount = node->getNumberOfSons();
				if( sonCount > 2 ) {
					errStr = "A network node has more than two children";
					return false;
				}
				for( int i=0; i<sonCount; i++ ) {
					MySpeciesNode *son = node->getSon( i );
					if( son->hasName() && son->getName().find('#')!=std::string::npos && son->getInfos().primaryFather ==NULL){    
						son->getInfos().primaryFather= son->getFather();
					}	
					if(son->getInfos().primaryFather !=NULL){
						if( node->getId() != son->getInfos().primaryFather->getId() ) {
							son->getInfos().secondaryFather =node;
						}
					}
				}

			}
			
			
			// species tree must be binary, except for hybrid nodes (2 parents)
			if( !speciesNetwork->checkNetwork( errStr ) ) {
				cerr << "ERROR: Invalid species network: " << errStr << endl;
				exit(1);
			}

			if( !speciesNetwork->checkNetwork( errStr ) ) {
				cerr << "ERROR: second check failed. Invalid species network " << errStr << endl;
				exit(1);
			}
	
	
			
			allNodes = speciesNetwork->getNodes();
			
					BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
				int sonCount = node->getNumberOfSons();
			for( int i=0; i<sonCount; i++ ) {
					MySpeciesNode *son = node->getSon( i );
					if(son->getInfos().primaryFather !=NULL){
						if( node->getId() != son->getInfos().primaryFather->getId() ) {
							//cout << "reticulation 2\n";
						}
					}
				}

			}
			
			

			// map of taxaNames required for the function restrictTreeToASetOfTaxa
			boost::unordered_map<string, int> taxaNamesSpecies;
			vector<string> leafNames = speciesNetwork->getLeavesNames();
			BOOST_FOREACH( string leafName, leafNames ) {
				size_t pos = leafName.find( "#" ); 
				leafName = leafName.substr(0,pos);
				taxaNamesSpecies.insert( make_pair(leafName,1) );
			}


			//reading the gene trees 
			vector<MyGeneTree*> geneTrees = MyGeneTree::readMyGeneTrees(
							gStringParams.find("gene.file")->second.c_str(),
							errString );
			if( errString != "" ) {
				cerr << "Error reading gene trees: " << errString << endl;
				exit(1);
			}
			if( geneTrees.size() < 1 ) {
				cerr << "No gene trees found" << endl;
				exit(1);
			}
			
			int maxTS =0;
			int g=0;
			
			BOOST_FOREACH( MyGeneTree *geneTree, geneTrees ) {
				g++;
				cout << "Gene " << g << endl;


				
				if(g==1){
					maxTS =speciesNetwork->assignNetworkPostOrderIdsDFS() -1;
	        		speciesNetwork->compute_RealPostOrder();  //needed?
	        		//allNodes = speciesNetwork->getNodes();
					//BOOST_FOREACH( MySpeciesNode *node, allNodes ) { 
					// 	cout << node->getInfos().realPostOrder << " " << node->getId() << endl; 
					// }
				}


				//////////////////////////////////////////////////////////
				// Gene trees
				//////////////////////////////////////////////////////////
		
				if( !geneTree->restrictTreeToASetOfTaxa( 
					taxaNamesSpecies, gStringParams.find("char.sep")->second[0],
					gBoolParams.find("verbose")->second )  )
				{
					cerr << "ERROR: Gene tree has no valid leaves." << endl;
					exit(1);
				}

		// Is there a limit?
				// needed by the algoritm!!!	
				if( geneTree->getNumberOfLeaves()<3 ) {
					cerr << "ERROR: Gene tree has only " 
						 << geneTree->getNumberOfLeaves()
						 << " leaves. At least three are required." << endl;
					exit(1);
				}

				if( !geneTree->isBinary() ) {
					if( gBoolParams.find("force.binary")->second )
						geneTree->makeBinary();
					else {
						cerr << "ERROR: Gene tree is not binary." << endl;

						vector<MyGeneNode*> nodes = geneTree->getNodes();		
						for( size_t i=0; i<nodes.size(); i++ ) {
							if( !nodes[i]->isLeaf() && nodes[i]->getNumberOfSons() != 2 ) { 
								cout << nodes[i]->getNumberOfSons() << endl;
							}    
						}

						exit(1);
					}
				}

				// gene trees must have unique leaves
				string dupName;
				if( !geneTree->uniqueLeaves( dupName ) ) {
					cerr << "ERROR: A gene tree has duplicate leaf names: <" 
						 << dupName << ">" << endl;
					exit(1);
				}

				NetAlg netAlg( speciesNetwork, geneTree,
							   gStringParams.find("char.sep")->second[0],
							   gDoubleParams.find("dupli.cost")->second,
							   gDoubleParams.find("loss.cost")->second,
							   gDoubleParams.find("HGT.cost")->second );
					   

		
		
		    double cost =-1;
			double costMinSwitch=-1;
			MyNetwork * bestSwitching = NULL;
		    	
				if( gBoolParams.find("run.brute")->second ) {
					cout << "======== BRUTE ALG====== " << endl;
					int numLosses = 0;
					int numDupli = 0;
					int numTransfers= 0;

					cost = netAlg.runBrute( numLosses, numDupli, numTransfers);        
								
					cout << numLosses << " losses, " << numDupli
						 << " duplications and " << numTransfers << " transfers" << endl;
					cout << "cost = " << cost << endl;
				}

				if( gBoolParams.find("best.switch")->second ) {
					cout << "===============SWITCH ALG============" << endl;

					int numLosses = 0;
					int numDupli = 0;
					int numTransfers= 0;
				
					vector<std::pair <int,int> > edgesBestSwitchings;
							
				        costMinSwitch = netAlg.runMinSwitch( numLosses, numDupli , numTransfers, edgesBestSwitchings );

					/*cout << "=============BEST===========\n";  

					for (int bs=0; bs < edgesBestSwitchings.size(); bs++){
						cout << edgesBestSwitchings[bs].first << " "  << edgesBestSwitchings[bs].second << endl;
					}
					cout << "=============++++===========\n";  
					*/
				
				
				
					//MyNetwork * bestSwitching = speciesNetwork ;
			 
					bestSwitching = new MyNetwork( * speciesNetwork) ;


					vector<MySpeciesNode*> allNodesSwitching = bestSwitching->getNodes();
				

				
					//bestSwitching->setCorrespondance(allNodesSwitching);


					BOOST_FOREACH( MySpeciesNode *node, allNodesSwitching ) {
						if(node->getInfos().secondaryFather !=NULL){
							std::pair <int,int> edge = std::make_pair(node->getInfos().secondaryFather->getId(),node->getId());
							//cout << "reticulation " << node->getId() << " " << node->getInfos().primaryFather->getId() <<  " " << node->getInfos().secondaryFather->getId() << "\n";	

							//if secondary father is kept, we delete the primary edge
							if ( std::find(edgesBestSwitchings.begin(), edgesBestSwitchings.end(), edge) != edgesBestSwitchings.end() ){
								//cout << "qui\n";
								MySpeciesNode* father = node->getInfos().primaryFather;
								MySpeciesNode* fatherNew = node->getInfos().secondaryFather;
								MySpeciesNode* son = bestSwitching->getNodeById(edge.second);
								//cout << " " << father->getId() << " " << son->getId() << "\n";	
								if(father->getId() == son->getFather()->getId()){
									father->removeSon(son);
									son->removeFather();
									//fatherNew->addSon(son);
									son->setFather(fatherNew);
								}
								//cout << son->getFather()->getId();
								//cout << father->getNumberOfSons() << endl;
							}	
							else{ //otherwise (ie if primary father is kept or none is kept -- in this case we keep the primary one as default )
								//cout << "qua\n";
								MySpeciesNode* father = node->getInfos().secondaryFather;
								MySpeciesNode* fatherNew = node->getInfos().primaryFather;
								MySpeciesNode* son = bestSwitching->getNodeById(edge.second);
								//cout << " " << father->getId() << " " << son->getId() << "\n";	
								if(father->getId() == son->getFather()->getId()){
									father->removeSon(son);
									son->removeFather();
									//fatherNew->addSon(son);
									son->setFather(fatherNew);
									//son->removeFather();							
									//son->setFather(fatherNew);
								}
								son->getInfos().secondaryFather=NULL;
								//cout << son->getFather()->getId();
								//son->getInfos().primaryFather=fatherNew;
							}
						
				
						
					
						
						}
					}
															
						//to finish, if we get the switching we should have the same min rec ??? think about it 
						bestSwitching->makeBinarySwitching(bestSwitching->getRootNode(), true);


						//int maxTS2 = 
						bestSwitching->assignNetworkPostOrderIdsDFS();
	        			bestSwitching->compute_RealPostOrder();  //needed?
	        			
	        			vector<std::pair <int,int> > edgesBestSwitchingsSwitching;
							
						numLosses =0 ;
						numDupli =0 ;
						numTransfers =0 ;
        				
        				NetAlg netAlgSwitching( bestSwitching, geneTree,
							   gStringParams.find("char.sep")->second[0],
							   gDoubleParams.find("dupli.cost")->second,
							   gDoubleParams.find("loss.cost")->second,
							   gDoubleParams.find("HGT.cost")->second );
					   
					   	double bestCostSwitching = netAlg.runMinSwitch( numLosses, numDupli , numTransfers, edgesBestSwitchingsSwitching );
						//cout << "bestCostSwitching = " << bestCostSwitching << endl; 


        				/*bestSwitching->printNetwork();


                    	vector<MySpeciesNode*>  allNodesSwitching = bestSwitching->getNodes();
				
						BOOST_FOREACH( MySpeciesNode *node, allNodesSwitching ) {
				
						if(node->hasFather())
							cout << node->getId() << " " << node->getFather()->getId() << endl;
					
						}
					
						//bestSwitching->makeBinarySwitching(bestSwitching->getRootNode(), true);

						vector<int> changedTimeSlicesSwitching;

						DTLMatrixNetwork *dtlMatrixSwitching = createAndRunMatrix( false, bestSwitching, cladesAndTripartitions, maxTS2, changedTimeSlicesSwitching, inEps );

						double bestCostSwitching = dtlMatrixSwitching->getBestCost();	
						cout << "bestCostSwitching = " << bestCostSwitching << endl; 
						*/

					

						bestSwitching->makeBinarySwitching();
						
						string pathName = gPathPrefix + 
						gStringParams.find("print.newick.best.switching.file")->second;
						bestSwitching->printNewick( pathName );
						cout << "Best switching printed in file: " << pathName << endl;
						
						if(bestCostSwitching!=costMinSwitch){
							cout << "ERRROR WRONG SWITCHING COST " << bestCostSwitching << " " << costMinSwitch << endl;
							return 0;
						}
						

                                        				
				        //bestSwitching->makeBinarySwitching(bestSwitching->getRootNode(), true);
                                        //bestSwitching->printNewick( "try" );
                                        

												
					/*allNodesSwitching = bestSwitching->getNodes();
				
					BOOST_FOREACH( MySpeciesNode *node, allNodesSwitching ) {
				
						if(node->getNumberOfSons()==1)
							cout << "problem binary A " << node->getId() << endl; // << " " << node->hasFather() << " " << node->getNumberOfSons() << endl;
						else if(node->getNumberOfSons()==0 && (!( node->hasName()) || node->getName().find('#')!=std::string::npos ))
							cout << "problem binary B " << node->getId() << endl;
					
					}*/
				
					cout << numLosses << " losses, " << numDupli
						 << " duplications and " << numTransfers << " transfers" << endl;
					cout << "Cost of a most parsimonious switching = " << costMinSwitch << endl;
					
					if(cost != -1 && costMinSwitch!=cost){
						cout << "ERRROR BRUTE AND SWITCHING COSTS ARE DIFFERENT" << costMinSwitch << " " << cost << endl;
						return 0;
					}
				}
		
				if( gBoolParams.find("min.recon")->second ) {
					cout << "=============== MIN RECON ALG============" << endl;
				
					allNodes = speciesNetwork->getNodes();
			
					//BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
					//if(node->getInfos().secondaryFather !=NULL){
					//		cout << "reticulation with number of son = " << node->getNumberOfSons() << endl;
					//	}
					//}


				
					CladesAndTripartitions *cladesAndTripartitions= new CladesAndTripartitions(gStringParams.find("char.sep")->second[0], *(geneTree) );
				
				
					// Create the matrix and print the best cost
					double inEps;
				
				
					vector<int> changedTimeSlices;
					DTLMatrix*dtlMatrix = createAndRunMatrix( false, speciesNetwork,
											cladesAndTripartitions, maxTS, 
											changedTimeSlices, inEps );
										
				
					double bestCost = dtlMatrix->getBestCost(
                            gBoolParams.find("gene.origination.species.root")
                                        ->second );
					if( bestCost == std::numeric_limits<double>::max() ) 
						cout << "Cost of a most parsimonious reconciliation: OVERFLOW" 
						 << endl;
					else
						cout << "Cost of a most parsimonious reconciliation: " 
						 << bestCost << endl;
					

					
					double costMinRec = netAlg.runMinRecon();
                    //netAlg.printRecon();
					
					if( gBoolParams.find("print.info")->second ) {
                        ofstream statFile;
                        vector<MyGeneTree*> emptyGeneTrees;
                        makeGraph( inEps, "", statFile, NULL, 
                                   emptyGeneTrees, cladesAndTripartitions, 
                                   dtlMatrix, false, speciesNetwork );
                    }
                    delete dtlMatrix;
                    

					//CS to test the new implementation, but I am not sure because we do not count losses exactly the same way....
					//cout << "cost old costMinRec= " << costMinRec << endl;
				    //if(gDoubleParams.find("HGT.cost")->second==0 && (bestCost!=costMinRec)){
                    //     cout << "ERRROR OLD AND NEW MIN COSTS ARE DIFFERENT " << bestCost << " " << costMinRec << endl;
                    //     return 0;
                    //}
        
	                                		

			

				}
				if(bestSwitching!=NULL)
					delete bestSwitching;	

			



				//BOOST_FOREACH( MyGeneTree *tree, geneTrees ) 
					delete geneTree;
			}
			delete speciesNetwork;    
			 
			
		   
			
		}

        if( gBoolParams.find("verbose")->second )
		    bpp::ApplicationTools::displayTime("Done:");

	}  catch(exception & e) {
		cerr << e.what() << endl;
		exit(-1);
	}
	
	return 0;
}








