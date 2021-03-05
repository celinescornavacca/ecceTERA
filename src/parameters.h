
/**

@file
@author Celine Scornavacca

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

In this file we read and check the consistency of the parameters passed to ecceTERA by the users.
We also specified the help message.
*/

#include <Bpp/Utils/AttributesTools.h>

string version = "1.2.4";

//name, type, default, description
std::map<string,bool> gBoolParams;
std::map<string,int> gIntParams;
std::map<string,double> gDoubleParams;
std::map<string,string> gStringParams; // string, char, and path 



// Ordered for the help message.
// Parameters with an empty description are not displayed in the help.
const int gParameterCount = 61; // This must matche the size of gParameters
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
 {"recPhyloXML.reconciliation", "bool", "false", "print the reconciliations using the RecPhyloXML format"},

 // amalgamation options
 {"ale", "bool", "false", "gene.file is an ALE file"},
 {"amalgamate", "bool", "false", "try all amalgamations"}, 
 {"fix.dtl.costs", "bool", "false", 
     "when iterating, only update weight, not costs"}, 
 {"max.iterations", "int", "1", 
     "max iterations to estimate costs (amalgamation)"}, 
 {"weight.amalgamation", "double", "0", 
     "weight multipler for clade costs (amalgamation)"}, 

 // ils options
 {"ils.cost", "double", "1", ""}, //"cost of a incomplete lineage sorting"},
 {"ils.cutoff", "double", "0", ""},
     //"branch length cutoff for incomplete lineage sorting (0 is disabled)"},
 {"ils.max.cluster.size", "int", "15", ""},
     //"maximum size of an ils cluster (abort if exceeded)"},

   
 // hidden
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

 // TO DEPRECATE
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



    if( gDoubleParams.find("ils.cutoff")->second > 0  
        && gIntParams.find("dated")->second != 2 ) 
    {
        cerr << "ILS not supported with undated or partially"
             << " dated species trees" << endl;
        exit(1);
    }
        
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

