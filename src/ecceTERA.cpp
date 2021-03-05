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


#include "Utils.h"

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
        int maxTS = processSpeciesTree( speciesTree, geneTrees);
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
                        counter, constructGraph, 
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
                    0, false, 
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








