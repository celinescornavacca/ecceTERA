#ifndef XMLUTILS_H_
#define XMLUTILS_H_

/*

This file contains XML util functions

Created the: 09-02-2015
by: Wandrille Duchemin

Last modified the: 09-02-2016
by: Wandrille Duchemin

*/
#include <iostream>
#include <string>
#include <map>

using namespace std;


string getLineBaliseName(string line);
string InterpretLineXML(string line, map <string, string> * properties, string * value);



#endif