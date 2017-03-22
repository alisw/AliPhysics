///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2011
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 228                       $: revision of last commit
// $Author:: butter             $: author of last commit
// $Date:: 2016-01-18 18:10:17 +#$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include "../include/inputParser.h"
#include <fstream>
#include <cstdlib>
#include <algorithm>

inputParser::inputParser()
{
}

inputParser::~inputParser()
{
}

int inputParser::parseFile(std::string filename)
{

    std::ifstream infile(filename.c_str());
    if ((!infile) || (!infile.good())) 
    {
      return -1;
    }
    
    int lineSize = 256;
    char tmp[lineSize];
    int nParameters = 0;
    while (!infile.getline(tmp, lineSize).eof())
    {

        std::string line(tmp);
	nParameters += parseString(line);
    }

    infile.close();
    return nParameters;

}
int inputParser::parseString(std::string str)
{

    std::string word;
    std::string name;
    std::string val;

    std::map<std::string, _parameter<int> >::iterator intIt;
    std::map<std::string, _parameter<unsigned int> >::iterator uIntIt;
    std::map<std::string, _parameter<float> >::iterator floatIt;
    std::map<std::string, _parameter<double> >::iterator doubleIt;
    std::map<std::string, _parameter<bool> >::iterator boolIt;
    std::map<std::string, _parameter<std::string> >::iterator stringIt;
    
    // Check if there is commented out stuff...
    size_t pos = str.find_first_of("#");

    // Cut the comment out of the str
    if (pos != str.npos) str.erase(pos, str.find_first_of('\n'));

    // Find the required equal sign and split the string into name and value
    size_t eqPos = str.find("=");
    std::string whitespaces (" \t\f\v\n\r");

    name = "";
    val = "";

    if (eqPos != str.npos)
    {
        name = str.substr(0, eqPos);
        name.erase(name.find_last_not_of(whitespaces)+1);
        val = str.substr(eqPos+1);
        val.erase(0, val.find_first_not_of(whitespaces));
        val.erase(val.find_last_not_of(whitespaces)+1);
    }

    if (name.length() > 0 && val.length() > 0)
    {
        intIt = _intParameters.find(name);
        if (intIt != _intParameters.end())
        {
            intIt->second._found = true;
            *(intIt->second._val) = atoi(val.c_str());
	    return true;
        }
        uIntIt = _uintParameters.find(name);
        if (uIntIt != _uintParameters.end())
        {
            uIntIt->second._found = true;
            *(uIntIt->second._val) = atoi(val.c_str());
	    return true;
        }
        floatIt = _floatParameters.find(name);
        if (floatIt != _floatParameters.end())
        {
            floatIt->second._found = true;
            *(floatIt->second._val) = atof(val.c_str());
	    return true;
        }
        doubleIt = _doubleParameters.find(name);
        if (doubleIt != _doubleParameters.end())
        {
            doubleIt->second._found = true;
            *(doubleIt->second._val) = atof(val.c_str());
	    return true;
        }
        boolIt = _boolParameters.find(name);
        if (boolIt != _boolParameters.end())
        {
            boolIt->second._found = true;
            *(boolIt->second._val) = atoi(val.c_str());
	    return true;
        }
        stringIt = _stringParameters.find(name);
        if (stringIt != _stringParameters.end())
        {
            stringIt->second._found = true;
            *(stringIt->second._val) = val;
	    return true;
        }

    }
    return false;
}
void inputParser::addIntParameter(std::string name, int *var, bool required)
{
    _parameter<int> par(name, var, required);
    _intParameters.insert(std::pair<std::string, _parameter<int> >(name, par));
}

void inputParser::addUintParameter(std::string name, unsigned int *var, bool required)
{
    _parameter<unsigned int> par(name, var, required);
    _uintParameters.insert(std::pair<std::string, _parameter<unsigned int> >(name, par));
}

void inputParser::addFloatParameter(std::string name, float *var, bool required)
{
    _parameter<float> par(name, var, required);
    _floatParameters.insert(std::pair<std::string, _parameter<float> >(name, par));
}

void inputParser::addDoubleParameter(std::string name, double *var, bool required)
{
    _parameter<double> par(name, var, required);
    _doubleParameters.insert(std::pair<std::string, _parameter<double> >(name, par));
}

void inputParser::addBoolParameter(std::string name, bool *var, bool required)
{
    _parameter<bool> par(name, var, required);
    _boolParameters.insert(std::pair<std::string, _parameter<bool> >(name, par));
}

void inputParser::addStringParameter(std::string name, std::string *var, bool required)
{
    _parameter<std::string> par(name, var, required);
    _stringParameters.insert(std::pair<std::string, _parameter<std::string> >(name, par));
}

void inputParser::printParameterInfo(std::ostream &out)
{

    std::map<std::string, _parameter<int> >::iterator intIt;
    std::map<std::string, _parameter<unsigned int> >::iterator uIntIt;
    std::map<std::string, _parameter<float> >::iterator floatIt;
    std::map<std::string, _parameter<double> >::iterator doubleIt;
    std::map<std::string, _parameter<bool> >::iterator boolIt;
    std::map<std::string, _parameter<std::string> >::iterator stringIt;
    out << "#########################################" << std::endl;
    out << "PARAMETER:\t\tVALUE:" << std::endl;
    out << "#########################################" << std::endl;
    out << "-----------------------------------------" << std::endl;
    for (intIt = _intParameters.begin(); intIt != _intParameters.end(); ++intIt)
    {
        intIt->second.printParameterInfo();
        out << "-----------------------------------------" << std::endl;
    }
    for (uIntIt = _uintParameters.begin(); uIntIt != _uintParameters.end(); ++uIntIt)
    {
        uIntIt->second.printParameterInfo();
        out << "-----------------------------------------" << std::endl;
    }
    for (floatIt = _floatParameters.begin(); floatIt != _floatParameters.end(); ++floatIt)
    {
        floatIt->second.printParameterInfo();
        out << "-----------------------------------------" << std::endl;
    }
    for (doubleIt = _doubleParameters.begin(); doubleIt != _doubleParameters.end(); ++doubleIt)
    {
        doubleIt->second.printParameterInfo();
        out << "-----------------------------------------" << std::endl;
    }
    for (boolIt = _boolParameters.begin(); boolIt != _boolParameters.end(); ++boolIt)
    {
        boolIt->second.printParameterInfo();
        out << "-----------------------------------------" << std::endl;
    }
    for (stringIt = _stringParameters.begin(); stringIt != _stringParameters.end(); ++stringIt)
    {
     	stringIt->second.printParameterInfo();
        out << "-----------------------------------------" << std::endl;
    }
    out << "#########################################" << std::endl;
}

bool inputParser::validateParameters(std::ostream& errOut)
{

    int nNonCriticalMissing = 0;
    int nCriticalMissing = 0;

    std::map<std::string, _parameter<int> >::iterator intIt;
    std::map<std::string, _parameter<float> >::iterator floatIt;
    std::map<std::string, _parameter<double> >::iterator doubleIt;
    std::map<std::string, _parameter<bool> >::iterator boolIt;

    for (intIt = _intParameters.begin(); intIt != _intParameters.end(); ++intIt)
    {
        if (!intIt->second._found)
        {
            if (intIt->second._required)
            {
                errOut << "Could not find parameter: " << intIt->second._name << " which is required. Please specify this parameter in the config file!" << std::endl;
                nCriticalMissing++;
            }
            else
            {
                nNonCriticalMissing++;
            }
        }
    }
    for (floatIt = _floatParameters.begin(); floatIt != _floatParameters.end(); ++floatIt)
    {
        if (!floatIt->second._found)
        {
            if (floatIt->second._required)
            {
                errOut << "Could not find parameter: " << floatIt->second._name << " which is required. Please specify this parameter in the config file!" << std::endl;
                nCriticalMissing++;
            }
            else
            {
                nNonCriticalMissing++;
            }
        }
    }
    for (doubleIt = _doubleParameters.begin(); doubleIt != _doubleParameters.end(); ++doubleIt)
    {
        if (!doubleIt->second._found)
        {
            if (doubleIt->second._required)
            {
                errOut << "Could not find parameter: " << doubleIt->second._name << " which is required. Please specify this parameter in the config file!" << std::endl;
                nCriticalMissing++;
            }
            else
            {
                nNonCriticalMissing++;
            }
        }
    }
    for (boolIt = _boolParameters.begin(); boolIt != _boolParameters.end(); ++boolIt)
    {
        if (!boolIt->second._found)
        {
            if (boolIt->second._required)
            {
                errOut << "Could not find parameter: " << boolIt->second._name << " which is required. Please specify this parameter in the config file!" << std::endl;
                nCriticalMissing++;
            }
            else
            {
                nNonCriticalMissing++;
            }
        }
    }
    if(nCriticalMissing > 0) return false;
    return true;
}

