///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
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
// $Rev:: 272                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-07-08 18:47:36 +0200 #$: date of last commit
//
// Description:
//      some simple streams for reporting plus some stream operators
//      for common STL classes
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef REPORTINGUTILS_H
#define REPORTINGUTILS_H


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>


//////////////////////////////////////////////////////////////////////////////
// macros for printing errors, warnings, and infos

// cuts out block "className::methodName" from __PRETTY_FUNCTION__ output
inline
std::string
getClassMethod__(std::string prettyFunction)
{
	size_t pos = prettyFunction.find("(");
	if (pos == std::string::npos)
		return prettyFunction;           // something is not right
	prettyFunction.erase(pos);         // cut away signature
	pos = prettyFunction.rfind(" ");
	if (pos == std::string::npos)
		return prettyFunction;           // something is not right
	prettyFunction.erase(0, pos + 1);  // cut away return type
	return prettyFunction;
}

#define printErr  std::cerr << "!!! " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: ERROR: "   << std::flush
#define printWarn std::cerr << ">>> " << __PRETTY_FUNCTION__ << " [" << __FILE__ << ":" << __LINE__ << "]: "<<std::endl<<"Warning: " << std::flush 
#define printInfo std::cout << ">>> " << getClassMethod__(__PRETTY_FUNCTION__) << "(): Info: "  << std::flush


//////////////////////////////////////////////////////////////////////////////
// functions to print version and compilation info

#ifndef SVN_VERSION  // SVN_VERSION set by Makefile
#define SVN_VERSION "undefined"
#endif
inline std::string svnVersion() { return SVN_VERSION; }

inline
void
printSvnVersion()
{
	const std::string ver = svnVersion();
	if (ver == "")
		printInfo << "subversion repository revision is unknown." << std::endl;
	else
	  printInfo << "subversion repository revision is '" << ver << "'" << std::endl;
}


#ifndef CMAKE_SOURCE_DIR  // CMAKE_SOURCE_DIR set by Makefile
#define CMAKE_SOURCE_DIR "undefined"
#endif
inline std::string compileDir() { return CMAKE_SOURCE_DIR; }

inline
void
printCompilerInfo()
{
	const std::string date = __DATE__;
	const std::string time = __TIME__;
	const std::string ver  = __VERSION__;
	const std::string dir  = compileDir();
	printInfo << "this executable was compiled in ";
	if (dir != "")
		std::cout << "'" << dir << "'";
	else
		std::cout << "unknown directory";
	std::cout << " on " << date << " " << time << " by compiler " << ver << std::endl;
}


//////////////////////////////////////////////////////////////////////////////
// simple stream operators for some STL classes

template<typename T1, typename T2>
inline
std::ostream&
operator << (std::ostream&            out,
             const std::pair<T1, T2>& pair)
{
	return out << "(" << pair.first << ", " << pair.second << ")";
}
	
	
template<typename T>
inline
std::ostream&
operator << (std::ostream&         out,
             const std::vector<T>& vec)
{
	out << "{";
	for (unsigned int i = 0; i < (vec.size() - 1); ++i)
		out << "[" << i << "] = " << vec[i] << ", ";
	return out << "[" << vec.size() - 1 << "] = " << vec[vec.size() - 1] << "}";
}


//////////////////////////////////////////////////////////////////////////////
// indicates progess by printing relative or absolute progress in regular intervals
inline
std::ostream&
progressIndicator(const unsigned int currentPos,
                  const unsigned int nmbTotal,
                  const bool         absolute   = false,
                  const unsigned int fieldWidth = 3,
                  const unsigned int nmbSteps   = 10,
                  std::ostream&      out        = std::cout)
{
	const double step = nmbTotal / (double)nmbSteps;
	if ((int)(currentPos / step) - (int)((currentPos - 1) / step) != 0) {
		if (absolute)
			out << "    " << std::setw(fieldWidth) << currentPos << " of " << nmbTotal << std::endl;
		else
			out << "    " << std::setw(fieldWidth) << (int)(currentPos / step) * nmbSteps << " %" << std::endl;
	}
	return out;
} 


// converts bool to "true"/"false" string
inline
std::string trueFalse(const bool val)
{
	if (val)
		return "true";
	else
		return "false";
}

// converts bool to "yes"/"no" string
inline
std::string yesNo(const bool val)
{
	if (val)
		return "yes";
	else
		return "no";
}

// converts bool to "on"/"off" string
inline
std::string onOff(const bool val)
{
	if (val)
		return "on";
	else
		return "off";
}

// converts bool to "enabled"/"disabled" string
inline
std::string enDisabled(const bool val)
{
	if (val)
		return "enabled";
	else
		return "disabled";
}


#endif  // REPORTINGUTILS_H
