#ifndef HEPMC_VERSION_H
#define HEPMC_VERSION_H
// ----------------------------------------------------------------------
//
// Version.h
// Author:  Lynn Garren
//
//  for now, these are free functions
//
// ----------------------------------------------------------------------

#include <string>
#include <iostream>
#include "HepMC/HepMCDefs.h"

namespace HepMC {

void version( std::ostream & os = std::cout );			//!< print HepMC version
void writeVersion( std::ostream & os );	//!< write HepMC version to os
std::string versionName( );	//!< return HepMC version

inline std::string versionName( )
{
    return HEPMC_VERSION;
}

inline void version( std::ostream & os )
{
    os << " --------------- HepMC Version " << versionName()
       << " --------------- " << std::endl;
}

inline void writeVersion( std::ostream & os )
{
    os << "             HepMC Version: " << versionName() << std::endl;
}

}	// HepMC

#endif // HEPMC_VERSION_H
