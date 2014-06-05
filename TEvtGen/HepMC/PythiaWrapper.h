//--------------------------------------------------------------------------
#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

//////////////////////////////////////////////////////////////////////////
// Matt.Dobbs@Cern.CH, November 2000
// Just a link to whichever pythia version is current.
//////////////////////////////////////////////////////////////////////////

// This pre-compiler directive is included (2002-01-16) to allow compatibility
// with MS Visual C++, which interfaces to fortran in a different manner.
// For it to work you need to define the _WIN32 variable when compiling.
#ifdef _WIN32 // Platform: Windows MS Visual C++
#include "HepMC/PythiaWrapper6_4_WIN32.h"

#else // Generic version, tested on Linux ecgs/gcc
#include "HepMC/PythiaWrapper6_4.h"

#endif // Platform

#include <cmath>

#include "HepMC/GenCrossSection.h"

namespace HepMC {

/// calculate the Pythia cross section and statistical error
inline GenCrossSection getPythiaCrossSection() {

  GenCrossSection xsec;
  // xsec(0,2) contains the sum of differential cross sections in mb
  // ngen(0,2) contains the combined number of generated events
  // convert to pb (HepMC convention)
  double xsecval = pyint5.xsec[2][0] * 1.0e9;
  // statistical error
  double xsecerr = xsecval / std::sqrt( (double)pyint5.ngen[2][0] );
  // set and return cross section information
  xsec.set_cross_section(xsecval, xsecerr);
  return xsec;
}



} // HepMC

#endif  // PYTHIA_WRAPPER_H
//--------------------------------------------------------------------------
