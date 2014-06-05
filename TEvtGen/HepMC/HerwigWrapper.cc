//--------------------------------------------------------------------------
//
// HerwigWrapper.cc
// Author:  Lynn Garren
//
// ----------------------------------------------------------------------

#ifdef _WIN32 // Platform: Windows MS Visual C++

//  Sorry, there is NO version currently available for Vusual C++.

#else 

#include <cmath>

#include "HepMC/HerwigWrapper.h"
#include "HepMC/GenCrossSection.h"

// declare the struct here to keep the shared library happy
struct hwgev hwevnt_;

namespace HepMC {

GenCrossSection getHerwigCrossSection(int ngen) {

  HepMC::GenCrossSection xsec;
  // set cross section information and convert to pb (HepMC convention)
  double xsecval = hwevnt.AVWGT * 1000.0;
  // statistical error
  // Herwig has a better calculation of the error, 
  // but that information does not appear to be saved anywhere
  double xsecerr = xsecval / std::sqrt((double)ngen);
  // set and return cross section information
  xsec.set_cross_section(xsecval, xsecerr);
  return xsec;
}

} // HepMC

#endif //Platform
