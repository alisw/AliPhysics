// PartonSystems.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// PartonSystem and PartonSystems classes.

#include "PartonSystems.h"

namespace Pythia8 {

//==========================================================================

// The PartonSystems class.

//--------------------------------------------------------------------------

// Replace the index of an incoming or outgoing parton by a new index.

void PartonSystems::replace(int iSys, int iPosOld, int iPosNew) {

  if (systems[iSys].iInA == iPosOld) {
    systems[iSys].iInA = iPosNew;
    return;
  }
  if (systems[iSys].iInB == iPosOld) {
    systems[iSys].iInB = iPosNew;
    return;
  }
  for (int i = 0; i < sizeOut(iSys); ++i) 
  if (systems[iSys].iOut[i] == iPosOld) {
    systems[iSys].iOut[i] = iPosNew;
    return;
  }
  
}

//--------------------------------------------------------------------------

// Return index of any parton in system, list starting with beam remnants.

int PartonSystems::getAll(int iSys, int iMem) const {

  if (hasInAB(iSys)) {
    if (iMem == 0) return systems[iSys].iInA;
    if (iMem == 1) return systems[iSys].iInB;
    return systems[iSys].iOut[iMem - 2];
  } else return systems[iSys].iOut[iMem];

}

//--------------------------------------------------------------------------

// Find system of given outgoing parton, optionally also incoming one.

int PartonSystems::getSystemOf(int iPos, bool alsoIn) const {

  // Loop over systems and over final-state members in each system.
  for (int iSys = 0; iSys < sizeSys(); ++iSys) { 
    if (alsoIn) {
      if (systems[iSys].iInA == iPos) return iSys;
      if (systems[iSys].iInB == iPos) return iSys;
    }
    for (int iMem = 0; iMem < sizeOut(iSys); ++iMem) 
      if (systems[iSys].iOut[iMem] == iPos) return iSys;
  }

  // Failure signalled by return value -1.
  return -1; 

}

//--------------------------------------------------------------------------

// Get the iMem index of iOut for an index into the event record

int PartonSystems::getIndexOfOut(int iSys, int iPos) const {
  for (int iMem = 0; iMem < sizeOut(iSys); ++iMem)
    if (systems[iSys].iOut[iMem] == iPos) return iMem;

  // Failure signalled by return value -1.
  return -1;
}


//--------------------------------------------------------------------------

// Print members in systems; for debug mainly.

void PartonSystems::list(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA Parton Systems Listing  -------------------" 
     << "--------------------------------- "
     << "\n \n  no  inA  inB  out members  \n";
  
  // Loop over system list and over members in each system.
  for (int iSys = 0; iSys < sizeSys(); ++iSys) {
    os << " " << setw(3) << iSys << " " << setw(4) << systems[iSys].iInA
       << " " << setw(4) << systems[iSys].iInB;
    for (int iMem = 0; iMem < sizeOut(iSys); ++iMem) {
      if (iMem%16 == 0 && iMem > 0) os << "\n              ";
      os << " " << setw(4) << systems[iSys].iOut[iMem];
    }
    os << "\n";
  }

  // Alternative if no systems. Done.
  if (sizeSys() == 0) os << "    no systems defined \n";
  os << "\n --------  End PYTHIA Parton Systems Listing  ---------------"
     << "---------------------------------" << endl;

}

//==========================================================================

} // end namespace Pythia8
