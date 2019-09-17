// PartonSystems.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// PartonSystem and PartonSystems classes.

#include "Pythia8/PartonSystems.h"

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
  if (systems[iSys].iInRes == iPosOld) {
    systems[iSys].iInRes = iPosNew;
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
  } else if (hasInRes(iSys)) {
    if (iMem == 0) return systems[iSys].iInRes;
    return systems[iSys].iOut[iMem - 1];
  }
  return systems[iSys].iOut[iMem];

}

//--------------------------------------------------------------------------

// Find system of given outgoing parton, optionally also incoming one.
// If the parton is outgoing in one system and incoming in another (eg a
// decaying resonance), the system in which it is incoming will be returned if
// alsoIn == true, else the system in which it is outgoing will be returned.

int PartonSystems::getSystemOf(int iPos, bool alsoIn) const {

  // If (alsoIn), first check if this parton appears as incoming in any system.
  if (alsoIn) {
    for (int iSys = 0; iSys < sizeSys(); ++iSys) {
      if (systems[iSys].iInA == iPos) return iSys;
      if (systems[iSys].iInB == iPos) return iSys;
      if (systems[iSys].iInRes == iPos) return iSys;
    }
  }

  // Then check if it appears as outgoing in any system.
  for (int iSys = 0; iSys < sizeSys(); ++iSys) {
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

void PartonSystems::list() const {

  // Header.
  cout << "\n --------  PYTHIA Parton Systems Listing  -------------------"
       << "--------------------------------- "
       << "\n \n  no  inA  inB  out members  \n";

  // Loop over system list and over members in each system.
  for (int iSys = 0; iSys < sizeSys(); ++iSys) {
    cout << " " << setw(3) << iSys << " ";
    if (hasInAB(iSys)) {
      cout << setw(4) << systems[iSys].iInA
           << " " << setw(4) << systems[iSys].iInB;
    } else if (hasInRes(iSys)) {
      cout << "  (" << setw(4) << systems[iSys].iInRes
           << ") ";
    } else cout<< setw(9) <<" "<<endl;
    for (int iMem = 0; iMem < sizeOut(iSys); ++iMem) {
      if (iMem%16 == 0 && iMem > 0) cout << "\n              ";
      cout << " " << setw(4) << systems[iSys].iOut[iMem];
    }
    cout << "\n";
  }

  // Alternative if no systems. Done.
  if (sizeSys() == 0) cout << "    no systems defined \n";
  cout << "\n --------  End PYTHIA Parton Systems Listing  ---------------"
       << "---------------------------------" << endl;

}

//==========================================================================

} // end namespace Pythia8
