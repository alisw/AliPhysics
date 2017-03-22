// PartonSystems.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains auxiliary classes for the parton-level processes.
// PartonSystem contains info on a single partonic subcollision.
// PartonSystems describes the set of subcollisions in the whole event.

#ifndef Pythia8_PartonSystems_H
#define Pythia8_PartonSystems_H

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// The PartonSystem class contains info on an individual singlet.
// Only to be used inside PartonSystems, so no private members.

class PartonSystem {

public:

  // Constructors.
  PartonSystem() : iInA(0), iInB(0), sHat(0.) {iOut.reserve(10);}

  // Stored quantities.
  int         iInA, iInB;
  vector<int> iOut;
  double      sHat, pTHat;

};

//==========================================================================

// The PartonSystems class describes the whole set of subcollisions.

class PartonSystems {

public:

  // Constructor.
  PartonSystems() {systems.resize(0);}

  // Reset system list to empty.
  void clear() {systems.resize(0);}

  // Add new subsystem to list; return its index. Number of subsystems.
  int addSys() {systems.push_back(PartonSystem());
    return systems.size() - 1;}
  int sizeSys() const {return systems.size();}

  // Set, add or replace info to one system.
  void setInA(int iSys, int iPos) {systems[iSys].iInA = iPos;}
  void setInB(int iSys, int iPos) {systems[iSys].iInB = iPos;}
  void addOut(int iSys, int iPos) {systems[iSys].iOut.push_back(iPos);}
  void setOut(int iSys, int iMem, int iPos) {systems[iSys].iOut[iMem] = iPos;}
  void replace(int iSys, int iPosOld, int iPosNew);
  void setSHat(int iSys, double sHatIn) {systems[iSys].sHat = sHatIn;}
  void setPTHat(int iSys, double pTHatIn) {systems[iSys].pTHat = pTHatIn;}
  void setSizeSys(int iSize) {systems.resize(iSize);}

  // Get info on one system.
  bool hasInAB(int iSys)         const {return ( (systems[iSys].iInA > 0)
                                        || (systems[iSys].iInB > 0) ) ;}
  int getInA(int iSys)           const {return systems[iSys].iInA;}
  int getInB(int iSys)           const {return systems[iSys].iInB;}
  int sizeOut(int iSys)          const {return systems[iSys].iOut.size();}
  int getOut(int iSys, int iMem) const {return systems[iSys].iOut[iMem];}
  int sizeAll(int iSys)          const {return (hasInAB(iSys))
    ? systems[iSys].iOut.size() + 2 : systems[iSys].iOut.size();}
  int getAll(int iSys, int iMem) const;
  double getSHat(int iSys)       const {return systems[iSys].sHat;}
  double getPTHat(int iSys)      const {return systems[iSys].pTHat;}

  // Find system of given outgoing parton, optionally also incoming one.
  int getSystemOf(int iPos, bool alsoIn = false) const;

  // Find iOut index of given system and event record index
  int getIndexOfOut(int iSys, int iPos) const;

  // List all current systems.
  void list(ostream& os = cout) const;

private:

  // List of all separate partonic subsystems.
  vector<PartonSystem> systems;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PartonSystems_H
