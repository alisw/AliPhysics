// PartonSystems.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
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
  PartonSystem() : hard(false), iInA(0), iInB(0), iInRes(0), sHat(0.) {
    iOut.reserve(10);}

  // Stored quantities.
  bool        hard;
  int         iInA, iInB, iInRes;
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
  void setHard(int iSys, bool hard) {systems[iSys].hard = hard;}
  void setInA(int iSys, int iPos) {systems[iSys].iInA = iPos;}
  void setInB(int iSys, int iPos) {systems[iSys].iInB = iPos;}
  void setInRes(int iSys, int iPos) {systems[iSys].iInRes = iPos;}
  void addOut(int iSys, int iPos) {systems[iSys].iOut.push_back(iPos);}
  void popBackOut(int iSys) {systems[iSys].iOut.pop_back();}
  void setOut(int iSys, int iMem, int iPos) {systems[iSys].iOut[iMem] = iPos;}
  void replace(int iSys, int iPosOld, int iPosNew);
  void setSHat(int iSys, double sHatIn) {systems[iSys].sHat = sHatIn;}
  void setPTHat(int iSys, double pTHatIn) {systems[iSys].pTHat = pTHatIn;}
  void setSizeSys(int iSize) {systems.resize(iSize);}

  // Get info on one system.
  bool hasInAB(int iSys)         const {return ( (systems[iSys].iInA > 0)
                                        && (systems[iSys].iInB > 0) ) ;}
  bool hasInRes(int iSys)        const {return (systems[iSys].iInRes > 0);}
  bool getHard(int iSys)         const {return systems[iSys].hard;}
  int getInA(int iSys)           const {return systems[iSys].iInA;}
  int getInB(int iSys)           const {return systems[iSys].iInB;}
  int getInRes(int iSys)         const {return systems[iSys].iInRes;}
  int sizeOut(int iSys)          const {return systems[iSys].iOut.size();}
  int getOut(int iSys, int iMem) const {return systems[iSys].iOut[iMem];}
  int sizeAll(int iSys)          const {return (systems[iSys].iOut.size()
      + (hasInAB(iSys) ? 2 : 0) + (hasInRes(iSys) ? 1 : 0));}
  int getAll(int iSys, int iMem) const;
  double getSHat(int iSys)       const {return systems[iSys].sHat;}
  double getPTHat(int iSys)      const {return systems[iSys].pTHat;}

  // Find system of given outgoing parton, optionally also incoming one.
  int getSystemOf(int iPos, bool alsoIn = false) const;

  // Find iOut index of given system and event record index
  int getIndexOfOut(int iSys, int iPos) const;

  // List all current systems.
  void list() const;

  // Remove the last system.
  void popBack() { systems.pop_back(); }

private:

  // List of all separate partonic subsystems.
  vector<PartonSystem> systems;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_PartonSystems_H
