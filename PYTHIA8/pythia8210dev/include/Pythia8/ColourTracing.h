// ColourTracing.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the class ColurTracing.
// ColourTracing traces colour lines in the event record.


#ifndef Pythia8_ColourTracing_H
#define Pythia8_ColourTracing_H

#include "Pythia8/Event.h"
#include "Pythia8/Info.h"

namespace Pythia8 {

//==========================================================================

// ColourTracing class. It is used to trace colours within the event record.

class ColourTracing {

public:

  void init( Info* infoPtrIn) {infoPtr = infoPtrIn;}

  // Setup the colour lists.
  bool setupColList(Event& event);

  // Trace a colour line, from a colour, from an anticolour, or in loop.
  bool traceFromAcol(int indxCol, Event& event, int iJun, int iCol,
    vector<int>& iParton);
  bool traceFromCol(int indxCol, Event& event, int iJun, int iCol,
    vector<int>& iParton);
  bool traceInLoop(Event& event, vector<int>& iParton);

  bool finished() { return (int(iColAndAcol.size()) == 0);}
  bool colFinished() { return (int(iColEnd.size()) == 0);}

  // Get junction chains where the junctions are directly connected.
  vector<vector<int > > getJunChains(Event& event);

private:

   vector<int> iColEnd, iAcolEnd, iColAndAcol;

  // Pointer to various information on the generation.
  Info*          infoPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ColourTracing_H
