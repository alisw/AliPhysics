// main93LinkDef.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header used to generate a ROOT dictionary for the PYTHIA classes.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ namespace Pythia8;
#pragma link C++ class Pythia8::Particle+;
#pragma link C++ class RootTrack+;
#pragma link C++ class RootEvent+;
#endif
