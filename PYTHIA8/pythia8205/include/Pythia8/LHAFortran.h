// LHAFortran.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Fortran Les Houches Accord user process information.
// LHAupFortran: derived class with the HEPRUP and HEPEUP Fortran info.
// You are expected to supply the fillHepRup and fillHepEup methods.

#ifndef Pythia8_LHAFortran_H
#define Pythia8_LHAFortran_H

#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Give access to the HEPRUP and HEPEUP Fortran commonblocks.

#ifdef _WIN32
  #define heprup_ HEPRUP
  #define hepeup_ HEPEUP
#endif

extern "C" {

  extern struct {
    int idbmup[2];
    double ebmup[2];
    int pdfgup[2], pdfsup[2], idwtup, nprup;
    double xsecup[100], xerrup[100], xmaxup[100];
    int lprup[100];
  } heprup_;

  extern struct {
    int nup, idprup;
    double xwgtup, scalup, aqedup, aqcdup;
    int idup[500], istup[500], mothup[500][2], icolup[500][2];
    double pup[500][5], vtimup[500],spinup[500];
  } hepeup_;

}

//==========================================================================

// A derived class with initialization information from the HEPRUP
// Fortran commonblock and event information from the HEPEUP one.

class LHAupFortran : public LHAup {

public:

  // Constructor.
  LHAupFortran() {}

  // Routine for doing the job of setting initialization info.
  bool setInit() {
    // Call the routine that does the job.
    if (!fillHepRup()) return false;
    // Store beam and strategy info.
    setBeamA(heprup_.idbmup[0], heprup_.ebmup[0], heprup_.pdfgup[0],
      heprup_.pdfsup[0]);
    setBeamB(heprup_.idbmup[1], heprup_.ebmup[1], heprup_.pdfgup[1],
      heprup_.pdfsup[1]);
    setStrategy(heprup_.idwtup);
    // Store process info. Protect against vanishing cross section.
    for (int ip = 0; ip < heprup_.nprup; ++ip) {
      double xsec = max( 1e-10, heprup_.xsecup[ip]);
      addProcess( heprup_.lprup[ip], xsec, heprup_.xerrup[ip],
        heprup_.xmaxup[ip] );
    }
    // Store the beam energies to calculate x values later.
    eBeamA = heprup_.ebmup[0];
    eBeamB = heprup_.ebmup[1];
    // Done.
    return true;
  }

  // Routine for doing the job of setting info on next event.
  bool setEvent(int idProcIn = 0) {
    // In some strategies the type of the next event has been set.
    hepeup_.idprup = idProcIn;
    // Call the routine that does the job.
    if (!fillHepEup()) return false;
    // Store process info.
    setProcess(hepeup_.idprup, hepeup_.xwgtup, hepeup_.scalup,
      hepeup_.aqedup, hepeup_.aqcdup);
    // Store particle info.
    for (int ip = 0; ip < hepeup_.nup; ++ip) addParticle(hepeup_.idup[ip],
      hepeup_.istup[ip], hepeup_.mothup[ip][0], hepeup_.mothup[ip][1],
      hepeup_.icolup[ip][0], hepeup_.icolup[ip][1], hepeup_.pup[ip][0],
      hepeup_.pup[ip][1], hepeup_.pup[ip][2], hepeup_.pup[ip][3],
      hepeup_.pup[ip][4], hepeup_.vtimup[ip], hepeup_.spinup[ip]) ;
    // Store x values (here E = pup[ip][3]), but note incomplete info.
    setPdf( hepeup_.idup[0], hepeup_.idup[1], hepeup_.pup[0][3]/eBeamA,
      hepeup_.pup[1][3]/eBeamB, 0., 0., 0., false);
    // Done.
    return true;
  }

private:

  // Save beam energies to calculate x values.
  double eBeamA, eBeamB;

  // User-written routine that does the intialization and fills heprup.
  bool fillHepRup();

  // User-written routine that does the event generation and fills hepeup.
  bool fillHepEup();

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_LHAFortran_H
