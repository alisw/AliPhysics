// RHadrons.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains a class for the production and decay
// of long-lived heavy coloured particles, for now the gluino.

#ifndef Pythia8_RHadrons_H
#define Pythia8_RHadrons_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/FragmentationFlavZpT.h"
#include "Pythia8/FragmentationSystems.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

// The RHadrons class contains the routines for the production and decay
// of long-lived heavy coloured particles.

class RHadrons {

public:

  // Constructor.
  RHadrons() : allowRH(), allowRSb(), allowRSt(), allowRGo(), allowSomeR(),
    setMassesRH(), idRSb(), idRSt(), idRGo(), maxWidthRH(), probGluinoballRH(),
    mOffsetCloudRH(), mCollapseRH(), diquarkSpin1RH(), m0Sb(), m0St(), m0Go(),
    nRHad(0), iRHad(), iBef(), iSys(), systemPtr(), infoPtr(),
    particleDataPtr(), rndmPtr(), flavSelPtr(), zSelPtr() {}

  // Initialization of R-hadron handling.
  bool init( Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn);

  // Pointers to flavours and z sent from HadronLevel.
  void fragPtrs( StringFlav* flavSelPtrIn, StringZ* zSelPtrIn)
    { flavSelPtr = flavSelPtrIn; zSelPtr = zSelPtrIn;}

  // Produce R-hadrons.
  bool produce( ColConfig& colConfig, Event& event);

  // Decay R-hadrons.
  bool decay( Event& event);

  // Tell whether a given particle is supposed to form R-hadrons.
  bool givesRHadron(int id);

  // Tell whether any R-hadrons have been formed.
  bool exist() {return (nRHad > 0);}

  // Tell whether a R-hadron production+decay happened, and trace down.
  int trace(int i) {
    for (int iR = 0; iR < nRHad; ++iR)
      if (iBefRHad[iR] == i || iCreRHad[iR] == i) return iAftRHad[iR];
    return 0;}

private:

  // Constants: could only be changed in the code itself.
  static const int    IDRHADSB[14], IDRHADST[14], IDRHADGO[38], NTRYMAX;
  static const double MSAFETY, EGBORROWMAX;

  // Initialization data, mainly read from Settings.
  bool   allowRH, allowRSb, allowRSt, allowRGo, allowSomeR, setMassesRH;
  int    idRSb, idRSt, idRGo;
  double maxWidthRH, probGluinoballRH, mOffsetCloudRH, mCollapseRH,
         diquarkSpin1RH, m0Sb, m0St, m0Go;

  // Current event properties.
  vector<int>  iBefRHad, iCreRHad, iRHadron, iAftRHad;
  vector<bool> isTriplet;
  int          nRHad, iRHad, iBef, iSys;
  ColSinglet*  systemPtr;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to classes for flavour and z generation.
  StringFlav*    flavSelPtr;
  StringZ*       zSelPtr;

  // Split a system that contains both a sparticle and a junction.
  bool splitOffJunction( ColConfig& colConfig, Event& event);

  // Open up a closed gluon/gluino loop.
  bool openClosedLoop( ColConfig& colConfig, Event& event);

  // Split a single colour singlet that contains two sparticles.
  bool splitSystem( ColConfig& colConfig, Event& event);

  // Produce a R-hadron from a squark.
  bool produceSquark( ColConfig& colConfig, Event& event);

  // Produce a R-hadron from a gluino.
  bool produceGluino( ColConfig& colConfig, Event& event);

  // Construct R-hadron code from squark and (di)quark codes.
  int toIdWithSquark( int id1, int id2);

  // Construct squark and (di)quark codes from R-hadron code.
  pair<int,int> fromIdWithSquark( int idRHad);

  // Construct R-hadron code from endpoints and a gluino.
  int toIdWithGluino( int id1, int id2);

  // Construct endpoint codes from R-hadron code with a gluino.
  pair<int,int> fromIdWithGluino( int idRHad);

  // Construct modified four-vectors to match modified masses.
  bool newKin( Vec4 pOld1, Vec4 pOld2, double mNew1, double mNew2,
    Vec4& pNew1, Vec4& pNew2, bool checkMargin = true);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_RHadrons_H
