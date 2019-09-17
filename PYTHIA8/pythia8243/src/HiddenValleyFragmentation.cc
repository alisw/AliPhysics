// HiddenValleyFragmentation.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// HiddenValleyFragmentation class and its helper classes.

#include "Pythia8/HiddenValleyFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The HVStringFlav class is used to select HV-quark and HV-hadron flavours.

//--------------------------------------------------------------------------

// Initialize data members of the flavour generation.

void HVStringFlav::init(Settings& settings, ParticleData* particleDataPtrIn,
  Rndm* rndmPtrIn, Info* infoPtrIn) {

  // Save pointers.
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  infoPtr         = infoPtrIn;

  // Read in data from Settings.
  nFlav           = settings.mode("HiddenValley:nFlav");
  probVector      = settings.parm("HiddenValley:probVector");
  thermalModel    = false;
  useWidthPre     = false;
  closePacking    = false;
  mT2suppression  = false;

}

//--------------------------------------------------------------------------

// Pick a new HV-flavour given an incoming one.

FlavContainer HVStringFlav::pick(FlavContainer& flavOld, double, double) {

  // Initial values for new flavour.
  FlavContainer flavNew;
  flavNew.rank = flavOld.rank + 1;

  // Pick new HV-flavour at random; keep track of sign.
  flavNew.id = 4900100 + min( 1 + int(nFlav * rndmPtr->flat()), nFlav);
  if (flavOld.id > 0) flavNew.id = -flavNew.id;

  // Done.
  return flavNew;

}

//--------------------------------------------------------------------------

// Combine two HV-flavours to produce an HV-hadron.
// This is simplified procedure, assuming only two HV mesons defined.

int HVStringFlav::combine(FlavContainer& flav1, FlavContainer& flav2) {

  // Positive and negative flavour. Note that with kinetic mixing
  // the Fv are really intended to represent qv, so remap.
  int idMeson = 0;
  int idPos =  max( flav1.id, flav2.id) - 4900000;
  int idNeg = -min( flav1.id, flav2.id) - 4900000;
  if (idPos < 20) idPos = 101;
  if (idNeg < 20) idNeg = 101;

  // Pick HV-meson code, spin either 0 or 1.
  if (idNeg == idPos)     idMeson =  4900111;
  else if (idPos > idNeg) idMeson =  4900211;
  else                    idMeson = -4900211;
  if (rndmPtr->flat() < probVector) idMeson += ((idMeson > 0) ? 2 : -2);

  // Done.
  return idMeson;

}

//==========================================================================

// The HVStringPT class is used to select pT in HV fragmentation.

//--------------------------------------------------------------------------

// Initialize data members of the string pT selection.

void HVStringPT::init(Settings& settings, ParticleData* particleDataPtrIn,
  Rndm* rndmPtrIn, Info* infoPtrIn) {

  // Save pointer.
  particleDataPtr  = particleDataPtrIn;
  rndmPtr          = rndmPtrIn;
  infoPtr          = infoPtrIn;

  // Parameter of the pT width. No enhancement, since this is finetuning.
  double sigmamqv  = settings.parm("HiddenValley:sigmamqv");
  double sigma     = sigmamqv * particleDataPtrIn->m0( 4900101);
  sigmaQ           = sigma / sqrt(2.);
  enhancedFraction = 0.;
  enhancedWidth    = 0.;

  // Parameter for pT suppression in MiniStringFragmentation.
  sigma2Had        = 2. * pow2( max( SIGMAMIN, sigma) );
  thermalModel     = false;
  useWidthPre      = false;
  closePacking     = false;

}

//==========================================================================

// The HVStringZ class is used to select z in HV fragmentation.

//--------------------------------------------------------------------------

// Initialize data members of the string z selection.

void HVStringZ::init(Settings& settings, ParticleData& particleData,
  Rndm* rndmPtrIn, Info* infoPtrIn) {

  // Save pointer.
  rndmPtr  = rndmPtrIn;
  infoPtr  = infoPtrIn;

  // Paramaters of Lund/Bowler symmetric fragmentation function.
  aLund    = settings.parm("HiddenValley:aLund");
  bmqv2    = settings.parm("HiddenValley:bmqv2");
  rFactqv  = settings.parm("HiddenValley:rFactqv");

  // Use qv mass to set scale of bEff = b * m^2;
  mqv2     = pow2( particleData.m0( 4900101) );
  bLund    = bmqv2 / mqv2;

  // Mass of qv meson used to set stop scale for fragmentation iteration.
  mhvMeson = particleData.m0( 4900111);

}

//--------------------------------------------------------------------------

// Generate the fraction z that the next hadron will take using Lund/Bowler.

double HVStringZ::zFrag( int , int , double mT2) {

  // Shape parameters of Lund symmetric fragmentation function.
  double bShape = bLund * mT2;
  double cShape = 1. + rFactqv * bmqv2;
  return zLund( aLund, bShape, cShape);

}

//==========================================================================

// The HiddenValleyFragmentation class.

//--------------------------------------------------------------------------

// Initialize and save pointers.

bool HiddenValleyFragmentation::init(Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;

  // Check whether Hidden Valley fragmentation switched on, and SU(N).
  doHVfrag = settings.flag("HiddenValley:fragment");
  if (settings.mode("HiddenValley:Ngauge") < 2) doHVfrag = false;
  if (!doHVfrag) return false;

  // Several copies of qv may be needed. Taken to have same mass.
  nFlav = settings.mode("HiddenValley:nFlav");
  if (nFlav > 1) {
    int spinType = particleDataPtr->spinType(4900101);
    double m0    = particleDataPtr->m0(4900101);
    for (int iFlav = 2; iFlav <= nFlav; ++iFlav)
      particleDataPtr->addParticle( 4900100 + iFlav, "qv", "qvbar",
      spinType, 0, 0, m0);
  }

  // Hidden Valley meson mass used to choose hadronization mode.
  mhvMeson = particleDataPtr->m0(4900111);

  // Initialize the hvEvent instance of an event record.
  hvEvent.init( "(Hidden Valley fragmentation)", particleDataPtr);

  // Create HVStringFlav instance for HV-flavour selection.
  hvFlavSelPtr = new HVStringFlav();
  hvFlavSelPtr->init( settings, particleDataPtr, rndmPtr, infoPtr);

  // Create HVStringPT instance for pT selection in HV fragmentation.
  hvPTSelPtr = new HVStringPT();
  hvPTSelPtr->init( settings, particleDataPtr, rndmPtr, infoPtr);

  // Create HVStringZ instance for z selection in HV fragmentation.
  hvZSelPtr = new HVStringZ();
  hvZSelPtr->init( settings, *particleDataPtr, rndmPtr, infoPtr);

  // Initialize auxiliary administrative class.
  hvColConfig.init(infoPtr, settings, hvFlavSelPtr);

  // Initialize HV-string and HV-ministring fragmentation.
  hvStringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr,
    hvFlavSelPtr, hvPTSelPtr, hvZSelPtr);
  hvMinistringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr,
    hvFlavSelPtr, hvPTSelPtr, hvZSelPtr);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Perform the fragmentation.

bool HiddenValleyFragmentation::fragment(Event& event) {

  // Reset containers for next event.
  hvEvent.reset();
  hvColConfig.clear();
  ihvParton.resize(0);

  // Extract HV-particles from event to hvEvent. Assign HV-colours.
  // Done if no HV-particles found.
  if (!extractHVevent(event)) return true;

  // Store found string system. Analyze its properties.
  if (!hvColConfig.insert(ihvParton, hvEvent)) return false;

  // Collect sequentially all partons in the HV subsystem.
  // Copy also if already in order, or else history tracing may fail.
  hvColConfig.collect(0, hvEvent, false);

  // Mass used to decide how to fragment system.
  mSys = hvColConfig[0].mass;

  // HV-string fragmentation when enough mass to produce >= 3 HV-mesons.
  if (mSys > 3.5 * mhvMeson) {
    if (!hvStringFrag.fragment( 0, hvColConfig, hvEvent)) return false;

  // HV-ministring fragmentation when enough mass to produce 2 HV-mesons.
  } else if (mSys > 2.1 * mhvMeson) {
    if (!hvMinistringFrag.fragment( 0, hvColConfig, hvEvent, true))
    return false;

  // If only enough mass for one HV-meson assume HV-glueballs emitted.
  } else if (!collapseToMeson()) return false;

  // Insert HV particles from hvEvent to event.
  insertHVevent(event);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Extract HV-particles from event to hvEvent. Assign HV-colours.

bool HiddenValleyFragmentation::extractHVevent(Event& event) {

  // Copy Hidden-Valley particles to special event record.
  for (int i = 0; i < event.size(); ++i) {
    int idAbs = event[i].idAbs();
    bool isHV = (idAbs > 4900000 && idAbs < 4900007)
             || (idAbs > 4900010 && idAbs < 4900017)
             || idAbs == 4900021
             || (idAbs > 4900100 && idAbs < 4900109);
    if (isHV) {
      int iHV = hvEvent.append( event[i]);
      // Convert HV-gluons into normal ones so as to use normal machinery.
      if (event[i].id() ==  4900021) hvEvent[iHV].id(21);
      // Second mother points back to position in complete event;
      // otherwise construct the HV history inside hvEvent.
      hvEvent[iHV].mothers( 0, i);
      hvEvent[iHV].daughters( 0, 0);
      int iMother = event[i].mother1();
      for (int iHVM = 1; iHVM < hvEvent.size(); ++iHVM)
      if (hvEvent[iHVM].mother2() == iMother) {
        hvEvent[iHV].mother1( iHVM);
        if (hvEvent[iHVM].daughter1() == 0) hvEvent[iHVM].daughter1(iHV);
        else                                hvEvent[iHVM].daughter2(iHV);
      }
    }
  }

  // Done if no HV particles found.
  hvOldSize = hvEvent.size();
  if (hvOldSize == 1) return false;

  // Initial colour - anticolour parton pair.
  int colBeg = hvEvent.nextColTag();
  for (int iHV = 1; iHV < hvOldSize; ++iHV)
  if (hvEvent[iHV].mother1() == 0) {
    if (hvEvent[iHV].id() > 0) hvEvent[iHV].col( colBeg);
    else                       hvEvent[iHV].acol( colBeg);
  }

  // Then trace colours down to daughters; new colour if two daughters.
  for (int iHV = 1; iHV < hvOldSize; ++iHV) {
    int dau1 = hvEvent[iHV].daughter1();
    int dau2 = hvEvent[iHV].daughter2();
    if (dau1 > 0 && dau2 == 0)
      hvEvent[dau1].cols( hvEvent[iHV].col(), hvEvent[iHV].acol());
    else if (dau2 > 0) {
      int colHV  = hvEvent[iHV].col();
      int acolHV = hvEvent[iHV].acol();
      int colNew = hvEvent.nextColTag();
      if (acolHV == 0) {
        hvEvent[dau1].cols( colNew, 0);
        hvEvent[dau2].cols( colHV, colNew);
      } else if (colHV == 0) {
        hvEvent[dau1].cols( 0, colNew);
        hvEvent[dau2].cols( colNew, acolHV);
      // Temporary: should seek recoiling dipole end!??
      } else if (rndmPtr->flat() > 0.5) {
        hvEvent[dau1].cols( colHV, colNew);
        hvEvent[dau2].cols( colNew, acolHV);
      } else {
        hvEvent[dau1].cols( colNew, acolHV);
        hvEvent[dau2].cols( colHV, colNew);
      }
    }
  }

  // Pick up the colour end.
  int colNow = 0;
  for (int iHV = 1; iHV < hvOldSize; ++iHV)
  if (hvEvent[iHV].isFinal() && hvEvent[iHV].acol() == 0) {
    ihvParton.push_back( iHV);
    colNow = hvEvent[iHV].col();
  }

  // Trace colour by colour until reached anticolour end.
  while (colNow > 0) {
    for (int iHV = 1; iHV < hvOldSize; ++iHV)
    if (hvEvent[iHV].isFinal() && hvEvent[iHV].acol() == colNow) {
      ihvParton.push_back( iHV);
      colNow = hvEvent[iHV].col();
      break;
    }
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Collapse of light system to one HV-meson, by the emission of HV-glueballs.

bool HiddenValleyFragmentation::collapseToMeson() {

  // If too low mass then cannot do anything. Should not happen.
  if (mSys < 1.001 * mhvMeson) {
    infoPtr->errorMsg("Error in HiddenValleyFragmentation::collapseToMeson:"
      " too low mass to do anything");
    return false;
  }

  // Choose mass of collective HV-glueball states flat between limits.
  double mhvGlue = (0.001 + 0.998 * rndmPtr->flat()) * (mSys - mhvMeson);

  // Find momentum in rest frame, with isotropic "decay" angles.
  double pAbs = 0.5 * sqrtpos( pow2(mSys*mSys - mhvMeson*mhvMeson
    - mhvGlue*mhvGlue) - pow2(2. * mhvMeson * mhvGlue) ) / mSys;
  double pz   = (2 * rndmPtr->flat() - 1.) * pAbs;
  double pT   = sqrtpos( pAbs*pAbs - pz*pz);
  double phi  = 2. * M_PI * rndmPtr->flat();
  double px   = pT * cos(phi);
  double py   = pT * sin(phi);

  // Construct four-vectors and boost them to event frame.
  Vec4 phvMeson( px, py, pz, sqrt(mhvMeson*mhvMeson + pAbs*pAbs) );
  Vec4 phvGlue( -px, -py, -pz, sqrt(mhvGlue*mhvGlue + pAbs*pAbs) );
  phvMeson.bst( hvColConfig[0].pSum );
  phvGlue.bst(  hvColConfig[0].pSum );

  // Add produced particles to the event record.
  vector<int> iParton = hvColConfig[0].iParton;
  int iFirst = hvEvent.append( 4900111, 82,  iParton.front(),
    iParton.back(), 0, 0, 0, 0, phvMeson, mhvMeson);
  int iLast  = hvEvent.append( 4900991, 82,  iParton.front(),
    iParton.back(), 0, 0, 0, 0, phvGlue, mhvGlue);

  // Mark original partons as hadronized and set their daughter range.
  for (int i = 0; i < int(iParton.size()); ++i) {
    hvEvent[ iParton[i] ].statusNeg();
    hvEvent[ iParton[i] ].daughters(iFirst, iLast);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Insert HV-particles from hvEvent to event.

bool HiddenValleyFragmentation::insertHVevent(Event& event) {

  // Offset for mother/daughter indices.
  hvNewSize = hvEvent.size();
  int nOffset = event.size() - hvOldSize;

  // Copy back HV-particles.
  int iNew, iMot1, iMot2, iDau1, iDau2;
  for (int iHV = hvOldSize; iHV < hvNewSize; ++iHV) {
    iNew = event.append( hvEvent[iHV]);

    // Restore HV-gluon codes. Do not keep HV-colours, to avoid confusion.
    if (hvEvent[iHV].id() == 21) event[iNew].id(4900021);
    event[iNew].cols( 0, 0);

    // Begin history construction.
    iMot1 = hvEvent[iHV].mother1();
    iMot2 = hvEvent[iHV].mother2();
    iDau1 = hvEvent[iHV].daughter1();
    iDau2 = hvEvent[iHV].daughter2();
    // Special mother for partons copied from event, else simple offset.
    // Also set daughters of mothers in original record.
    if (iMot1 > 0 && iMot1 < hvOldSize) {
      iMot1 = hvEvent[iMot1].mother2();
      event[iMot1].statusNeg();
      event[iMot1].daughter1(iNew);
    } else if (iMot1 > 0) iMot1 += nOffset;
    if (iMot2 > 0 && iMot2 < hvOldSize) {
      iMot2 = hvEvent[iMot2].mother2();
      event[iMot2].statusNeg();
      if (event[iMot2].daughter1() == 0) event[iMot2].daughter1(iNew);
      else                               event[iMot2].daughter2(iNew);
    } else if (iMot2 > 0) iMot2 += nOffset;
    if (iDau1 > 0) iDau1 += nOffset;
    if (iDau2 > 0) iDau2 += nOffset;
    event[iNew].mothers( iMot1, iMot2);
    event[iNew].daughters( iDau1, iDau2);
  }

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
