// BeamRemnants.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// BeamRemnants class.

#include "Pythia8/BeamRemnants.h"

namespace Pythia8 {

//==========================================================================

// The BeamRemnants class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// If same (anti)colour appears twice in final state, repair or reject.
const bool   BeamRemnants::ALLOWCOLOURTWICE = true;

// Maximum number of tries to match colours and kinematics in the event.
const int    BeamRemnants::NTRYCOLMATCH     = 10;
const int    BeamRemnants::NTRYKINMATCH     = 10;

// Overall correction step for energy-momentum conservation; only
// becomes relevant in rescattering scenarios when FSR dipole emissions
// and primordial kT is added. Should hopefully not be needed.
const bool   BeamRemnants::CORRECTMISMATCH  = false;

//--------------------------------------------------------------------------

// Initialization.

bool BeamRemnants::init( Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
  BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
  PartonSystems* partonSystemsPtrIn, PartonVertex* partonVertexPtrIn,
  ParticleData* particleDataPtrIn,
  ColourReconnection* colourReconnectionPtrIn) {

  // Save pointers.
  infoPtr               = infoPtrIn;
  rndmPtr               = rndmPtrIn;
  beamAPtr              = beamAPtrIn;
  beamBPtr              = beamBPtrIn;
  partonSystemsPtr      = partonSystemsPtrIn;
  partonVertexPtr       = partonVertexPtrIn;
  colourReconnectionPtr = colourReconnectionPtrIn;
  particleDataPtr       = particleDataPtrIn;

  // Width of primordial kT distribution.
  doPrimordialKT      = settings.flag("BeamRemnants:primordialKT");
  primordialKTsoft    = settings.parm("BeamRemnants:primordialKTsoft");
  primordialKThard    = settings.parm("BeamRemnants:primordialKThard");
  primordialKTremnant = settings.parm("BeamRemnants:primordialKTremnant");
  halfScaleForKT      = settings.parm("BeamRemnants:halfScaleForKT");
  halfMassForKT       = settings.parm("BeamRemnants:halfMassForKT");
  reducedKTatHighY    = settings.parm("BeamRemnants:reducedKTatHighY");

  // Handling of rescattering kinematics uncertainties from primodial kT.
  allowRescatter    = settings.flag("MultipartonInteractions:allowRescatter");
  doRescatterRestoreY = settings.flag("BeamRemnants:rescatterRestoreY");

  // Choice of beam remnant and colour reconnection scenarios.
  remnantMode         = settings.mode("BeamRemnants:remnantMode");
  doReconnect         = settings.flag("ColourReconnection:reconnect");
  reconnectMode       = settings.mode("ColourReconnection:mode");

  // Do multiparton interactions.
  doMPI               = settings.flag("PartonLevel:MPI");

  // Check that remnant model and colour reconnection model work together.
  if (remnantMode == 1 && reconnectMode == 0) {
    infoPtr->errorMsg("Abort from BeamRemnants::init: The remnant model"
      " and colour reconnection model does not work together");
    return false;
  }

  // Total and squared CM energy at nominal energy.
  eCM                 = infoPtr->eCM();
  sCM                 = eCM * eCM;

  // Initialize junction splitting class.
  junctionSplitting.init(infoPtr, settings, rndmPtr, particleDataPtrIn);

  // Possibility to set parton vertex information.
  doPartonVertex      = settings.flag("PartonVertex:setVertex")
                     && (partonVertexPtr != 0);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select the flavours/kinematics/colours of the two beam remnants.
// Notation: iPar = all partons, iSys = matched systems of two beams,
//           iRem = additional partons in remnants.

bool BeamRemnants::add( Event& event, int iFirst, bool doDiffCR) {

  // Update to current CM energy.
  eCM     = infoPtr->eCM();
  sCM     = eCM * eCM;

  // Check that flavour bookkept in event and in beam particles agree.
  for (int i = 0; i < beamAPtr->size(); ++i) {
    int j = (*beamAPtr)[i].iPos();
    if ((*beamAPtr)[i].id() != event[j].id()) {
      infoPtr->errorMsg("Error in BeamRemnants::add: "
        "event and beam flavours do not match");
      return false;
    }
  }
  for (int i = 0; i < beamBPtr->size(); ++i) {
    int j =  (*beamBPtr)[i].iPos();
    if ((*beamBPtr)[i].id() != event[j].id()) {
      infoPtr->errorMsg("Error in BeamRemnants::add: "
        "event and beam flavours do not match");
      return false;
    }
  }

  // Deeply inelastic scattering needs special remnant handling.
  // Here getGammaMode separates from photoproduction.
  isDIS = (beamAPtr->isLepton() && !beamBPtr->isLepton()
           && beamAPtr->getGammaMode() == 0)
       || (beamBPtr->isLepton() && !beamAPtr->isLepton()
           && beamBPtr->getGammaMode() == 0);

  // Number of scattering subsystems. Size of event record before treatment.
  nSys    = partonSystemsPtr->sizeSys();
  oldSize = event.size();

  // Store event as it was before adding anything.
  Event eventSave = event;
  BeamParticle beamAsave = (*beamAPtr);
  BeamParticle beamBsave = (*beamBPtr);
  PartonSystems partonSystemsSave = (*partonSystemsPtr);

  // Two different methods to add the beam remnants.
  if (remnantMode == 0) {
    if (!addOld(event)) return false;
  } else
    if (!addNew(event)) return false;

  if (isDIS) return true;

  // Store event before doing colour reconnections.
  Event eventTmpSave = event;
  bool colCorrect = false;
  for (int i = 0; i < 10; ++i) {
    if (doReconnect && doDiffCR
    && (reconnectMode == 1 || reconnectMode == 2)) {
      colourReconnectionPtr->next(event, iFirst);

      // Check that the new colour structure is physical.
      if (!junctionSplitting.checkColours(event))
        event = eventTmpSave;
      else {
        colCorrect = true;
        break;
      }
      // If no colour reconnection, just check the colour configuration once.
    } else {
      if (junctionSplitting.checkColours(event))
        colCorrect = true;
      break;
    }
  }

  // Possibility to add vertex information to beam particles.
  if (doPartonVertex) {
    BeamParticle* beamPtr = beamAPtr;
    // Add vertex information for both beams.
    for (int beam = 0; beam < 2; ++beam) {
      for (int i = 0; i < beamPtr->size(); ++i) {
        int j = (*beamPtr)[i].iPos();
        // We might have daughters.
        vector<int> dList = event[j].daughterList();
        // First the beam remnant particle itself.
        partonVertexPtr->vertexBeam(j, beam, event);
        // Then possible daughters.
        for(int k = 0, N = dList.size(); k < N; ++k )
                partonVertexPtr->vertexBeam(dList[k],beam,event);
      }
      // Switch beam.
      beamPtr = beamBPtr;
    }
  }

  // Restore event and return false if colour reconnection failed.
  if (!colCorrect) {
    event = eventSave;
    (*beamAPtr) = beamAsave;
    (*beamBPtr) = beamBsave;
    (*partonSystemsPtr) = partonSystemsSave;
    infoPtr->errorMsg("Error in BeamRemnants::add: "
      "failed to find physical colour state after colour reconnection");
    return false;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Old function for adding beam remnant.

bool BeamRemnants::addOld( Event& event) {

  // Add required extra remnant flavour content.
  // Start all over if fails (in option where junctions not allowed).
  if ( !beamAPtr->remnantFlavours(event, isDIS)
    || !beamBPtr->remnantFlavours(event, isDIS) ) {
    infoPtr->errorMsg("Error in BeamRemnants::add:"
      " remnant flavour setup failed");
    return false;
  }

  // Do the kinematics of the collision subsystems and two beam remnants.
  if (!setKinematics(event)) return false;

  // Allow colour reconnections.
  if (doReconnect && reconnectMode == 0 && !isDIS)
    colourReconnectionPtr->next(event, oldSize);

  // Save current modifiable colour configuration for fast restoration.
  vector<int> colSave;
  vector<int> acolSave;
  for (int i = oldSize; i < event.size(); ++i) {
    colSave.push_back( event[i].col() );
    acolSave.push_back( event[i].acol() );
  }
  event.saveJunctionSize();

  // Allow several tries to match colours of initiators and remnants.
  // Frequent "failures" since shortcutting colours separately on
  // the two event sides may give "colour singlet gluons" etc.
  bool physical = true;
  for (int iTry = 0; iTry < NTRYCOLMATCH; ++iTry) {
    physical = true;

    // Reset list of colour "collapses" (transformations).
    colFrom.resize(0);
    colTo.resize(0);

    // First process each set of beam colours on its own.
    if (!beamAPtr->remnantColours(event, colFrom, colTo))
      physical = false;
    if (!beamBPtr->remnantColours(event, colFrom, colTo))
      physical = false;

    // Then check that colours and anticolours are matched in whole event.
    if ( physical && !checkColours(event) ) physical = false;

    // If no problems then done, else restore and loop.
    if (physical) break;
    for (int i = oldSize; i < event.size(); ++i)
      event[i].cols( colSave[i - oldSize], acolSave[i - oldSize] );
    event.restoreJunctionSize();
    infoPtr->errorMsg("Warning in BeamRemnants::add:"
      " colour tracing failed; will try again");
  }

  // If no solution after several tries then failed.
  if (!physical) {
    infoPtr->errorMsg("Error in BeamRemnants::add:"
      " colour tracing failed after several attempts");
    return false;
  }

  // Possibility to add vertex information to beam particles.
  if (doPartonVertex) {
    BeamParticle* beamPtr = beamAPtr;
    // Add vertex information for both beams.
    for (int beam = 0; beam < 2; ++beam) {
      for (int i = 0; i < beamPtr->size(); ++i) {
        int j = (*beamPtr)[i].iPos();
        // We might have daughters.
        vector<int> dList = event[j].daughterList();
        // First the beam remnant particle itself.
        partonVertexPtr->vertexBeam(j, beam, event);
        // Then possible daughters.
        for(int k = 0, N = dList.size(); k < N; ++k )
                partonVertexPtr->vertexBeam(dList[k],beam,event);
      }
      // Switch beam.
      beamPtr = beamBPtr;
    }
  }
  // Done.
  return true;
}

//--------------------------------------------------------------------------

// New function for adding beam remnant.

bool BeamRemnants::addNew( Event& event) {

   // Start by saving a copy of the event, if the beam remnant fails.
  Event eventSave = event;
  BeamParticle beamAsave = (*beamAPtr);
  BeamParticle beamBsave = (*beamBPtr);
  PartonSystems partonSystemsSave = (*partonSystemsPtr);

  // Do several tries in case an unphysical colour contruction is made.
  bool beamRemnantFound = false;
  int nMaxTries = 10;

  for (int iTries = 0;iTries < nMaxTries; ++iTries) {

    // Set the initial colours.
    beamAPtr->setInitialCol(event);
    beamBPtr->setInitialCol(event);

    // Find colour state of outgoing partons and reconnect colours to match it.
    beamAPtr->findColSetup(event);
    beamBPtr->updateCol(beamAPtr->getColUpdates());

    beamBPtr->findColSetup(event);
    beamAPtr->updateCol(beamBPtr->getColUpdates());

    // Add beam remnants.
    beamAPtr->remnantFlavoursNew(event);
    beamBPtr->remnantFlavoursNew(event);

    // It is possible junctions were added, so update list.
    event.saveJunctionSize();

    // Do the kinematics of the collision subsystems and two beam remnants.
    if (!setKinematics(event)) {
      // If it does not work, try parton level again.
      event = eventSave;
      (*beamAPtr) = beamAsave;
      (*beamBPtr) = beamBsave;
      (*partonSystemsPtr) = partonSystemsSave;
      return false;
    }

    // Update the colour changes for all final state particles.
    updateColEvent(event, beamAPtr->getColUpdates());
    updateColEvent(event, beamBPtr->getColUpdates());

    // Check whether the new colour structure can be accepted.
    if (junctionSplitting.checkColours(event)) {
      beamRemnantFound = true;
      break;
    }

    // If failed, restore earlier configuration and try to find new
    // colour structure.
    else {
      event = eventSave;
      (*beamAPtr) = beamAsave;
      (*beamBPtr) = beamBsave;
      (*partonSystemsPtr) = partonSystemsSave;
    }
  }

  // Return if it was not possible to find physical colour structure.
  if (!beamRemnantFound) {
    infoPtr->errorMsg("Error in BeamRemnants::add: "
        "failed to find physical colour structure");
    // Restore event to previous state.
    event = eventSave;
    (*beamAPtr) = beamAsave;
    (*beamBPtr) = beamBsave;
    (*partonSystemsPtr) = partonSystemsSave;
    return false;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Set up trial transverse and longitudinal kinematics for each beam
// separately. Final decisions involve comparing the two beams.

bool BeamRemnants::setKinematics( Event& event) {

  // References to beams to simplify indexing.
  BeamParticle& beamA = *beamAPtr;
  BeamParticle& beamB = *beamBPtr;

  // Nothing to do for lepton-lepton scattering with all energy already used.
  if ( beamA.isUnresolvedLepton() && beamB.isUnresolvedLepton() )
    return true;

  // Photon from photon beam is unresolved.
  if ( beamA.isGamma() && beamA[0].id() == 22 ) {
    beamA.resolvedGamma(false);
  }
  if ( beamB.isGamma() && beamB[0].id() == 22 )  {
    beamB.resolvedGamma(false);
  }

  // Check if photon beams need remnants.
  bool beamAisGamma     = beamA.isGamma();
  bool beamBisGamma     = beamB.isGamma();
  bool gammaAResolved   = beamA.resolvedGamma();
  bool gammaBResolved   = beamB.resolvedGamma();
  bool gammaOneResolved = false;

  // Do nothing if ISR have found the original beam photons and only one MPI.
  if ( ( beamAisGamma && !gammaAResolved && beamBisGamma && !gammaBResolved )
    && infoPtr->nMPI() == 1 )
    return true;

  // Check that has not already used up beams. Unless direct photon in photon.
  if ( beamA.isGamma() && beamA[0].id() == 22 ) { }
  else if ( beamB.isGamma() && beamB[0].id() == 22 )  { }
  else if ( ( !(beamA.isLepton() || (beamAisGamma && !gammaAResolved) )
         && beamA.xMax(-1) <= 0.) ||
       ( !(beamB.isLepton() || (beamBisGamma && !gammaBResolved) )
         && beamB.xMax(-1) <= 0.) ) {
    infoPtr->errorMsg("Error in BeamRemnants::setKinematics:"
      " no momentum left for beam remnants");
    return false;
  }

  // Check if exactly one remnant for photon-photon collisions.
  if ( (beamAisGamma && beamBisGamma) &&
      ( (gammaAResolved && !gammaBResolved) ||
        (!gammaAResolved && gammaBResolved) ) )
    gammaOneResolved = true;

  // Unresolved photon + hadron.
  if ( ( (beamAisGamma && !gammaAResolved)
      || (beamBisGamma && !gammaBResolved) )
      && !(beamAisGamma && beamBisGamma) )
    gammaOneResolved = true;

  // Unresolved photon from lepton + hadron.
  if ( ( beamA.getGammaMode() == 2 && beamB.isHadron() )
    || ( beamB.getGammaMode() == 2 && beamA.isHadron() ) )
    gammaOneResolved = true;

  // Special kinematics setup for one-remnant systems (DIS).
  if ( (gammaOneResolved && infoPtr->nMPI() == 1) || isDIS )
    return setOneRemnKinematics(event, iDS);

  // Last beam-status particles. Offset relative to normal beam locations.
  int nBeams   = 3;
  for (int i = 3; i < event.size(); ++i)
    if (event[i].statusAbs() < 20) nBeams = i + 1;
  int nOffset  = nBeams - 3;

  // If extra photons in event fix the offset.
  if ( !(beamA.isHadron() || beamB.isHadron()) ) {
    if (beamA.hasResGamma()) --nOffset;
    if (beamB.hasResGamma()) --nOffset;
  }

  // Reserve space for extra information on the systems and beams.
  int nMaxBeam = max( beamA.size(), beamB.size() );
  vector<double> sHatSys(nMaxBeam);
  vector<double> kTwidth(nMaxBeam);
  vector<double> kTcomp(nMaxBeam);
  vector<RotBstMatrix> Msys(nSys);

  // Loop over all subsystems. Default values. Find invariant mass.
  double kTcompSumA   = 0.;
  double kTcompSumB   = 0.;
  for (int iSys = 0; iSys < nSys; ++iSys) {
    double kTwidthNow = 0.;
    double mHatDamp   = 1.;
    int iInA          = partonSystemsPtr->getInA(iSys);
    int iInB          = partonSystemsPtr->getInB(iSys);
    double sHatNow    = (event[iInA].p() + event[iInB].p()).m2Calc();

    // Allow primordial kT reduction for small-mass and small-pT systems
    // (for hardest interaction pT -> renormalization scale so also 2 -> 1).
    if (doPrimordialKT) {
      double mHat     = sqrt(sHatNow);
      double yDamp    = pow( (event[iInA].e() + event[iInB].e()) / mHat,
                        reducedKTatHighY );
      mHatDamp        = mHat / (mHat + halfMassForKT * yDamp);
      double scale    = (iSys == 0) ? infoPtr->QRen(iDS)
                      : partonSystemsPtr->getPTHat(iSys);
      kTwidthNow      = ( (halfScaleForKT * primordialKTsoft
      + scale * primordialKThard) / (halfScaleForKT + scale) ) * mHatDamp;
    }

    // Store properties of compensation systems and total compensation power.
    // Rescattered partons do not compensate, but may be massive.
    sHatSys[iSys] = sHatNow;
    kTwidth[iSys] = kTwidthNow ;
    kTcomp[iSys]  = mHatDamp;
    if (beamA[iSys].isFromBeam()) kTcompSumA += mHatDamp;
    else beamA[iSys].m( event[iInA].m() );
    if (beamB[iSys].isFromBeam()) kTcompSumB += mHatDamp;
    else beamB[iSys].m( event[iInB].m() );
  }

  // Primordial kT and compensation power among remnants.
  double kTwidthNow = (doPrimordialKT) ? primordialKTremnant : 0.;
  for (int iRem = nSys; iRem < nMaxBeam; ++iRem) {
    sHatSys[iRem] = 0.;
    kTwidth[iRem] = kTwidthNow ;
    kTcomp[iRem]  = 1.;
  }
  kTcompSumA += beamA.size() - nSys;
  kTcompSumB += beamB.size() - nSys;

  // Allow ten tries to construct kinematics (but normally works first).
  bool physical;
  double xSum[2], xInvM[2], w2Beam[2], wPosRem, wNegRem, w2Rem;
  for (int iTry = 0; iTry < NTRYKINMATCH; ++iTry) {
    physical = true;

    // Loop over the two beams.
    for (int iBeam = 0; iBeam < 2; ++iBeam) {
      BeamParticle& beam = (iBeam == 0) ? beamA : beamB;
      int nPar = beam.size();

      // Generate Gaussian pT for initiator partons inside hadrons.
      // Store/restore rescattered parton momenta before primordial kT.
      if ( (beam.isHadron() || beam.isGamma()) && doPrimordialKT ) {
        double pxSum = 0.;
        double pySum = 0.;
        for (int iPar = 0; iPar < nPar; ++iPar) {
          // Do nothing if initiator parton from gamma->qqbar
          // where gamma is the original beam particle or for
          // the other valence arising from the same splitting.
          if ( iPar == beam.gamVal() || ( beam.gamVal() >= 0 &&
              ( iPar != beam.gamVal() && beam[iPar].isValence()) ) ) ;
          else if ( beam[iPar].isFromBeam() ) {
            pair<double, double> gauss2 = rndmPtr->gauss2();
            double px = kTwidth[iPar] * gauss2.first;
            double py = kTwidth[iPar] * gauss2.second;
            beam[iPar].px(px);
            beam[iPar].py(py);
            pxSum += px;
            pySum += py;
          } else {
            int iInAB = (iBeam == 0) ? partonSystemsPtr->getInA(iPar)
                                     : partonSystemsPtr->getInB(iPar);
            beam[iPar].p( event[iInAB].p() );
          }
        }

        // Share recoil between all initiator partons, rescatterers excluded.
        double kTcompSum = (iBeam == 0) ? kTcompSumA : kTcompSumB;

        // If parton from gamma->qqbar in ISR, set pT accordingly.
        if ( beam.gamVal() >= 0 ) {

          // Get the pT2 value for the given gamma->qqbar splitting.
          double pT2corr = beam.pT2gamma2qqbar();

          // Sample the direction and set px and py.
          double phi = 2*M_PI*rndmPtr->flat();
          double px = cos(phi) * sqrt(pT2corr);
          double py = sin(phi) * sqrt(pT2corr);
          beam[beam.gamVal()].px(px);
          beam[beam.gamVal()].py(py);
          pxSum += px;
          pySum += py;

          // Do not change the momentum of these partons when sharing recoil.
          kTcompSum -= kTcomp[beam.gamVal()];

          // Other valence parton takes (most of) the recoil but this will
          // still be part of pT balancing.
          for (int iPar = 0; iPar < nPar; ++iPar)
            if ( iPar != beam.gamVal() && beam[iPar].isValence() ) {
              beam[iPar].px(-px);
              beam[iPar].py(-py);
              pxSum += -px;
              pySum += -py;
              break;
            }
        }

        for (int iPar = 0; iPar < nPar; ++iPar) {
          if (beam[iPar].isFromBeam() && iPar != beam.gamVal() ) {
            beam[iPar].px( beam[iPar].px() - pxSum * kTcomp[iPar] / kTcompSum);
            beam[iPar].py( beam[iPar].py() - pySum * kTcomp[iPar] / kTcompSum);
          }
        }

      // Without primordial kT: still need to store rescattered partons.
      } else if ( beam.isHadron() || beam.isGamma() ) {
        for (int iPar = 0; iPar < nPar; ++iPar) {
          if ( !beam[iPar].isFromBeam() ) {
            int iInAB = (iBeam == 0) ? partonSystemsPtr->getInA(iPar)
                                     : partonSystemsPtr->getInB(iPar);
            beam[iPar].p( event[iInAB].p() );
          }
        }

        // And set the pT in case gamma->qqbar splitting with MPIs.
        if ( beam.gamVal() >= 0 ) {
          double phi = 2*M_PI*rndmPtr->flat();
          double px = cos(phi) * sqrt(beam.pT2gamma2qqbar());
          double py = sin(phi) * sqrt(beam.pT2gamma2qqbar());
          beam[beam.gamVal()].px(px);
          beam[beam.gamVal()].py(py);

          // Other valence parton takes the recoil but this will still be
          // part of pT balancing.
          for (int iPar = 0; iPar < nPar; ++iPar) {
            if ( iPar != beam.gamVal() && beam[iPar].isValence() ) {
              beam[iPar].px(-px);
              beam[iPar].py(-py);
              break;
            }
          }
        }

      }

      // Pick unrescaled x values for remnants. Sum up (unscaled) p+ and p-.
      xSum[iBeam]  = 0.;
      xInvM[iBeam] = 0.;
      for (int iRem = nSys; iRem < nPar; ++iRem) {
        double xPrel = beam.xRemnant( iRem);
        beam[iRem].x(xPrel);
        xSum[iBeam]  += xPrel;
        xInvM[iBeam] += beam[iRem].mT2()/xPrel;
      }

      // Squared transverse mass for each beam, using lightcone x.
      w2Beam[iBeam] = xSum[iBeam] * xInvM[iBeam];

    // End separate treatment of the two beams.
    }

    // Recalculate kinematics of initiator systems with primordial kT.
    wPosRem = eCM;
    wNegRem = eCM;
    for (int iSys = 0; iSys < nSys; ++iSys) {
      int iA          = beamA[iSys].iPos();
      int iB          = beamB[iSys].iPos();
      double sHat     = sHatSys[iSys];

      // Classify system: rescattering on either or both sides?
      bool normalA    = beamA[iSys].isFromBeam();
      bool normalB    = beamB[iSys].isFromBeam();
      bool normalSys  = normalA && normalB;
      bool doubleRes  = !normalA && !normalB;

      // Check whether final-state system momentum matches initial-state one.
      if (allowRescatter && CORRECTMISMATCH) {
        Vec4 pInitial = event[iA].p() + event[iB].p();
        Vec4 pFinal;
        for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem) {
          int iAB      = partonSystemsPtr->getOut(iSys, iMem);
          if (event[iAB].isFinal()) pFinal += event[iAB].p();
        }

        // Scale down primordial kT from side A if p+ increased.
        if (normalA && pFinal.pPos() > pInitial.pPos())
          beamA[iSys].scalePT( pInitial.pPos() / pFinal.pPos() );

        // Scale down primordial kT from side B if p- increased.
        if (normalB && pFinal.pNeg() > pInitial.pNeg())
          beamB[iSys].scalePT( pInitial.pNeg() / pFinal.pNeg() );
      }

      // Rescatter: possible change in sign of lightcone momentum of a
      //            rescattered parton. If this happens, try to pick
      //            new primordial kT values
      if (allowRescatter
         && (event[iA].pPos() / beamA[iSys].pPos() < 0
         ||  event[iB].pNeg() / beamB[iSys].pNeg() < 0) ) {
        infoPtr->errorMsg("Warning in BeamRemnants::setKinematics:"
          " change in lightcone momentum sign; retrying kinematics");
        physical = false;
        break;
      }

      // Begin kinematics of partons after primordial kT has been added.
      double sHatTAft = sHat + pow2( beamA[iSys].px() + beamB[iSys].px())
                             + pow2( beamA[iSys].py() + beamB[iSys].py());
      double w2A      = beamA[iSys].mT2();
      double w2B      = beamB[iSys].mT2();
      double w2Diff   = sHatTAft - w2A - w2B;
      double lambda   = pow2(w2Diff) - 4. * w2A * w2B;

      // Too large transverse momenta means that kinematics will not work.
      if (lambda <= 0.) {
        physical      = false;
        break;
      }
      double lamRoot  = sqrtpos( lambda );

      // Mirror solution if the two incoming have reverse rapidity ordering.
      if (allowRescatter && doubleRes && (event[iA].pPos() * event[iB].pNeg()
        < event[iA].pNeg() * event[iB].pPos()) ) lamRoot = -lamRoot;

      // Two procedures, which agree for normal scattering, separate here.
      // First option keeps rapidity (and mass) of system unchanged by
      // primordial kT, by modification of rescattered parton.
      if (normalSys || doRescatterRestoreY || doubleRes) {

        // Find kinematics of initial system: transverse mass, and
        // longitudinal momentum carried by all or rescattered partons.
        double sHatTBef = sHat;
        double wPosBef, wNegBef, wPosBefRes, wNegBefRes;
        // Normal scattering.
        if (normalSys) {
          wPosBef       = beamA[iSys].x() * eCM;
          wNegBef       = beamB[iSys].x() * eCM;
          wPosBefRes    = 0.;
          wNegBefRes    = 0.;
        // Rescattering on side A.
        } else if (normalB) {
          sHatTBef     += event[iA].pT2();
          wPosBef       = event[iA].pPos();
          wNegBef       = event[iA].pNeg() + beamB[iSys].x() * eCM;
          wPosBefRes    = beamA[iSys].pPos();
          wNegBefRes    = beamA[iSys].pNeg();
        // Rescattering on side B.
        } else if (normalA) {
          sHatTBef     += event[iB].pT2();
          wPosBef       = beamA[iSys].x() * eCM + event[iB].pPos();
          wNegBef       = event[iB].pNeg();
          wPosBefRes    = beamB[iSys].pPos();
          wNegBefRes    = beamB[iSys].pNeg();
        // Rescattering on both sides.
        } else {
          sHatTBef     += pow2( event[iA].px() + event[iB].px())
                        + pow2( event[iA].py() + event[iB].py());
          wPosBef       = event[iA].pPos() + event[iB].pPos();
          wNegBef       = event[iA].pNeg() + event[iB].pNeg();
          wPosBefRes    = beamA[iSys].pPos() + beamB[iSys].pPos();
          wNegBefRes    = beamA[iSys].pNeg() + beamB[iSys].pNeg();
        }

        // Rescale outgoing momenta to keep same mass and rapidity of system
        // as originally, and solve for kinematics.
        double rescale  = sqrt(sHatTAft / sHatTBef);
        double wPosAft  = rescale * wPosBef;
        double wNegAft  = rescale * wNegBef;
        wPosRem        -= wPosAft - wPosBefRes;
        wNegRem        -= wNegAft - wNegBefRes;
        double wPosA    = 0.5 * (sHatTAft + w2A - w2B + lamRoot) / wNegAft;
        double wNegB    = 0.5 * (sHatTAft + w2B - w2A + lamRoot) / wPosAft;

        // Store modified beam parton momenta.
        beamA[iSys].e(  0.5 * (wPosA + w2A / wPosA) );
        beamA[iSys].pz( 0.5 * (wPosA - w2A / wPosA) );
        beamB[iSys].e(  0.5 * (w2B / wNegB + wNegB) );
        beamB[iSys].pz( 0.5 * (w2B / wNegB - wNegB) );

      // Second option keeps rescattered parton (and system mass) unchanged
      // by primordial kT, by modification of system rapidity.
      } else {

        // Rescattering on side A: preserve already scattered parton.
        if (normalB) {
          double wPosA  = beamA[iSys].pPos();
          double wNegB  = 0.5 * (w2Diff + lamRoot) / wPosA;
          beamB[iSys].e(  0.5 * (w2B / wNegB + wNegB) );
          beamB[iSys].pz( 0.5 * (w2B / wNegB - wNegB) );
          wPosRem      -= w2B / wNegB;
          wNegRem      -= wNegB;


        // Rescattering on side B: preserve already scattered parton.
        } else if (normalA) {
          double wNegB  = beamB[iSys].pNeg();
          double wPosA  = 0.5 * (w2Diff + lamRoot) / wNegB;
          beamA[iSys].e(  0.5 * (wPosA + w2A / wPosA) );
          beamA[iSys].pz( 0.5 * (wPosA - w2A / wPosA) );
          wPosRem      -= wPosA;
          wNegRem      -= w2A / wPosA;

        // Primordial kT in double rescattering does change the mass of
        // the system without any possible fix in this framework, which
        // leads to incorrect boosts. Current choice is therefore always
        // to handle it with the first procedure, where mass is retained.
        } else {
        }
      }

      // Construct system rotation and boost caused by primordial kT.
      Msys[iSys].reset();
      Msys[iSys].toCMframe( event[iA].p(), event[iB].p() );
      Msys[iSys].fromCMframe( beamA[iSys].p(), beamB[iSys].p() );

      // Boost rescattered partons in subsequent beam A list.
      for (int iSys2 = iSys + 1; iSys2 < nSys; ++iSys2) {
        if (!beamA[iSys2].isFromBeam()) {
          int iBefResc = event[ beamA[iSys2].iPos() ].mother1();
          for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem)
          if (partonSystemsPtr->getOut(iSys, iMem) == iBefResc) {
            Vec4 pTemp = event[iBefResc].p();
            pTemp.rotbst( Msys[iSys] );
            beamA[iSys2].p( pTemp );
          }
        }

        // Boost rescattered partons in subsequent beam B list.
        if (!beamB[iSys2].isFromBeam()) {
          int iBefResc = event[ beamB[iSys2].iPos() ].mother1();
          for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem)
          if (partonSystemsPtr->getOut(iSys, iMem) == iBefResc) {
            Vec4 pTemp = event[iBefResc].p();
            pTemp.rotbst( Msys[iSys] );
            beamB[iSys2].p( pTemp );
          }
        }
      }
    }

    // Check that remaining momentum is enough for remnants.
    if (wPosRem < 0. || wNegRem < 0.) physical = false;
    w2Rem = wPosRem * wNegRem;
    if (sqrtpos(w2Rem) < sqrt(w2Beam[0]) + sqrt(w2Beam[1]))
      physical = false;

    // End of loop over ten tries. Do not loop when solution found.
    if (physical) break;
  }

  // If no solution after ten tries then failed.
  if (!physical) {
    infoPtr->errorMsg("Error in BeamRemnants::setKinematics:"
      " kinematics construction failed");
    return false;
  }

  // For successful initiator kinematics process whole systems.
  Vec4 pSumOut;
  for (int iSys = 0; iSys < nSys; ++iSys) {

    // Copy initiators and their systems and boost them accordingly.
    // Update subsystem and beams info on new positions of partons.
    // Update daughter info of mothers, i.e. of beams, for hardest interaction.
    if (beamA[iSys].isFromBeam()) {
      int iA       = beamA[iSys].iPos();
      int iAcopy   = event.copy(iA, -61);
      event[iAcopy].rotbst(Msys[iSys], false);
      partonSystemsPtr->setInA(iSys, iAcopy);
      beamA[iSys].iPos( iAcopy);
      if (iSys == 0) {
        int mother = event[iAcopy].mother1();
        event[mother].daughter1(iAcopy);
      }
    }
    if (beamB[iSys].isFromBeam()) {
      int iB       = beamB[iSys].iPos();
      int iBcopy   = event.copy(iB, -61);
      event[iBcopy].rotbst(Msys[iSys], false);
      partonSystemsPtr->setInB(iSys, iBcopy);
      beamB[iSys].iPos( iBcopy);
      if (iSys == 0) {
        int mother = event[iBcopy].mother1();
        event[mother].daughter1(iBcopy);
      }
    }

    for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem) {
      int iAB      = partonSystemsPtr->getOut(iSys, iMem);
      if (event[iAB].isFinal()) {
        int iABcopy = event.copy(iAB, 62);
        event[iABcopy].rotbst(Msys[iSys], false);
        partonSystemsPtr->setOut(iSys, iMem, iABcopy);
        pSumOut   += event[iABcopy].p();
      }
    }

  }

  // Colour dipoles spanning systems gives mismatch between FSR recoils
  // and primordial kT boosts.
  if (allowRescatter && CORRECTMISMATCH) {

    // Find summed pT of beam remnants = - wanted pT of systems.
    double pxBeams = 0.;
    double pyBeams = 0.;
    for (int iRem = nSys; iRem < beamA.size(); ++iRem) {
      pxBeams     += beamA[iRem].px();
      pyBeams     += beamA[iRem].py();
    }
    for (int iRem = nSys; iRem < beamB.size(); ++iRem) {
      pxBeams     += beamB[iRem].px();
      pyBeams     += beamB[iRem].py();
    }

    // Boost all final partons in systems transversely, and also their sum.
    Vec4 pSumTo( -pxBeams, -pyBeams, pSumOut.pz(), sqrt( pow2(pxBeams)
      + pow2(pyBeams) + pow2(pSumOut.pz()) + pSumOut.m2Calc() ) );
    RotBstMatrix Mmismatch;
    Mmismatch.bst( pSumOut, pSumTo);
    for (int iSys = 0; iSys < nSys; ++iSys)
    for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem) {
      int iAB = partonSystemsPtr->getOut(iSys, iMem);
      if (event[iAB].isFinal()) event[iAB].rotbst(Mmismatch, false);
    }
    pSumOut.rotbst(Mmismatch);

   // Reset energy and momentum sum, to be compensated by beam remnants.
    wPosRem = eCM - (pSumOut.e() + pSumOut.pz());
    wNegRem = eCM - (pSumOut.e() - pSumOut.pz());
    w2Rem    = wPosRem * wNegRem;
    if ( wPosRem < 0. || wNegRem < 0.
      || sqrtpos(w2Rem) < sqrt(w2Beam[0]) + sqrt(w2Beam[1])) {
      infoPtr->errorMsg("Error in BeamRemnants::setKinematics:"
        " kinematics construction failed owing to recoil mismatch");
      return false;
    }
  }

  // Construct x rescaling factors for the two remants.
  double lambdaRoot = sqrtpos( pow2(w2Rem - w2Beam[0] - w2Beam[1])
    - 4. * w2Beam[0] * w2Beam[1] );
  double rescaleA   = (w2Rem + w2Beam[0] - w2Beam[1] + lambdaRoot)
    / (2. * w2Rem * xSum[0]) ;
  double rescaleB   = (w2Rem + w2Beam[1] - w2Beam[0] + lambdaRoot)
    / (2. * w2Rem * xSum[1]) ;

  // Construct energy and pz for remnants in first beam.
  for (int iRem = nSys; iRem < beamA.size(); ++iRem) {
    double pPos = rescaleA * beamA[iRem].x() * wPosRem;
    double pNeg = beamA[iRem].mT2() / pPos;
    beamA[iRem].e( 0.5 * (pPos + pNeg) );
    beamA[iRem].pz( 0.5 * (pPos - pNeg) );

    // Add these partons to the normal event record.
    int iNew = event.append( beamA[iRem].id(), 63, 1 + nOffset, 0, 0, 0,
      beamA[iRem].col(), beamA[iRem].acol(), beamA[iRem].p(),
      beamA[iRem].m() );
    beamA[iRem].iPos( iNew);
  }

  // Construct energy and pz for remnants in second beam.
  for (int iRem = nSys; iRem < beamB.size(); ++iRem) {
    double pNeg = rescaleB * beamB[iRem].x() * wNegRem;
    double pPos = beamB[iRem].mT2() / pNeg;
    beamB[iRem].e( 0.5 * (pPos + pNeg) );
    beamB[iRem].pz( 0.5 * (pPos - pNeg) );

    // Add these partons to the normal event record.
    int iNew = event.append( beamB[iRem].id(), 63, 2 + nOffset, 0, 0, 0,
      beamB[iRem].col(), beamB[iRem].acol(), beamB[iRem].p(),
      beamB[iRem].m() );
    beamB[iRem].iPos( iNew);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Special beam remnant kinematics for gamma-gamma collision with one
// remnant system, other created by ISR, and for Deeply Inelastic Scattering.
// Currently assumes unresolved lepton.

bool BeamRemnants::setOneRemnKinematics( Event& event, int beamOffset) {

  // Identify beams with and without remnant.
  int iBeamHad;
  if ( beamAPtr->isGamma() && beamBPtr->isGamma() )
    iBeamHad = beamAPtr->resolvedGamma() ? 1 : 2;
  else iBeamHad = beamAPtr->isHadron() ? 1 : 2;
  BeamParticle& beamHad   = (iBeamHad == 1) ? *beamAPtr : *beamBPtr;
  BeamParticle& beamOther = (iBeamHad == 2) ? *beamAPtr : *beamBPtr;

  // Beam offsets in case of gamma+gamma.
  int iBeamA = 1 + beamOffset;
  int iBeamB = 2 + beamOffset;

  // Identify remnant-side hadronic four-momentum and scattered lepton if DIS.
  int iLepScat = isDIS ? (beamOther[0].iPos() + 2) : -1;

  // Definition of scattered lepton more intricate in the presence of QED
  // radiation. For now, define highest energy lepton as scattered.
  int iLepOut = particleDataPtr->isLepton(beamAPtr->id())
              ? 1 + beamOffset + 2 + 2 : 2 + beamOffset + 2 + 2;
  if ( !beamOther.isGamma()
    && (!event[iLepOut].isLepton()
    || !event[iLepScat].isLepton()
    || !event[iLepScat].isFinal()
    || iLepScat > event.size()-1)) {
    double eMax = -1.0;
    for (int i = event.size()-1; i > 0 ; --i) {
      if ( event[i].isFinal()
        && event[iLepOut].id() == event[i].id()
        && event[i].e() > eMax) {
        iLepScat = i;
        eMax = event[i].e();
        break;
      }
    }
  }

  // In the presence of QED radiation from the beam particles, we should now
  // find beam particles and transform to CM frame of the final beam
  // configuration.
  if (beamAPtr->isLepton() && beamBPtr->isHadron())
    { iBeamA = (*beamAPtr)[0].iPos(); }
  if (beamBPtr->isLepton() && beamAPtr->isHadron())
    { iBeamB = (*beamBPtr)[0].iPos(); }
  double etot = (event[iBeamA].p() + event[iBeamB].p()).mCalc();
  RotBstMatrix MtoRest;
  MtoRest.bst( event[iBeamA].p() + event[iBeamB].p(), Vec4(0.,0.,0.,etot));
  if ( abs(event[iBeamA].pz() + event[iBeamB].pz()) > 1e-5)
    event.rotbst(MtoRest);

  Vec4 pHadScat;
  for (int i = 5 + beamOffset; i < event.size(); ++i)
    if (event[i].isFinal() && i != iLepScat) pHadScat += event[i].p();

  // Set scattered lepton momentum if DIS and find the hadronic system.
  Vec4 pLepScat = isDIS ? event[iLepScat].p() : Vec4();
  Vec4 pHadTot  = event[iBeamA].p() + event[iBeamB].p();
  if ( isDIS ) pHadTot -= pLepScat;
  Vec4 pRemnant = pHadTot - pHadScat;
  double w2Tot  = pHadTot.m2Calc();
  double w2Scat = pHadScat.m2Calc();

  // Find a boost to the hadronic rest frame.
  RotBstMatrix MtoHadRest;
  // Do not flip the direction with boosts.
  if (iBeamHad == 1) MtoHadRest.toCMframe( pRemnant, pHadScat);
  else MtoHadRest.toCMframe( pHadScat, pRemnant);

  // Boost the event to hadronic rest frame.
  event.rotbst( MtoHadRest);
  pHadScat.rotbst( MtoHadRest);

  // Set the variables for kT-generation.
  int nBeam = beamHad.size();
  vector<double> kTwidth(nBeam);
  vector<double> kTcomp(nBeam);
  // Currently no MPI's for one-remnant systems.
  int iSys         = 0;
  double kTcompSum = 0;
  if (doPrimordialKT) {

    // Calculate the kT width from hard process.
    int iInA         = partonSystemsPtr->getInA(iSys);
    int iInB         = partonSystemsPtr->getInB(iSys);
    double sHatNow   = (event[iInA].p() + event[iInB].p()).m2Calc();
    double mHat      = sqrt(sHatNow);
    double mHatDamp  = mHat / (mHat + halfMassForKT);
    double scale     = (iSys == 0) ? infoPtr->QRen(iDS)
      : partonSystemsPtr->getPTHat(iSys);
    kTwidth[0]       = ( (halfScaleForKT * primordialKTsoft
      + scale * primordialKThard) / (halfScaleForKT + scale) ) * mHatDamp;
    kTcomp[0]        = mHatDamp;
    kTcompSum = mHatDamp + nBeam - 1.;

    // Calculate the kT width for remnant partons.
    for ( int iRem = 1; iRem < nBeam; ++iRem) {
      kTwidth[iRem] = primordialKTremnant;
      kTcomp[iRem]  = 1.;
    }

  }

  // Allow ten tries to construct kinematics (but normally works first).
  bool isPhysical = true;
  double xSum, xInvM, w2Remn, wDiff;
  for (int iTry = 0; iTry < NTRYKINMATCH; ++iTry) {
    isPhysical = true;

    // Add primordial kT for one-remnant system.
    if (doPrimordialKT) {

      // Generate Gaussian pT for initiator partons inside hadron.
      double pxSum = 0.;
      double pySum = 0.;
      for (int iPar = 0; iPar < nBeam; ++iPar) {
        pair<double, double> gauss2 = rndmPtr->gauss2();
        double px = kTwidth[iPar] * gauss2.first;
        double py = kTwidth[iPar] * gauss2.second;
        beamHad[iPar].px(px);
        beamHad[iPar].py(py);
        pxSum += px;
        pySum += py;
      }

      // Share recoil between all initiator partons.
      for (int iPar = 0; iPar < nBeam; ++iPar) {
        beamHad[iPar].px( beamHad[iPar].px() - pxSum*kTcomp[iPar]/kTcompSum );
        beamHad[iPar].py( beamHad[iPar].py() - pySum*kTcomp[iPar]/kTcompSum );
      }
    }

    // Pick unrescaled x values for remnants. Sum up (unscaled) p+ and p-.
    xSum  = 0.;
    xInvM = 0.;
    for (int iRem = 1; iRem < beamHad.size(); ++iRem) {
      double xPrel = beamHad.xRemnant( iRem);
      beamHad[iRem].x(xPrel);
      xSum  += xPrel;
      xInvM += beamHad[iRem].mT2() / xPrel;
    }

    // Include the primordial kT for scattered system if any.
    double mT2Scat = w2Scat;
    if (doPrimordialKT) mT2Scat += pow2( beamHad[0].pT());

    // Squared transverse mass for remnant, may give failure.
    w2Remn = xSum * xInvM;
    wDiff  = sqrt(w2Tot) - sqrtpos(mT2Scat) - sqrtpos(w2Remn);
    if (wDiff < 0.) isPhysical = false;
    if (isPhysical) break;
  }

  if (!isPhysical) {
    // Boost back event.
    MtoHadRest.invert();
    event.rotbst( MtoHadRest);
    infoPtr->errorMsg("Error in BeamRemnants::setOneRemnKinematics:"
      " too big beam remnant invariant mass");
    return false;
  }

  // Calculate kinematics for the scattered and remnant system with kT.
  if (doPrimordialKT) {

    int iInA  = partonSystemsPtr->getInA(iSys);
    int iInB  = partonSystemsPtr->getInB(iSys);
    int iInit = (iBeamHad == 1) ? iInA : iInB;

    // Beam momentum for boosts.
    Vec4 pBeam = beamOther.p();

    // Calculate the energy and pz of initiator with added kT.
    double w2        = (beamOther.p() + event[iInit].p()).m2Calc();
    double m2        = pow2(beamHad.m());
    double mT2       = pow2(beamHad[iSys].pT()) + m2;
    double eInitNew  = beamOther.e()*mT2/(w2-m2) + 0.25*(w2-m2)/beamOther.e();
    double pzInitNew = beamOther.e()*mT2/(w2-m2) - 0.25*(w2-m2)/beamOther.e();
    if(iBeamHad == 1) pzInitNew *= -1.0;
    beamHad[iSys].pz( pzInitNew);
    beamHad[iSys].e( eInitNew);

    // Define the boost from event rest frame to frame where initiator
    // has some kT.
    RotBstMatrix MtoInitWithkT;

    // Find the rest frame of beam particle and initiator.
    if (iBeamHad == 1) {
      MtoInitWithkT.toCMframe( beamHad[iSys].p(), beamOther.p());
    } else {
      MtoInitWithkT.toCMframe( beamOther.p(), beamHad[iSys].p());
    }

    // Boost to frame where beam particle has its original momentum.
    pBeam.rotbst(MtoInitWithkT);
    MtoInitWithkT.bst( pBeam, beamOther.p());

    // Invert the boost.
    MtoInitWithkT.invert();

    // Add the initiator with kT to the event record.
    int iBeam = beamHad[iSys].iPos();
    int iCopy = event.copy(iBeam, -61);
    event[iCopy].rotbst(MtoInitWithkT);
    if (iBeamHad == 1) partonSystemsPtr->setInA(iSys, iCopy);
    else partonSystemsPtr->setInB(iSys, iCopy);
    beamHad[iSys].iPos( iCopy);
    if (iSys == 0) {
      int mother = event[iCopy].mother1();
      event[mother].daughter1(iCopy);
    }

    // Calculate transverse mass for scattered system with kT.
    w2Scat += pow2( beamHad[iSys].pT());
  }

  // Boost of scattered system to compensate for remnant mass.
  double lambda   = pow2( w2Tot - w2Scat - w2Remn) - 4. * w2Scat * w2Remn;
  double pzNew    = 0.5 * sqrt( lambda / w2Tot);
  // Keep the direction of momentum to ensure stable boost.
  if(iBeamHad == 1) pzNew *= -1.0;
  double eNewScat = 0.5 * (w2Tot + w2Scat - w2Remn) / sqrt(w2Tot);
  Vec4 pNewScat( beamHad[iSys].px(), beamHad[iSys].py(), pzNew, eNewScat);
  RotBstMatrix MforScat;
  MforScat.bst( pHadScat, pNewScat);
  int sizeSave = event.size();
  for (int i = 5 + beamOffset; i < sizeSave; ++i)
    if ( i != iLepScat && event[i].isFinal() ) {
      int iNew = event.copy( i, 62);
      event[iNew].rotbst( MforScat, false);
  }

  // Calculate kinematics of remnants and insert into event record.
  double eNewRemn = 0.5 * (w2Tot + w2Remn - w2Scat) / sqrt(w2Tot);
  double wNewRemn = eNewRemn + pzNew;
  for (int iRem = 1; iRem < beamHad.size(); ++iRem) {
    double wNegNow = wNewRemn * beamHad[iRem].x() / xSum;
    double wPosNow = beamHad[iRem].mT2() / wNegNow;
    beamHad[iRem].pz(-0.5 * (wNegNow - wPosNow));
    beamHad[iRem].e ( 0.5 * (wPosNow + wNegNow));
    int iNew = event.append( beamHad[iRem].id(), 63, iBeamHad + beamOffset,
      0, 0, 0, beamHad[iRem].col(), beamHad[iRem].acol(), beamHad[iRem].p(),
      beamHad[iRem].m() );
    beamHad[iRem].iPos( iNew);
  }

  // Boost back event.
  MtoHadRest.invert();
  event.rotbst( MtoHadRest);

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Collapse colours and check that they are consistent.

bool BeamRemnants::checkColours( Event& event) {

  // No colours in lepton beams so no need to do anything.
  if (beamAPtr->isLepton() && beamBPtr->isLepton()) return true;

  // Remove ambiguities when one colour collapses two ways.
  // Resolve chains where one colour is mapped to another.
  for (int iCol = 1; iCol < int(colFrom.size()); ++iCol)
  for (int iColRef = 0; iColRef < iCol; ++iColRef) {
    if (colFrom[iCol] == colFrom[iColRef]) {
      colFrom[iCol] = colTo[iCol];
      colTo[iCol] = colTo[iColRef];
    }
    if (colTo[iCol] == colFrom[iColRef]) colTo[iCol] = colTo[iColRef];
  }

  // Transform event record colours from beam remnant colour collapses.
  for (int i = oldSize; i < event.size(); ++i) {
    int col = event[i].col();
    int acol = event[i].acol();
    for (int iCol = 0; iCol < int(colFrom.size()); ++iCol) {
      if (col == colFrom[iCol]) {col = colTo[iCol]; event[i].col(col);}
      if (acol == colFrom[iCol]) {acol = colTo[iCol]; event[i].acol(acol);}
      // Sextets have extra, negative, tags.
      if (col == -colFrom[iCol]) {col = -colTo[iCol]; event[i].col(col);}
      if (acol == -colFrom[iCol]) {acol = -colTo[iCol]; event[i].acol(acol);}
    }
  }

  // Transform junction colours from beam remnant colour collapses.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)
    for (int leg = 0; leg < 3; ++leg) {
      int col = event.colJunction(iJun, leg);
      for (int iCol = 0; iCol < int(colFrom.size()); ++iCol) {
        if (col == colFrom[iCol]) {
          col = colTo[iCol];
          event.colJunction(iJun, leg, col);
        }
      }
    }

  // Arrays for current colours and anticolours, and for singlet gluons.
  vector<int> colList;
  vector<int> acolList;
  vector<int> iSingletGluon;

  // Find current colours and anticolours in the event record.
  for (int i = oldSize; i < event.size(); ++i)
  if (event[i].isFinal()) {
    int id   = event[i].id();
    int col  = event[i].col();
    int acol = event[i].acol();
    int colType = event[i].colType();

    // Quarks must have colour set, antiquarks anticolour, gluons both.
    if ( (id > 0 && id < 9 && (col <= 0 || acol != 0) )
      || (id < 0 && id > -9 && (col != 0 || acol <= 0) )
      || (id == 21 && (col <= 0 || acol <= 0) ) ) {
      infoPtr->errorMsg("Error in BeamRemnants::checkColours: "
        "q/qbar/g has wrong colour slots set");
      return false;
    }

    // Sextets must have one positive and one negative tag
    if ( (colType == 3 && (col <= 0 || acol >= 0))
         || (colType == -3 && (col >= 0 || acol <= 0)) ) {
      infoPtr->errorMsg("Error in BeamRemnants::checkColours: "
                        "sextet has wrong colours");
    }

    // Save colours/anticolours, and position of colour singlet gluons.
    if ( col > 0)  colList.push_back(  col );
    if (acol > 0) acolList.push_back( acol );
    if (col > 0 && acol == col) iSingletGluon.push_back(i);
    // Colour sextets
    if ( col < 0) acolList.push_back( -col );
    if (acol < 0) colList.push_back( -acol );
  }

  // Run through list of singlet gluons and put them on final-state dipole
  // (i,j) that offers smallest (p_g p_i) * (p_g p_j) / (p_i p_j).
  for (int iS = 0; iS < int(iSingletGluon.size()); ++iS) {
    int    iGlu      = iSingletGluon[iS];
    int    iAcolDip  = -1;
    double pT2DipMin = sCM;
    for (int iC = oldSize; iC < event.size(); ++iC)
      if (iC != iGlu && event[iC].isFinal()) {
      int colDip = event[iC].col();
      if (colDip > 0 && event[iC].acol() !=colDip)
      for (int iA = oldSize; iA < event.size(); ++iA)
        if (iA != iGlu && iA != iC && event[iA].isFinal()
        && event[iA].acol() == colDip && event[iA].col() !=colDip) {
        double pT2Dip = (event[iGlu].p() * event[iC].p())
          * (event[iGlu].p() * event[iA].p())
          / (event[iC].p() * event[iA].p());
        if (pT2Dip < pT2DipMin) {
          iAcolDip  = iA;
          pT2DipMin = pT2Dip;
        }
      }
    }

    // Fail if no dipole. Else insert singlet gluon onto relevant dipole.
    if (iAcolDip == -1)  return false;
    event[iGlu].acol( event[iAcolDip].acol() );
    event[iAcolDip].acol( event[iGlu].col() );

    // Update any junction legs that match reconnected dipole.
    for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {

      // Only junctions need to be updated, not antijunctions.
      if (event.kindJunction(iJun) % 2 == 0) continue;
      for (int leg = 0; leg < 3; ++leg) {
        int col = event.colJunction(iJun, leg);
        if (col == event[iGlu].acol())
          event.colJunction(iJun, leg, event[iGlu].col());
      }
    }

  }

  // Check that not the same colour or anticolour appears twice.
  for (int iCol = 0; iCol < int(colList.size()) - 1; ++iCol) {
    int col = colList[iCol];
    for (int iCol2 = iCol + 1; iCol2 < int(colList.size()); ++iCol2)
    if (colList[iCol2] == col) {
      infoPtr->errorMsg("Warning in BeamRemnants::checkColours:"
        " colour appears twice");
      if (!ALLOWCOLOURTWICE) return false;
    }
  }
  for (int iAcol = 0; iAcol < int(acolList.size()) - 1; ++iAcol) {
    int acol = acolList[iAcol];
    for (int iAcol2 = iAcol + 1; iAcol2 < int(acolList.size()); ++iAcol2)
    if (acolList[iAcol2] == acol) {
      infoPtr->errorMsg("Warning in BeamRemnants::checkColours:"
        " anticolour appears twice");
      if (!ALLOWCOLOURTWICE) return false;
    }
  }

  // Remove all matching colour-anticolour pairs.
  bool foundPair = true;
  while (foundPair && colList.size() > 0 && acolList.size() > 0) {
    foundPair = false;
    for (int iCol = 0; iCol < int(colList.size()); ++iCol) {
      for (int iAcol = 0; iAcol < int(acolList.size()); ++iAcol) {
        if (acolList[iAcol] == colList[iCol]) {
          colList[iCol] = colList.back(); colList.pop_back();
          acolList[iAcol] = acolList.back(); acolList.pop_back();
          foundPair = true;
          break;
        }
      }
      if (foundPair) break;
    }
  }

  // Check that remaining (anti)colours are accounted for by junctions.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
    int kindJun = event.kindJunction(iJun);
    for (int leg = 0; leg < 3; ++leg) {
      int colEnd = event.colJunction(iJun, leg);

      // Junction connected to three colours.
      if (kindJun == 1) {
        for (int iCol = 0; iCol < int(colList.size()); ++iCol)
        if (colList[iCol] == colEnd) {
          // Found colour match: remove and exit.
          colList[iCol] = colList.back();
          colList.pop_back();
          break;
        }
      }

      // Junction connected to three anticolours.
      else if (kindJun == 2) {
        for (int iAcol = 0; iAcol < int(acolList.size()); ++iAcol)
        if (acolList[iAcol] == colEnd) {
          // Found colour match: remove and exit.
          acolList[iAcol] = acolList.back();
          acolList.pop_back();
          break;
        }
      }

      // Other junction types
      else if ( kindJun == 3 || kindJun == 5) {
        for (int iCol = 0; iCol < int(colList.size()); ++iCol)
        if (colList[iCol] == colEnd) {
          // Found colour match: remove and exit.
          colList[iCol] = colList.back();
          colList.pop_back();
          break;
        }
      }

      // Other antijunction types
      else if ( kindJun == 4 || kindJun == 6) {
        for (int iAcol = 0; iAcol < int(acolList.size()); ++iAcol)
        if (acolList[iAcol] == colEnd) {
          // Found colour match: remove and exit.
          acolList[iAcol] = acolList.back();
          acolList.pop_back();
          break;
        }
      }

      // End junction check.
    }
  }


  // Repair step - sometimes needed when rescattering allowed.
  if (colList.size() > 0 || acolList.size() > 0) {
    infoPtr->errorMsg("Warning in BeamRemnants::checkColours:"
                      " need to repair unmatched colours");
  }
  while (colList.size() > 0 && acolList.size() > 0) {

    // Replace one colour and one anticolour index by a new common one.
    int  colMatch =  colList.back();
    int acolMatch = acolList.back();
    int  colNew   = event.nextColTag();
    colList.pop_back();
    acolList.pop_back();
    for (int i = oldSize; i < event.size(); ++i) {
      if (event[i].isFinal() && event[i].col() == colMatch) {
        event[i].col( colNew);
        break;
      }
      else if (event[i].isFinal() && event[i].acol() == -colMatch) {
        event[i].acol( -colNew);
        break;
      }
    }
    for (int i = oldSize; i < event.size(); ++i) {
      if (event[i].isFinal() && event[i].acol() == acolMatch) {
        event[i].acol( colNew);
        break;
      }
      if (event[i].isFinal() && event[i].col() == -acolMatch) {
        event[i].col( -colNew);
        break;
      }
    }
  }

  // Done.
  return (colList.size() == 0 && acolList.size() == 0);

}

//--------------------------------------------------------------------------

// Update colours of outgoing particles in the event record.

void BeamRemnants::updateColEvent( Event& event,
  vector<pair <int,int> > colChanges) {

  for (int iCol = 0; iCol < int(colChanges.size()); ++iCol) {

    int oldCol = colChanges[iCol].first;
    int newCol = colChanges[iCol].second;
    if (oldCol == newCol)
      continue;

    // Add a copy of final particles with old colour and change the colour.
    for (int j = 0; j < event.size(); ++j) {
      if (event[j].isFinal() && event[j].col() == oldCol)
        event[event.copy(j, 64)].col(newCol);
      if (event[j].isFinal() && event[j].acol() == -oldCol)
        event[event.copy(j, 64)].acol(-newCol);

      if (event[j].isFinal() && event[j].acol() == oldCol)
        event[event.copy(j,64)].acol(newCol);
      if (event[j].isFinal() && event[j].col() == -oldCol)
        event[event.copy(j,64)].col(-newCol);
    }

    // Update junction.
    for (int j = 0;j < event.sizeJunction(); ++j)
      for (int jCol = 0; jCol < 3; ++jCol)
        if (event.colJunction(j,jCol) == oldCol)
          event.colJunction(j,jCol,newCol);
  }

}

//==========================================================================

} // end namespace Pythia8
