// BeamRemnants.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// BeamRemnants class.

#include "BeamRemnants.h"

namespace Pythia8 {

//==========================================================================

// The BeamDipole class is purely internal to reconnectColours.

class BeamDipole {

public:

  // Constructor.
  BeamDipole( int colIn = 0, int iColIn = 0, int iAcolIn = 0) 
    : col(colIn), iCol(iColIn), iAcol(iAcolIn) {}  

  // Members.
  int    col, iCol, iAcol;
  double p1p2;
 
};

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
  PartonSystems* partonSystemsPtrIn) {

  // Save pointers.
  infoPtr             = infoPtrIn;
  rndmPtr             = rndmPtrIn;
  beamAPtr            = beamAPtrIn;
  beamBPtr            = beamBPtrIn;
  partonSystemsPtr    = partonSystemsPtrIn;

  // Width of primordial kT distribution.
  doPrimordialKT      = settings.flag("BeamRemnants:primordialKT");
  primordialKTsoft    = settings.parm("BeamRemnants:primordialKTsoft");
  primordialKThard    = settings.parm("BeamRemnants:primordialKThard");
  primordialKTremnant = settings.parm("BeamRemnants:primordialKTremnant");
  halfScaleForKT      = settings.parm("BeamRemnants:halfScaleForKT");
  halfMassForKT       = settings.parm("BeamRemnants:halfMassForKT");

  // Handling of rescattering kinematics uncertainties from primodial kT.
  allowRescatter    = settings.flag("MultipartonInteractions:allowRescatter");
  doRescatterRestoreY = settings.flag("BeamRemnants:rescatterRestoreY");

  // Parameters for colour reconnection scenario, partly borrowed from
  // multiparton interactions not to introduce too many new ones.
  doReconnect         = settings.flag("BeamRemnants:reconnectColours");
  reconnectRange      = settings.parm("BeamRemnants:reconnectRange");
  pT0Ref              = settings.parm("MultipartonInteractions:pT0Ref");
  ecmRef              = settings.parm("MultipartonInteractions:ecmRef");
  ecmPow              = settings.parm("MultipartonInteractions:ecmPow");

  // Total and squared CM energy at nominal energy.
  eCM                 = infoPtr->eCM();
  sCM                 = eCM * eCM;

  // The MPI pT0 smoothening scale and its reconnection-strength combination.
  pT0                 = pT0Ref * pow(eCM / ecmRef, ecmPow);
  pT20Rec             = pow2(reconnectRange * pT0); 
  
  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Select the flavours/kinematics/colours of the two beam remnants. 
// Notation: iPar = all partons, iSys = matched systems of two beams,
//           iRem = additional partons in remnants.

bool BeamRemnants::add( Event& event) {

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

  // Number of scattering subsystems. Size of event record before treatment.
  nSys    = partonSystemsPtr->sizeSys();
  oldSize = event.size();

  // Add required extra remnant flavour content. 
  // Start all over if fails (in option where junctions not allowed).
  if ( !beamAPtr->remnantFlavours(event) 
    || !beamBPtr->remnantFlavours(event) ) {
    infoPtr->errorMsg("Error in BeamRemnants::add:"
      " remnant flavour setup failed"); 
    return false;
  }

  // Do the kinematics of the collision subsystems and two beam remnants.
  if (!setKinematics(event)) return false;

  // Allow colour reconnections.
  if (doReconnect) reconnectColours(event);

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

  // Check that has not already used up beams.
  if ( (!beamA.isLepton() && beamA.xMax(-1) <= 0.)
    || (!beamB.isLepton() && beamB.xMax(-1) <= 0.) ) {
    infoPtr->errorMsg("Error in BeamRemnants::setKinematics:"
      " no momentum left for beam remnants"); 
    return false;
  }

  // Last beam-status particles. Offset relative to normal beam locations.
  int nBeams   = 3;
  for (int i = 3; i < event.size(); ++i) 
    if (event[i].statusAbs() < 20) nBeams = i + 1; 
  int nOffset  = nBeams - 3; 

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
      mHatDamp        = mHat / (mHat + halfMassForKT);
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
      if (beam.isHadron() && doPrimordialKT) { 
        double pxSum = 0.;
        double pySum = 0.;
        for (int iPar = 0; iPar < nPar; ++iPar) { 
          if ( beam[iPar].isFromBeam() ) {
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
        for (int iPar = 0; iPar < nPar; ++iPar) 
        if (beam[iPar].isFromBeam() ) {
          beam[iPar].px( beam[iPar].px() - pxSum * kTcomp[iPar] / kTcompSum );
          beam[iPar].py( beam[iPar].py() - pySum * kTcomp[iPar] / kTcompSum );
        }

      // Without primordial kT: still need to store rescattered partons.
      } else if (beam.isHadron()) {
        for (int iPar = 0; iPar < nPar; ++iPar)  
        if ( !beam[iPar].isFromBeam() ) {
          int iInAB = (iBeam == 0) ? partonSystemsPtr->getInA(iPar)    
                                   : partonSystemsPtr->getInB(iPar);       
          beam[iPar].p( event[iInAB].p() );
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
      event[iAcopy].rotbst(Msys[iSys]);
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
      event[iBcopy].rotbst(Msys[iSys]);
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
        event[iABcopy].rotbst(Msys[iSys]); 
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
      if (event[iAB].isFinal()) event[iAB].rotbst(Mmismatch);      
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

// Allow colour reconnections by mergings of collision subsystems.
// iRec is system that may be reconnected, by moving its gluons to iSys,   
// where minimal pT (or equivalently Lambda) is used to pick location.
// Therefore all dipoles in iSys have to be found, and all gluons in iRec.
// Matching q-qbar pairs are treated by analogy with gluons.
// Note: owing to rescatterings some outgoing partons must be skipped.

bool BeamRemnants::reconnectColours( Event&  event) {

  // References to beams to simplify indexing.
  BeamParticle& beamA = *beamAPtr;  
  BeamParticle& beamB = *beamBPtr;  

  // Prepare record of which systems should be merged onto another.
  // The iSys system must have colour in final state to attach to it.
  vector<int>  iMerge(nSys);
  vector<bool> hasColour(nSys);
  for (int iSys = 0; iSys < nSys; ++iSys) {
    iMerge[iSys] = iSys;
    bool hasCol = false;
    for (int iMem = 0; iMem < partonSystemsPtr->sizeOut(iSys); ++iMem) {
      int iNow = partonSystemsPtr->getOut( iSys, iMem);
      if (event[iNow].isFinal() && (event[iNow].col() > 0 
        || event[iNow].acol() > 0) ) {
        hasCol = true;
        break;
      }    
    }
    hasColour[iSys] = hasCol;
  }

  // Loop over systems to decide which should be reconnected. 
  for (int iRec = nSys - 1; iRec > 0; --iRec) {

    // Determine reconnection strength from pT scale of system.
    double pT2Rec  = pow2( partonSystemsPtr->getPTHat(iRec) );
    double probRec = pT20Rec / (pT20Rec + pT2Rec); 

    // Loop over other systems iSys at higher pT scale and 
    // decide whether to reconnect the iRec gluons onto one of them.
    for (int iSys = iRec - 1; iSys >= 0; --iSys)
    if (hasColour[iSys] && probRec > rndmPtr->flat()) { 

      // The iRec system and all merged with it to be merged with iSys.
      iMerge[iRec] = iSys;       
      for (int iRec2 = iRec + 1; iRec2 < nSys; ++iRec2)
      if (iMerge[iRec2] == iRec) iMerge[iRec2] = iSys;    

      // Once a system has been merged do not test it anymore.
      break;
    }
  }

  // Loop over systems. Check whether other systems to be merged with it.
  for (int iSys = 0; iSys < nSys; ++iSys) {
    int nMerge = 0;
    for (int iRec = iSys + 1; iRec < nSys; ++iRec)
    if (iMerge[iRec] == iSys) ++nMerge;
    if (nMerge == 0) continue; 

    // Incoming partons not counted if rescattered.
    int  iInASys = partonSystemsPtr->getInA(iSys);
    bool hasInA  = (beamA[iSys].isFromBeam());   
    int  iInBSys = partonSystemsPtr->getInB(iSys);    
    bool hasInB  = (beamB[iSys].isFromBeam()); 

    // Begin find dipoles in iSys system.
    vector<BeamDipole> dipoles;
    int sizeOut = partonSystemsPtr->sizeOut(iSys);
    for (int iMem = 0; iMem < sizeOut; ++iMem) {

      // Find colour dipoles to beam remnant.
      int iNow = partonSystemsPtr->getOut( iSys, iMem);
      if (!event[iNow].isFinal()) continue;
      int col = event[iNow].col();  
      if (col > 0) {
        if      (hasInA && event[iInASys].col() == col)
          dipoles.push_back( BeamDipole( col, iNow, iInASys ) );
        else if (hasInB && event[iInBSys].col() == col)
          dipoles.push_back( BeamDipole( col, iNow, iInBSys ) );
 
        // Find colour dipole between final-state partons.
        else for (int iMem2 = 0; iMem2 < sizeOut; ++iMem2) 
        if (iMem2 != iMem) {
          int iNow2 = partonSystemsPtr->getOut( iSys, iMem2); 
          if (!event[iNow2].isFinal()) continue;
          if (event[iNow2].acol() == col) {
            dipoles.push_back( BeamDipole( col, iNow, iNow2) );
            break;
          }
        }
      }

      // Find anticolour dipoles to beam remnant.
      int acol = event[iNow].acol();  
      if (acol > 0) {
        if      (hasInA && event[iInASys].acol() == acol)
          dipoles.push_back( BeamDipole( acol, iInASys, iNow ) );
        else if (hasInB && event[iInBSys].acol() == acol)
          dipoles.push_back( BeamDipole( acol, iInBSys, iNow ) ); 
      }
    }
   
    // Skip mergings if no dipoles found.
    if (dipoles.size() == 0) continue; 

    // Find dipole sizes.
    for (int iDip = 0; iDip < int(dipoles.size()); ++iDip) 
      dipoles[iDip].p1p2 = event[dipoles[iDip].iCol].p() 
                         * event[dipoles[iDip].iAcol].p();
    
    // Loop over systems iRec to be merged with iSys.
    for (int iRec = iSys + 1; iRec < nSys; ++iRec) {
      if (iMerge[iRec] != iSys) continue;

      // Information on iRec. Vectors for gluons and anything else.
      int sizeRec = partonSystemsPtr->sizeOut(iRec);
      int iInARec = partonSystemsPtr->getInA(iRec);
      int iInBRec = partonSystemsPtr->getInB(iRec);   
      int nGluRec = 0;
      vector<int>    iGluRec;
      vector<double> pT2GluRec;
      int nAnyRec = 0;
      vector<int>    iAnyRec;
      vector<bool>   freeAnyRec;

      // Copy of gluon positions in descending order. 
      for (int iMem = 0; iMem < sizeRec; ++iMem) {
        int iNow = partonSystemsPtr->getOut( iRec, iMem);
        if (!event[iNow].isFinal()) continue;
        if (event[iNow].isGluon()) {
          ++nGluRec;
          iGluRec.push_back( iNow );  
          pT2GluRec.push_back( event[iNow].pT2() );
          for (int i = nGluRec - 1; i > 1; --i) {
            if (pT2GluRec[i - 1] > pT2GluRec[i]) break;
            swap(   iGluRec[i - 1],   iGluRec[i] );    
            swap( pT2GluRec[i - 1], pT2GluRec[i] ); 
          }  
        // Copy of anything else, mainly quarks, in no particular order. 
        } else {
          ++nAnyRec;
          iAnyRec.push_back( iNow ); 
          freeAnyRec.push_back( true ); 
        }
      }

      // For each gluon in iRec now find the dipole that gives the smallest
      // (pGlu * pI) (pGlu * pJ) / (pI * pJ), i.e. minimal pT (and Lambda). 
      for (int iGRec = 0; iGRec < nGluRec; ++iGRec) {
        int    iGlu      = iGluRec[iGRec];
        Vec4   pGlu      = event[iGlu].p();
        int    iDipMin   = 0;
        double pT2DipMin = sCM;
        for (int iDip = 0; iDip < int(dipoles.size()); ++iDip) {
          double pT2Dip = (pGlu * event[dipoles[iDip].iCol].p())
            * (pGlu * event[dipoles[iDip].iAcol].p()) / dipoles[iDip].p1p2;
          if (pT2Dip < pT2DipMin) {
            iDipMin   = iDip;
            pT2DipMin = pT2Dip;
          }
        }  

        // Attach the gluon to the dipole, i.e. split the dipole in two.
        int colGlu   = event[iGlu].col();
        int acolGlu  = event[iGlu].acol();
        int colDip   = dipoles[iDipMin].col;
        int iColDip  = dipoles[iDipMin].iCol;
        int iAcolDip = dipoles[iDipMin].iAcol;
        event[iGlu].acol( colDip );
        if (event[iAcolDip].acol() == colDip) 
             event[iAcolDip].acol( colGlu );
        else event[iAcolDip].col(  colGlu ); 
        dipoles[iDipMin].iAcol = iGlu;
        dipoles[iDipMin].p1p2 = event[iColDip].p() * pGlu;
        dipoles.push_back( BeamDipole( colGlu, iGlu, iAcolDip ) );
        dipoles.back().p1p2 = pGlu * event[iAcolDip].p();
     
        // Remove gluon from old system: reconnect colours.
        for (int i = oldSize; i < event.size(); ++i)
        if (i != iGlu && i != iAcolDip) { 
          if (event[i].isFinal()) {     
            if (event[i].acol() == colGlu) event[i].acol( acolGlu ); 
          } else {      
              if (event[i].col()  == colGlu) event[i].col( acolGlu );  
          }       
        }

	// Update any junction legs that match reconnected dipole.
	for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {

	  // Only junctions need to be updated, not antijunctions.
	  if (event.kindJunction(iJun) % 2 == 0) continue;
	  for (int leg = 0; leg < 3; ++leg) {
	    int col = event.colJunction(iJun, leg); 
	    if (col == colDip) 
	      event.colJunction(iJun, leg, colGlu);
	  }	
	}
	
      }

      // See if any matching quark-antiquark pairs among the rest.
      for (int iQRec = 0; iQRec < nAnyRec; ++iQRec) {
        int iQ  = iAnyRec[iQRec];
        int idQ = event[iQ].id();
        if (freeAnyRec[iQRec] && idQ > 0 && idQ < 6) 
        for (int iQbarRec = 0; iQbarRec < nAnyRec; ++iQbarRec) {
          int iQbar  = iAnyRec[iQbarRec];
          if (freeAnyRec[iQbarRec] && event[iQbar].id() == -idQ) {

            // Check that these can be traced back to same gluon splitting.
            // For now also avoid qqbar pairs produced in rescatterings.??
            int iTopQ    = event.iTopCopyId(iQ);
            int iTopQbar = event.iTopCopyId(iQbar);
            int iMother  = event[iTopQ].mother1();
            if (event[iTopQbar].mother1() == iMother
              && event[iMother].isGluon() && event[iMother].status() != -34
              && event[iMother + 1].status() != -34 ) {

              // Now find the dipole that gives the smallest
              // ((pQ + pQbar) * pI) ((pQ + pQbar) * pJ) / (pI * pJ). 
              Vec4   pGlu      = event[iQ].p() + event[iQbar].p();
              int    iDipMin   = 0;
              double pT2DipMin = sCM;
              for (int iDip = 0; iDip < int(dipoles.size()); ++iDip) {
                double pT2Dip = (pGlu * event[dipoles[iDip].iCol].p())
                  * (pGlu * event[dipoles[iDip].iAcol].p()) 
                  / dipoles[iDip].p1p2;
                if (pT2Dip < pT2DipMin) {
                  iDipMin   = iDip;
                  pT2DipMin = pT2Dip;
                }
              }  

              // Attach the q-qbar pair to the dipole, i.e. split the dipole.
              int colGlu   = event[iQ].col();
              int acolGlu  = event[iQbar].acol();
              int colDip   = dipoles[iDipMin].col;
              int iColDip  = dipoles[iDipMin].iCol;
              int iAcolDip = dipoles[iDipMin].iAcol;
              event[iQbar].acol( colDip );
              if (event[iAcolDip].acol() == colDip) 
                   event[iAcolDip].acol( colGlu );
              else event[iAcolDip].col(  colGlu ); 
              dipoles[iDipMin].iAcol = iQbar;
              dipoles[iDipMin].p1p2 = event[iColDip].p() * event[iQbar].p();
              dipoles.push_back( BeamDipole( colGlu, iQ, iAcolDip ) );
              dipoles.back().p1p2 = event[iQ].p() * event[iAcolDip].p();
     
              // Remove q-qbar pair from old system: reconnect colours. 
              freeAnyRec[iQRec]    = false;
              freeAnyRec[iQbarRec] = false;
              for (int i = oldSize; i < event.size(); ++i)
              if (i != iQRec && i != iQbarRec && i != iColDip 
                && i != iAcolDip) { 
                if (event[i].isFinal()) {     
                  if (event[i].acol() == colGlu) event[i].acol( acolGlu ); 
                } else {      
                    if (event[i].col()  == colGlu) event[i].col( acolGlu );  
                }       
              }
               
	      // Update any junction legs that match reconnected dipole.
	      for (int iJun = 0; iJun < event.sizeJunction(); ++iJun) {
		
		// Only junctions need to be updated, not antijunctions.
		if (event.kindJunction(iJun) % 2 == 0) continue;
		for (int leg = 0; leg < 3; ++leg) {
		  int col = event.colJunction(iJun, leg); 
		  if (col == colDip) 
		    event.colJunction(iJun, leg, colGlu);
		}	
	      }
	      
	      // Done with processing of q-qbar pairs.
            }
          }
        }
      }

      // If only two beam gluons left of system, set their colour = anticolour.
      // Used by BeamParticle::remnantColours to skip irrelevant gluons.
      if ( event[iInARec].isGluon() && !event[iInARec].isRescatteredIncoming()
        && event[iInBRec].isGluon() && !event[iInBRec].isRescatteredIncoming()
        && event[iInARec].col() == event[iInBRec].acol() 
        && event[iInARec].acol() == event[iInBRec].col() ) { 
          event[iInARec].acol( event[iInARec].col() );
          event[iInBRec].acol( event[iInBRec].col() );
      }

    // End of loops over iRec and iSys systems.
    }
  }

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

  // Run though list of singlet gluons and put them on final-state dipole
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
    for (int i = oldSize; i < event.size(); ++i) 
    if (event[i].isFinal() && event[i].col() == colMatch) {
      event[i].col( colNew);
      break;
    }
    for (int i = oldSize; i < event.size(); ++i) 
    if (event[i].isFinal() && event[i].acol() == acolMatch) {
      event[i].acol( colNew);
      break;
    }
  }

  // Done.
  return (colList.size() == 0 && acolList.size() == 0);

}

//==========================================================================

} // end namespace Pythia8
