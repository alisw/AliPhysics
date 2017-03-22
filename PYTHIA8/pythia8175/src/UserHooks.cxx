// UserHooks.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the UserHooks class.

// Note: compilation crashes if PhaseSpace.h is moved to UserHooks.h.
#include "PhaseSpace.h"
#include "UserHooks.h"

namespace Pythia8 {
 
//==========================================================================

// The UserHooks class.

//--------------------------------------------------------------------------

// multiplySigmaBy allows the user to introduce a multiplicative factor 
// that modifies the cross section of a hard process. Since it is called
// from before the event record is generated in full, the normal analysis
// does not work. The code here provides a rather extensive summary of
// which methods actually do work. It is a convenient starting point for 
// writing your own derived routine.

double UserHooks::multiplySigmaBy( const SigmaProcess* sigmaProcessPtr, 
  const PhaseSpace* phaseSpacePtr, bool inEvent) {

  // Process code, necessary when some to be treated differently.
  //int code       = sigmaProcessPtr->code();

  // Final multiplicity, i.e. whether 2 -> 1 or 2 -> 2.
  //int nFinal     = sigmaProcessPtr->nFinal();

  // Incoming x1 and x2 to the hard collision, and factorization scale.
  //double x1      = phaseSpacePtr->x1();
  //double x2      = phaseSpacePtr->x2();
  //double Q2Fac   = sigmaProcessPtr->Q2Fac();

  // Renormalization scale and assumed alpha_strong and alpha_EM.
  //double Q2Ren   = sigmaProcessPtr->Q2Ren();
  //double alphaS  = sigmaProcessPtr->alphaSRen();
  //double alphaEM = sigmaProcessPtr->alphaEMRen();
  
  // Subprocess mass-square.
  //double sHat = phaseSpacePtr->sHat();

  // Now methods only relevant for 2 -> 2.
  //if (nFinal == 2) {
    
    // Mandelstam variables and hard-process pT.
    //double tHat  = phaseSpacePtr->tHat();
    //double uHat  = phaseSpacePtr->uHat();
    //double pTHat = phaseSpacePtr->pTHat();
  
    // Masses of the final-state particles. (Here 0 for light quarks.)
    //double m3    = sigmaProcessPtr->m(3);
    //double m4    = sigmaProcessPtr->m(4);
  //}

  // Dummy statement to avoid compiler warnings.
  return ((inEvent && sigmaProcessPtr->code() == 0 
    && phaseSpacePtr->sHat() < 0.) ? 0. : 1.);

}

//--------------------------------------------------------------------------

// biasSelectionBy allows the user to introduce a multiplicative factor 
// that modifies the cross section of a hard process. The event is assigned
// a wegith that is the inverse of the selection bias, such that the
// cross section is unchanged. Since it is called from before the 
// event record is generated in full, the normal analysis does not work. 
// The code here provides a rather extensive summary of which methods 
// actually do work. It is a convenient starting point for writing 
// your own derived routine.

double UserHooks::biasSelectionBy( const SigmaProcess* sigmaProcessPtr, 
  const PhaseSpace* phaseSpacePtr, bool inEvent) {

  // Process code, necessary when some to be treated differently.
  //int code       = sigmaProcessPtr->code();

  // Final multiplicity, i.e. whether 2 -> 1 or 2 -> 2.
  //int nFinal     = sigmaProcessPtr->nFinal();

  // Incoming x1 and x2 to the hard collision, and factorization scale.
  //double x1      = phaseSpacePtr->x1();
  //double x2      = phaseSpacePtr->x2();
  //double Q2Fac   = sigmaProcessPtr->Q2Fac();

  // Renormalization scale and assumed alpha_strong and alpha_EM.
  //double Q2Ren   = sigmaProcessPtr->Q2Ren();
  //double alphaS  = sigmaProcessPtr->alphaSRen();
  //double alphaEM = sigmaProcessPtr->alphaEMRen();
  
  // Subprocess mass-square.
  //double sHat = phaseSpacePtr->sHat();

  // Now methods only relevant for 2 -> 2.
  //if (nFinal == 2) {
    
    // Mandelstam variables and hard-process pT.
    //double tHat  = phaseSpacePtr->tHat();
    //double uHat  = phaseSpacePtr->uHat();
    //double pTHat = phaseSpacePtr->pTHat();
  
    // Masses of the final-state particles. (Here 0 for light quarks.)
    //double m3    = sigmaProcessPtr->m(3);
    //double m4    = sigmaProcessPtr->m(4);
  //}

  // Insert here your calculation of the selection bias. 
  // Here illustrated by a weighting up of events at high pT.
  //selBias = pow4(phaseSpacePtr->pTHat()); 

  // Return the selBias weight. 
  // Warning: if you use another variable than selBias
  // the compensating weight will not be set correctly.
  //return selBias;

  // Dummy statement to avoid compiler warnings.
  return ((inEvent && sigmaProcessPtr->code() == 0 
    && phaseSpacePtr->sHat() < 0.) ? 0. : 1.);
}

//--------------------------------------------------------------------------

// omitResonanceDecays omits resonance decay chains from process record.

void UserHooks::omitResonanceDecays(const Event& process, bool finalOnly) {

  // Reset work event to be empty
  workEvent.clear(); 

  // Loop through all partons. Beam particles should be copied.
  for (int i = 0; i < process.size(); ++i) {
    bool doCopy  = false;
    bool isFinal = false;
    if (i < 3) doCopy = true;

    // Daughters of beams should normally be copied.
    else {
      int iMother = process[i].mother1();
      if (iMother == 1 || iMother == 2) doCopy = true;
       
      // Granddaughters of beams should normally be copied and are final.
      else if (iMother > 2) {
        int iGrandMother =  process[iMother].mother1(); 
        if (iGrandMother == 1 || iGrandMother == 2) {
          doCopy  = true;
          isFinal = true;
        }  
      }
    }

    // Optionally non-final are not copied.
    if (finalOnly && !isFinal) doCopy = false;
   
    // Do copying and modify status/daughters of final.
    if (doCopy) {
      int iNew = workEvent.append( process[i]);
      if (isFinal) {
        workEvent[iNew].statusPos(); 
        workEvent[iNew].daughters( 0, 0);
        // When final only : no mothers; position in full event as daughters. 
        if (finalOnly) {  
          workEvent[iNew].mothers( 0, 0);
          workEvent[iNew].daughters( i, i);
        }
      }
    }
  }

}

//--------------------------------------------------------------------------

// subEvent extracts currently resolved partons in the hard process.

void UserHooks::subEvent(const Event& event, bool isHardest) {

  // Reset work event to be empty. 
  workEvent.clear();  

  // At the PartonLevel final partons are bookkept by subsystem.
  if (partonSystemsPtr->sizeSys() > 0) {

    // Find which subsystem to study.
    int iSys = 0;
    if (!isHardest) iSys = partonSystemsPtr->sizeSys() - 1;

    // Loop through all the final partons of the given subsystem.
    for (int i = 0; i < partonSystemsPtr->sizeOut(iSys); ++i) {
      int iOld = partonSystemsPtr->getOut( iSys, i);

      // Copy partons to work event.
      int iNew = workEvent.append( event[iOld]); 

      // No mothers. Position in full event as daughters.  
      workEvent[iNew].mothers( 0, 0);
      workEvent[iNew].daughters( iOld, iOld);
    }

  // At the ProcessLevel no subsystems have been defined.
  } else {

    // Loop through all partons, and copy all final ones.
    for (int iOld = 0; iOld < event.size(); ++iOld) 
    if (event[iOld].isFinal()) {
      int iNew = workEvent.append( event[iOld]); 

      // No mothers. Position in full event as daughters.  
      workEvent[iNew].mothers( 0, 0);
      workEvent[iNew].daughters( iOld, iOld);
    }
  }
 
}
 
//==========================================================================

// The SuppressSmallPT class, derived from UserHooks.

//--------------------------------------------------------------------------

// Modify event weight at the trial level, before selection.

double SuppressSmallPT::multiplySigmaBy( const SigmaProcess* sigmaProcessPtr, 
  const PhaseSpace* phaseSpacePtr, bool ) {

  // Need to initialize first time this method is called.
  if (!isInit) {
    
    // Calculate pT0 as for multiparton interactions.
    // Fudge factor allows offset relative to MPI framework.
    double eCM    = phaseSpacePtr->ecm();
    double pT0Ref = settingsPtr->parm("MultipartonInteractions:pT0Ref");
    double ecmRef = settingsPtr->parm("MultipartonInteractions:ecmRef");
    double ecmPow = settingsPtr->parm("MultipartonInteractions:ecmPow");
    double pT0    = pT0timesMPI * pT0Ref * pow(eCM / ecmRef, ecmPow);
    pT20          = pT0 * pT0;
  
    // Initialize alpha_strong object as for multiparton interactions,
    // alternatively as for hard processes.
    double alphaSvalue;
    int    alphaSorder;    
    if (useSameAlphaSasMPI) {
      alphaSvalue = settingsPtr->parm("MultipartonInteractions:alphaSvalue");
      alphaSorder = settingsPtr->mode("MultipartonInteractions:alphaSorder");
    } else {
      alphaSvalue = settingsPtr->parm("SigmaProcess:alphaSvalue");
      alphaSorder = settingsPtr->mode("SigmaProcess:alphaSorder");
    }
    alphaS.init( alphaSvalue, alphaSorder); 

    // Initialization finished.
    isInit = true;
  }
        
  // Only modify 2 -> 2 processes.
  int nFinal = sigmaProcessPtr->nFinal();
  if (nFinal != 2) return 1.;

  // pT scale of process. Weight pT^4 / (pT^2 + pT0^2)^2 
  double pTHat     = phaseSpacePtr->pTHat();
  double pT2       = pTHat * pTHat;
  double wt        = pow2( pT2 / (pT20 + pT2) );

  if (numberAlphaS > 0) {
    // Renormalization scale and assumed alpha_strong.
    double Q2RenOld  = sigmaProcessPtr->Q2Ren();
    double alphaSOld = sigmaProcessPtr->alphaSRen();

    // Reweight to new alpha_strong at new scale.
    double Q2RenNew  = pT20 + Q2RenOld;
    double alphaSNew = alphaS.alphaS(Q2RenNew);
    wt              *= pow( alphaSNew / alphaSOld, numberAlphaS);
  }

  // End weight calculation.
  return wt;

}

 
//==========================================================================

} // end namespace Pythia8
