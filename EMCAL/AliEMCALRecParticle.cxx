/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */
//_________________________________________________________________________
//  A Reconstructed Particle in EMCAL    
//  To become a general class of AliRoot ?       
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                 
//       
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---


// --- AliRoot header files ---
#include "AliEMCALRecParticle.h"
#include "AliEMCALGetter.h" 
#include "TParticle.h"

ClassImp(AliEMCALRecParticle)


//____________________________________________________________________________
AliEMCALRecParticle::AliEMCALRecParticle(const AliEMCALRecParticle & rp)
  : AliEMCALFastRecParticle(rp)
{
  // copy ctor

  fEMCALTrackSegment = rp.fEMCALTrackSegment ; 
  fDebug            = kFALSE ; 
  fType             = rp.fType ; 
  fIndexInList      = rp.fIndexInList ;

  fPdgCode     = rp.fPdgCode;
  fStatusCode  = rp.fStatusCode;
  fMother[0]   = rp.fMother[0];
  fMother[1]   = rp.fMother[1];
  fDaughter[0] = rp.fDaughter[0];
  fDaughter[1] = rp.fDaughter[1];
  fWeight      = rp.fWeight;
  fCalcMass    = rp.fCalcMass;
  fPx          = rp.fPx;
  fPy          = rp.fPy;
  fPz          = rp.fPz;
  fE           = rp.fE;
  fVx          = rp.fVx;
  fVy          = rp.fVy;
  fVz          = rp.fVz;
  fVt          = rp.fVt;
  fPolarTheta  = rp.fPolarTheta;
  fPolarPhi    = rp.fPolarPhi;
  fParticlePDG = rp.fParticlePDG; 
  
}

//____________________________________________________________________________
const Int_t AliEMCALRecParticle::GetNPrimaries() const  
{   
  return -1;
}

//____________________________________________________________________________
const Int_t AliEMCALRecParticle::GetNPrimariesToRecParticles() const  
{ 
  Int_t rv = 0 ;
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  Int_t ecaRPindex = dynamic_cast<AliEMCALTrackSegment*>(gime->TrackSegments()->At(GetEMCALTSIndex()))->GetECAIndex();
  dynamic_cast<AliEMCALTowerRecPoint*>(gime->ECARecPoints()->At(ecaRPindex))->GetPrimaries(rv) ; 
  return rv ; 
}

//____________________________________________________________________________
const TParticle * AliEMCALRecParticle::GetPrimary(Int_t index) const  
{
  if ( index > GetNPrimariesToRecParticles() ) { 
    if (fDebug) 
      Warning("GetPrimary", "AliEMCALRecParticle::GetPrimary -> %d is larger that the number of primaries %d", 
	      index, GetNPrimaries()) ;
    return 0 ; 
  } 
  Int_t dummy ; 
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 

  Int_t ecaRPindex = dynamic_cast<AliEMCALTrackSegment*>(gime->TrackSegments()->At(GetEMCALTSIndex()))->GetECAIndex();
  Int_t primaryindex = dynamic_cast<AliEMCALTowerRecPoint*>(gime->ECARecPoints()->At(ecaRPindex))->GetPrimaries(dummy)[index] ; 
  return gime->Primary(primaryindex) ;
}

