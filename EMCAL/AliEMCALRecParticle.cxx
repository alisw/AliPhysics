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
  AliEMCALRecParticle::AliEMCALRecParticle(): fEMCALRecPoint(0), fDebug(kFALSE)
{
  // ctor
  const Int_t nSPECIES = AliESDtrack::kSPECIESN;
  for(Int_t i = 0; i<nSPECIES ; i++)
    fPID[i]=0.;
}
//____________________________________________________________________________
AliEMCALRecParticle::AliEMCALRecParticle(const AliEMCALRecParticle & rp)
  : AliEMCALFastRecParticle(rp)
{
  // copy ctor

  fEMCALRecPoint    = rp.fEMCALRecPoint ; 
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
  const Int_t nSPECIES = AliESDtrack::kSPECIESN;
  for(Int_t i = 0; i<nSPECIES ; i++)
    fPID[i]=rp.fPID[i];  
}

//____________________________________________________________________________
Int_t AliEMCALRecParticle::GetNPrimaries() const  
{   
  return -1;
}

//____________________________________________________________________________
Int_t AliEMCALRecParticle::GetNPrimariesToRecParticles() const  
{ 
  // Returns the number of primaries at the origine of a RecParticle
  Int_t rv = 0 ;
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  dynamic_cast<AliEMCALRecPoint*>(gime->ECARecPoints()->At(GetEMCALRPIndex()))->GetPrimaries(rv) ; 

  return rv ; 
}

//____________________________________________________________________________
const TParticle * AliEMCALRecParticle::GetPrimary(Int_t index) const  
{
  // Getts the list of primary particles at the origine of the RecParticle
  if ( index > GetNPrimariesToRecParticles() ) { 
    if (fDebug) 
      Warning("GetPrimary", "AliEMCALRecParticle::GetPrimary -> %d is larger that the number of primaries %d", 
	      index, GetNPrimaries()) ;
    return 0 ; 
  } 
  Int_t dummy ; 
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 

  Int_t primaryindex = dynamic_cast<AliEMCALRecPoint*>(gime->ECARecPoints()->At(GetEMCALRPIndex()))->GetPrimaries(dummy)[index] ; 

  return gime->Primary(primaryindex) ;
}

//____________________________________________________________________________
const Double_t * AliEMCALRecParticle::GetPID()
{
  // Get the probability densities that this reconstructed particle
  // has a type of i:
  // i       particle types
  // ----------------------
  // 0       electron
  // 1       muon
  // 2       pi+-
  // 3       K+-
  // 4       p/pbar
  // 5       photon
  // 6       pi0 at high pt
  // 7       neutron
  // 8       K0L

 
  if (IsElectron()     ) fPID[0] = 1.0;
  if (IsChargedHadron()) {
    fPID[1] = 0.25;
    fPID[2] = 0.25;
    fPID[3] = 0.25;
    fPID[4] = 0.25;
  }
  if (IsFastChargedHadron()) {
    fPID[1] = 0.33;
    fPID[2] = 0.33;
    fPID[3] = 0.33;
    fPID[4] = 0.00;
  }
  if (IsSlowChargedHadron()) {
    fPID[1] = 0.00;
    fPID[2] = 0.00;
    fPID[3] = 0.00;
    fPID[4] = 1.00;
  }

  if (IsPhoton() || IsHardPhoton()) fPID[5] = 1.0;
  if (IsHardPi0())                  fPID[6] = 1.0;
  if (IsFastNeutralHadron())        fPID[7] = 1.0;
  if (IsSlowNeutralHadron())        fPID[8] = 1.0;

  if (IsEleCon()) fPID[9] = 1.0;
  return fPID;
}
