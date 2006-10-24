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
//  A Reconstructed Particle in PHOS    
//  To become a general class of AliRoot ?       
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                 
//       
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---


// --- AliRoot header files ---
#include "AliPHOSRecParticle.h"
#include "AliPHOSGetter.h" 
#include "AliPHOSGeometry.h" 
#include "AliLog.h"

//____________________________________________________________________________
AliPHOSRecParticle::AliPHOSRecParticle(): 
  fPHOSTrackSegment(0),
  fDebug(kFALSE),
  fPos()
{
  // ctor
  const Int_t nSPECIES = AliPID::kSPECIESN;
  for(Int_t i = 0; i<nSPECIES ; i++)
    fPID[i]=0.;
}


//____________________________________________________________________________
AliPHOSRecParticle::AliPHOSRecParticle(const AliPHOSRecParticle & rp):
  AliPHOSFastRecParticle(rp),
  fPHOSTrackSegment(rp.fPHOSTrackSegment),
  fDebug(kFALSE),
  fPos()
{
  // copy ctor
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
  const Int_t nSPECIES = AliPID::kSPECIESN;
  for(Int_t i = 0; i<nSPECIES ; i++)
    fPID[i]=rp.fPID[i];
}

//____________________________________________________________________________
Int_t AliPHOSRecParticle::GetNPrimaries() const  
{ 
  return -1;
}

//____________________________________________________________________________
Int_t AliPHOSRecParticle::GetNPrimariesToRecParticles() const  
{ 
  // Get the number of primaries at the origine of the RecParticle
  Int_t rv = 0 ;
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  Int_t emcRPindex = dynamic_cast<AliPHOSTrackSegment*>(gime->TrackSegments()->At(GetPHOSTSIndex()))->GetEmcIndex();
  dynamic_cast<AliPHOSEmcRecPoint*>(gime->EmcRecPoints()->At(emcRPindex))->GetPrimaries(rv) ; 
  return rv ; 
}

//____________________________________________________________________________
const TParticle * AliPHOSRecParticle::GetPrimary() const  
{
  // Get the primary particle at the origine of the RecParticle and 
  // which has deposited the largest energy in SDigits
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  if (!gime) 
    Error("GetPrimary", "Getter not yet instantiated") ; 
  gime->Event(gime->EventNumber(), "SRTPX") ; 
  if(GetNPrimaries() == 0)
    return 0 ;
  if(GetNPrimaries() == 1)
    return GetPrimary(0) ;
  Int_t AbsId = 0;
  Int_t module ;
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ;
   Double_t x,z ;
  geom->ImpactOnEmc(static_cast<double>(Theta()),static_cast<double>(Phi()), module,z,x);
  Int_t amp = 0 ;
  Int_t iPrim=-1 ;
  if(module != 0){
    geom->RelPosToAbsId(module,x,z,AbsId) ;
   TClonesArray * sdigits = gime->SDigits() ;
   AliPHOSDigit * sdig ;
    
   for(Int_t i = 0 ; i < sdigits->GetEntriesFast() ; i++){
     sdig = static_cast<AliPHOSDigit *>(sdigits->At(i)) ;
     if((sdig->GetId() == AbsId)){
       if((sdig->GetAmp() > amp) && (sdig->GetNprimary())){
	 amp = sdig->GetAmp() ;
	 iPrim = sdig->GetPrimary(1) ;
       } 
       // do not scan rest of list
       if(sdig->GetId() > AbsId)
	 continue ; 
     }
   }
  }
  if(iPrim >= 0)
    return gime->Primary(iPrim) ;
  else
    return 0 ;
} 
  
//____________________________________________________________________________
Int_t AliPHOSRecParticle::GetPrimaryIndex() const  
{
  // Get the primary track index in TreeK which deposits the most energy
  // in Digits which forms EmcRecPoint, which produces TrackSegment,
  // which the RecParticle is created from


  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  if (!gime) 
    AliError(Form("Getter not yet instantiated")) ; 
  //PH  gime->Event(gime->EventNumber(), "DRTX") ; 
  gime->Event(gime->EventNumber(), "DRT") ; 
  
  // Get TrackSegment corresponding to this RecParticle
  AliPHOSTrackSegment *ts          = gime->TrackSegment(fPHOSTrackSegment);

  // Get EmcRecPoint corresponding to this TrackSegment
  Int_t emcRecPointIndex = ts->GetEmcIndex();

  AliPHOSEmcRecPoint  *emcRecPoint = gime->EmcRecPoint(emcRecPointIndex);

  // Get the list of digits forming this EmcRecParticle
  Int_t  nDigits   = emcRecPoint->GetDigitsMultiplicity();
  Int_t *digitList = emcRecPoint->GetDigitsList();

  // Find the digit with maximum amplitude
  AliPHOSDigit *digit = 0;
  Int_t maxAmp = 0;
  Int_t bestDigitIndex = -1;
  for (Int_t iDigit=0; iDigit<nDigits; iDigit++) {
    digit = gime->Digit(digitList[iDigit]);
    if (digit->GetAmp() > maxAmp) {
      maxAmp = digit->GetAmp();
      bestDigitIndex = iDigit;
    }
  }
  if (bestDigitIndex>-1)
    digit = gime->Digit(digitList[bestDigitIndex]);
  if (digit==0) return -12345;

  
  // Get the list of primary tracks producing this digit
  // and find which track has more track energy.
  //Int_t nPrimary = digit->GetNprimary();
  //TParticle *track = 0;
  //Double_t energyEM     = 0;
  //Double_t energyHadron = 0;
  //Int_t    trackEM      = -1;
  //Int_t    trackHadron  = -1;
  //for (Int_t iPrim=0; iPrim<nPrimary; iPrim++) {
  //  Int_t iPrimary = digit->GetPrimary(iPrim+1);
  //  if (iPrimary == -1 || TMath::Abs(iPrimary)>10000000)
  //    continue ;  //PH 10000000 is the shift of the track 
  //                //PH indexes in the underlying event
  //  track = gime->Primary(iPrimary);
  //  Int_t pdgCode   = track->GetPdgCode();
  //  Double_t energy = track->Energy();
  //  if (pdgCode==22 || TMath::Abs(pdgCode)==11) {
  //    if (energy > energyEM) {
  //	energyEM = energy;
  //	trackEM = iPrimary;
  //      }
  //   }
  //  else {
  //     if (energy > energyHadron) {
  //	energyHadron = energy;
  //	trackHadron = iPrimary;
	//    }
  //  }
  //}
  // Preferences are given to electromagnetic tracks
  //if (trackEM     != -1) return trackEM;     // track is gamma or e+-
  //if (trackHadron != -1) return trackHadron; // track is hadron
  //return -12345;                             // no track found :(


  // Get the list of hits producing this digit,
  // find which hit has deposited more energy 
  // and find the primary track.

  AliPHOSHit *hit = 0;
  TClonesArray *hits = gime->Hits();
  if (hits==0) return -12345;

  Double_t maxedep  =  0;
  Int_t    maxtrack = -1;
  Int_t    nHits    = hits ->GetEntries();
  Int_t    id       = digit->GetId();

  for (Int_t iHit=0; iHit<nHits; iHit++) {
    hit = gime->Hit(iHit);
    if(hit->GetId() == id){
      Double_t edep  = hit->GetEnergy();
      Int_t    track = hit->GetPrimary();
      if(edep > maxedep){
	maxedep  = edep;
	maxtrack = track;
      }
    }
  }

  if (maxtrack != -1) return maxtrack; // track is hadron
  return -12345;                       // no track found :(
}

//____________________________________________________________________________
const TParticle * AliPHOSRecParticle::GetPrimary(Int_t index) const  
{
  // Get one of the primary particles at the origine of the RecParticle
  if ( index > GetNPrimariesToRecParticles() ) { 
    if (fDebug) 
      Warning("GetPrimary", "AliPHOSRecParticle::GetPrimary -> %d is larger that the number of primaries %d", 
	      index, GetNPrimaries()) ;
    return 0 ; 
  } 
  else { 
    Int_t dummy ; 
    AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 

    Int_t emcRPindex = dynamic_cast<AliPHOSTrackSegment*>(gime->TrackSegments()->At(GetPHOSTSIndex()))->GetEmcIndex();
    Int_t primaryindex = dynamic_cast<AliPHOSEmcRecPoint*>(gime->EmcRecPoints()->At(emcRPindex))->GetPrimaries(dummy)[index] ; 
    return gime->Primary(primaryindex) ;
   } 
  //  return 0 ; 
}

//____________________________________________________________________________
void AliPHOSRecParticle::SetPID(Int_t type, Float_t weight)
{
  // Set the probability densities that this reconstructed particle
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
  // 9       Conversion electron
  
  fPID[type] = weight ; 
}
