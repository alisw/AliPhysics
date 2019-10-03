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

//--------------------------------------------------
// Filling of CalTrkTrack objects in the MC reader task
//
// Author: magali.estienne@subatech.in2p3.fr
//         alexandre.shabetai@cern.ch
//-------------------------------------------------

// --- ROOT system ---
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TPDGCode.h>

// --- AliRoot system ---
#include "AliJetFillCalTrkTrackMC.h"
#include "AliAODMCParticle.h"
#include "AliJetKineReaderHeader.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"

using std::cout;
using std::endl;
ClassImp(AliJetFillCalTrkTrackMC)

/////////////////////////////////////////////////////////////////////////

AliJetFillCalTrkTrackMC::AliJetFillCalTrkTrackMC(): 
  AliJetFillCalTrkEvent(),
  fHadCorr(0),
  fApplyMIPCorrection(kTRUE),
  fVEvt(0x0),
  fMCEvent(0x0)
{
  // constructor
}

//-----------------------------------------------------------------------
AliJetFillCalTrkTrackMC::AliJetFillCalTrkTrackMC(AliVEvent* evt): 
  AliJetFillCalTrkEvent(),
  fHadCorr(0),
  fApplyMIPCorrection(kTRUE),
  fVEvt(evt),
  fMCEvent(0x0)
{
  // constructor
}

//-----------------------------------------------------------------------
AliJetFillCalTrkTrackMC::AliJetFillCalTrkTrackMC(const AliJetFillCalTrkTrackMC &det): 
  AliJetFillCalTrkEvent(det),
  fHadCorr(det.fHadCorr),
  fApplyMIPCorrection(det.fApplyMIPCorrection),
  fVEvt(det.fVEvt),
  fMCEvent(det.fMCEvent)
{
  // Copy constructor
}

//-----------------------------------------------------------------------
AliJetFillCalTrkTrackMC& AliJetFillCalTrkTrackMC::operator=(const AliJetFillCalTrkTrackMC& other)
{
  // Assignment
  if (this != &other) { 
   fHadCorr = other.fHadCorr;
   fApplyMIPCorrection = other.fApplyMIPCorrection;
   fVEvt = other.fVEvt;
   fMCEvent = other.fMCEvent;
  }
 
  return (*this);

}

//-----------------------------------------------------------------------
AliJetFillCalTrkTrackMC::~AliJetFillCalTrkTrackMC()
{
  // destructor
}

//-----------------------------------------------------------------------
void AliJetFillCalTrkTrackMC::Exec(Option_t const * /*option*/)
{
  // Main method.
  // Fill AliJetFillCalTrkTrackMC the with the charged particle information

  fDebug = fReaderHeader->GetDebug();
  fOpt = fReaderHeader->GetDetector();
  TString datatype = fReaderHeader->GetDataType();
  datatype.ToUpper();

  Int_t type=0;

  if      (datatype.Contains("AODMC2B")) { type=kTrackAODMCChargedAcceptance; }
  else if (datatype.Contains("AODMC2"))  { type=kTrackAODMCCharged; }
  else if (datatype.Contains("AODMC"))   { type=kTrackAODMCAll; }
 
  else if (!datatype.Contains("AOD") && datatype.Contains("MC2")) {type=kTrackKineCharged;} 
  else if (!datatype.Contains("AOD") && datatype.Contains("MC"))  {type=kTrackKineAll;}
  
  else { cout<< "Unknown Data type !" << endl; }

  // temporary storage of signal and pt cut flag
  Bool_t sflag = 0;
  Bool_t cflag = 0;

  // get cuts set by user
  Float_t ptMin  = fReaderHeader->GetPtCut();
  Float_t etaMin = fReaderHeader->GetFiducialEtaMin();
  Float_t etaMax = fReaderHeader->GetFiducialEtaMax();
  Float_t phiMin = fReaderHeader->GetFiducialPhiMin();
  Float_t phiMax = fReaderHeader->GetFiducialPhiMax();

  if(fDebug>2) Printf("%s:%d Selecting tracks with %d",(char*)__FILE__,__LINE__,type);

  Int_t goodTrackStd = 0;
  Int_t goodTrackNonStd = 0;

  if (type == kTrackKineAll || type == kTrackKineCharged){

    FillKine();
  } // Kine

  else if (type == kTrackAODMCCharged || type == kTrackAODMCAll || type == kTrackAODMCChargedAcceptance) {
    if(!fVEvt) return; 
    TClonesArray *mcArray = dynamic_cast<TClonesArray*>(fVEvt->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!mcArray) {
      Printf("%s:%d No MC particle branch found",(char*)__FILE__,__LINE__);
      return; 
    }
    for(int it = 0;it < mcArray->GetEntriesFast();++it){
      cflag=sflag=0;
      AliAODMCParticle *part = (AliAODMCParticle*)(mcArray->At(it));
      if(!part->IsPhysicalPrimary())continue;
      Int_t   pdg     = TMath::Abs(part->GetPdgCode());
      // exclude neutrinos anyway
      if((pdg == 12 || pdg == 14 || pdg == 16)) continue;
      cflag = ( part->Pt()>ptMin ) ? 1 : 0;
      sflag = ( TMath::Abs(part->GetLabel()) < 10000 ) ? 1 : 0;
      if(type == kTrackAODMCAll){
        fCalTrkEvent->AddCalTrkTrack(part,cflag,sflag);
        goodTrackStd++;
      }
      else if (type == kTrackAODMCCharged || type == kTrackAODMCChargedAcceptance ){
        if(part->Charge()==0)continue;
        if(kTrackAODMCCharged){
          fCalTrkEvent->AddCalTrkTrack(part,cflag,sflag);
          goodTrackStd++;
        }
        else {
          if((part->Eta()>etaMax) || (part->Eta()<etaMin)) continue;
          Float_t phi = ( (part->Phi()) < 0) ? (part->Phi()) + 2. * TMath::Pi() : (part->Phi());
          if((phi>phiMax) || (phi<phiMin)) continue;
          fCalTrkEvent->AddCalTrkTrack(part,cflag,sflag);
          goodTrackStd++;
        }
      }
    }
   if(fDebug>0) printf("Number of good tracks selected: %5d \n", goodTrackStd+goodTrackNonStd); 
  } // AODMCparticle

}

//-----------------------------------------------------------------------
Bool_t AliJetFillCalTrkTrackMC::FillKine()
{
  // Fill event
  Int_t goodTrack = 0;
  
  // Get the stack
  if(!fMCEvent) {cout<<"could not get MCEvent!"<<endl; return kFALSE;}
  
  Int_t nt = fMCEvent->GetNumberOfTracks();
  
  // Get cuts set by user and header
  Double_t ptMin = ((AliJetKineReaderHeader*) fReaderHeader)->GetPtCut();
  Float_t etaMin = fReaderHeader->GetFiducialEtaMin();
  Float_t etaMax = fReaderHeader->GetFiducialEtaMax();
  
  TLorentzVector p4;
  // Loop over particles
  for (Int_t it = 0; it < nt; it++) {
    if(!fMCEvent->IsPhysicalPrimary(it)) continue;
    AliVParticle* part = fMCEvent->GetTrack(it);
    Int_t   pdg    = TMath::Abs(part->PdgCode());
    Float_t pt     = part->Pt();
    
    if( (pdg == 12 || pdg == 14 || pdg == 16)) continue;
    
    Float_t p      = part->P();
    Float_t p0     = p;
    Float_t eta    = part->Eta();
    Float_t phi    = part->Phi();
    Float_t charge = part->Charge();
    
    if (((AliJetKineReaderHeader*)fReaderHeader)->ChargedOnly()) {
      if (charge == 0) continue;
    } // End charged only
    
    // Fast simulation of EMCAL if requested
    if (((AliJetKineReaderHeader*)fReaderHeader)->FastSimEMCAL()) {
      // Charged particles only
      if (charge != 0){
	// Simulate efficiency
	if (!Efficiency(p0, eta, phi)) continue;
	// Simulate resolution
	p = SmearMomentum(4, p0);
      } // end "if" charged particles
      // Neutral particles (exclude K0L, n, nbar)
      if (pdg == kNeutron || pdg == kK0Long) continue;
    } // End fast EMCAL
    
    // Fast simulation of TPC if requested
    if (((AliJetKineReaderHeader*)fReaderHeader)->FastSimTPC()) {
      // Charged particles only
      if (charge == 0)               continue;
      // Simulate efficiency
      if (!Efficiency(p0, eta, phi)) continue;
      // Simulate resolution
      p = SmearMomentum(4, p0);
    } // End fast TPC

    // Fill momentum array
    Float_t r  = p/p0;
    Float_t px = r * part->Px();
    Float_t py = r * part->Py();
    Float_t pz = r * part->Pz();
    Float_t m =  part->M();
    
    Float_t e  = TMath::Sqrt(px * px + py * py + pz * pz + m * m);
    p4 = TLorentzVector(px, py, pz, e);
    if ((p4.Eta()>etaMax) || (p4.Eta()<etaMin)) continue;

    // Flag used to store the r factor
    Float_t ptReso = r;
    Bool_t cflag = kFALSE, sflag = kTRUE;
    if (pt > ptMin) cflag = kTRUE; // track surviving pt cut

    fCalTrkEvent->AddCalTrkTrackKine(part,cflag,sflag,ptReso);

    goodTrack++;

  } // track loop

  if(fDebug>0) printf(" Number of good tracks %d \n", goodTrack);

  return kTRUE;

}

//-----------------------------------------------------------------------
Float_t AliJetFillCalTrkTrackMC::SmearMomentum(Int_t ind, Float_t p)
{
  //  The relative momentum error, i.e.
  //  (delta p)/p = sqrt (a**2 + (b*p)**2) * 10**-2,
  //  where typically a = 0.75 and b = 0.16 - 0.24 depending on multiplicity
  //  (the lower value is for dn/d(eta) about 2000, and the higher one for 8000)
  //
  //  If we include information from TRD, b will be a factor 2/3 smaller.
  //
  //  ind = 1: high multiplicity
  //  ind = 2: low  multiplicity
  //  ind = 3: high multiplicity + TRD
  //  ind = 4: low  multiplicity + TRD

  Float_t pSmeared;
  Float_t a = 0.75;
  Float_t b = 0.12;

  if (ind == 1) b = 0.12;
  if (ind == 2) b = 0.08;
  if (ind == 3) b = 0.12;
  if (ind == 4) b = 0.08;

  Float_t sigma =  p * TMath::Sqrt(a * a + b * b * p * p)*0.01;
  pSmeared = p + gRandom->Gaus(0., sigma);
  return pSmeared;

}

//-----------------------------------------------------------------------
Bool_t AliJetFillCalTrkTrackMC::Efficiency(Float_t p, Float_t /*eta*/, Float_t phi)
{
  // Fast simulation of geometrical acceptance and tracking efficiency
  
  //  Tracking
  Float_t eff = 0.99;
  if (p < 0.5) eff -= (0.5-p)*0.2/0.3;
  // Geometry
  if (p > 0.5) {
    phi *= 180. / TMath::Pi();
    // Sector number 0 - 17
    Int_t isec  = Int_t(phi / 20.);
    // Sector centre
    Float_t phi0 = isec * 20. + 10.;
    Float_t phir = TMath::Abs(phi-phi0);
    // 2 deg of dead space
    if (phir > 9.) eff = 0.;
  }
  
  if (gRandom->Rndm() > eff) {
    return kFALSE;
  } else {
    return kTRUE;
  }
  
}


