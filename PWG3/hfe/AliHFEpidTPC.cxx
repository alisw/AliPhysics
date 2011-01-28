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
//
// Class for TPC PID
// Implements the abstract base class AliHFEpidBase
// 
// Class contains TPC specific cuts and QA histograms
// Two selection strategies are offered: Selection of certain value
// regions in the TPC dE/dx (by IsSelected), and likelihoods
//
// Authors: 
//
//   Markus Fasel <M.Fasel@gsi.de> 
//   Markus Heide <mheide@uni-muenster.de> 
//  
#include <TF1.h>
#include <TMath.h>

#include "AliAODpidUtil.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliESDpid.h"

#include "AliHFEpidTPC.h"
#include "AliHFEpidQAmanager.h"

ClassImp(AliHFEpidTPC)

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC() :
  // add a list here
  AliHFEpidBase()
  , fLineCrossingsEnabled(0)
  , fUpperSigmaCut(NULL)
  , fLowerSigmaCut(NULL)
  , fElectronMeanCorrection(NULL)
  , fNsigmaTPC(3)
  , fRejectionEnabled(0)
  , fPID(NULL)
{
  //
  // default  constructor
  // 
  memset(fRejection, 0, sizeof(Float_t) * 4 * AliPID::kSPECIES);
  memset(fLineCrossingSigma, 0, sizeof(Double_t) * AliPID::kSPECIES);
  memset(fPAsigCut, 0, sizeof(Float_t) * 2);
  memset(fNAsigmaTPC, 0, sizeof(Float_t) * 2);

}

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC(const char* name) :
  // add a list here
  AliHFEpidBase(name)
  , fLineCrossingsEnabled(0)
  , fUpperSigmaCut(NULL)
  , fLowerSigmaCut(NULL)
  , fElectronMeanCorrection(NULL)
  , fNsigmaTPC(3)
  , fRejectionEnabled(0)
  , fPID(NULL)
{
  //
  // default  constructor
  // 
  //
  memset(fRejection, 0, sizeof(Float_t) * 4 * AliPID::kSPECIES);
  memset(fLineCrossingSigma, 0, sizeof(Double_t) * AliPID::kSPECIES);
  memset(fPAsigCut, 0, sizeof(Float_t) * 2);
  memset(fNAsigmaTPC, 0, sizeof(Float_t) * 2);
  fPID = new AliPID;
}

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC(const AliHFEpidTPC &ref) :
  AliHFEpidBase("")
  , fLineCrossingsEnabled(0)
  , fUpperSigmaCut(NULL)
  , fLowerSigmaCut(NULL)
  , fElectronMeanCorrection(NULL)
  , fNsigmaTPC(2)
  , fRejectionEnabled(0)
  , fPID(NULL)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidTPC &AliHFEpidTPC::operator=(const AliHFEpidTPC &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  } 
  return *this;
}
//___________________________________________________________________
void AliHFEpidTPC::Copy(TObject &o) const{
  //
  // Copy function 
  // called in copy constructor and assigment operator
  //
  AliHFEpidTPC &target = dynamic_cast<AliHFEpidTPC &>(o);

  target.fLineCrossingsEnabled = fLineCrossingsEnabled;
  target.fUpperSigmaCut = fUpperSigmaCut;
  target.fLowerSigmaCut = fLowerSigmaCut;
  target.fElectronMeanCorrection = fElectronMeanCorrection;
  target.fNsigmaTPC = fNsigmaTPC;
  target.fRejectionEnabled = fRejectionEnabled;
  target.fPID = new AliPID(*fPID);
  memcpy(target.fLineCrossingSigma, fLineCrossingSigma, sizeof(Double_t) * AliPID::kSPECIES);
  memcpy(target.fPAsigCut, fPAsigCut, sizeof(Float_t) * 2);
  memcpy(target.fNAsigmaTPC, fNAsigmaTPC, sizeof(Float_t) * 2);
 
  AliHFEpidBase::Copy(target);
}

//___________________________________________________________________
AliHFEpidTPC::~AliHFEpidTPC(){
  //
  // Destructor
  //
  if(fPID) delete fPID;
}

//___________________________________________________________________
Bool_t AliHFEpidTPC::InitializePID(){
  //
  // Add TPC dE/dx Line crossings
  //
  //AddTPCdEdxLineCrossing(AliPID::kKaon, 0.3, 0.018);
  //AddTPCdEdxLineCrossing(AliPID::kProton, 0.9, 0.054);
  return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidTPC::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const
{
  //
  // For the TPC pid we use the 2-sigma band around the bethe bloch curve
  // for electrons
  // exclusion of the crossing points
  //

  if((!fESDpid && track->IsESDanalysis()) || (!fAODpid && track->IsAODanalysis())) return 0;
  
  // QA before selection
  if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTPCpid, AliHFEdetPIDqa::kBeforePID);
  AliDebug(1, "Doing TPC PID based on n-Sigma cut approach");
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t nsigma = NumberOfSigmas(track->GetRecTrack(), AliPID::kElectron, anatype);
  AliDebug(1, Form("TPC NSigma: %f", nsigma));
  // exclude crossing points:
  // Determine the bethe values for each particle species
  Bool_t isLineCrossing = kFALSE;
  for(Int_t ispecies = 0; ispecies < AliPID::kSPECIES; ispecies++){
    if(ispecies == AliPID::kElectron) continue;
    if(!(fLineCrossingsEnabled & 1 << ispecies)) continue;
    if(TMath::Abs(NumberOfSigmas(track->GetRecTrack(), (AliPID::EParticleType)ispecies, anatype)) < fLineCrossingSigma[ispecies] && TMath::Abs(nsigma) < fNsigmaTPC){
      // Point in a line crossing region, no PID possible, but !PID still possible ;-)
      isLineCrossing = kTRUE;      
      break;
    }
  }
  if(isLineCrossing) return 0;

  // Check particle rejection
  if(HasParticleRejection()){
    Int_t reject = Reject(track->GetRecTrack(), anatype);
    if(reject != 0) return reject;
  }

  // Check if we have an asymmetric sigma model set
  Int_t pdg = 0;
  if(fUpperSigmaCut || fLowerSigmaCut){
    pdg = CutSigmaModel(track->GetRecTrack(), anatype) ? 11 : 0;
  } else { 
    // Perform Asymmetric n-sigma cut if required, else perform symmetric TPC sigma cut
    Float_t p = 0.;
    if(HasAsymmetricSigmaCut() && (p = track->GetRecTrack()->P()) >= fPAsigCut[0] && p <= fPAsigCut[1]){ 
      if(nsigma >= fNAsigmaTPC[0] && nsigma <= fNAsigmaTPC[1]) pdg = 11; 
    } else {
      if(TMath::Abs(nsigma) < fNsigmaTPC ) pdg = 11;
    }
  }
  if(pidqa && pdg != 0) pidqa->ProcessTrack(track, AliHFEpid::kTPCpid, AliHFEdetPIDqa::kAfterPID);
  return pdg;

}

//___________________________________________________________________
Bool_t AliHFEpidTPC::CutSigmaModel(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const {
  //
  // N SigmaCut using parametrization of the cuts
  //
  Bool_t isSelected = kTRUE;
  Float_t nsigma = NumberOfSigmas(track, AliPID::kElectron, anaType);
  Double_t p = GetP(track, anaType);
  if(fUpperSigmaCut && nsigma > fUpperSigmaCut->Eval(p)) isSelected = kFALSE;
  if(fLowerSigmaCut && nsigma < fLowerSigmaCut->Eval(p)) isSelected = kFALSE;
  return isSelected;
}

//___________________________________________________________________
Int_t AliHFEpidTPC::Reject(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const{
  //
  // reject particles based on asymmetric sigma cut
  //
  Int_t pdc[AliPID::kSPECIES] = {11,13,211,321,2212};
  Double_t p = GetP(track, anaType);
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    if(!TESTBIT(fRejectionEnabled, ispec)) continue;
    // Particle rejection enabled
    if(p < fRejection[4*ispec] || p > fRejection[4*ispec+2]) continue;
    Double_t sigma = NumberOfSigmas(track, static_cast<AliPID::EParticleType>(ispec), anaType);
    if(sigma >= fRejection[4*ispec+1] && sigma <= fRejection[4*ispec+3]) return pdc[ispec] * track->Charge();
  }
  return 0;
}

//___________________________________________________________________
void AliHFEpidTPC::AddTPCdEdxLineCrossing(Int_t species, Double_t sigma){
  //
  // Add exclusion point for the TPC PID where a dEdx line crosses the electron line
  // Stores line center and line sigma
  //
  if(species >= AliPID::kSPECIES){
    AliError("Species doesn't exist");
    return;
  }
  fLineCrossingsEnabled |= 1 << species;
  fLineCrossingSigma[species] = sigma;
}

//___________________________________________________________________
Double_t AliHFEpidTPC::NumberOfSigmas(const AliVParticle *track, AliPID::EParticleType species, AliHFEpidObject::AnalysisType_t anaType) const {
  //    
  // Get the number of sigmas
  //
  Double_t nSigmas = 100;
  if(anaType == AliHFEpidObject::kESDanalysis){
    // ESD analysis
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack && fESDpid) nSigmas = fESDpid->NumberOfSigmasTPC(esdtrack, species);
  } else {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    if(aodtrack && fAODpid) nSigmas = fAODpid->NumberOfSigmasTPC(aodtrack, species);
  }
  // Correct for the mean o
  if(fElectronMeanCorrection)
    nSigmas -= fElectronMeanCorrection->Eval(GetP(track, anaType));   
  return nSigmas;
}

//___________________________________________________________________
Double_t AliHFEpidTPC::GetP(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anatype) const {
  //
  // Get the momentum at the inner wall of the TPC
  //
  Double_t p = -1;
  if(anatype == AliHFEpidObject::kESDanalysis){
    // ESD analysis: Use Inner Params for the momentum estimate
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack) p = esdtrack->GetInnerParam() ? esdtrack->GetInnerParam()->GetP() : esdtrack->P();
  } else { 
    // AOD analysis: Use TPC momentum stored in the AliAODpid object
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    if(aodtrack) p = aodtrack->GetDetPid() ? aodtrack->GetDetPid()->GetTPCmomentum() : aodtrack->P();
  }
  return p;
}
