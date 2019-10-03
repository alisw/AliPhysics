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
// Class for TOF PID
// Implements the abstract base class AliHFEpidBase
// IsInitialized() does the PID decision
// 
// Authors:
//   Markus Fasel  <M.Fasel@gsi.de>
//   Matus Kalisky <matus.kalisky@cern.ch>  (contact)
//

#include <TMath.h>

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"

#include "AliHFEdetPIDqa.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidQAmanager.h"


ClassImp(AliHFEpidTOF)

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF():
  AliHFEpidBase()
  , fNsigmaTOF(3)
  , fUseOnlyIfAvailable(kFALSE)
  , fRejectMismatch(kFALSE)
  , fGenerateTOFmismatch(kFALSE)
  , fNmismatchTracks(0)
{
  //
  // Constructor
  //
  
  memset(fSigmaBordersTOFLower, 0, sizeof(Float_t) * 12);
  memset(fSigmaBordersTOFUpper, 0, sizeof(Float_t) * 12);
} 

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const Char_t *name):
  AliHFEpidBase(name)
  , fNsigmaTOF(3)
  , fUseOnlyIfAvailable(kFALSE)
  , fRejectMismatch(kFALSE)
  , fGenerateTOFmismatch(kFALSE)
  , fNmismatchTracks(0)
{
  //
  // Constructor
  //
  memset(fSigmaBordersTOFLower, 0, sizeof(Float_t) * 12);
  memset(fSigmaBordersTOFUpper, 0, sizeof(Float_t) * 12);

}

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const AliHFEpidTOF &c):
  AliHFEpidBase("")
  , fNsigmaTOF(3)
  , fUseOnlyIfAvailable(kFALSE)
  , fRejectMismatch(kFALSE)
  , fGenerateTOFmismatch(c.fGenerateTOFmismatch)
  , fNmismatchTracks(c.fNmismatchTracks)
{  
  // 
  // Copy operator
  //

  c.Copy(*this);
}
//___________________________________________________________________
AliHFEpidTOF &AliHFEpidTOF::operator=(const AliHFEpidTOF &ref){
  //
  // Assignment operator
  //

  if(this != &ref){
    ref.Copy(*this);
  }

  return *this;
}
//___________________________________________________________________
AliHFEpidTOF::~AliHFEpidTOF(){
  //
  // Destructor
  //
}
//___________________________________________________________________
void AliHFEpidTOF::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidTOF &target = dynamic_cast<AliHFEpidTOF &>(ref);

  target.fNsigmaTOF = fNsigmaTOF;
  target.fUseOnlyIfAvailable = fUseOnlyIfAvailable;
  target.fRejectMismatch = fRejectMismatch;
  target.fGenerateTOFmismatch = fGenerateTOFmismatch;
  target.fNmismatchTracks = fNmismatchTracks;
  memcpy(target.fSigmaBordersTOFLower, fSigmaBordersTOFLower, sizeof(Float_t) * 12);
  memcpy(target.fSigmaBordersTOFUpper, fSigmaBordersTOFUpper, sizeof(Float_t) * 12);

  AliHFEpidBase::Copy(ref);
}
//___________________________________________________________________
Bool_t AliHFEpidTOF::InitializePID(Int_t /*run*/){
  //
  // InitializePID: TOF experts have to implement code here
  //
  return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidTOF::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const
{
  //
  // TOF PID based on n-Sigma cut
  // Selects Protons and Kaons via n-sigma cut up to 3 GeV/c
  // In addition histos for n-sigma before (all species) and after (only closest species) are filled
  //
  if(!fkPIDResponse) return 0;
  AliDebug(2, "PID object available");

  const AliVTrack *vtrack = dynamic_cast<const AliVTrack *>(track->GetRecTrack());
  if(!vtrack) return 0;
  //Bool_t hasTOFpid = vtrack->GetStatus() & AliESDtrack::kTOFpid;
  AliPIDResponse::EDetPidStatus statuspidtof = fkPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,vtrack);
  Bool_t hasTOFpid = kFALSE;
  if(statuspidtof==AliPIDResponse::kDetPidOk) hasTOFpid = kTRUE;
  if(fUseOnlyIfAvailable && !hasTOFpid){
    AliDebug(2, "No TOF PID, but PID required only if available");
    return 11;   
  } else if(!hasTOFpid){
    AliDebug(2, "No TOF PID, and TOF PID is required always");
    return 0;
  }
  AliDebug(2, "Track Has TOF PID");
  
  if(fRejectMismatch){
    if(IsMismatch(vtrack)) return 0;
  }

  if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTOFpid, AliHFEdetPIDqa::kBeforePID);

  // Fill before selection
  Double_t sigEle = fkPIDResponse->NumberOfSigmasTOF(track->GetRecTrack(), AliPID::kElectron);
  AliDebug(2, Form("Number of sigmas in TOF: %f", sigEle));
  Int_t pdg = 0;
  if(TestBit(kSigmaBand)){
    Int_t centrality = track->IsPbPb() ? track->GetCentrality() + 1 : 0;
    AliDebug(2, Form("Centrality: %d\n", centrality));
    if(centrality > 11) return kFALSE;
    if(sigEle > fSigmaBordersTOFLower[centrality] && sigEle < fSigmaBordersTOFUpper[centrality]) pdg = 11;
  } else {
    // Fixed nsigma cut
    if(TMath::Abs(sigEle) < fNsigmaTOF) pdg = 11;
  }

  if(pdg == 11 && pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTOFpid, AliHFEdetPIDqa::kAfterPID);
 
  return pdg;
}
//___________________________________________________________________
void AliHFEpidTOF::SetTOFnSigmaBand(Float_t lower, Float_t upper)
{
  //
  // Lower and higher cut independant of the centrality
  //

  for(Int_t k=0; k < 12; k++) {
    fSigmaBordersTOFLower[k] = lower;
    fSigmaBordersTOFUpper[k] = upper;
  }

  SetBit(kSigmaBand, kTRUE);

}
//___________________________________________________________________
void AliHFEpidTOF::SetTOFnSigmaBandCentrality(Float_t lower, Float_t upper, Int_t centralityBin)
{
  //
  // Lower and higher cut as function of centrality
  //

  if(centralityBin < 11) {
    fSigmaBordersTOFLower[centralityBin+1] = lower;
    fSigmaBordersTOFUpper[centralityBin+1] = upper;
  }

  SetBit(kSigmaBand, kTRUE);

}

//___________________________________________________________________
Bool_t AliHFEpidTOF::IsMismatch(const AliVTrack * const track) const {
  //
  // Check for mismatch
  //
  if(!fkPIDResponse) return kFALSE;
  Double_t probs[AliPID::kSPECIESC];
  AliPIDResponse::EDetPidStatus status = fkPIDResponse->ComputeTOFProbability(track, AliPID::kSPECIESC, probs);
  return status == AliPIDResponse::kDetMismatch;
}

//___________________________________________________________________
void AliHFEpidTOF::GenerateTOFmismatch(const AliVTrack * const trk, int ntrk, TArrayD &sigmaEl){
  //
  // Function generate randomised TOF mismatch hits for a given input track. The number of generated
  // mismatch tracks is steered by the parameter ntrk. For all mismatch tracks the number of sigmas
  // to the electron time-of-flight hypothesis is calculated, and the resulting numbers of sigmas
  // are stored in the array sigmaEl for further processing
  //
  if(sigmaEl.GetSize() < ntrk) sigmaEl.Set(ntrk);
  sigmaEl.Reset();

  // work on copy
  AliVTrack *copytrk(NULL);
  Bool_t isAOD = kFALSE;
  if(dynamic_cast<const AliESDtrack *>(trk)){
    copytrk = new AliESDtrack(*(static_cast<const AliESDtrack *>(trk)));
  } else {
    copytrk = new AliAODTrack(*(static_cast<const AliAODTrack *>(trk)));
    isAOD = kTRUE;
  }

  // Generate mismatch values for number of sigmas to the electron hypothesis and store then in the 
  // output array
  for(int itrk = 0; itrk < ntrk; itrk++){
    Double_t tofsignal = AliTOFPIDResponse::GetMismatchRandomValue(copytrk->Eta());
    if(isAOD){
      AliAODTrack *aodtrk = static_cast<AliAODTrack *>(copytrk);
      AliAODPid *aodpid = aodtrk->GetDetPid();
      if(aodpid) aodpid->SetTOFsignal(tofsignal);
    } else {
      AliESDtrack *esdtrk = static_cast<AliESDtrack *>(copytrk);
      esdtrk->SetTOFsignal(tofsignal);
    }
    sigmaEl[itrk] = fkPIDResponse->NumberOfSigmasTOF(copytrk, AliPID::kElectron);
  }
  delete copytrk;
}
