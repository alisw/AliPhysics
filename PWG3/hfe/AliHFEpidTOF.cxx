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

#include "AliAODpidUtil.h"
#include "AliAODTrack.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliPID.h"

#include "AliHFEdetPIDqa.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidQAmanager.h"


ClassImp(AliHFEpidTOF)

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF():
  AliHFEpidBase()
  , fPID(NULL)
  , fNsigmaTOF(3)
{
  //
  // Constructor
  //
} 

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const Char_t *name):
  AliHFEpidBase(name)
  , fPID(NULL)
  , fNsigmaTOF(3)
{
  //
  // Constructor
  //
}

//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const AliHFEpidTOF &c):
  AliHFEpidBase("")
  , fPID(NULL)
  , fNsigmaTOF(3)
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
  if(fPID) delete fPID;
}
//___________________________________________________________________
void AliHFEpidTOF::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidTOF &target = dynamic_cast<AliHFEpidTOF &>(ref);

  target.fPID = fPID;          
  target.fNsigmaTOF = fNsigmaTOF;

  AliHFEpidBase::Copy(ref);
}
//___________________________________________________________________
Bool_t AliHFEpidTOF::InitializePID(){
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
  if((!fESDpid && track->IsESDanalysis()) || (!fAODpid && track->IsAODanalysis())) return 0;
  AliDebug(2, "PID object available");

  AliHFEpidObject::AnalysisType_t anaType = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  const AliVTrack *vtrack = dynamic_cast<const AliVTrack *>(track->GetRecTrack());
  if(!(vtrack && (vtrack->GetStatus() & AliESDtrack::kTOFpid))) return 0;
  AliDebug(2, "Track Has TOF PID");

  if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTOFpid, AliHFEdetPIDqa::kBeforePID);

  // Fill before selection
  Double_t sigEle = NumberOfSigmas(track->GetRecTrack(), AliPID::kElectron, anaType);
  AliDebug(2, Form("Number of sigmas in TOF: %f", sigEle));
  Int_t pdg = 0;
  if(TMath::Abs(sigEle) < fNsigmaTOF){
    pdg = 11;
    if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTOFpid, AliHFEdetPIDqa::kAfterPID);
  }
 
  return pdg;
}

/////////////////////////////////////////////////////////
//
// Wrappers for functions which are not virtual but
// only available for ESD and AOD case separately
//
/////////////////////////////////////////////////////////

//___________________________________________________________________
Double_t AliHFEpidTOF::NumberOfSigmas(const AliVParticle *track, AliPID::EParticleType species, AliHFEpidObject::AnalysisType_t anaType) const {
  //    
  // Get the number of sigmas
  //
  Double_t nSigmas = 100;
  if(anaType == AliHFEpidObject::kESDanalysis){
    // ESD analysis
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack && fESDpid) nSigmas = fESDpid->NumberOfSigmasTOF(esdtrack, species, 0.);
  } else {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    if(aodtrack && fAODpid) nSigmas = fAODpid->NumberOfSigmasTOF(aodtrack, species);
  }
  return nSigmas;
}

//___________________________________________________________________
Double_t AliHFEpidTOF::GetTOFsignal(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anatype) const{
  //
  // Get the TOF signal
  //
  Double_t tofSignal = 0.;
  if(anatype == AliHFEpidObject::kESDanalysis){
    // ESD analysis
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack) tofSignal = esdtrack->GetTOFsignal();
  } else {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    if(aodtrack) tofSignal = aodtrack->GetDetPid() ? aodtrack->GetDetPid()->GetTOFsignal() : 0.;
  }
  return tofSignal;
}

//___________________________________________________________________
Double_t AliHFEpidTOF::GetTime0(AliHFEpidObject::AnalysisType_t anatype) const{
  //
  // Get Time0
  //
  Double_t time0 = 0.;
  if(anatype == AliHFEpidObject::kESDanalysis){
    time0 = fESDpid->GetTOFResponse().GetTimeZero();
  }else {
    time0 = fAODpid->GetTOFResponse().GetTimeZero();
  }
  return time0;
}

//___________________________________________________________________
void AliHFEpidTOF::GetIntegratedTimes(const AliVParticle *track, Double_t *times, AliHFEpidObject::AnalysisType_t anatype) const{
  //
  // Return integrated times
  //
  if(anatype == AliHFEpidObject::kESDanalysis){
    // ESD analysis
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack) esdtrack->GetIntegratedTimes(times);
  } else {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    if(aodtrack) aodtrack->GetDetPid()->GetIntegratedTimes(times);
  }

}

