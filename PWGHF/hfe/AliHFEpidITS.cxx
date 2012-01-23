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
// ITS PID class
// checks ITS PID based on ITS dE/dx truncated mean
//
// Authors: Matus Kalisky <matus.kalisky@cern.ch>
//          Markus Fasel <M.Fasel@gsi.de>
//
#include <TClass.h>
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>

//#include "AliAODTrack.h"
//#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliVParticle.h"

#include "AliHFEpidITS.h"

ClassImp(AliHFEpidITS)

//___________________________________________________________________
AliHFEpidITS::AliHFEpidITS(const Char_t *name):
    AliHFEpidBase(name)
{
  //
  // Default constructor
  //
}

//___________________________________________________________________
AliHFEpidITS::AliHFEpidITS(const AliHFEpidITS &ref):
    AliHFEpidBase("")
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidITS &AliHFEpidITS::operator=(const AliHFEpidITS &ref){
  //
  // Assignment operator
  //
  if(this != &ref) ref.Copy(*this);
  return *this;
}

//___________________________________________________________________
AliHFEpidITS::~AliHFEpidITS(){
  //
  // Destructor
  //
}

//___________________________________________________________________
void AliHFEpidITS::Copy(TObject &o) const {
  //
  // Copy function
  // Provides a deep copy
  //
  AliHFEpidBase::Copy(o);
}

//___________________________________________________________________
Bool_t AliHFEpidITS::InitializePID(Int_t /*run*/){
  //
  // ITS PID initialization
  //
  return kTRUE;
}


//___________________________________________________________________
Int_t AliHFEpidITS::IsSelected(const AliHFEpidObject* /*track*/, AliHFEpidQAmanager* /*pidqa*/) const {
  //
  // Does PID decision for ITS
  // 
  return 11;  // @TODO: Implement ITS PID decision
}

//___________________________________________________________________
Double_t AliHFEpidITS::GetITSSignalV1(AliVParticle *vtrack){
  //
  // Calculate the ITS signal according to the mean charge of the clusters
  //
  if(!TString(vtrack->IsA()->GetName()).CompareTo("AliAODTrack")){
    AliError("PID for AODs not implemented yet");
    return 0.;
  }
  AliESDtrack *track = dynamic_cast<AliESDtrack *>(vtrack);
  if(!track) return 0.;
  Double_t signal = 0.;
#ifdef TRUNK
  Double_t dedx[4];
  track->GetITSdEdxSamples(dedx);
  signal = TMath::Mean(4, dedx);
#else
  signal = track->GetITSsignal();
#endif
  Double_t p = track->GetTPCInnerParam() ? track->GetTPCInnerParam()->P() : track->P();
  AliDebug(1, Form("Momentum: %f, ITS Signal: %f", p, signal));
  return signal;
}

//___________________________________________________________________
Double_t AliHFEpidITS::GetITSSignalV2(AliVParticle *vtrack){
  //
  // Calculates the ITS signal. Truncated mean is used.
  //
  if(!TString(vtrack->IsA()->GetName()).CompareTo("AliAODTrack")){
    AliError("PID for AODs not implemented yet");
    return 0.;
  }
  AliESDtrack *track = dynamic_cast<AliESDtrack *>(vtrack);
  if(!track) return 0.;
  Double_t dedx[4], tmp[4];
  Int_t indices[4];
  track->GetITSdEdxSamples(tmp);
  TMath::Sort(4, tmp, indices);
  for(Int_t ien = 0; ien < 4; ien++) dedx[ien] = tmp[indices[ien]];
  Double_t signal = TMath::Mean(3, dedx); 
  Double_t p = track->GetTPCInnerParam() ? track->GetTPCInnerParam()->P() : track->P();
  AliDebug(1, Form("Momentum: %f, ITS Signal: %f", p, signal));
  return signal;
}

