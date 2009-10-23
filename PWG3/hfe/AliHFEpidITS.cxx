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

//___________________________________________________________________
AliHFEpidITS::AliHFEpidITS(const Char_t *name):
    AliHFEpidBase(name)
  , fQAlist(0x0)
{
  //
  // Default constructor
  //
}

//___________________________________________________________________
AliHFEpidITS::AliHFEpidITS(const AliHFEpidITS &ref):
    AliHFEpidBase("")
  , fQAlist(0x0)
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
  if(fQAlist){
    fQAlist->Clear();
    delete fQAlist;
  }
}

//___________________________________________________________________
void AliHFEpidITS::Copy(TObject &o) const {
  //
  // Copy function
  // Provides a deep copy
  //
  AliHFEpidITS &target = dynamic_cast<AliHFEpidITS &>(o);

  target.fQAlist = dynamic_cast<TList *>(fQAlist->Clone());
  AliHFEpidBase::Copy(target);
}

//___________________________________________________________________
Bool_t AliHFEpidITS::InitializePID(){
  //
  // ITS PID initialization
  //
  return kTRUE;
}


//___________________________________________________________________
Int_t AliHFEpidITS::IsSelected(AliHFEpidObject* /*track*/){
  //
  // Does PID decision for ITS
  // 
  return 11;  // @TODO: Implement ITS PID decision
}

//___________________________________________________________________
void AliHFEpidITS::AddQAhistograms(TList *l){
  //
  // Adding QA histograms for ITS PID
  //
  // QA histograms are:
  // - all Particles ITS signal vs. p (both methods)
  // - single particle ITS signal vs. p (both methods)
  // 
  if(!fQAlist) fQAlist = new TList;
  fQAlist->SetName("fITSqaHistograms");

  // prepare axis
  const Int_t kMomentumBins = 41;
  const Double_t kPtMin = 0.1;
  const Double_t kPtMax = 10.;
  const Int_t kSigBins = 300;
  Double_t momentumBins[kMomentumBins];
  for(Int_t ibin = 0; ibin < kMomentumBins; ibin++)
    momentumBins[ibin] = static_cast<Double_t>(TMath::Power(10,TMath::Log10(kPtMin) + (TMath::Log10(kPtMax)-TMath::Log10(kPtMin))/(kMomentumBins-1)*static_cast<Double_t>(ibin)));

  TH2 *histo = NULL;
  fQAlist->AddAt((histo = new TH2F("fITSsigV1all", "ITS signal vs. p (all species, Method1):p / GeV/c:ITS signal / a.u.", kMomentumBins - 1, momentumBins, kSigBins, 0., kSigBins)), kITSsigV1);
  fQAlist->AddAt((histo = new TH2F("fITSsigV2all", "ITS signal vs. p (all species, Method2):p / GeV/c:ITS signal / a.u.", kMomentumBins - 1, momentumBins, kSigBins, 0., kSigBins)), kITSsigV2);
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    fQAlist->AddAt((histo = new TH2F(Form("fITSsigV1%s", AliPID::ParticleName(ispec)), Form("ITS signal vs. p (%s, Method1):p / GeV/c:ITS signal / a.u.", AliPID::ParticleName(ispec)), kMomentumBins - 1, momentumBins, kSigBins, 0., kSigBins)), 2 * ispec + kHistosSigAll);
    fQAlist->AddAt((histo = new TH2F(Form("fITSsigV2%s", AliPID::ParticleName(ispec)), Form("ITS signal vs. p (%s, Method2):p / GeV/c:ITS signal / a.u.", AliPID::ParticleName(ispec)), kMomentumBins - 1, momentumBins, kSigBins, 0., kSigBins)), 2 * ispec + 1 + kHistosSigAll);
  }
  l->Add(fQAlist);
}

//___________________________________________________________________
Double_t AliHFEpidITS::GetITSSignalV1(AliVParticle *vtrack, Int_t mcPID){
  //
  // Calculate the ITS signal according to the mean charge of the clusters
  //
  if(!TString(vtrack->IsA()->GetName()).CompareTo("AliAODtrack")){
    AliError("PID for AODs not implemented yet");
    return 0.;
  }
  AliESDtrack *track = dynamic_cast<AliESDtrack *>(vtrack);
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
  if(IsQAon()) FillHistogramsSignalV1(p, signal, mcPID);
  return signal;
}

//___________________________________________________________________
Double_t AliHFEpidITS::GetITSSignalV2(AliVParticle *vtrack, Int_t mcPID){
  //
  // Calculates the ITS signal. Truncated mean is used.
  //
  if(!TString(vtrack->IsA()->GetName()).CompareTo("AliAODtrack")){
    AliError("PID for AODs not implemented yet");
    return 0.;
  }
  AliESDtrack *track = dynamic_cast<AliESDtrack *>(vtrack);
  Double_t dedx[4], tmp[4];
  Int_t indices[4];
  track->GetITSdEdxSamples(tmp);
  TMath::Sort(4, tmp, indices);
  for(Int_t ien = 0; ien < 4; ien++) dedx[ien] = tmp[indices[ien]];
  Double_t signal = TMath::Mean(3, dedx); 
  Double_t p = track->GetTPCInnerParam() ? track->GetTPCInnerParam()->P() : track->P();
  AliDebug(1, Form("Momentum: %f, ITS Signal: %f", p, signal));
  if(IsQAon()) FillHistogramsSignalV2(p, signal, mcPID);
  return signal;
}

//___________________________________________________________________
void AliHFEpidITS::FillHistogramsSignalV1(Double_t p, Double_t signal, Int_t species){
  (dynamic_cast<TH2 *>(fQAlist->At(kITSsigV1)))->Fill(p, signal);
  if(species >= 0 && species < AliPID::kSPECIES)
    (dynamic_cast<TH2 *>(fQAlist->At(kHistosSigAll + 2 * species)))->Fill(p, signal);
}

//___________________________________________________________________
void AliHFEpidITS::FillHistogramsSignalV2(Double_t p, Double_t signal, Int_t species){
  (dynamic_cast<TH2 *>(fQAlist->At(kITSsigV2)))->Fill(p, signal);
  if(species >= 0 && species < AliPID::kSPECIES)
    (dynamic_cast<TH2 *>(fQAlist->At(kHistosSigAll + 2 * species + 1)))->Fill(p, signal);
}

