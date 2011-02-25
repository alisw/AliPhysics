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
// PID Response class for the TRD detector
// Based on 1D Likelihood approach
// Calculation of probabilities using Bayesian approach
// Attention: This method is only used to separate electrons from pions
//
// Authors:
//  Markus Fasel <M.Fasel@gsi.de>
//  Anton Andronic <A.Andronic@gsi.de>
//
#include <TClass.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObjArray.h>
#include <TString.h>
#include <TROOT.h> 
#include <TSystem.h>
#include <TDirectory.h>

#include "AliLog.h"

#include "AliTRDPIDResponse.h"

ClassImp(AliTRDPIDResponse)

const Double_t AliTRDPIDResponse::fgkPBins[kNPBins] = {1, 2, 3, 4, 5, 6};

//____________________________________________________________
AliTRDPIDResponse::AliTRDPIDResponse():
  TObject()
  ,fReferences(NULL)
  ,fGainNormalisationFactor(1.)
  ,fPIDmethod(kLQ1D)
{
  //
  // Default constructor
  //
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++)
    for(Int_t ipbin = 0; ipbin < kNPBins; ipbin++)
      fMapRefHists[ispec][ipbin] = -1;
}

//____________________________________________________________
AliTRDPIDResponse::AliTRDPIDResponse(const AliTRDPIDResponse &ref):
  TObject(ref)
  ,fReferences(ref.fReferences)
  ,fGainNormalisationFactor(ref.fGainNormalisationFactor)
  ,fPIDmethod(ref.fPIDmethod)
{
  //
  // Copy constructor
  // Flat copy of the reference histos
  // For creating a deep copy call SetOwner
  //
  Int_t size  = (AliPID::kSPECIES)*(kNPBins);
  memcpy(fMapRefHists, ref.fMapRefHists, sizeof(Double_t) * size);
}

//____________________________________________________________
AliTRDPIDResponse &AliTRDPIDResponse::operator=(const AliTRDPIDResponse &ref){
  //
  // Assignment operator
  // Performs a flat copy of the reference histos
  //
  if(this == &ref) return *this;
  
  // Make copy
  TObject::operator=(ref);
  if(IsOwner()){
    if(fReferences){
      fReferences->Delete();
      delete fReferences;
    }
    SetBit(kIsOwner, kFALSE);
  } else if(fReferences) delete fReferences;
  fReferences = ref.fReferences;
  Int_t size  = (AliPID::kSPECIES)*(kNPBins);
  memcpy(fMapRefHists, ref.fMapRefHists, sizeof(Double_t) * size);
  fPIDmethod = ref.fPIDmethod;
  
  return *this;
}

//____________________________________________________________
AliTRDPIDResponse::~AliTRDPIDResponse(){
  //
  // Destructor
  //
  if(IsOwner()){
    // Destroy histos
    if(fReferences) fReferences->Delete();
  }
  if(fReferences) delete fReferences;
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::Load(const Char_t * filename){
  //
  // Load References into the toolkit
  //
  AliDebug(1, "Loading reference histos from root file");
  TDirectory *owd = gDirectory;// store old working directory
  
  if(!filename)
    filename = Form("%s/STEER/LQ1dRef_v1.root",gSystem->ExpandPathName("$ALICE_ROOT"));
  TFile *in = TFile::Open(filename);
  if(!in){
    AliError("Ref file not available.");
    return kFALSE;
  }
  
  gROOT->cd();
  fReferences = new TObjArray(AliPID::kSPECIES*kNPBins);
  fReferences->SetOwner();
  TIter iter(in->GetListOfKeys());
  TKey *histkey = NULL;
  TObject *tmp = NULL;
  Int_t arrayPos = 0, pbin = 0; 
  AliPID::EParticleType species;
  TString histname;
  TH1 *hnew = NULL;
  while((histkey = dynamic_cast<TKey *>(iter()))){
    tmp = histkey->ReadObj();
    histname = tmp->GetName();
    if(histname.BeginsWith("fHQel")){ // Electron histogram
      histname.ReplaceAll("fHQel_p","");
      species = AliPID::kElectron;
    } else if(histname.BeginsWith("fHQmu")){ // Muon histogram
      histname.ReplaceAll("fHQmu_p","");
      species = AliPID::kMuon;
    } else if(histname.BeginsWith("fHQpi")){ // Pion histogram
      histname.ReplaceAll("fHQpi_p","");
      species = AliPID::kPion;
    }else if(histname.BeginsWith("fHQka")){ // Kaon histogram
      histname.ReplaceAll("fHQka_p","");
      species = AliPID::kKaon;
    }else if(histname.BeginsWith("fHQpr")){ // Proton histogram
      histname.ReplaceAll("fHQpr_p","");
      species = AliPID::kProton;
    } else continue;
    pbin = histname.Atoi() - 1;
    AliDebug(1, Form("Species %d, Histname %s, Pbin %d, Position in container %d", species, histname.Data(), pbin, arrayPos));
    fMapRefHists[species][pbin] = arrayPos;
    fReferences->AddAt((hnew =new TH1F(*dynamic_cast<TH1F *>(tmp))), arrayPos);
    hnew->SetDirectory(gROOT);
    arrayPos++;
  } 
  in->Close();
  owd->cd();
  AliDebug(2, Form("Successfully loaded %d Reference Histograms", arrayPos));
  SetBit(kIsOwner, kTRUE);
  delete in;
  return kTRUE;
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::Load(const TObjArray *histos){
  //
  // Load Reference into the PID Response (from OCDB/OADB)
  //
  AliDebug(1, "Setting references (from the database)");
  if(fReferences) fReferences->Clear();
  else
    fReferences = new TObjArray(AliPID::kSPECIES*kNPBins);
  fReferences->SetOwner();
  TIter iter(histos);
  Int_t arrayPos = 0, pbin = 0; 
  AliPID::EParticleType species;
  TObject *tmp = NULL;
  TString histname;
  TH1 *hnew = NULL;
  while((tmp = iter())){
    histname = tmp->GetName();
    if(histname.BeginsWith("fHQel")){ // Electron histogram
      histname.ReplaceAll("fHQel_p","");
      species = AliPID::kElectron;
    } else if(histname.BeginsWith("fHQmu")){ // Muon histogram
      histname.ReplaceAll("fHQmu_p","");
      species = AliPID::kMuon;
    } else if(histname.BeginsWith("fHQpi")){ // Pion histogram
      histname.ReplaceAll("fHQpi_p","");
      species = AliPID::kPion;
    }else if(histname.BeginsWith("fHQka")){ // Kaon histogram
      histname.ReplaceAll("fHQka_p","");
      species = AliPID::kKaon;
    }else if(histname.BeginsWith("fHQpr")){ // Proton histogram
      histname.ReplaceAll("fHQpr_p","");
      species = AliPID::kProton;
    } else continue;
    pbin = histname.Atoi() - 1;
    AliDebug(1, Form("Species %d, Histname %s, Pbin %d, Position in container %d", species, histname.Data(), pbin, arrayPos));
    fMapRefHists[species][pbin] = arrayPos;
    fReferences->AddAt((hnew = new TH1F(*dynamic_cast<TH1F *>(tmp))), arrayPos);
    hnew->SetDirectory(gROOT);
    arrayPos++;
  }
  if(!arrayPos){
    AliError("Failed loading References");
    return kFALSE;
  }
  AliDebug(2, Form("Successfully loaded %d Reference Histograms", arrayPos));
  SetBit(kIsOwner, kTRUE);
  return kTRUE;
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::GetResponse(Int_t n, Double_t *dedx, Float_t *p, Double_t prob[AliPID::kSPECIES], Bool_t kNorm) const
{
  //
  // Calculate TRD likelihood values for the track based on dedx and 
  // momentum values. The likelihoods are calculated by query the 
  // reference data depending on the PID method selected
  //
  // Input parameter :
  //   n - number of dedx slices/chamber
  //   dedx - array of dedx slices organized layer wise
  //   p - array of momentum measurements organized layer wise
  // 
  // Return parameters
  //   prob - probabilities allocated by TRD for particle specis
  //   kNorm - switch to normalize probabilities to 1. By default true. If false return not normalized prob.
  // 
  // Return value
  //   true if calculation success
  // 

  if(!fReferences){
    AliWarning("Missing reference data. PID calculation not possible.");
    return kFALSE;
  }

  for(Int_t is(AliPID::kSPECIES); is--;) prob[is]=.2;
  Double_t prLayer[AliPID::kSPECIES];
  Double_t DE[10], s(0.);
  for(Int_t il(kNlayer); il--;){
    memset(prLayer, 0, AliPID::kSPECIES*sizeof(Double_t));
    if(!CookdEdx(n, &dedx[il*n], &DE[0])) continue;

    s=0.;
    for(Int_t is(AliPID::kSPECIES); is--;){
      if((DE[0] > 0.) && (p[il] > 0.)) prLayer[is] = GetProbabilitySingleLayer(is, p[il], DE[0]);
      AliDebug(3, Form("Probability for Species %d in Layer %d: %f", is, il, prLayer[is]));
      s+=prLayer[is];
    }
    if(s<1.e-30){
      AliDebug(2, Form("Null all species prob layer[%d].", il));
      continue;
    }
    for(Int_t is(AliPID::kSPECIES); is--;){
      if(kNorm) prLayer[is] /= s;
      prob[is] *= prLayer[is];
    }
  }
  if(!kNorm) return kTRUE;

  s=0.;
  for(Int_t is(AliPID::kSPECIES); is--;) s+=prob[is];
  if(s<1.e-30){
    AliDebug(2, "Null total prob.");
    return kFALSE;
  }
  for(Int_t is(AliPID::kSPECIES); is--;) prob[is]/=s;
  return kTRUE;
}

//____________________________________________________________
Double_t AliTRDPIDResponse::GetProbabilitySingleLayer(Int_t species, Double_t plocal, Double_t dEdx) const {
  //
  // Get the non-normalized probability for a certain particle species as coming
  // from the reference histogram
  // Interpolation between momentum bins
  //
  AliDebug(1, Form("Make Probability for Species %d with Momentum %f", species, plocal));
  Int_t pbin = GetLowerMomentumBin(plocal);  
  AliDebug(1, Form("Bin %d", pbin));
  Double_t pLayer = 0.;
  // Do Interpolation exept for underflow and overflow
  if(pbin >= 0 && pbin < kNPBins-1){
    TH1 *refHistUpper = NULL, *refHistLower = NULL;
    if(fMapRefHists[species][pbin] >= 0)
      refHistLower = dynamic_cast<TH1 *>(fReferences->UncheckedAt(fMapRefHists[species][pbin]));
    if(fMapRefHists[species][pbin+1] >= 0)
      refHistUpper = dynamic_cast<TH1 *>(fReferences->UncheckedAt(fMapRefHists[species][pbin+1]));
    AliDebug(1, Form("Reference Histos (Upper/Lower): [%p|%p]", refHistUpper, refHistLower));
  
    if (refHistLower && refHistUpper ) {
      Double_t pLower = refHistLower->GetBinContent(refHistLower->GetXaxis()->FindBin(dEdx));
      Double_t pUpper = refHistUpper->GetBinContent(refHistUpper->GetXaxis()->FindBin(dEdx));
  
      pLayer = pLower + (pUpper - pLower)/(fgkPBins[pbin+1]-fgkPBins[pbin]) * (plocal - fgkPBins[pbin]); 
    }
  }
  else{
    TH1 *refHist = NULL;
    if(pbin < 0){
      // underflow
      if(fMapRefHists[species][0] >= 0){
        refHist = dynamic_cast<TH1 *>(fReferences->UncheckedAt(fMapRefHists[species][0])); 
      }
    } else {
      // overflow
      if(fMapRefHists[species][kNPBins-1] > 0)
        refHist = dynamic_cast<TH1 *>(fReferences->UncheckedAt(fMapRefHists[species][kNPBins-1])); 
    }
    if (refHist)
      pLayer = refHist->GetBinContent(refHist->GetXaxis()->FindBin(dEdx));
  }
  AliDebug(1, Form("Probability %f", pLayer));
  return pLayer;
}

//____________________________________________________________
Int_t AliTRDPIDResponse::GetLowerMomentumBin(Double_t p) const {
  //
  // Get the momentum bin for a given momentum value
  //
  Int_t bin = -1;
  if(p > fgkPBins[kNPBins-1]) return kNPBins-1;
  for(Int_t ibin = 0; ibin < kNPBins - 1; ibin++){
    if(p >= fgkPBins[ibin] && p < fgkPBins[ibin+1]){
      bin = ibin;
      break;
    }
  }
  return bin;
}

//____________________________________________________________
void AliTRDPIDResponse::SetOwner(){
  //
  // Make Deep Copy of the Reference Histograms
  //
  if(!fReferences || IsOwner()) return;
  TObjArray *ctmp = new TObjArray();
  for(Int_t ien = 0; ien < fReferences->GetEntriesFast(); ien++)
    ctmp->AddAt(fReferences->UncheckedAt(ien)->Clone(), ien);
  delete fReferences;
  fReferences = ctmp;
  SetBit(kIsOwner, kTRUE);
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::CookdEdx(Int_t nSlice, Double_t *in, Double_t *out) const
{
  switch(fPIDmethod){
  case kNN: // NN 
    break;
  case kLQ2D: // 2D LQ 
    break;
  case kLQ1D: // 1D LQ 
    out[0]= 0.;
    for(Int_t islice = 0; islice < nSlice; islice++) out[0] += in[islice] * fGainNormalisationFactor;
    if(out[0] < 1e-6) return kFALSE;
    break;
  default:
    return kFALSE;
  }
  return kTRUE;
}

