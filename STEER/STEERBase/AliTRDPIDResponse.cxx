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
#include <TAxis.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TROOT.h> 
#include <TString.h>
#include <TSystem.h>
#include <TVectorT.h>

#include "AliLog.h"

#include "AliTRDPIDReference.h"
#include "AliTRDPIDResponse.h"

ClassImp(AliTRDPIDResponse)

//____________________________________________________________
AliTRDPIDResponse::AliTRDPIDResponse():
  TObject()
  ,fkPIDReference(NULL)
  ,fkPIDParams(NULL)
  ,fGainNormalisationFactor(1.)
  ,fPIDmethod(kLQ1D)
{
  //
  // Default constructor
  //
}

//____________________________________________________________
AliTRDPIDResponse::AliTRDPIDResponse(const AliTRDPIDResponse &ref):
  TObject(ref)
  ,fkPIDReference(ref.fkPIDReference)
  ,fkPIDParams(ref.fkPIDParams)
  ,fGainNormalisationFactor(ref.fGainNormalisationFactor)
  ,fPIDmethod(ref.fPIDmethod)
{
  //
  // Copy constructor
  //
}

//____________________________________________________________
AliTRDPIDResponse &AliTRDPIDResponse::operator=(const AliTRDPIDResponse &ref){
  //
  // Assignment operator
  //
  if(this == &ref) return *this;
  
  // Make copy
  TObject::operator=(ref);
  fGainNormalisationFactor = ref.fGainNormalisationFactor;
  fkPIDReference = ref.fkPIDReference;
  fkPIDParams = ref.fkPIDParams;
  fPIDmethod = ref.fPIDmethod;
  
  return *this;
}

//____________________________________________________________
AliTRDPIDResponse::~AliTRDPIDResponse(){
  //
  // Destructor
  //
  if(IsOwner()) delete fkPIDReference;
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::Load(const Char_t * filename, const Char_t *refName){
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
  fkPIDReference = dynamic_cast<const AliTRDPIDReference *>(in->Get(refName)->Clone());
  in->Close(); delete in;
  owd->cd();
  SetBit(kIsOwner, kTRUE);
  AliDebug(2, Form("Successfully loaded References for %d Momentum bins", fkPIDReference->GetNumberOfMomentumBins()));
  return kTRUE;
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::GetResponse(Int_t n, const Double_t * const dedx, const Float_t * const p, Double_t prob[AliPID::kSPECIES], Bool_t kNorm) const
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

  if(!fkPIDReference){
    AliWarning("Missing reference data. PID calculation not possible.");
    return kFALSE;
  }

  for(Int_t is(AliPID::kSPECIES); is--;) prob[is]=.2;
  Double_t prLayer[AliPID::kSPECIES];
  Double_t dE[10], s(0.);
  for(Int_t il(kNlayer); il--;){
    memset(prLayer, 0, AliPID::kSPECIES*sizeof(Double_t));
    if(!CookdEdx(n, &dedx[il*n], &dE[0])) continue;

    s=0.;
    for(Int_t is(AliPID::kSPECIES); is--;){
      if((dE[0] > 0.) && (p[il] > 0.)) prLayer[is] = GetProbabilitySingleLayer(is, p[il], dE[0]);
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
  Float_t pLower, pUpper;
  Double_t probLayer = 0.;
  TH1 *refUpper = dynamic_cast<TH1 *>(fkPIDReference->GetUpperReference((AliPID::EParticleType)species, plocal, pUpper)),
      *refLower = dynamic_cast<TH1 *>(fkPIDReference->GetLowerReference((AliPID::EParticleType)species, plocal, pLower));
  // Do Interpolation exept for underflow and overflow
  if(refLower && refUpper){
    Double_t probLower = refLower->GetBinContent(refLower->GetXaxis()->FindBin(dEdx));
    Double_t probUpper = refUpper->GetBinContent(refUpper->GetXaxis()->FindBin(dEdx));
  
    probLayer = probLower + (probUpper - probLower)/(pUpper-pLower) * (plocal - pLower);
  } else if(refLower){
    // underflow
    probLayer = refLower->GetBinContent(refLower->GetXaxis()->FindBin(dEdx));
  } else if(refUpper){
    // overflow
    probLayer = refUpper->GetBinContent(refUpper->GetXaxis()->FindBin(dEdx));
  } else {
    AliError("No references available");
  }
  AliDebug(1, Form("Probability %f", probLayer));
  return probLayer;
}

//____________________________________________________________
void AliTRDPIDResponse::SetOwner(){
  //
  // Make Deep Copy of the Reference Histograms
  //
  if(!fkPIDReference || IsOwner()) return;
  const AliTRDPIDReference *tmp = fkPIDReference;
  fkPIDReference = dynamic_cast<const AliTRDPIDReference *>(tmp->Clone());
  SetBit(kIsOwner, kTRUE);
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::CookdEdx(Int_t nSlice, const Double_t * const in, Double_t *out) const
{
	//
	// Recalculate dE/dx
	//
  switch(fPIDmethod){
  case kNN: // NN 
    break;
  case kLQ2D: // 2D LQ 
    break;
  case kLQ1D: // 1D LQ 
    out[0]= 0.;
    for(Int_t islice = 0; islice < nSlice; islice++) 
      if(in[islice] > 0) out[0] += in[islice] * fGainNormalisationFactor;   // Protect against negative values for slices having no dE/dx information
    if(out[0] < 1e-6) return kFALSE;
    break;
  default:
    return kFALSE;
  }
  return kTRUE;
}

//____________________________________________________________
Bool_t AliTRDPIDResponse::IdentifiedAsElectron(Int_t nTracklets, const Double_t *like, Double_t p, Double_t level) const {
  //
  // Check whether particle is identified as electron assuming a certain electron efficiency level
  // Only electron and pion hypothesis is taken into account
  //
  // Inputs:
  //         Number of tracklets
  //         Likelihood values
  //         Momentum
  //         Electron efficiency level
  //
  // If the function fails when the params are not accessible, the function returns true
  //
  if(!fkPIDParams){
    AliError("No PID Param object available");
    return kTRUE;
  } 
  Double_t probEle = like[AliPID::kElectron]/(like[AliPID::kElectron] + like[AliPID::kPion]);
  const TVectorD *params = GetParams(nTracklets, level);
  if(!params){
    AliError("No Params found for the given configuration");
    return kTRUE;
  }
  Double_t threshold = 1. - (*params)[0] - (*params)[1] * p - (*params)[2] * TMath::Exp(-(*params)[3] * p);
  if(probEle > TMath::Max(TMath::Min(threshold, 0.99), 0.2)) return kTRUE; // truncate the threshold upperwards to 0.999 and lowerwards to 0.2 and exclude unphysical values
  return kFALSE;
}

//____________________________________________________________
const TVectorD* AliTRDPIDResponse::GetParams(Int_t ntracklets, Double_t level) const {
  //
  // returns the threshold for a given number of tracklets and a given efficiency level
  //tby definition the lower of step is given.
  //
  if(ntracklets > 6 || ntracklets <=0) return NULL;
  TObjArray * entry = dynamic_cast<TObjArray *>(fkPIDParams->At(ntracklets - 1));
  if(!entry) return NULL;
  
  TObjArray*cut = NULL;
  TVectorF *effLevel = NULL; const TVectorD *parameters = NULL;
  Float_t currentLower = 0.;
  TIter cutIter(entry);
  while((cut = dynamic_cast<TObjArray *>(cutIter()))){
    effLevel = static_cast<TVectorF *>(cut->At(0));
    if((*effLevel)[0] > currentLower && (*effLevel)[0] <= level){
      // New Lower entry found
      parameters = static_cast<const TVectorD *>(cut->At(1));
      currentLower = (*effLevel)[0];
    }
  }  
  AliDebug(2, Form("Taking params for %d tracklets and %f electron Efficiency\n", ntracklets, currentLower));

  return parameters;
}
