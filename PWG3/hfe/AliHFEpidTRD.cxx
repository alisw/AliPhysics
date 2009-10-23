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
/************************************************************************
 *                                                                      *
 * Class for TRD PID                                                    *
 * Implements the abstract base class AliHFEpidbase                     *
 * Make PID does the PID decision                                       *
 * Class further contains TRD specific cuts and QA histograms           *
 *                                                                      *
 * Authors:                                                             *
 *   Markus Fasel <M.Fasel@gsi.de>                                      *
 *                                                                      *
 ************************************************************************/
#include <TAxis.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TIterator.h>
#include <TKey.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>

#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliPID.h"

#include "AliHFEpidTRD.h"

ClassImp(AliHFEpidTRD)

const Double_t AliHFEpidTRD::fgkVerySmall = 1e-12;

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const char* name) :
    AliHFEpidBase(name)
  , fPIDMethod(kNN)
  , fContainer(0x0)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
}

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const AliHFEpidTRD &ref):
    AliHFEpidBase("")
  , fPIDMethod(kLQ)
  , fContainer(0x0)
{
  //
  // Copy constructor
  //
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
  ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidTRD &AliHFEpidTRD::operator=(const AliHFEpidTRD &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  }
  return *this;
}

//___________________________________________________________________
void AliHFEpidTRD::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidTRD &target = dynamic_cast<AliHFEpidTRD &>(ref);

  target.fPIDMethod = fPIDMethod;
  memcpy(target.fThreshParams, fThreshParams, sizeof(Double_t) * kThreshParams);
  if(fContainer) target.fContainer = dynamic_cast<TList *>(fContainer->Clone());
  AliHFEpidBase::Copy(ref);
}

//___________________________________________________________________
AliHFEpidTRD::~AliHFEpidTRD(){
  //
  // Destructor
  //
  if(fContainer){
    fContainer->Clear();
    delete fContainer;
  }
}

//______________________________________________________
Bool_t AliHFEpidTRD::InitializePID(){
  //
  // InitializePID: Load TRD thresholds and create the electron efficiency axis
  // to navigate 
  //
  InitParameters();
  return kTRUE;
}

//______________________________________________________
Int_t AliHFEpidTRD::IsSelected(AliHFEpidObject *track){
  //
  // Does PID for TRD alone:
  // PID thresholds based on 90% Electron Efficiency level approximated by a linear 
  // step function
  //
  if(track->fAnalysisType == AliHFEpidObject::kESDanalysis){
      AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(track->fRecTrack);
      if(!esdTrack) return 0;
      AliMCParticle *esdmc = dynamic_cast<AliMCParticle *>(track->fMCtrack);
      return MakePIDesd(esdTrack, esdmc);
   } else {
     AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(track->fRecTrack);
     if(!aodTrack) return 0;
     AliAODMCParticle *aodmc = dynamic_cast<AliAODMCParticle *>(track->fMCtrack);
     return MakePIDaod(aodTrack, aodmc);
   }
}

//______________________________________________________
Int_t AliHFEpidTRD::MakePIDaod(AliAODTrack * /*aodTrack*/, AliAODMCParticle * /*aodmc*/){
  AliError("AOD PID not yet implemented");
  return 0;
}

//______________________________________________________
Int_t AliHFEpidTRD::MakePIDesd(AliESDtrack *esdTrack, AliMCParticle * /*mcTrack*/){
  Double_t p = esdTrack->GetOuterParam() ? esdTrack->GetOuterParam()->P() : esdTrack->P();
  if(p < 2.0) return 0;

  Double_t pidProbs[AliPID::kSPECIES];
  esdTrack->GetTRDpid(pidProbs);
  if(IsQAon())FillHistogramsLikelihood(0, p, pidProbs[AliPID::kElectron]);
  Double_t threshold = GetTRDthresholds(0.91, p);
  AliDebug(1, Form("Threshold: %f\n", threshold));
  if(IsQAon()) (dynamic_cast<TH2F *>(fContainer->At(kHistTRDthresholds)))->Fill(p, threshold);
  if(pidProbs[AliPID::kElectron] > threshold){
    if(IsQAon()) FillHistogramsLikelihood(1, p, pidProbs[AliPID::kElectron]);
    return 11;
  }
  return 211;
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetTRDthresholds(Double_t electronEff, Double_t p){ 
  //
  // Return momentum dependent and electron efficiency dependent TRD thresholds
  // 
  Double_t params[4];
  GetParameters(electronEff, params);
  Double_t threshold = 1. - params[0] - params[1] * p - params[2] * TMath::Exp(-params[3] * p);
  return TMath::Max(TMath::Min(threshold, 0.99), 0.2); // truncate the threshold upperwards to 0.999 and lowerwards to 0.2 and exclude unphysical values
}

//___________________________________________________________________
void AliHFEpidTRD::InitParameters(){
  //
  // Fill the Parameters into an array
  //

  // Parameters for 6 Layers
  fThreshParams[0] = -0.001839; // 0.7 electron eff
  fThreshParams[1] = 0.000276;
  fThreshParams[2] = 0.044902; 
  fThreshParams[3] = 1.726751;
  fThreshParams[4] = -0.002405; // 0.75 electron eff
  fThreshParams[5] = 0.000372;
  fThreshParams[6] = 0.061775;
  fThreshParams[7] = 1.739371;
  fThreshParams[8] = -0.003178; // 0.8 electron eff
  fThreshParams[9] = 0.000521;
  fThreshParams[10] = 0.087585;
  fThreshParams[11] = 1.749154;
  fThreshParams[12] = -0.004058; // 0.85 electron eff
  fThreshParams[13] = 0.000748;
  fThreshParams[14] = 0.129583;
  fThreshParams[15] = 1.782323;
  fThreshParams[16] = -0.004967; // 0.9 electron eff
  fThreshParams[17] = 0.001216;
  fThreshParams[18] = 0.210128;
  fThreshParams[19] = 1.807665;
  fThreshParams[20] = -0.000996; // 0.95 electron eff
  fThreshParams[21] = 0.002627;
  fThreshParams[22] = 0.409099;
  fThreshParams[23] = 1.787076;
}

//___________________________________________________________________
void AliHFEpidTRD::GetParameters(Double_t electronEff, Double_t *parameters){
  //
  // return parameter set for the given efficiency bin
  //
  Int_t effbin = static_cast<Int_t>((electronEff - 0.7)/0.05);
  memcpy(parameters, fThreshParams + effbin * 4, sizeof(Double_t) * 4);
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetTRDSignalV1(AliESDtrack *track, Int_t mcPID){
  //
  // Calculation of the TRD Signal via truncated mean
  // Method 1: Take all Slices available
  // cut out 0s
  // Order them in increasing order
  // Cut out the upper third
  // Calculate mean over the last 2/3 slices
  //
  const Int_t kNSlices = 48;
  AliDebug(2, Form("Number of Tracklets: %d\n", track->GetTRDpidQuality()));
  Double_t trdSlices[kNSlices], tmp[kNSlices];
  Int_t indices[48];
  Int_t icnt = 0;
  for(Int_t idet = 0; idet < 6; idet++)
    for(Int_t islice = 0; islice < 8; islice++){
      AliDebug(2, Form("Chamber[%d], Slice[%d]: TRDSlice = %f", idet, islice, track->GetTRDslice(idet, islice)));
      if(TMath::Abs(track->GetTRDslice(idet, islice)) < fgkVerySmall) continue;;
      trdSlices[icnt++] = track->GetTRDslice(idet, islice);
    }
  AliDebug(1, Form("Number of Slices: %d\n", icnt));
  if(icnt < 6) return 0.;   // We need at least 6 Slices for the truncated mean
  TMath::Sort(icnt, trdSlices, indices, kFALSE);
  memcpy(tmp, trdSlices, sizeof(Double_t) * icnt);
  for(Int_t ien = 0; ien < icnt; ien++)
    trdSlices[ien] = tmp[indices[ien]];
  Double_t trdSignal = TMath::Mean(static_cast<Int_t>(static_cast<Double_t>(icnt) * 2./3.), trdSlices);
  Double_t mom = track->GetOuterParam() ? track->GetOuterParam()->P() : -1;
  AliDebug(1, Form("PID Meth. 1: p[%f], TRDSignal[%f]", mom, trdSignal));
  if(IsQAon() && mom > 0.) FillHistogramsTRDSignalV1(trdSignal, mom, mcPID);
  return trdSignal;
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetTRDSignalV2(AliESDtrack *track, Int_t mcPID){
  //
  // Calculation of the TRD Signal via truncated mean
  // Method 2: Take only first 5 slices per chamber
  // Order them in increasing order
  // Cut out upper half 
  // Now do mean with the reamining 3 slices per chamber
  //
  Double_t trdSlicesLowTime[30], trdSlicesRemaining[32];
  Int_t indices[30];
  Int_t cntLowTime=0, cntRemaining = 0;
  for(Int_t idet = 0; idet < 6; idet++)
    for(Int_t islice = 0; islice < 8; islice++){
      if(TMath::Abs(track->GetTRDslice(idet, islice)) < fgkVerySmall) continue;;
      if(islice < 5){
        AliDebug(2, Form("Part 1, Det[%d], Slice[%d], TRDSlice: %f", idet, islice, track->GetTRDslice(idet, islice)));
        trdSlicesLowTime[cntLowTime++] = track->GetTRDslice(idet, islice);
      } else{
        AliDebug(2, Form("Part 1, Det[%d], Slice[%d], TRDSlice: %f", idet, islice, track->GetTRDslice(idet, islice)));
        trdSlicesRemaining[cntRemaining++] = track->GetTRDslice(idet, islice);
      }
    }
  if(cntLowTime < 4 || cntRemaining < 2) return 0.; // Min. Number of Slices at high time is 2 (matches with 1 layer), for the truncated mean we need at least 4 Slices
  TMath::Sort(cntLowTime, trdSlicesLowTime, indices, kFALSE);
  // Fill the second array with the lower half of the first time bins
  for(Int_t ien = 0; ien < static_cast<Int_t>(static_cast<Double_t>(cntLowTime) * 0.5); ien++)
    trdSlicesRemaining[cntRemaining++] = trdSlicesLowTime[indices[ien]];
  Double_t trdSignal = TMath::Mean(cntRemaining, trdSlicesRemaining);
  Double_t mom = track->GetOuterParam() ? track->GetOuterParam()->P() : -1;
  AliDebug(1, Form("PID Meth. 2: p[%f], TRDSignal[%f]", mom, trdSignal));
  if(IsQAon() && mom > 0.) FillHistogramsTRDSignalV2(trdSignal, mom, mcPID);
  return trdSignal;
}

//___________________________________________________________________
void AliHFEpidTRD::FillHistogramsTRDSignalV1(Double_t signal, Double_t mom, Int_t species){
  //
  // Fill histograms for TRD Signal from Method 1 vs. p for different particle species
  //
  (dynamic_cast<TH2F *>(fContainer->At(kHistTRDSigV1)))->Fill(mom, signal);
  if(species < 0 || species >= AliPID::kSPECIES) return;
  (dynamic_cast<TH2F *>(fContainer->At(kHistOverallSpecies + species)))->Fill(mom, signal);
}

//___________________________________________________________________
void AliHFEpidTRD::FillHistogramsTRDSignalV2(Double_t signal, Double_t mom, Int_t species){
  //
  // Fill histograms for TRD Signal from Method 2 vs. p for different particle species
  //
  (dynamic_cast<TH2F *>(fContainer->At(kHistTRDSigV2)))->Fill(mom, signal);  
  if(species < 0 || species >= AliPID::kSPECIES) return;
  (dynamic_cast<TH2F *>(fContainer->At(kHistOverallSpecies + AliPID::kSPECIES + species)))->Fill(mom, signal);
}

//___________________________________________________________________
void AliHFEpidTRD::AddQAhistograms(TList *l){
  //
  // Adding QA histograms for the TRD PID
  // QA histograms are:
  // + TRD Signal from Meth. 1 vs p for all species
  // + TRD Signal from Meth. 2 vs p for all species
  // + For each species
  //    - TRD Signal from Meth. 1 vs p
  //    - TRD Signal from Meth. 2 vs p
  //
  const Int_t kMomentumBins = 41;
  const Double_t kPtMin = 0.1;
  const Double_t kPtMax = 10.;
  const Int_t kSigBinsMeth1 = 100; 
  const Int_t kSigBinsMeth2 = 100;
  const Double_t kMinSig = 0.;
  const Double_t kMaxSigMeth1 = 10000.;
  const Double_t kMaxSigMeth2 = 10000.;
  
  if(!fContainer) fContainer = new TList;
  fContainer->SetName("fTRDqaHistograms");

  Double_t momentumBins[kMomentumBins];
  for(Int_t ibin = 0; ibin < kMomentumBins; ibin++)
    momentumBins[ibin] = static_cast<Double_t>(TMath::Power(10,TMath::Log10(kPtMin) + (TMath::Log10(kPtMax)-TMath::Log10(kPtMin))/(kMomentumBins-1)*static_cast<Double_t>(ibin)));
  // Likelihood Histograms
  fContainer->AddAt(new TH2F("fTRDlikeBefore", "TRD Electron Likelihood before cut", kMomentumBins - 1, momentumBins, 1000, 0, 1), kHistTRDlikeBefore);
  fContainer->AddAt(new TH2F("fTRDlikeAfter", "TRD Electron Likelihood after cut", kMomentumBins - 1, momentumBins, 1000, 0, 1), kHistTRDlikeAfter);
  fContainer->AddAt(new TH2F("fTRDthesholds", "TRD Electron thresholds", kMomentumBins - 1, momentumBins, 1000, 0, 1), kHistTRDthresholds);
  // Signal Histograms
  fContainer->AddAt(new TH2F("fTRDSigV1all", "TRD Signal (all particles, Method 1)", kMomentumBins - 1, momentumBins, kSigBinsMeth1, kMinSig, kMaxSigMeth1), kHistTRDSigV1);
  fContainer->AddAt(new TH2F("fTRDSigV2all", "TRD Signal (all particles, Method 2)", kMomentumBins - 1, momentumBins, kSigBinsMeth2, kMinSig, kMaxSigMeth2), kHistTRDSigV2);
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    fContainer->AddAt(new TH2F(Form("fTRDSigV1%s", AliPID::ParticleName(ispec)), Form("TRD Signal (%s, Method 1)", AliPID::ParticleName(ispec)), kMomentumBins - 1, momentumBins, kSigBinsMeth1, kMinSig, kMaxSigMeth1), kHistOverallSpecies + ispec);
    fContainer->AddAt(new TH2F(Form("fTRDSigV2%s", AliPID::ParticleName(ispec)), Form("TRD Signal (%s, Method 2)", AliPID::ParticleName(ispec)), kMomentumBins - 1, momentumBins, kSigBinsMeth1, kMinSig, kMaxSigMeth2), kHistOverallSpecies + AliPID::kSPECIES + ispec);
  }
  l->AddLast(fContainer);
}

//___________________________________________________________________
void AliHFEpidTRD::FillHistogramsLikelihood(Int_t whenFilled, Float_t p, Float_t elProb){
  //
  // Fill Likelihood Histogram before respectively after decision
  //
  TH2F *histo = NULL;
  if(whenFilled)
    histo = dynamic_cast<TH2F *>(fContainer->At(kHistTRDlikeAfter));
  else
    histo = dynamic_cast<TH2F *>(fContainer->At(kHistTRDlikeBefore));
  if(!histo){
    AliError("QA histograms not found");
    return;
  }
  histo->Fill(p, elProb);
}
