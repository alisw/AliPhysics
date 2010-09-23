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
// Class for TRD PID
// Implements the abstract base class AliHFEpidbase 
// Make PID does the PID decision 
// Class further contains TRD specific cuts and QA histograms 
//  
// Authors: 
//   Markus Fasel <M.Fasel@gsi.de>  
// 
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TString.h>

#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliPID.h"

#include "AliHFEcollection.h"
#include "AliHFEpidTRD.h"

ClassImp(AliHFEpidTRD)

const Double_t AliHFEpidTRD::fgkVerySmall = 1e-12;

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD() :
    AliHFEpidBase()
  , fMinP(1.)
  , fElectronEfficiency(0.91)
  , fPIDMethod(kNN)
  , fContainer(NULL)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
}

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const char* name) :
    AliHFEpidBase(name)
  , fMinP(1.)
  , fElectronEfficiency(0.91)
  , fPIDMethod(kNN)
  , fContainer(NULL)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
}

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const AliHFEpidTRD &ref):
    AliHFEpidBase("")
  , fMinP(ref.fMinP)
  , fElectronEfficiency(ref.fElectronEfficiency)
  , fPIDMethod(ref.fPIDMethod)
  , fContainer(NULL)
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

  target.fMinP = fMinP;
  target.fPIDMethod = fPIDMethod;
  target.fElectronEfficiency = fElectronEfficiency;
  memcpy(target.fThreshParams, fThreshParams, sizeof(Double_t) * kThreshParams);
  if(fContainer) target.fContainer = dynamic_cast<AliHFEcollection *>(fContainer->Clone());
  AliHFEpidBase::Copy(ref);
}

//___________________________________________________________________
AliHFEpidTRD::~AliHFEpidTRD(){
  //
  // Destructor
  //
  if(fContainer)
    delete fContainer;
}

//______________________________________________________
Bool_t AliHFEpidTRD::InitializePID(){
  //
  // InitializePID: Load TRD thresholds and create the electron efficiency axis
  // to navigate 
  //
  if(fPIDMethod == kLQ)
    InitParameters1DLQ();
  else
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
  //
  // Does PID decision for AOD tracks as discribed above
  //
  AliError("AOD PID not yet implemented");
  return 0;
}

//______________________________________________________
Int_t AliHFEpidTRD::MakePIDesd(AliESDtrack *esdTrack, AliMCParticle * /*mcTrack*/){
  //
  // Does PID decision for ESD tracks as discribed above
  //
  Double_t p = esdTrack->GetOuterParam() ? esdTrack->GetOuterParam()->P() : esdTrack->P();
  if(IsQAon())
    FillStandardQA(0, esdTrack);
  if(p < fMinP) return 0;

  Double_t pidProbs[AliPID::kSPECIES];
  esdTrack->GetTRDpid(pidProbs);
  Double_t threshold = GetTRDthresholds(fElectronEfficiency, p);
  AliDebug(1, Form("Threshold: %f\n", threshold));
  if(IsQAon()) fContainer->Fill("fTRDthresholds", p, threshold);
  if(pidProbs[AliPID::kElectron] > threshold){
    if(IsQAon()) 
      FillStandardQA(1, esdTrack);
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

  AliDebug(2, "Loading threshold Parameter");
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
void AliHFEpidTRD::InitParameters1DLQ(){
  //
  // Init Parameters for 1DLQ PID (M. Fasel, Sept. 6th, 2010)
  //

  // Parameters for 6 Layers
  AliDebug(2, Form("Loading threshold parameter for Method 1DLQ"));
  fThreshParams[0] = -0.02241; // 0.7 electron eff
  fThreshParams[1] = 0.05043;
  fThreshParams[2] = 0.7925; 
  fThreshParams[3] = 2.625;
  fThreshParams[4] = 0.07438; // 0.75 electron eff
  fThreshParams[5] = 0.05158;
  fThreshParams[6] = 2.864;
  fThreshParams[7] = 4.356;
  fThreshParams[8] = 0.1977; // 0.8 electron eff
  fThreshParams[9] = 0.05956;
  fThreshParams[10] = 2.853;
  fThreshParams[11] = 3.713;
  fThreshParams[12] = 0.5206; // 0.85 electron eff
  fThreshParams[13] = 0.03077;
  fThreshParams[14] = 2.966;
  fThreshParams[15] = 4.07;
  fThreshParams[16] = 0.8808; // 0.9 electron eff
  fThreshParams[17] = 0.002092;
  fThreshParams[18] = 1.17;
  fThreshParams[19] = 4.506;
  fThreshParams[20] = 1.; // 0.95 electron eff
  fThreshParams[21] = 0.;
  fThreshParams[22] = 0.;
  fThreshParams[23] = 0.;

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
  AliDebug(3, Form("Number of Tracklets: %d\n", track->GetTRDpidQuality()));
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
  AliDebug(3, Form("PID Meth. 1: p[%f], TRDSignal[%f]", mom, trdSignal));
  if(IsQAon() && mom > 0.) FillHistogramsTRDSignal(trdSignal, mom, mcPID, 1);
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
        AliDebug(3, Form("Part 1, Det[%d], Slice[%d], TRDSlice: %f", idet, islice, track->GetTRDslice(idet, islice)));
        trdSlicesLowTime[cntLowTime++] = track->GetTRDslice(idet, islice);
      } else{
        AliDebug(3, Form("Part 1, Det[%d], Slice[%d], TRDSlice: %f", idet, islice, track->GetTRDslice(idet, islice)));
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
  AliDebug(3, Form("PID Meth. 2: p[%f], TRDSignal[%f]", mom, trdSignal));
  if(IsQAon() && mom > 0.) FillHistogramsTRDSignal(trdSignal, mom, mcPID, 2);
  return trdSignal;
}

//___________________________________________________________________
void AliHFEpidTRD::FillHistogramsTRDSignal(Double_t signal, Double_t p, Int_t species, UInt_t version){
  //
  // Fill histograms for TRD Signal from Method 2 vs. p for different particle species
  // Non-standard QA content
  //
  if(version == 0 || version > 2) return;
  if(species >= AliPID::kSPECIES || species < -1) species = -1;
  Double_t content[3] = {species, p, signal};
  Char_t histname[256]; sprintf(histname, "fTRDsignalV%d", version);
  THnSparseF *hsig = dynamic_cast<THnSparseF *>(fContainer->Get(histname));
  if(hsig) hsig->Fill(content);
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
  const Int_t kMomentumBins = 100;
  const Double_t kPtMin = 0.1;
  const Double_t kPtMax = 20.;
  
  if(!fContainer) fContainer = new AliHFEcollection("fQAhistosTRD", "TRD QA histos");
  fContainer->CreateTH2F("fTRDlikeBefore", "TRD Electron Likelihood before cut; p [GeV/c]; TRD Electron Likelihood", kMomentumBins, kPtMin, kPtMax, 1000, 0, 1);
  fContainer->CreateTH2F("fTRDlikeAfter", "TRD Electron Likelihood after cut; p [GeV/c]; TRD Electron Likelihood", kMomentumBins, kPtMin, kPtMax, 1000, 0, 1);
  fContainer->CreateTH2F("fTRDthresholds", "TRD Electron Thresholds; p [GeV/c]; Thresholds", kMomentumBins, kPtMin, kPtMax, 1000, 0., 1.);
  fContainer->CreateTH2F("fTPCsignalBefore", "TPC dEdx Spectrum before TRD PID; p [GeV/c]; TPC Signal [a.u.]", kMomentumBins, kPtMin, kPtMax, 100, 0., 100.);
  fContainer->CreateTH2F("fTPCsignalAfter", "TPC dEdx Spectrum before TRD PID; p [GeV/c]; TPC Signal [a.u.]", kMomentumBins, kPtMin, kPtMax, 100, 0., 100.);
  fContainer->CreateTH2F("fTPCsigmaBefore", "TPC dEdx before TRD PID; p [GeV/c]; Normalized TPC distance to the electron line [n#sigma]", kMomentumBins, kPtMin, kPtMax, 100, -10., 10.);
  fContainer->CreateTH2F("fTPCsigmaAfter", "TPC dEdx Spectrum before TRD PID; p [GeV/c]; Normalized TPC distance to the electron line [n#sigma]", kMomentumBins, kPtMin, kPtMax, 100, -10., 10.);

  // Monitor THnSparse for TRD Signal
  const Int_t kBinsTRDsignal = 3; 
  Int_t nBins[kBinsTRDsignal] = {AliPID::kSPECIES +1, kMomentumBins, 100};
  Double_t binMin[kBinsTRDsignal] = {-1, kPtMin, 0};
  Double_t binMax[kBinsTRDsignal] = {AliPID::kSPECIES, kPtMax, 1000};
  fContainer->CreateTHnSparse("fTRDsignalV1", "TRD Signal V1", kBinsTRDsignal, nBins, binMin, binMax);
  fContainer->CreateTHnSparse("fTRDsignalV2", "TRD Signal V2", kBinsTRDsignal, nBins, binMin, binMax);

  // Make logatithmic binning
  TString hnames[7] = {"fTRDlikeBefore", "fTRDlikeAfter", "fTRDthresholds", "fTRDsignalBefore", "fTRDsignalAfter", "fTPCsigmaBefore", "fTPCsigmaAfter"};
  for(Int_t ihist = 0; ihist < 7; ihist++)
    fContainer->BinLogAxis(hnames[ihist].Data(), 0);
  fContainer->BinLogAxis("fTRDsignalV1", 1);
  fContainer->BinLogAxis("fTRDsignalV2", 1);

  l->Add(fContainer->GetList());
}

//___________________________________________________________________
void AliHFEpidTRD::FillStandardQA(Int_t whenFilled, AliESDtrack *esdTrack){
  //
  // Fill Likelihood Histogram before respectively after decision
  //
  Double_t p =  esdTrack->P();
  Double_t like[AliPID::kSPECIES];
  esdTrack->GetTRDpid(like);
  TString step;
  if(whenFilled)
    step = "After";
  else
    step = "Before";
  TString histos[3] = {"fTRDlike", "fTPCsignal", "fTPCsigma"};
  for(Int_t ihist = 0; ihist < 3; ihist++) histos[ihist] += step;
  fContainer->Fill(histos[0].Data(), esdTrack->P(), like[AliPID::kElectron]);
  fContainer->Fill(histos[1].Data(), p, esdTrack->GetTPCsignal());
  const Double_t sigmaShift = 1.;
  if(fESDpid)
    fContainer->Fill(histos[2].Data(), p, fESDpid->NumberOfSigmasTPC(esdTrack, AliPID::kElectron) - sigmaShift);
}
