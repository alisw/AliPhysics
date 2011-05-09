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
// Class further contains TRD specific cuts and QA histograms 
//  
// Authors:
//   Jens Wiechula <Jens.Wiechula@cern.ch>
//   Markus Fasel <M.Fasel@gsi.de>  
// 
#include <TList.h>
#include <TString.h>

#include <AliVParticle.h>
#include <AliESDtrack.h>
#include <AliLog.h>
#include <AliPID.h>

#include "AliDielectronTRDpidCut.h"

ClassImp(AliDielectronTRDpidCut)

const Double_t AliDielectronTRDpidCut::fgkVerySmall = 1e-12;

//___________________________________________________________________
AliDielectronTRDpidCut::AliDielectronTRDpidCut() :
    AliAnalysisCuts()
  , fMinP(1.)
  , fElectronEfficiency(0.91)
  , fPIDMethod(kNN)
  , fRequirePIDbit(0)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
  SetUseDefaultParameters();
}

//___________________________________________________________________
AliDielectronTRDpidCut::AliDielectronTRDpidCut(const char* name) :
    AliAnalysisCuts(name,name)
  , fMinP(1.)
  , fElectronEfficiency(0.91)
  , fPIDMethod(kNN)
  , fRequirePIDbit(0)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
  SetUseDefaultParameters();
}

//___________________________________________________________________
AliDielectronTRDpidCut::AliDielectronTRDpidCut(const AliDielectronTRDpidCut &ref):
    AliAnalysisCuts(ref.GetName(),ref.GetTitle())
  , fMinP(ref.fMinP)
  , fElectronEfficiency(ref.fElectronEfficiency)
  , fPIDMethod(ref.fPIDMethod)
  , fRequirePIDbit(ref.fRequirePIDbit)
{
  //
  // Copy constructor
  //
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
  ref.Copy(*this);
}

//___________________________________________________________________
AliDielectronTRDpidCut &AliDielectronTRDpidCut::operator=(const AliDielectronTRDpidCut &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  }
  return *this;
}

//___________________________________________________________________
void AliDielectronTRDpidCut::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliDielectronTRDpidCut &target = dynamic_cast<AliDielectronTRDpidCut &>(ref);
  
  Bool_t defaultParameters = UseDefaultParameters();
  target.SetUseDefaultParameters(defaultParameters);
  target.fMinP = fMinP;
  target.fPIDMethod = fPIDMethod;
  target.fRequirePIDbit = fRequirePIDbit;
  target.fElectronEfficiency = fElectronEfficiency;
  memcpy(target.fThreshParams, fThreshParams, sizeof(Double_t) * kThreshParams);
  AliAnalysisCuts::Copy(ref);
}

//___________________________________________________________________
AliDielectronTRDpidCut::~AliDielectronTRDpidCut(){
  //
  // Destructor
  //
}

//______________________________________________________
Bool_t AliDielectronTRDpidCut::InitializePID(){
  //
  // InitializePID: Load TRD thresholds and create the electron efficiency axis
  // to navigate 
  //
  if(UseDefaultParameters()){
    if(fPIDMethod == kLQ)
      InitParameters1DLQ();
    else
      InitParameters();
  }
  return kTRUE;
}

//______________________________________________________
Bool_t AliDielectronTRDpidCut::IsSelected(TObject *track) {
  //
  // Does PID for TRD alone:
  // PID thresholds based on 90% Electron Efficiency level approximated by a linear 
  // step function
  //
  const AliESDtrack *part = dynamic_cast<const AliESDtrack *>(track);
  if (!part) return kFALSE;
  
  if (fRequirePIDbit==AliDielectronTRDpidCut::kRequire&&!(part->GetStatus()&AliESDtrack::kTRDpid)) return kFALSE;
  if (fRequirePIDbit==AliDielectronTRDpidCut::kIfAvailable&&!(part->GetStatus()&AliESDtrack::kTRDpid)) return kTRUE;
  
  Double_t p = GetP((AliVParticle*)track);
  if(p < fMinP){ 
    AliDebug(1, Form("Track momentum below %f", fMinP));
    if(fRequirePIDbit==AliDielectronTRDpidCut::kRequire) return kFALSE;
	if(fRequirePIDbit==AliDielectronTRDpidCut::kIfAvailable) return kTRUE;
  }

  Double_t electronLike = GetElectronLikelihood((AliVParticle*)track);
  Double_t threshold = GetTRDthresholds(fElectronEfficiency, p);
  AliDebug(1, Form("Threshold: %f\n", threshold));
  if(electronLike > threshold){
    return kTRUE;
  }
  return kFALSE;

}

//___________________________________________________________________
Double_t AliDielectronTRDpidCut::GetTRDthresholds(Double_t electronEff, Double_t p) const {
  //
  // Return momentum dependent and electron efficiency dependent TRD thresholds
  // 
  Double_t params[4];
  GetParameters(electronEff, params);
  Double_t threshold = 1. - params[0] - params[1] * p - params[2] * TMath::Exp(-params[3] * p);
  return TMath::Max(TMath::Min(threshold, 0.99), 0.2); // truncate the threshold upperwards to 0.999 and lowerwards to 0.2 and exclude unphysical values
}

//___________________________________________________________________
void AliDielectronTRDpidCut::SetThresholdParameters(Double_t electronEff, Double_t *params){
  //
  // Set threshold parameters for the given bin
  //
  if(electronEff >= 1. || electronEff < 0.7) return;
  Int_t effbin = static_cast<Int_t>((electronEff - 0.7)/0.05); 
  memcpy(&fThreshParams[effbin * 4], params, sizeof(Double_t) * 4); 
  SetUseDefaultParameters(kFALSE);
}

//___________________________________________________________________
void AliDielectronTRDpidCut::InitParameters(){
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
void AliDielectronTRDpidCut::InitParameters1DLQ(){
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
void AliDielectronTRDpidCut::RenormalizeElPi(const Double_t *likein, Double_t *likeout) const {
  //
  // Renormalize likelihoods for electrons and pions neglecting the 
  // likelihoods for protons, kaons and muons
  //
  memset(likeout, 0, sizeof(Double_t) * AliPID::kSPECIES);
  Double_t norm = likein[AliPID::kElectron] + likein[AliPID::kPion];
  if(norm == 0.) norm = 1.;   // Safety
  likeout[AliPID::kElectron] = likein[AliPID::kElectron] / norm;
  likeout[AliPID::kPion] = likein[AliPID::kPion] / norm;
}

//___________________________________________________________________
void AliDielectronTRDpidCut::GetParameters(Double_t electronEff, Double_t *parameters) const {
  //
  // return parameter set for the given efficiency bin
  //
  Int_t effbin = static_cast<Int_t>((electronEff - 0.7)/0.05);
  memcpy(parameters, fThreshParams + effbin * 4, sizeof(Double_t) * 4);
}

//___________________________________________________________________
Double_t AliDielectronTRDpidCut::GetElectronLikelihood(const AliVParticle *track) const {
  //
  // Get TRD likelihoods for ESD respectively AOD tracks
  //
  Double_t pidProbs[AliPID::kSPECIES]; memset(pidProbs, 0, sizeof(Double_t) * AliPID::kSPECIES);
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  if(esdtrack) esdtrack->GetTRDpid(pidProbs);
  if(!IsRenormalizeElPi()) return pidProbs[AliPID::kElectron];
  Double_t probsNew[AliPID::kSPECIES];
  RenormalizeElPi(pidProbs, probsNew);
  return probsNew[AliPID::kElectron];
}

//___________________________________________________________________
Double_t AliDielectronTRDpidCut::GetP(const AliVParticle *track) const {
  //
  // Get the Momentum in the TRD
  //
  Double_t p = 0.;
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  if(esdtrack) p = esdtrack->GetOuterParam() ? esdtrack->GetOuterParam()->P() : esdtrack->P();
  return p;
}

//___________________________________________________________________
Double_t AliDielectronTRDpidCut::GetChargeLayer(const AliVParticle *track, UInt_t layer) const {
  //
  // Get the Charge in a single TRD layer
  //
  if(layer >= 6) return 0.;
  Double_t charge = 0.;
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  if(esdtrack)
    for(Int_t islice = 0; islice < esdtrack->GetNumberOfTRDslices(); islice++)
      charge += esdtrack->GetTRDslice(static_cast<UInt_t>(layer), islice);
  return charge;
}

//___________________________________________________________________
void AliDielectronTRDpidCut::GetTRDmomenta(const AliVParticle *track, Double_t *mom) const {
  //
  // Fill Array with momentum information at the TRD tracklet
  //
  const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
  if(esdtrack)
    for(Int_t itl = 0; itl < 6; itl++)
      mom[itl] = esdtrack->GetTRDmomentum(itl);
}

//___________________________________________________________________
Double_t AliDielectronTRDpidCut::GetTRDSignalV1(const AliESDtrack *track, Float_t truncation) const {
  //
  // Calculation of the TRD Signal via truncated mean
  // Method 1: Take all Slices available
  // cut out 0s
  // Order them in increasing order
  // Cut out the upper third
  // Calculate mean over the last 2/3 slices
  //
  const Int_t kNSlices = 48;
  const Int_t kSlicePerLayer = 7;
  // Weight the slice to equalize the MPV of the dQ/dl-distribution per slice to the one in the first slice
  // Pions are used as reference for the equalization
  const Double_t kWeightSlice[8] = {1., 2.122, 1.8, 1.635, 1.595, 1.614, 1.16, 7.0};
  AliDebug(3, Form("Number of Tracklets: %d\n", track->GetTRDntrackletsPID()));
  Double_t trdSlices[kNSlices], tmp[kNSlices];
  Int_t indices[48];
  Int_t icnt = 0;
  for(Int_t idet = 0; idet < 6; idet++)
    for(Int_t islice = 0; islice < kSlicePerLayer; islice++){
      AliDebug(2, Form("Chamber[%d], Slice[%d]: TRDSlice = %f", idet, islice, track->GetTRDslice(idet, islice)));
      if(TMath::Abs(track->GetTRDslice(idet, islice)) < fgkVerySmall) continue;;
      trdSlices[icnt++] = track->GetTRDslice(idet, islice) * kWeightSlice[islice];
    }
  AliDebug(1, Form("Number of Slices: %d\n", icnt));
  if(icnt < 6) return 0.;   // We need at least 6 Slices for the truncated mean
  TMath::Sort(icnt, trdSlices, indices, kFALSE);
  memcpy(tmp, trdSlices, sizeof(Double_t) * icnt);
  for(Int_t ien = 0; ien < icnt; ien++)
    trdSlices[ien] = tmp[indices[ien]];
  Double_t trdSignal = TMath::Mean(static_cast<Int_t>(static_cast<Float_t>(icnt) * truncation), trdSlices);
  Double_t mom = track->GetOuterParam() ? track->GetOuterParam()->P() : -1;
  AliDebug(3, Form("PID Meth. 1: p[%f], TRDSignal[%f]", mom, trdSignal));
  return trdSignal;
}

//___________________________________________________________________
Double_t AliDielectronTRDpidCut::GetTRDSignalV2(const AliESDtrack *track, Float_t truncation) const {
  //
  // Calculation of the TRD Signal via truncated mean
  // Method 2: Take only first 5 slices per chamber
  // Order them in increasing order
  // Cut out upper half 
  // Now do mean with the reamining 3 slices per chamber
  //
  const Double_t kWeightSlice[8] = {1., 2.122, 1.8, 1.635, 1.595, 1.614, 1.16, 7.0};
  const Int_t kLayers = 6;
  const Int_t kSlicesLow = 6;
  const Int_t kSlicesHigh = 1;
  Double_t trdSlicesLowTime[kLayers*kSlicesLow], trdSlicesRemaining[kLayers*(kSlicesHigh + kSlicesLow)];
  Int_t indices[kLayers*kSlicesLow];
  Int_t cntLowTime=0, cntRemaining = 0;
  for(Int_t idet = 0; idet < 6; idet++)
    for(Int_t islice = 0; islice < kSlicesLow+kSlicesHigh; islice++){
      if(TMath::Abs(track->GetTRDslice(idet, islice)) < fgkVerySmall) continue;;
      if(islice < kSlicesLow){
        AliDebug(3, Form("Part 1, Det[%d], Slice[%d], TRDSlice: %f", idet, islice, track->GetTRDslice(idet, islice)));
        trdSlicesLowTime[cntLowTime++] = track->GetTRDslice(idet, islice) * kWeightSlice[islice];
      } else{
        AliDebug(3, Form("Part 1, Det[%d], Slice[%d], TRDSlice: %f", idet, islice, track->GetTRDslice(idet, islice)));
        trdSlicesRemaining[cntRemaining++] = track->GetTRDslice(idet, islice) * kWeightSlice[islice];
      }
    }
  if(cntLowTime < 4 || cntRemaining < 2) return 0.; // Min. Number of Slices at high time is 2 (matches with 1 layer), for the truncated mean we need at least 4 Slices
  TMath::Sort(cntLowTime, trdSlicesLowTime, indices, kFALSE);
  // Fill the second array with the lower half of the first time bins
  for(Int_t ien = 0; ien < static_cast<Int_t>(static_cast<Float_t>(cntLowTime) * truncation); ien++)
    trdSlicesRemaining[cntRemaining++] = trdSlicesLowTime[indices[ien]];
  Double_t trdSignal = TMath::Mean(cntRemaining, trdSlicesRemaining);
  Double_t mom = track->GetOuterParam() ? track->GetOuterParam()->P() : -1;
  AliDebug(3, Form("PID Meth. 2: p[%f], TRDSignal[%f]", mom, trdSignal));
  return trdSignal;
}
