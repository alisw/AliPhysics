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

#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliOADBContainer.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliHFEOADBThresholdsTRD.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEpidTRD.h"

ClassImp(AliHFEpidTRD)

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD() :
    AliHFEpidBase()
  , fOADBThresholds(NULL)
  , fMinP(0.5)
  , fNTracklets(6)
  , fCutNTracklets(0)
  , fRunNumber(0)
  , fElectronEfficiency(0.90)
  , fTotalChargeInSlice0(kFALSE)
  , fTRD2DPID(kFALSE)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
}

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const char* name) :
    AliHFEpidBase(name)
  , fOADBThresholds(NULL)
  , fMinP(0.5)
  , fNTracklets(6)
  , fCutNTracklets(0)
  , fRunNumber(0)
  , fElectronEfficiency(0.91)
  , fTotalChargeInSlice0(kFALSE)
  , fTRD2DPID(kFALSE)
{
  //
  // default  constructor
  // 
  memset(fThreshParams, 0, sizeof(Double_t) * kThreshParams);
}

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const AliHFEpidTRD &ref):
    AliHFEpidBase("")
  , fOADBThresholds(NULL)
  , fMinP(ref.fMinP)
  , fNTracklets(ref.fNTracklets)
  , fCutNTracklets(ref.fCutNTracklets)
  , fRunNumber(ref.fRunNumber)
  , fElectronEfficiency(ref.fElectronEfficiency)
  , fTotalChargeInSlice0(ref.fTotalChargeInSlice0)
  , fTRD2DPID(ref.fTRD2DPID)
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
  target.fNTracklets = fNTracklets;
  target.fCutNTracklets = fCutNTracklets;
  target.fRunNumber = fRunNumber;
  target.fTotalChargeInSlice0 = fTotalChargeInSlice0;
  target.fTRD2DPID = fTRD2DPID;
  target.fElectronEfficiency = fElectronEfficiency;
  memcpy(target.fThreshParams, fThreshParams, sizeof(Double_t) * kThreshParams);
  AliHFEpidBase::Copy(ref);
}

//___________________________________________________________________
AliHFEpidTRD::~AliHFEpidTRD(){
  //
  // Destructor
  //
}

//______________________________________________________
Bool_t AliHFEpidTRD::InitializePID(Int_t run){
  //
  // InitializePID call different init function depending on TRD PID method
  //
  //

  //if(fTRD2DPID) return Initialize2D(run);
  if(fTRD2DPID) return kTRUE;
  else return Initialize1D(run);



}

//______________________________________________________
Bool_t AliHFEpidTRD::Initialize1D(Int_t run){
  //
  // InitializePID: Load TRD thresholds and create the electron efficiency axis
  // to navigate 
  //
  AliDebug(1, Form("Initializing TRD PID for run %d", run));
  if(InitParamsFromOADB(run)){
    SetBit(kThresholdsInitialized);
    return kTRUE;
  }
  AliDebug(1, Form("Threshold Parameters for %d tracklets and an electron efficiency %f loaded:", fNTracklets, fElectronEfficiency));
  AliDebug(1, Form("Params: [%f|%f|%f|%f]", fThreshParams[0], fThreshParams[1], fThreshParams[2], fThreshParams[3]));
  fRunNumber = run;
  return kFALSE;
}

/*
//______________________________________________________
Bool_t AliHFEpidTRD::Initialize2D(Int_t run){
  //
  // Initialize2DimPID
  //
  //

    return kTRUE;

}
*/

//______________________________________________________
Int_t AliHFEpidTRD::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const {
  //
  // Does PID for TRD alone:
  // PID algorithm selected according to flag
  //
  //

    if(fTRD2DPID) return IsSelected2D(track, pidqa);
    else return IsSelected1D(track, pidqa);


}

//______________________________________________________
Int_t AliHFEpidTRD::IsSelected1D(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const {
  //
  // Does PID for TRD alone:
  // PID thresholds based on 90% Electron Efficiency level approximated by a linear 
  // step function
  //
  if(!TestBit(kThresholdsInitialized)) {
    AliDebug(1,"Threshold Parameters not available");
    return 0;
  }
  AliDebug(2, "Applying TRD PID");
  if(!fkPIDResponse){
    AliDebug(2, "Cannot process track");
    return 0;
  }


/*
  const AliESDtrack *esdt = dynamic_cast<const AliESDtrack *>(track->GetRecTrack());
  printf("checking IdentifiedAsElectronTRD, number of Tracklets: %d\n", esdt->GetTRDntrackletsPID());
  if(fkPIDResponse->IdentifiedAsElectronTRD(dynamic_cast<const AliVTrack *>(track->GetRecTrack()), 0.8)) printf("Track identified as electron\n");
  else printf("Track rejected\n");
*/
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis: AliHFEpidObject::kAODanalysis;
  Double_t p = GetP(track->GetRecTrack(), anatype);
  if(p < fMinP){ 
    AliDebug(2, Form("Track momentum below %f", fMinP));
    return 0;
  }

  if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTRDpid, AliHFEdetPIDqa::kBeforePID); 

  if(fCutNTracklets > 0){
    AliDebug(1, Form("Number of tracklets cut applied: %d\n", fCutNTracklets));
    Int_t ntracklets = track->GetRecTrack() ? track->GetRecTrack()->GetTRDntrackletsPID() : 0;
    if(TestBit(kExactTrackletCut)){
      AliDebug(1, Form("Exact cut applied: %d tracklets found\n", ntracklets));
      if(ntracklets != fCutNTracklets) return 0;
    } else {
      AliDebug(1, Form("Greater Equal cut applied: %d tracklets found\n", ntracklets));
      if(ntracklets < fCutNTracklets) return 0;
    }
  }
  AliDebug(1,"Track selected\n");

  Double_t electronLike = GetElectronLikelihood(track->GetRecTrack(), anatype);
  Double_t threshold;
  if(TestBit(kSelectCutOnTheFly)){ 
    threshold = GetTRDthresholds(p, track->GetRecTrack()->GetTRDntrackletsPID());
  } else {
    threshold = GetTRDthresholds(p);
  }
  AliDebug(2, Form("Threshold: %f\n", threshold));
  if(electronLike > threshold){
    if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTRDpid, AliHFEdetPIDqa::kAfterPID);
    return 11;
  }
  return 211;

}

//______________________________________________________
Int_t AliHFEpidTRD::IsSelected2D(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const {
  //
  // 2D TRD PID
  // 
  // 
  //
  AliDebug(2, "Applying TRD PID");
  if(!fkPIDResponse){
    AliDebug(2, "Cannot process track");
    return 0;
  }

  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis: AliHFEpidObject::kAODanalysis;
  Double_t p = GetP(track->GetRecTrack(), anatype);
  if(p < fMinP){ 
    AliDebug(2, Form("Track momentum below %f", fMinP));
    return 0;
  }

  if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTRDpid, AliHFEdetPIDqa::kBeforePID); 

  if(fCutNTracklets > 0){
    AliDebug(1, Form("Number of tracklets cut applied: %d\n", fCutNTracklets));
    Int_t ntracklets = track->GetRecTrack() ? track->GetRecTrack()->GetTRDntrackletsPID() : 0;
    if(TestBit(kExactTrackletCut)){
      AliDebug(1, Form("Exact cut applied: %d tracklets found\n", ntracklets));
      if(ntracklets != fCutNTracklets) return 0;
    } else {
      AliDebug(1, Form("Greater Equal cut applied: %d tracklets found\n", ntracklets));
      if(ntracklets < fCutNTracklets) return 0;
    }
  }
  AliDebug(1,"Track selected\n");

  Int_t centralitybin = track->IsPbPb() ? track->GetCentrality() + 1 : 0;
  Float_t fCentralityLimitsdefault[12]= {0.,5.,10., 20., 30., 40., 50., 60.,70.,80., 90., 100.};
  Float_t centrality=-1;
  if(centralitybin>=0) centrality=fCentralityLimitsdefault[centralitybin]+1;


  if(fkPIDResponse->IdentifiedAsElectronTRD(track->GetRecTrack(),fElectronEfficiency,centrality,AliTRDPIDResponse::kLQ2D)){
      AliDebug(2, Form("Electron effi: %f\n", fElectronEfficiency));
      return 11;
  } else return 211;



}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetTRDthresholds(Double_t p, UInt_t nTracklets) const { 
  //
  // Return momentum dependent and electron efficiency dependent TRD thresholds
  // Determine threshold based on the number of tracklets on the fly, electron efficiency not modified
  // 
  Double_t threshParams[4];
  AliDebug(1, Form("Select cut for %d tracklets\n", nTracklets));
  // Get threshold paramters for the given number of tracklets from OADB container
  AliHFEOADBThresholdsTRD *thresholds = dynamic_cast<AliHFEOADBThresholdsTRD *>(fOADBThresholds->GetObject(fRunNumber));
  if(!thresholds){
    AliDebug(1, Form("Thresholds for run %d not in the OADB", fRunNumber));
    return 0.;
  }
  if(!thresholds->GetThresholdParameters(nTracklets, fElectronEfficiency, threshParams)){
    AliDebug(1, "loading thresholds failed\n");
    return 0.;
  }
  Double_t threshold = 1. - threshParams[0] - threshParams[1] * p - threshParams[2] * TMath::Exp(-threshParams[3] * p);
  return TMath::Max(TMath::Min(threshold, 0.99), 0.2); // truncate the threshold upperwards to 0.999 and lowerwards to 0.2 and exclude unphysical values
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetTRDthresholds(Double_t p) const { 
  //
  // Return momentum dependent and electron efficiency dependent TRD thresholds
  // 
  Double_t threshold = 1. - fThreshParams[0] - fThreshParams[1] * p - fThreshParams[2] * TMath::Exp(-fThreshParams[3] * p);
  return TMath::Max(TMath::Min(threshold, 0.99), 0.2); // truncate the threshold upperwards to 0.999 and lowerwards to 0.2 and exclude unphysical values
}


//___________________________________________________________________
Bool_t AliHFEpidTRD::InitParamsFromOADB(Int_t run){
  //
  // The name of the function says it all
  //
  AliHFEOADBThresholdsTRD *thresholds = dynamic_cast<AliHFEOADBThresholdsTRD *>(fOADBThresholds->GetObject(run));
  if(!thresholds){
    AliDebug(1, Form("Thresholds for run %d not in the OADB", run));
    return kFALSE;
  }
  return thresholds->GetThresholdParameters(fNTracklets, fElectronEfficiency, fThreshParams);
}

//___________________________________________________________________
void AliHFEpidTRD::RenormalizeElPi(const Double_t * const likein, Double_t * const likeout) const {
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
Double_t AliHFEpidTRD::GetElectronLikelihood(const AliVTrack *track, AliHFEpidObject::AnalysisType_t anaType) const {
  //
  // Get TRD likelihoods for ESD respectively AOD tracks
  //
  Double_t pidProbs[AliPID::kSPECIES]; memset(pidProbs, 0, sizeof(Double_t) * AliPID::kSPECIES);
  if(anaType == AliHFEpidObject::kESDanalysis){
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack) esdtrack->GetTRDpid(pidProbs);
  } else {
    fkPIDResponse->ComputeTRDProbability(track, AliPID::kSPECIES, pidProbs);
  }
  if(!IsRenormalizeElPi()) return pidProbs[AliPID::kElectron];
  Double_t probsNew[AliPID::kSPECIES];
  RenormalizeElPi(pidProbs, probsNew);
  return probsNew[AliPID::kElectron];
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetP(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const {
  //
  // Get the Momentum in the TRD
  //
  Double_t p = 0.;
  if(anaType == AliHFEpidObject::kESDanalysis){
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack) p = esdtrack->GetOuterParam() ? esdtrack->GetOuterParam()->P() : esdtrack->P();
  } else {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    if(aodtrack) p = aodtrack->P();
  }
  return p;
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetChargeLayer(const AliVParticle *track, UInt_t layer, AliHFEpidObject::AnalysisType_t anaType) const {
  //
  // Get the Charge in a single TRD layer
  //
  if(layer >= 6) return 0.;
  Double_t charge = 0.;
  if(anaType == AliHFEpidObject::kESDanalysis){
    const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
    if(esdtrack){
      // Distinction between old and new reconstruction: in the new reconstruction, the total charge is stored in slice 0, slices 1 to 8 are used for the slices for 
      // the neural network. 
      if(fTotalChargeInSlice0)
        charge = esdtrack->GetTRDslice(static_cast<UInt_t>(layer), 0);
      else
       for(Int_t islice = 0; islice < esdtrack->GetNumberOfTRDslices(); islice++) charge += esdtrack->GetTRDslice(static_cast<UInt_t>(layer), islice);
    }
  } else {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    AliAODPid *aoddetpid = aodtrack ? aodtrack->GetDetPid() : NULL;
    if(aoddetpid){
      if(fTotalChargeInSlice0)
        charge = aoddetpid->GetTRDslices()[layer * aoddetpid->GetTRDnSlices()];
      else
       for(Int_t islice = 0; islice < aoddetpid->GetTRDnSlices(); islice++) charge += aoddetpid->GetTRDslices()[layer * aoddetpid->GetTRDnSlices() + islice];
    }
  }
  return charge;
}

//___________________________________________________________________
void AliHFEpidTRD::GetTRDmomenta(const AliVTrack *track, Double_t *mom) const {
  //
  // Fill Array with momentum information at the TRD tracklet
  //
  for(Int_t itl = 0; itl < 6; itl++) 
    mom[itl] = track->GetTRDmomentum(itl);
}

//___________________________________________________________________
Double_t AliHFEpidTRD::GetTRDSignalV1(const AliESDtrack *track, Float_t truncation) const {
  //
  // Calculation of the TRD Signal via truncated mean
  // Method 1: Take all Slices available
  // cut out 0s
  // Order them in increasing order
  // Cut out the upper third
  // Calculate mean over the last 2/3 slices
  //
  const Int_t kNSlices = 48;
  const Int_t kLastSlice = 6; // Slice 7 is taken out from the truncated mean calculation
  const Double_t kVerySmall = 1e-12;
  // Weight the slice to equalize the MPV of the dQ/dl-distribution per slice to the one in the first slice
  // Pions are used as reference for the equalization
  const Double_t kWeightSlice[8] = {1., 2.122, 1.8, 1.635, 1.595, 1.614, 1.16, 7.0};
  const Double_t kWeightSliceNo0[8] = {1., 1., 1.271, 1.451, 1.531, 1.543, 1.553, 2.163};  // Weighting factors in case slice 0 stores the total charge
  const Double_t *kWeightFactor = fTotalChargeInSlice0 ? kWeightSliceNo0 : kWeightSlice;
  AliDebug(3, Form("Number of Tracklets: %d\n", track->GetTRDntrackletsPID()));
  Double_t trdSlices[kNSlices], tmp[kNSlices];
  Int_t indices[48];
  Int_t icnt = 0;
  for(Int_t idet = 0; idet < 6; idet++)
    for(Int_t islice = fTotalChargeInSlice0 ? 1 : 0 ; islice <= kLastSlice; islice++){
      AliDebug(2, Form("Chamber[%d], Slice[%d]: TRDSlice = %f", idet, islice, track->GetTRDslice(idet, islice)));
      if(TMath::Abs(track->GetTRDslice(idet, islice)) < kVerySmall) continue;;
      trdSlices[icnt++] = track->GetTRDslice(idet, islice) * kWeightFactor[islice];
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
Double_t AliHFEpidTRD::GetTRDSignalV2(const AliESDtrack *track, Float_t truncation) const {
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
  const Double_t kVerySmall = 1e-12;
  Double_t trdSlicesLowTime[kLayers*kSlicesLow], trdSlicesRemaining[kLayers*(kSlicesHigh + kSlicesLow)];
  Int_t indices[kLayers*kSlicesLow];
  Int_t cntLowTime=0, cntRemaining = 0;
  for(Int_t idet = 0; idet < 6; idet++)
    for(Int_t islice = fTotalChargeInSlice0 ? 1 : 0; islice < kSlicesLow+kSlicesHigh; islice++){
      if(TMath::Abs(track->GetTRDslice(idet, islice)) < kVerySmall) continue;;
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
