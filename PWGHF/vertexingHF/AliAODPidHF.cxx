/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 * *************************************************************************/

/* $Id$ */

///***********************************************************
/// \class Class AliAODPidHF
/// \brief class for PID with AliAODRecoDecayHF
/// \author Authors: D. Caffarri caffarri@pd.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de P. Antonioli pietro.antonioli@bo.infn.it, J. van der Maarel j.vandermaarel@cern.ch
///***********************************************************
#include <TCanvas.h>
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>

#include "AliAODPidHF.h"
#include "AliAODPid.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"
#include "AliESDtrack.h"

/// \cond CLASSIMP
ClassImp(AliAODPidHF);
/// \endcond

//------------------------------
AliAODPidHF::AliAODPidHF():
TObject(),
fnNSigma(5),
fnSigma(0),
fTOFSigma(160.),
fCutTOFmismatch(0.01),
fMinNClustersTPCPID(0),
fnPriors(5),
fPriors(0),
fnPLimit(2),
fPLimit(0),
fAsym(kFALSE),
fTPC(kFALSE),
fTOF(kFALSE),
fITS(kFALSE),
fTRD(kFALSE),
fMatch(0),
fForceTOFforKaons(kFALSE),
fCompat(kFALSE),
fPCompatTOF(1.5),
fUseAsymTOF(kFALSE),
fLownSigmaTOF(-3.),
fUpnSigmaTOF(3.),
fLownSigmaCompatTOF(-3.),
fUpnSigmaCompatTOF(3.),
fnNSigmaCompat(2),
fnSigmaCompat(0),
fMC(kFALSE),
fOnePad(kFALSE),
fMCLowEn2011(kFALSE),
fppLowEn2011(kFALSE),
fPbPb(kFALSE),
fTOFdecide(kFALSE),
fOldPid(kFALSE),
fPtThresholdTPC(999999.),
fMaxTrackMomForCombinedPID(999999.),
fPidResponse(0),
fPidCombined(new AliPIDCombined()),
fTPCResponse(new AliTPCPIDResponse()),
fPriorsH(),
fCombDetectors(kTPCTOF),
fUseCombined(kFALSE),
fDefaultPriors(kTRUE),
fApplyNsigmaTPCDataCorr(kFALSE),
fMeanNsigmaTPCPionData{},
fMeanNsigmaTPCKaonData{},
fMeanNsigmaTPCProtonData{},
fSigmaNsigmaTPCPionData{},
fSigmaNsigmaTPCKaonData{},
fSigmaNsigmaTPCProtonData{},
fPlimitsNsigmaTPCDataCorr{},
fNPbinsNsigmaTPCDataCorr(0),
fEtalimitsNsigmaTPCDataCorr{},
fNEtabinsNsigmaTPCDataCorr(0)
{
  ///
  /// Default constructor
  ///
  fPLimit=new Double_t[fnPLimit];
  fnSigma=new Double_t[fnNSigma];
  fPriors=new Double_t[fnPriors];
  fnSigmaCompat=new Double_t[fnNSigmaCompat];
  
  for(Int_t i=0;i<fnNSigma;i++){
    fnSigma[i]=0.;
  }
  for(Int_t i=0;i<fnPriors;i++){
    fPriors[i]=0.;
  }
  for(Int_t i=0;i<fnPLimit;i++){
    fPLimit[i]=0.;
  }
  for(Int_t i=0;i<fnNSigmaCompat;i++){
    fnSigmaCompat[i]=3.;
  }
  for(Int_t i=0; i<3; i++){ // pi, K, proton
    fMaxnSigmaCombined[i]=3.;
    fMinnSigmaTPC[i]=-3;
    fMaxnSigmaTPC[i]=3;
    fMinnSigmaTOF[i]=-3;
    fMaxnSigmaTOF[i]=3;
  }
  for (Int_t s=0;s<AliPID::kSPECIES;s++) {
    for (Int_t d=0;d<4;d++) {
      fIdBandMin[s][d] = NULL;
      fIdBandMax[s][d] = NULL;
      fCompBandMin[s][d] = NULL;
      fCompBandMax[s][d] = NULL;
    }
  }

  for(int iP=0; iP<=kMaxPBins; iP++) {
    fPlimitsNsigmaTPCDataCorr[iP] = 0.;
  }
  for(int iEta=0; iEta<=kMaxEtaBins; iEta++) {
    fEtalimitsNsigmaTPCDataCorr[iEta]=0;
  }
}
//----------------------
AliAODPidHF::~AliAODPidHF()
{
  /// destructor
  if(fPLimit) delete [] fPLimit;
  if(fnSigma) delete [] fnSigma;
  if(fPriors) delete [] fPriors;
  if(fnSigmaCompat) delete [] fnSigmaCompat;
  delete fPidCombined;
  
  delete fTPCResponse;
  for (Int_t ispecies=0;ispecies<AliPID::kSPECIES;++ispecies) {
    delete fPriorsH[ispecies];
  }

  for (Int_t s=0;s<AliPID::kSPECIES;s++) {
    for (Int_t d=0;d<4;d++) {
      delete fIdBandMin[s][d];
      delete fIdBandMax[s][d];
      delete fCompBandMin[s][d];
      delete fCompBandMax[s][d];
    }
  }
}
//------------------------
AliAODPidHF::AliAODPidHF(const AliAODPidHF& pid) :
TObject(),
fnNSigma(pid.fnNSigma),
fnSigma(0),
fTOFSigma(pid.fTOFSigma),
fCutTOFmismatch(pid.fCutTOFmismatch),
fMinNClustersTPCPID(pid.fMinNClustersTPCPID),
fnPriors(pid.fnPriors),
fPriors(0),
fnPLimit(pid.fnPLimit),
fPLimit(0),
fAsym(pid.fAsym),
fTPC(pid.fTPC),
fTOF(pid.fTOF),
fITS(pid.fITS),
fTRD(pid.fTRD),
fMatch(pid.fMatch),
fForceTOFforKaons(pid.fForceTOFforKaons),
fCompat(pid.fCompat),
fPCompatTOF(pid.fPCompatTOF),
fUseAsymTOF(pid.fUseAsymTOF),
fLownSigmaTOF(pid.fLownSigmaTOF),
fUpnSigmaTOF(pid.fUpnSigmaTOF),
fLownSigmaCompatTOF(pid.fLownSigmaCompatTOF),
fUpnSigmaCompatTOF(pid.fUpnSigmaCompatTOF),
fnNSigmaCompat(pid.fnNSigmaCompat),
fnSigmaCompat(0x0),
fMC(pid.fMC),
fOnePad(pid.fOnePad),
fMCLowEn2011(pid.fMCLowEn2011),
fppLowEn2011(pid.fppLowEn2011),
fPbPb(pid.fPbPb),
fTOFdecide(pid.fTOFdecide),
fOldPid(pid.fOldPid),
fPtThresholdTPC(pid.fPtThresholdTPC),
fMaxTrackMomForCombinedPID(pid.fMaxTrackMomForCombinedPID),
fPidResponse(0x0),
fPidCombined(0x0),
fTPCResponse(0x0),
fCombDetectors(pid.fCombDetectors),
fUseCombined(pid.fUseCombined),
fDefaultPriors(pid.fDefaultPriors),
fApplyNsigmaTPCDataCorr(pid.fApplyNsigmaTPCDataCorr),
fNPbinsNsigmaTPCDataCorr(pid.fNPbinsNsigmaTPCDataCorr),
fNEtabinsNsigmaTPCDataCorr(pid.fNEtabinsNsigmaTPCDataCorr)
{
  
  fnSigmaCompat=new Double_t[fnNSigmaCompat];
  for(Int_t i=0;i<fnNSigmaCompat;i++){
    fnSigmaCompat[i]=pid.fnSigmaCompat[i];
  }
  fnSigma = new Double_t[fnNSigma];
  for(Int_t i=0;i<fnNSigma;i++){
    fnSigma[i]=pid.fnSigma[i];
  }
  fPriors = new Double_t[fnPriors];
  for(Int_t i=0;i<fnPriors;i++){
    fPriors[i]=pid.fPriors[i];
  }
  fPLimit = new Double_t[fnPLimit];
  for(Int_t i=0;i<fnPLimit;i++){
    fPLimit[i]=pid.fPLimit[i];
  }
  for(Int_t i=0;i<AliPID::kSPECIES;i++){
    fPriorsH[i] = pid.fPriorsH[i] ? new TH1F(*pid.fPriorsH[i]) : NULL;
  }
  for(Int_t i=0; i<3; i++){ // pi, K, proton
    fMaxnSigmaCombined[i]=pid.fMaxnSigmaCombined[i];
    fMinnSigmaTPC[i]=pid.fMinnSigmaTPC[i];
    fMaxnSigmaTPC[i]=pid.fMaxnSigmaTPC[i];
    fMinnSigmaTOF[i]=pid.fMinnSigmaTOF[i];
    fMaxnSigmaTOF[i]=pid.fMaxnSigmaTOF[i];
  }
  
  //  if(pid.fTPCResponse) fTPCResponse = new AliTPCPIDResponse(*(pid.fTPCResponse));
  fTPCResponse = new AliTPCPIDResponse();
  SetBetheBloch();
  fPidCombined = new AliPIDCombined();
  //fPidResponse = new AliPIDResponse(*(pid.fPidResponse));
  //fPidCombined = new AliPIDCombined(*(pid.fPidCombined));

  //Copy bands
  for (Int_t s=0;s<AliPID::kSPECIES;s++) {
    for (Int_t d=0;d<4;d++) {
      fIdBandMin[s][d]   = pid.fIdBandMin[s][d]   ? new TF1(*pid.fIdBandMin[s][d])   : NULL;
      fIdBandMax[s][d]   = pid.fIdBandMax[s][d]   ? new TF1(*pid.fIdBandMax[s][d])   : NULL;
      fCompBandMin[s][d] = pid.fCompBandMin[s][d] ? new TF1(*pid.fCompBandMin[s][d]) : NULL;
      fCompBandMax[s][d] = pid.fCompBandMax[s][d] ? new TF1(*pid.fCompBandMax[s][d]) : NULL;
    }
  }
  
  for(Int_t iBin=0; iBin<fNPbinsNsigmaTPCDataCorr; iBin++) {
    for(Int_t iEta=0; iEta<fNEtabinsNsigmaTPCDataCorr; iEta++) {
      fMeanNsigmaTPCPionData[iEta][iBin] = pid.fMeanNsigmaTPCPionData[iEta][iBin];
      fMeanNsigmaTPCKaonData[iEta][iBin] = pid.fMeanNsigmaTPCKaonData[iEta][iBin];
      fMeanNsigmaTPCProtonData[iEta][iBin] = pid.fMeanNsigmaTPCProtonData[iEta][iBin];
      fSigmaNsigmaTPCPionData[iEta][iBin] = pid.fSigmaNsigmaTPCPionData[iEta][iBin];
      fSigmaNsigmaTPCKaonData[iEta][iBin] = pid.fSigmaNsigmaTPCKaonData[iEta][iBin];
      fSigmaNsigmaTPCProtonData[iEta][iBin] = pid.fSigmaNsigmaTPCProtonData[iEta][iBin];
    }
    fPlimitsNsigmaTPCDataCorr[iBin] = pid.fPlimitsNsigmaTPCDataCorr[iBin];
  }
  fPlimitsNsigmaTPCDataCorr[fNPbinsNsigmaTPCDataCorr] = pid.fPlimitsNsigmaTPCDataCorr[fNPbinsNsigmaTPCDataCorr];
  for(int iEta=0; iEta<=fNEtabinsNsigmaTPCDataCorr; iEta++) {
    fEtalimitsNsigmaTPCDataCorr[iEta] = pid.fEtalimitsNsigmaTPCDataCorr[iEta];
  }
}
//----------------------
Int_t AliAODPidHF::RawSignalPID(AliAODTrack *track, TString detector) const{
  /// raw PID for single detectors, returns the particle type with smaller sigma
  Int_t specie=-1;
  if(detector.Contains("ITS")) return ApplyPidITSRaw(track,specie);
  if(detector.Contains("TPC")) return ApplyPidTPCRaw(track,specie);
  if(detector.Contains("TOF")) return ApplyPidTOFRaw(track,specie);
  
  return specie;
  
}
//---------------------------
Bool_t AliAODPidHF::IsKaonRaw(AliAODTrack *track, TString detector) const{
  /// checks if the track can be a kaon, raw PID applied for single detectors
  Int_t specie=0;
  
  if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,3);
  if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,3);
  if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,3);
  
  if(specie==3) return kTRUE;
  return kFALSE;
}
//---------------------------
Bool_t AliAODPidHF::IsPionRaw (AliAODTrack *track, TString detector) const{
  /// checks if the track can be a pion, raw PID applied for single detectors
  
  Int_t specie=0;
  
  if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,2);
  if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,2);
  if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,2);
  
  if(specie==2) return kTRUE;
  return kFALSE;
}
//---------------------------
Bool_t AliAODPidHF::IsProtonRaw (AliAODTrack *track, TString detector) const{
  /// checks if the track can be a proton raw PID applied for single detectors
  
  Int_t specie=0;
  if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,4);
  if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,4);
  if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,4);
  
  if(specie==4) return kTRUE;
  
  return kFALSE;
}
//--------------------------
Bool_t AliAODPidHF::IsElectronRaw(AliAODTrack *track, TString detector) const{
  /// checks if the track can be an electron raw PID applied for single detectors
  
  Int_t specie=-1;
  if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,0);
  if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,0);
  if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,0);
  
  if(specie==0) return kTRUE;
  
  return kFALSE;
}
//--------------------------
Int_t AliAODPidHF::ApplyPidTPCRaw(AliAODTrack *track,Int_t specie) const{
  /// n-sigma cut, TPC PID
  
  Double_t nsigma=-999.;
  Int_t pid=-1;
  
  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
    Double_t nsigmaMin=999.;
    for(Int_t ipart=0;ipart<5;ipart++){
      if(GetnSigmaTPC(track,ipart,nsigma)==1){
        nsigma=TMath::Abs(nsigma);
        if((nsigma<nsigmaMin) && (nsigma<fnSigma[0])) {
          pid=ipart;
          nsigmaMin=nsigma;
        }
      }
    }
  }else{ // asks only for one particle specie
    if(GetnSigmaTPC(track,specie,nsigma)==1){
      nsigma=TMath::Abs(nsigma);
      if (nsigma>fnSigma[0]) pid=-1;
      else pid=specie;
    }
  }
  
  return pid;
}
//----------------------------
Int_t AliAODPidHF::ApplyPidITSRaw(AliAODTrack *track,Int_t specie) const{
  /// truncated mean, ITS PID
  
  Double_t nsigma=-999.;
  Int_t pid=-1;
  
  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
    Double_t nsigmaMin=999.;
    for(Int_t ipart=0;ipart<5;ipart++){
      if(GetnSigmaITS(track,ipart,nsigma)==1){
        nsigma=TMath::Abs(nsigma);
        if((nsigma<nsigmaMin) && (nsigma<fnSigma[4])) {
          pid=ipart;
          nsigmaMin=nsigma;
        }
      }
    }
  }else{ // asks only for one particle specie
    if(GetnSigmaITS(track,specie,nsigma)==1){
      nsigma=TMath::Abs(nsigma);
      if (nsigma>fnSigma[4]) pid=-1;
      else pid=specie;
    }
  }
  
  return pid;
}
//----------------------------
Int_t AliAODPidHF::ApplyPidTOFRaw(AliAODTrack *track,Int_t specie) const{
  /// n-sigma cut, TOF PID
  
  Double_t nsigma=-999.;
  Int_t pid=-1;
  
  if(specie<0){
    Double_t nsigmaMin=999.;
    for(Int_t ipart=0;ipart<5;ipart++){
      if(GetnSigmaTOF(track,ipart,nsigma)==1){
        nsigma=TMath::Abs(nsigma);
        if((nsigma<nsigmaMin)&& (nsigma<fnSigma[3])){
          pid=ipart;
          nsigmaMin=nsigma;
        }
      }
    }
  }else{ // asks only for one particle specie
    Double_t nSigmaMin,nSigmaMax;
    if(fUseAsymTOF){
      nSigmaMin=fLownSigmaTOF;
      nSigmaMax=fUpnSigmaTOF;
    }else{
      nSigmaMin=-fnSigma[3];
      nSigmaMax=fnSigma[3];
    }
    if(GetnSigmaTOF(track,specie,nsigma)==1){
      if(nsigma<nSigmaMin || nsigma>nSigmaMax) pid=-1;
      else pid=specie;
    }
  }
  return pid;
}
//----------------------------
Int_t AliAODPidHF::ApplyTOFCompatibilityBand(AliAODTrack *track,Int_t specie) const{
  /// n-sigma cut, TOF PID
  
  if(specie<0) return -1;
  Double_t nsigma=-999.;
  Int_t pid=-1;
  
  Double_t nSigmaMin,nSigmaMax;
  if(fUseAsymTOF){
    nSigmaMin=fLownSigmaCompatTOF;
    nSigmaMax=fUpnSigmaCompatTOF;
  }else{
    nSigmaMin=-fnSigmaCompat[1];
    nSigmaMax=fnSigmaCompat[1];
  }
  if(GetnSigmaTOF(track,specie,nsigma)==1){
    if(nsigma<nSigmaMin || nsigma>nSigmaMax) pid=-1;
    else pid=specie;
  }
  return pid;
}
//------------------------------
void AliAODPidHF::CombinedProbability(AliAODTrack *track,Bool_t *type) const{
  /// combined PID stored inside the AOD track
  
  const Double_t *pid=track->PID();
  Float_t max=0.;
  Int_t k=-1;
  for (Int_t i=0; i<10; i++) {
    if (pid[i]>max) {k=i; max=pid[i];}
  }
  
  if(k==2) type[0]=kTRUE;
  if(k==3) type[1]=kTRUE;
  if(k==4) type[2]=kTRUE;
  
  return;
}
//--------------------------------
Bool_t AliAODPidHF::CheckITSPIDStatus(AliAODTrack *track) const{
  /// Check if the track is good for ITS PID
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kITS,track);
  if (status != AliPIDResponse::kDetPidOk) return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckTPCPIDStatus(AliAODTrack *track) const{
  /// Check if the track is good for TPC PID
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTPC,track);
  if (status != AliPIDResponse::kDetPidOk) return kFALSE;
  UInt_t nclsTPCPID = track->GetTPCsignalN();
  if(nclsTPCPID<fMinNClustersTPCPID) return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckTOFPIDStatus(AliAODTrack *track) const{
  /// Check if the track is good for TOF PID
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
  if (status != AliPIDResponse::kDetPidOk) return kFALSE;
  Float_t probMis = fPidResponse->GetTOFMismatchProbability(track);
  if (probMis > fCutTOFmismatch) return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckTRDPIDStatus(AliAODTrack *track) const{
  /// Check if the track is good for TRD PID
  AliPIDResponse::EDetPidStatus status = fPidResponse->CheckPIDStatus(AliPIDResponse::kTRD,track);
  if (status != AliPIDResponse::kDetPidOk) return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckStatus(AliAODTrack *track,TString detectors) const{
  /// Quality cuts on the tracks, detector by detector
  if(detectors.Contains("ITS")) return CheckITSPIDStatus(track);
  else if(detectors.Contains("TPC")) return CheckTPCPIDStatus(track);
  else if(detectors.Contains("TOF")) return CheckTOFPIDStatus(track);
  else if(detectors.Contains("TRD")) return CheckTRDPIDStatus(track);
  else{
    AliError("Wrong detector name");
    return kFALSE;
  }
}
//--------------------------------------------
Bool_t AliAODPidHF::TPCRawAsym(AliAODTrack* track,Int_t specie) const{
  /// TPC nsigma cut PID, different sigmas in different p bins
  
  AliAODPid *pidObj = track->GetDetPid();
  Double_t mom = pidObj->GetTPCmomentum();
  if(mom>fPtThresholdTPC) return kTRUE;
  
  Double_t nsigma;
  if(GetnSigmaTPC(track,specie,nsigma)!=1) return kFALSE;
  nsigma=TMath::Abs(nsigma);
  
  
  if(mom<fPLimit[0] && nsigma<fnSigma[0]) return kTRUE;
  if(mom<fPLimit[1] && mom>fPLimit[0] && nsigma<fnSigma[1]) return kTRUE;
  if(mom>fPLimit[1] && nsigma<fnSigma[2]) return kTRUE;
  
  return kFALSE;
}
//------------------
Int_t AliAODPidHF::MatchTPCTOF(AliAODTrack *track, Int_t specie){
  /// combination of the PID info coming from TPC and TOF
  
  Double_t ptrack=track->P();
  if(ptrack>fMaxTrackMomForCombinedPID) return 1;
  
  Bool_t okTPC=CheckTPCPIDStatus(track);
  if(ptrack>fPtThresholdTPC) okTPC=kFALSE;
  Bool_t okTOF=CheckTOFPIDStatus(track);
  
  if(fMatch==1){
    //TOF || TPC (a la' Andrea R.)
    // convention:
    // for the single detectors: -1 = kFALSE, 1 = kTRUE, 0 = compatible
    // the method returns the sum of the response of the 2 detectors
    
    if(fTPC && fTOF) {
      if(!okTPC && !okTOF) return 0;
    }
    
    Int_t tTPCinfo=0;
    if(fTPC && okTPC){
      tTPCinfo=-1;
      if(fAsym) {
        if(TPCRawAsym(track,specie)) tTPCinfo=1;
      }else{
        if(ApplyPidTPCRaw(track,specie)==specie) tTPCinfo=1;
      }
      if(fCompat && tTPCinfo<0){
        Double_t sig0tmp=fnSigma[0];
        SetSigma(0,fnSigmaCompat[0]);
        if(ApplyPidTPCRaw(track,specie)==specie) tTPCinfo=0;
        SetSigma(0,sig0tmp);
      }
    }
    
    Int_t tTOFinfo=0;
    if(fTOF){
      if(!okTOF && fTPC) return tTPCinfo;
      tTOFinfo=-1;
      if(ApplyPidTOFRaw(track,specie)==specie) tTOFinfo=1;
      if(fCompat && tTOFinfo>0){
        if(ptrack>fPCompatTOF) {
          if(ApplyTOFCompatibilityBand(track,specie)==specie) tTOFinfo=0;
        }
      }
    }
    
    
    if(tTPCinfo+tTOFinfo==0 && fTOFdecide){
      if(!okTOF) return tTPCinfo;
      return tTOFinfo;
    }
    
    if(tTPCinfo+tTOFinfo==0 && fITS){
      if(!CheckITSPIDStatus(track)) return tTPCinfo+tTOFinfo;
      Int_t tITSinfo = -1;
      if(ApplyPidITSRaw(track,specie)==specie) tITSinfo=1;
      return tITSinfo;
    }
    return tTPCinfo+tTOFinfo;
  }
  
  if(fMatch==2){
    //TPC & TOF (a la' Yifei)
    // convention: -1 = kFALSE, 1 = kTRUE, 0 = not identified
    Int_t tTPCinfo=0;
    
    if(fTPC && okTPC) {
      tTPCinfo=1;
      if(fAsym){
        if(!TPCRawAsym(track,specie)) tTPCinfo=-1;
      }else{
        if(ApplyPidTPCRaw(track,specie)!=specie) tTPCinfo=-1;
      }
    }
    
    Int_t tTOFinfo=1;
    if(fTOF){
      if(fTPC && !okTOF) return tTPCinfo;
      if(ApplyPidTPCRaw(track,specie)!=specie) tTOFinfo=-1;
    }
    
    if(tTOFinfo==1 && tTPCinfo==1) return 1;
    
    if(tTPCinfo+tTOFinfo==0 && fITS){
      if(!CheckITSPIDStatus(track)) return tTPCinfo+tTOFinfo;
      Int_t tITSinfo = -1;
      if(ApplyPidITSRaw(track,specie)==specie) tITSinfo=1;
      return tITSinfo;
    }
    return -1;
  }
  
  if(fMatch==3){
    //TPC for p<fPLimit[0], TOF for p>=fPLimit[0] (a la' Andrea A.)
    // convention (temporary): -1 = kFALSE, 1 = kTRUE, 0 = not identified
    if(fTPC && fTOF) if(!okTPC && !okTOF) return 0;
    
    
    Int_t tTPCinfo=-1;
    if(ptrack>=fPLimit[0] && ptrack<fPLimit[1] && fTPC) {
      if(!okTPC) return 0;
      if(fAsym) {
        if(TPCRawAsym(track,specie)) tTPCinfo=1;
      }else{
        if(ApplyPidTPCRaw(track,specie)==specie) tTPCinfo=1;
      }
      return tTPCinfo;
    }
    
    Int_t tTOFinfo=-1;
    if(ptrack>=fPLimit[1] && fTOF){
      if(!okTOF) return 0;
      if(ApplyPidTOFRaw(track,specie)==specie) tTOFinfo=1;
      return tTOFinfo;
    }
    
    Int_t tITSinfo=-1;
    if(ptrack<fPLimit[0] && fITS){
      if(!CheckITSPIDStatus(track)) return 0;
      if(ApplyPidITSRaw(track,specie)==specie) tITSinfo=1;
      return tITSinfo;
    }
  }
  
  if(fMatch==4 || fMatch==5){
    
    // fMatch == 4 ---> "circular cut" in nSigmaTPC, nSimgaTOF plane
    //             ---> nsigmaTPC^2+nsigmaTOF^2 < cut^2
    // fMatch == 5 ---> "rectangular cut" in nSigmaTPC, nsigmaTOF plane
    //             ---> ns1<nSigmaTPC<NS1  && ns2<nSigmaTOF<NS2
    
    Double_t nSigmaTPC=0.;
    if(okTPC) {
      nSigmaTPC = fPidResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)specie);
      if(fApplyNsigmaTPCDataCorr && nSigmaTPC>-990.) { 
        Float_t mean=0., sigma=1.; 
        GetNsigmaTPCMeanSigmaData(mean, sigma, (AliPID::EParticleType)specie, track->GetTPCmomentum(),track->Eta());
        nSigmaTPC = (nSigmaTPC-mean)/sigma;
      }
      if(nSigmaTPC<-990.) nSigmaTPC=0.;
    }
    Double_t nSigmaTOF=0.;
    if(okTOF) {
      nSigmaTOF=fPidResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)specie);
    }
    Int_t iPart=specie-2; //species is 2 for pions,3 for kaons and 4 for protons
    if(iPart<0 || iPart>2) return -1;
    if(fMatch==4){
      Double_t nSigma2=nSigmaTPC*nSigmaTPC+nSigmaTOF*nSigmaTOF;
      if(nSigma2<fMaxnSigmaCombined[iPart]*fMaxnSigmaCombined[iPart]) return 1;
      else return -1;
    }
    else if(fMatch==5){
      if(fForceTOFforKaons && iPart==1 && !okTOF) return -1;
      if((nSigmaTPC>fMinnSigmaTPC[iPart] && nSigmaTPC<fMaxnSigmaTPC[iPart]) &&
         (nSigmaTOF>fMinnSigmaTOF[iPart] && nSigmaTOF<fMaxnSigmaTOF[iPart])) return 1;
      else return -1;
    }
  }
  
  //Asymmetric cuts using user defined bands
  if (fMatch == 10) {
    if (fTPC && fTOF && !okTPC && !okTOF) {
      return 0;
    }
    
    Int_t tTPCinfo = 0;
    if (fTPC && okTPC) {
      tTPCinfo = CheckBands((AliPID::EParticleType) specie, AliPIDResponse::kTPC, track);
    }
    
    Int_t tTOFinfo = 0;
    if (fTOF) {
      if (!okTOF && fTPC) {
        return tTPCinfo;
      }
      tTOFinfo = CheckBands((AliPID::EParticleType) specie, AliPIDResponse::kTOF, track);
    }
    
    
    if (tTPCinfo+tTOFinfo == 0 && fTOFdecide) {
      if (!okTOF) {
        return tTPCinfo;
      }
      return tTOFinfo;
    }
    
    if (tTPCinfo+tTOFinfo == 0 && fITS) {
      if (!CheckITSPIDStatus(track)) {
        return tTPCinfo+tTOFinfo;
      }
      Int_t tITSinfo = CheckBands((AliPID::EParticleType) specie, AliPIDResponse::kITS, track);
      return tITSinfo;
    }
    return tTPCinfo+tTOFinfo;
  }

  return -1;
  
}

//--------------------------------------------------------------
Int_t AliAODPidHF::MatchTPCTOFMin(AliAODTrack *track, Int_t specie){
  /// combination of the PID info coming from TPC and TOF

  Bool_t okTPC=CheckTPCPIDStatus(track);
  Bool_t okTOF=CheckTOFPIDStatus(track);

  if(fTPC && fTOF){
    if(!okTPC && !okTOF) return 0;
  }

  Int_t pid=-1;
  Double_t nsigmaTPC[5]={999.,999.,999.,999.,999.};
  Double_t nsigmaTOF[5]={999.,999.,999.,999.,999.};
  Double_t nsigmaMin=999.;
  Double_t nsigma[5]={999.,999.,999.,999.,999.};

  if(okTPC) {
    for(Int_t ipart=0;ipart<5;ipart++){
      if(GetnSigmaTPC(track,ipart,nsigmaTPC[ipart])<1) nsigmaTPC[ipart]=0.;
    }
  }else{
    for(Int_t ipart=0;ipart<5;ipart++){nsigmaTPC[ipart]=0.;}
  }

  if(okTOF){
    for(Int_t ipart=0;ipart<5;ipart++){
      if(GetnSigmaTOF(track,ipart,nsigmaTOF[ipart])<1) nsigmaTOF[ipart]=0.;
    }
  }else{
    for(Int_t ipart=0;ipart<5;ipart++){nsigmaTOF[ipart]=0.;}
  }

  for(Int_t ipart=0;ipart<5;ipart++){
    nsigma[ipart]=TMath::Sqrt(nsigmaTPC[ipart]*nsigmaTPC[ipart]+nsigmaTOF[ipart]*nsigmaTOF[ipart]);
    if(nsigma[ipart]<nsigmaMin) {nsigmaMin=nsigma[ipart];pid=ipart;}
  }

  if(pid==specie) return 1;

  return 0;
}


//----------------------------------
Int_t AliAODPidHF::MakeRawPid(AliAODTrack *track, Int_t specie){
  /// general method to compute PID
  if(fMatch>0){
    return MatchTPCTOF(track,specie);
  }else{
    if(fTPC && !fTOF && !fITS) {
      Int_t tTPCres=0;
      if(!fAsym){
        tTPCres=ApplyPidTPCRaw(track,specie);
        if(tTPCres==specie) return 1;
        else return tTPCres;
      }else{
        if(TPCRawAsym(track,specie)) tTPCres=1;
        else tTPCres=-1;
      }
      return tTPCres;
    }else if(fTOF && !fTPC && !fITS) {
      Int_t tTOFres=ApplyPidTOFRaw(track,specie);
      if(tTOFres==specie) return 1;
      else return tTOFres;
    }else if(fITS && !fTPC && !fTOF) {
      Int_t tITSres=ApplyPidITSRaw(track,specie);
      if(tITSres==specie) return 1;
      else return tITSres;
    }else{
      AliError("You should enable just one detector if you don't want to match");
      return 0;
    }
  }
}
//--------------------------------------------
void AliAODPidHF::GetTPCBetheBlochParams(Double_t alephParameters[5]) const {
  /// TPC bethe bloch parameters
  if(fMC) {  // MC
    
    if(fPbPb) { // PbPb MC
      
      alephParameters[0] = 1.44405/50.;
      alephParameters[1] = 2.35409e+01;
      alephParameters[2] = TMath::Exp(-2.90330e+01);
      alephParameters[3] = 2.10681e+00;
      alephParameters[4] = 4.62254e+00;
      
    } else {  // pp MC
      if(fMCLowEn2011){
        alephParameters[0]=0.0207667;
        alephParameters[1]=29.9936;
        alephParameters[2]=3.87866e-11;
        alephParameters[3]=2.17291;
        alephParameters[4]=7.1623;
      }else if(fOnePad){
        alephParameters[0]=0.029021;
        alephParameters[1]=25.4181;
        alephParameters[2]=4.66596e-08;
        alephParameters[3]=1.90008;
        alephParameters[4]=4.63783;
      }else{
        alephParameters[0] = 2.15898/50.;
        alephParameters[1] = 1.75295e+01;
        alephParameters[2] = 3.40030e-09;
        alephParameters[3] = 1.96178e+00;
        alephParameters[4] = 3.91720e+00;
      }
    }
    
  } else { // Real Data
    
    if(fOnePad) { // pp 1-pad (since LHC10d)
      
      alephParameters[0] =1.34490e+00/50.;
      alephParameters[1] = 2.69455e+01;
      alephParameters[2] = TMath::Exp(-2.97552e+01);
      alephParameters[3] = 2.35339e+00;
      alephParameters[4] = 5.98079e+00;
      
    } else if(fPbPb) { // PbPb
      
      // alephParameters[0] = 1.25202/50.;
      // alephParameters[1] = 2.74992e+01;
      // alephParameters[2] = TMath::Exp(-3.31517e+01);
      // alephParameters[3] = 2.46246;
      // alephParameters[4] = 6.78938;
      
      alephParameters[0] = 5.10207e+00/50.;
      alephParameters[1] = 7.94982e+00;
      alephParameters[2] = TMath::Exp(-9.07942e+00);
      alephParameters[3] = 2.38808e+00;
      alephParameters[4] = 1.68165e+00;
      
    } else if(fppLowEn2011){ // pp low energy
      
      alephParameters[0]=0.031642;
      alephParameters[1]=22.353;
      alephParameters[2]=4.16239e-12;
      alephParameters[3]=2.61952;
      alephParameters[4]=5.76086;
      
    } else {  // pp no 1-pad (LHC10bc)
      
      alephParameters[0] = 0.0283086/0.97;
      alephParameters[1] = 2.63394e+01;
      alephParameters[2] = 5.04114e-11;
      alephParameters[3] = 2.12543e+00;
      alephParameters[4] = 4.88663e+00;
      
    }
    
  }
  
}

//-----------------------
void AliAODPidHF::SetBetheBloch() {
  /// Set Bethe Bloch Parameters
  
  Double_t alephParameters[5];
  GetTPCBetheBlochParams(alephParameters);
  fTPCResponse->SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2],alephParameters[3],alephParameters[4]);
  
  return;
}


//--------------------------------------------------------------------------
Int_t AliAODPidHF::GetnSigmaITS(AliAODTrack *track,Int_t species, Double_t &nsigma) const{
  /// get n sigma for ITS
  
  
  if (!CheckITSPIDStatus(track)) return -1;
  
  Double_t nsigmaITS=-999;
  
  if (fOldPid) {
    Double_t mom=track->P();
    AliAODPid *pidObj = track->GetDetPid();
    Double_t dedx=pidObj->GetITSsignal();
    
    AliITSPIDResponse itsResponse;
    AliPID::EParticleType type=AliPID::EParticleType(species);
    nsigmaITS = itsResponse.GetNumberOfSigmas(mom,dedx,type);
    
  } // old pid
  else { // new pid
    
    AliPID::EParticleType type=AliPID::EParticleType(species);
    nsigmaITS = fPidResponse->NumberOfSigmasITS(track,type);
    
  } //new pid
  
  nsigma = nsigmaITS;
  
  return 1;
  
}
//--------------------------------------------------------------------------
Int_t AliAODPidHF::GetnSigmaTPC(AliAODTrack *track, Int_t species, Double_t &nsigma) const{
  /// get n sigma for TPC
  
  if(!CheckTPCPIDStatus(track)) return -1;
  
  Double_t nsigmaTPC=-999;
  
  if(fOldPid){
    AliAODPid *pidObj = track->GetDetPid();
    Double_t dedx=pidObj->GetTPCsignal();
    Double_t mom = pidObj->GetTPCmomentum();
    if(mom>fPtThresholdTPC) return -2;
    UShort_t nTPCClus=pidObj->GetTPCsignalN();
    if(nTPCClus==0) {nTPCClus=track->GetTPCNcls();}
    AliPID::EParticleType type=AliPID::EParticleType(species);
    nsigmaTPC = fTPCResponse->GetNumberOfSigmas(mom,dedx,nTPCClus,type);
    nsigma=nsigmaTPC;
  } else{
    if(!fPidResponse) return -1;
    AliPID::EParticleType type=AliPID::EParticleType(species);
    nsigmaTPC = fPidResponse->NumberOfSigmasTPC(track,type);
    if(fApplyNsigmaTPCDataCorr && nsigmaTPC>-990.) {
      Float_t mean=0., sigma=1.; 
      GetNsigmaTPCMeanSigmaData(mean, sigma, type, track->GetTPCmomentum(), track->Eta());
      nsigmaTPC = (nsigmaTPC-mean)/sigma;
    }
    nsigma=nsigmaTPC;
  }
  return 1;
}

//-----------------------------

Int_t AliAODPidHF::GetnSigmaTOF(AliAODTrack *track,Int_t species, Double_t &nsigma) const{
  /// get n sigma for TOF
  
  if(!CheckTOFPIDStatus(track)) return -1;
  
  if(fPidResponse){
    track->SetTOFsignalTunedOnData(100000); // force tune-on-data to have latest development of tail parametrisation in old AODs
    nsigma = fPidResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)species);
    return 1;
  }else{
    AliFatal("To use TOF PID you need to attach AliPIDResponseTask");
    nsigma=-999.;
    return -1;
  }
}

//-----------------------
Bool_t AliAODPidHF::IsExcluded(AliAODTrack *track, Int_t labelTrack, Double_t nsigmaCut, TString detectors) {
  /// Exclude a given hypothesis (labelTracks) in detector
  
  if (detectors.Contains("ITS")) {
    
    AliInfo("Nothing to be done");
    /*
     Double_t nsigma=0.;
     if (GetnSigmaITS(track,labelTrack,nsigma)==1){
     if(nsigma>nsigmaCut) return kTRUE;
     }
     */
    return kFALSE;
    
  } else if (detectors.Contains("TPC")) {
    
    Double_t nsigma=0.;
    if (GetnSigmaTPC(track,labelTrack,nsigma)==1){
      if(nsigma>nsigmaCut) return kTRUE;
    }
    return kFALSE;
    
  } else if (detectors.Contains("TOF")) {
    
    Double_t nsigma=0.;
    if (GetnSigmaTOF(track,labelTrack,nsigma)==1){
      if(nsigma>nsigmaCut) return kTRUE;
    }
    return kFALSE;
    
  }
  return kFALSE;
  
}
//-----------------------
Bool_t AliAODPidHF::IsTOFPiKexcluded(AliAODTrack *track,Double_t nsigmaK){
  /// TOF proton compatibility
  
  if(!CheckTOFPIDStatus(track)) return 0;
  
  Double_t nsigma;
  if(GetnSigmaTOF(track,3,nsigma)==1){
    if(nsigma>nsigmaK) return kTRUE;
  }
  return kFALSE;
  /*  Double_t time[AliPID::kSPECIESN];
   Double_t sigmaTOFPid[AliPID::kSPECIES];
   AliAODPid *pidObj = track->GetDetPid();
   pidObj->GetIntegratedTimes(time);
   Double_t sigTOF=pidObj->GetTOFsignal();
   
   AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();
   if (event) {
   AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
   if (tofH && fPidResponse) {
   AliTOFPIDResponse TOFres = (AliTOFPIDResponse)fPidResponse->GetTOFResponse();
   sigTOF -= TOFres.GetStartTime(track->P());
   sigmaTOFPid[3]=TOFres.GetExpectedSigma(track->P(),time[3],AliPID::ParticleMass(3));
   }
   else  pidObj->GetTOFpidResolution(sigmaTOFPid);
   } else  pidObj->GetTOFpidResolution(sigmaTOFPid);
   Double_t sigmaTOFtrack;
   if (sigmaTOFPid[3]>0) sigmaTOFtrack=sigmaTOFPid[3];
   else sigmaTOFtrack=fTOFSigma;  // backward compatibility for old AODs
   
   if((sigTOF-time[3])>nsigmaK*sigmaTOFtrack)return kTRUE;// K, Pi excluded (->LIKELY A PROTON)
   
   return kFALSE;
   */
}

//--------------------------------------------------------------------------
void AliAODPidHF::SetPriorDistribution(AliPID::EParticleType type,TH1F *prior){
  
	///
	/// method setting the prior distributions to the AliPIDCombined object of the AliAODPidHF data member
	/// all the checks are done directly in the AliPIDCombined object
	///
  
	GetPidCombined()->SetPriorDistribution(type,prior);
}
//--------------------------------------------------------------------------
void AliAODPidHF::DrawPrior(AliPID::EParticleType type){
  
	///
	/// Drawing prior distribution for type "type"
  
	new TCanvas();
	GetPidCombined()->GetPriorDistribution(type)->Draw();
}

//-----------------------------
void AliAODPidHF::SetPriors(Double_t *priors, Int_t npriors){
  /// Set the values for the priors
  if (fnPriors != npriors) {
    if (fPriors) delete fPriors;
    fPriors = new Double_t[npriors];
    fnPriors = npriors;
  }
  for(Int_t i = 0; i < fnPriors; i++) fPriors[i] = priors[i];
}

//-----------------------------
void AliAODPidHF::SetPLimit(Double_t *plim, Int_t npLim) {
  /// Set limits of momentum ranges where different PID selections are applied
  if (fnPLimit != npLim) {
    if (fPLimit) delete fPLimit;
    fPLimit = new Double_t[npLim];
    fnPLimit = npLim;
  }
  for(Int_t i = 0; i < fnPLimit; i++) fPLimit[i] = plim[i];
}
//-----------------------------
void AliAODPidHF::SetPriorsHistos(TString priorFileName){
  /// Set histograms with priors
  
  for (Int_t ispecies=0;ispecies<AliPID::kSPECIES;++ispecies) {
    if(fPriorsH[ispecies]) delete fPriorsH[ispecies];
    TString nt ="name";
    nt+="_prior_";
    nt+=AliPID::ParticleName(ispecies);
  }
  TDirectory *current = gDirectory;
  TFile *priorFile=TFile::Open(priorFileName);
  if (priorFile) {
    TH1F* h3=static_cast<TH1F*>(priorFile->Get("priors3step9"));
    TH1F* h2=static_cast<TH1F*>(priorFile->Get("priors2step9"));
    TH1F* h1=static_cast<TH1F*>(priorFile->Get("priors1step9"));
    current->cd();
    fPriorsH[AliPID::kProton] = new TH1F(*h3);
    fPriorsH[AliPID::kKaon  ] = new TH1F(*h2);
    fPriorsH[AliPID::kPion  ] = new TH1F(*h1);
    priorFile->Close();
    delete priorFile;
    TF1 *salt=new TF1("salt","1.e-10",0,10);
    fPriorsH[AliPID::kProton]->Add(salt);
    fPriorsH[AliPID::kKaon  ]->Add(salt);
    fPriorsH[AliPID::kPion  ]->Add(salt);
    delete salt;
  }
}
//----------------------------------
void AliAODPidHF::SetUpCombinedPID(){
  /// Configuration of combined Bayesian PID
  
  fPidCombined->SetSelectedSpecies(AliPID::kSPECIES);
  if(!fDefaultPriors){
  	for (Int_t ispecies=0;ispecies<AliPID::kSPECIES;++ispecies) {
    	fPidCombined->SetPriorDistribution(static_cast<AliPID::EParticleType>(ispecies),fPriorsH[ispecies]);
  	}
  }else{
  	fPidCombined->SetDefaultTPCPriors();
  }
  switch (fCombDetectors){
    case kTPCTOF:
      fPidCombined->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
      break;
    case kTPCITS:
      fPidCombined->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetITS);
      break;
    case kTPC:
      fPidCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      break;
    case kTOF:
      fPidCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
      break;
  }
}


//-----------------------------
void AliAODPidHF::PrintAll() const {
  /// print the configuration
  printf("Detectors used for PID: ");
  if(fITS) printf("ITS ");
  if(fTPC) printf("TPC ");
  if(fTRD) printf("TRD ");
  if(fTOF) printf("TOF ");
  printf("\n");
  printf("Minimum TPC PID clusters = %d\n",fMinNClustersTPCPID);
  printf("Maximum momentum for using TPC PID = %f\n",fPtThresholdTPC);
  printf("TOF Mismatch probablility cut = %f\n",fCutTOFmismatch);
  printf("Maximum momentum for combined PID TPC PID = %f\n",fMaxTrackMomForCombinedPID);
  if(fOldPid){
    printf("Use OLD PID");
    printf("  fMC = %d\n",fMC);
    printf("  fPbPb = %d\n",fPbPb);
    printf("  fOnePad = %d\n",fOnePad);
    printf("  fMCLowEn2011 = %d\n",fMCLowEn2011);
    printf("  fppLowEn2011 = %d\n",fppLowEn2011);
  }
  printf("--- Matching algorithm = %d ---\n",fMatch);
  if(fMatch==1){
    if(fITS) printf("nSigmaITS = %.2f\n",fnSigma[4]);
    if(fTOF){
      printf("nSigmaTOF = %.2f\n",fnSigma[3]);
      if(fCompat) printf("Compatibility band at nSigmaTOF=%.2f for p>%.2f\n",fnSigmaCompat[1],fPCompatTOF);
    }
    if(fTPC){
      if(fAsym){
        printf("nSigmaTPC:\n");
        printf("   pt<%.2f      \t nsigmaTPC= %.2f\n",fPLimit[0],fnSigma[0]);
        printf("   %.2f<pt<%.2f \t nsigmaTPC= %.2f\n",fPLimit[0],fPLimit[1],fnSigma[1]);
        printf("   pt>%.2f      \t nsigmaTPC= %.2f\n",fPLimit[1],fnSigma[2]);
      }else{
        printf("nSigmaTPC = %.2f\n",fnSigma[0]);
      }
      if(fCompat) printf("Compatibility band at nSigmaTPC=%.2f\n",fnSigmaCompat[0]);
    }
  }else if(fMatch==4){
    printf("Cuts on sqrt(nSigmaTPC^2+nSigmaTOF^2):\n");
    printf(" Pions:   nSigma = %.2f\n",fMaxnSigmaCombined[0]);
    printf(" Kaons:   nSigma = %.2f\n",fMaxnSigmaCombined[1]);
    printf(" Protons: nSigma = %.2f\n",fMaxnSigmaCombined[2]);
  }else if(fMatch==5){
    printf("nSigma ranges:\n");
    printf(" Pions:   %.2f<nSigmaTPC<%.2f   %.2f<nSigmaTOF<%.2f\n",
           fMinnSigmaTPC[0],fMaxnSigmaTPC[0],fMinnSigmaTOF[0],fMaxnSigmaTOF[0]);
    printf(" Kaons:   %.2f<nSigmaTPC<%.2f   %.2f<nSigmaTOF<%.2f\n",
           fMinnSigmaTPC[1],fMaxnSigmaTPC[1],fMinnSigmaTOF[1],fMaxnSigmaTOF[1]);
    printf(" Protons: %.2f<nSigmaTPC<%.2f   %.2f<nSigmaTOF<%.2f\n",
           fMinnSigmaTPC[2],fMaxnSigmaTPC[2],fMinnSigmaTOF[2],fMaxnSigmaTOF[2]);
  } else if (fMatch == 10) {
    printf("Asymmetric PID using identification/compatibility bands as a function of track momentum p\n");
    printf("The following bands are set:\n");
    TString species[] = {"electron", "muon", "pion", "kaon", "proton"};
    TString detectors[] = {"ITS", "TPC", "TRD", "TOF"};
    for (Int_t s=0;s<AliPID::kSPECIES;s++) {
      for (Int_t d=0;d<4;d++) {
        if (fIdBandMin[s][d] && fIdBandMax[s][d]) {
          printf("  Identification band %s %s\n", species[s].Data(), detectors[d].Data());
        }
        if (fCompBandMin[s][d] && fCompBandMax[s][d]) {
          printf("  Compatibility band %s %s\n", species[s].Data(), detectors[d].Data());
        }
      }
    }
  }
}

//------------------
void AliAODPidHF::SetIdBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TH1F *min, TH1F *max) {
  Int_t spe = (Int_t) specie;
  Int_t det = (Int_t) detector;

  if (spe >= AliPID::kSPECIES || det > 3 || !min || !max) {
    AliError("Identification band not set");
    return;
  }

  TAxis *axis;
  HistFunc *histFunc;

  axis = min->GetXaxis();
  histFunc = new HistFunc(min);
  TF1 *minFunc = new TF1(Form("IdMin_%d_%d", spe, det), *histFunc, axis->GetBinLowEdge(axis->GetFirst()), axis->GetBinUpEdge(axis->GetLast()), 0, "HistFunc");

  axis = max->GetXaxis();
  histFunc = new HistFunc(max);
  TF1 *maxFunc = new TF1(Form("IdMax_%d_%d", spe, det), *histFunc, axis->GetBinLowEdge(axis->GetFirst()), axis->GetBinUpEdge(axis->GetLast()), 0, "HistFunc");

  SetIdBand(specie, detector, minFunc, maxFunc);
}

//------------------
void AliAODPidHF::SetIdBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TF1 *min, TF1 *max) {
  Int_t spe = (Int_t) specie;
  Int_t det = (Int_t) detector;

  if (spe >= AliPID::kSPECIES || det > 3 || !min || !max) {
    AliError("Identification band not set");
    return;
  }

  if (fIdBandMin[spe][det]) {
    delete fIdBandMin[spe][det];
  }
  fIdBandMin[spe][det] = new TF1(*min);

  if (fIdBandMax[spe][det]) {
    delete fIdBandMax[spe][det];
  }
  fIdBandMax[spe][det] = new TF1(*max);
}

//------------------
void AliAODPidHF::SetCompBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TH1F *min, TH1F *max) {
  Int_t spe = (Int_t) specie;
  Int_t det = (Int_t) detector;

  if (spe >= AliPID::kSPECIES || det > 3 || !min || !max) {
    AliError("Compatibility band not set");
    return;
  }

  TAxis *axis;
  HistFunc *histFunc;

  axis = min->GetXaxis();
  histFunc = new HistFunc(min);
  TF1 *minFunc = new TF1(Form("CompMin_%d_%d", spe, det), *histFunc, axis->GetBinLowEdge(axis->GetFirst()), axis->GetBinUpEdge(axis->GetLast()), 0, "HistFunc");

  axis = max->GetXaxis();
  histFunc = new HistFunc(max);
  TF1 *maxFunc = new TF1(Form("CompMax_%d_%d", spe, det), *histFunc, axis->GetBinLowEdge(axis->GetFirst()), axis->GetBinUpEdge(axis->GetLast()), 0, "HistFunc");

  SetCompBand(specie, detector, minFunc, maxFunc);
}

//------------------
void AliAODPidHF::SetCompBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TF1 *min, TF1 *max) {
  Int_t spe = (Int_t) specie;
  Int_t det = (Int_t) detector;

  if (spe >= AliPID::kSPECIES || det > 3 || !min || !max) {
    AliError("Compatibility band not set");
    return;
  }

  if (fCompBandMin[spe][det]) {
    delete fCompBandMin[spe][det];
  }
  fCompBandMin[spe][det] = new TF1(*min);

  if (fCompBandMax[spe][det]) {
    delete fCompBandMax[spe][det];
  }
  fCompBandMax[spe][det] = new TF1(*max);
}

//------------------
Bool_t AliAODPidHF::CheckDetectorPIDStatus(AliPIDResponse::EDetector detector, AliAODTrack* track) {
  switch (detector) {
    case AliPIDResponse::kITS:
      return CheckITSPIDStatus(track);
      break;
    case AliPIDResponse::kTPC:
      return CheckTPCPIDStatus(track);
      break;
    case AliPIDResponse::kTRD:
      return CheckTRDPIDStatus(track);
      break;
    case AliPIDResponse::kTOF:
      return CheckTOFPIDStatus(track);
      break;
    default:
      return kFALSE;
      break;
  }
}

//------------------
Float_t AliAODPidHF::NumberOfSigmas(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, AliAODTrack *track) {
  switch (detector) {
    case AliPIDResponse::kITS:
    {
      return fPidResponse->NumberOfSigmasITS(track, specie);
      break;
    }
    case AliPIDResponse::kTPC:
    {
      Double_t nsigmaTPC = fPidResponse->NumberOfSigmasTPC(track, specie);
      if(fApplyNsigmaTPCDataCorr && nsigmaTPC>-990.) {
        Float_t mean=0., sigma=1.; 
        GetNsigmaTPCMeanSigmaData(mean, sigma, specie, track->GetTPCmomentum(), track->Eta());
        nsigmaTPC = (nsigmaTPC-mean)/sigma;
      }
      return nsigmaTPC;
      break;
    }
    case AliPIDResponse::kTOF:
    {
      return fPidResponse->NumberOfSigmasTOF(track, specie);
      break;
    }
    default:
    {
      return -999.;
      break;
    }
  }
}

//------------------
Int_t AliAODPidHF::CheckBands(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, AliAODTrack *track) {
  /// \return Return: -1 for no match, 0 for compatible, 1 for identified

  Int_t spe = (Int_t) specie;
  Int_t det = (Int_t) detector;

  if (!fPidResponse || spe >= AliPID::kSPECIES) {
    return -1;
  }

  if (!CheckDetectorPIDStatus(detector, track)) {
    return 0;
  }

  Double_t P = track->P();

  Float_t nSigma = NumberOfSigmas(specie, detector, track);
  Float_t minContent, maxContent;
  Bool_t hasAnyBand = kFALSE;

  //Check if within identification band, return 1
  TF1 *IdBandMin = fIdBandMin[spe][det];
  TF1 *IdBandMax = fIdBandMax[spe][det];

  if (IdBandMin && IdBandMax) {
    minContent = IdBandMin->IsInside(&P) ? IdBandMin->Eval(P) : 0;
    maxContent = IdBandMax->IsInside(&P) ? IdBandMax->Eval(P) : 0;
    if (minContent != 0 || maxContent != 0) {
      //At least one identification band is set at this momentum
      hasAnyBand = kTRUE;
      if ((minContent == 0 || nSigma >= minContent) && (maxContent == 0 || nSigma <= maxContent)) {
        return 1;
      }
    }
  }

  //Check if within compatibility band, return 0
  TF1 *CompBandMin = fCompBandMin[spe][det];
  TF1 *CompBandMax = fCompBandMax[spe][det];

  if (CompBandMin && CompBandMax) {
    minContent = CompBandMin->IsInside(&P) ? CompBandMin->Eval(P) : 0;
    maxContent = CompBandMax->IsInside(&P) ? CompBandMax->Eval(P) : 0;
    if (minContent != 0 || maxContent != 0) {
      //At least one compatibility band is set at this momentum
      hasAnyBand = kTRUE;
      if ((minContent == 0 || nSigma >= minContent) && (maxContent == 0 || nSigma <= maxContent)) {
        return 0;
      }
    }
  }

  if (!hasAnyBand) {
    //No bands
    return 0;
  }

  //Bands were set and checked, but no match
  return -1;
}

//------------------
void AliAODPidHF::SetShiftedAsymmetricPID() {
  SetMatch(10);
  SetTPC(kTRUE);
  SetTOF(kTRUE);

  //TPC K: shift by -0.2
  TF1 *TPCCompBandMinK = new TF1("TPCCompBandMinK", "[0]", 0, 24); TPCCompBandMinK->SetParameter(0, -3.2);
  TF1 *TPCCompBandMaxK = new TF1("TPCCompBandMaxK", "[0]", 0, 24); TPCCompBandMaxK->SetParameter(0, 2.8);
  SetCompBand(AliPID::kKaon, AliPIDResponse::kTPC, TPCCompBandMinK, TPCCompBandMaxK);

  TF1 *TPCIdBandMinK = new TF1("TPCIdBandMinK", "[0]", 0, 24); TPCIdBandMinK->SetParameter(0, -2.2);
  TF1 *TPCIdBandMaxK = new TF1("TPCIdBandMaxK", "[0]", 0, 24); TPCIdBandMaxK->SetParameter(0, 1.8);
  SetIdBand(AliPID::kKaon, AliPIDResponse::kTPC, TPCIdBandMinK, TPCIdBandMaxK);

  //TPC pi: shift by -0.14
  TF1 *TPCCompBandMinPi = new TF1("TPCCompBandMinPi", "[0]", 0, 24); TPCCompBandMinPi->SetParameter(0, -3.14);
  TF1 *TPCCompBandMaxPi = new TF1("TPCCompBandMaxPi", "[0]", 0, 24); TPCCompBandMaxPi->SetParameter(0, 2.86);
  SetCompBand(AliPID::kPion, AliPIDResponse::kTPC, TPCCompBandMinPi, TPCCompBandMaxPi);

  TF1 *TPCIdBandMinPi = new TF1("TPCIdBandMinPi", "[0]", 0, 24); TPCIdBandMinPi->SetParameter(0, -2.14);
  TF1 *TPCIdBandMaxPi = new TF1("TPCIdBandMaxPi", "[0]", 0, 24); TPCIdBandMaxPi->SetParameter(0, 1.86);
  SetIdBand(AliPID::kPion, AliPIDResponse::kTPC, TPCIdBandMinPi, TPCIdBandMaxPi);

  //TOF K: shift by -0.1
  TF1 *TOFCompBandMinK = new TF1("TOFCompBandMinK", "[0]", 2, 24); TOFCompBandMinK->SetParameter(0, -3.1);
  TF1 *TOFCompBandMaxK = new TF1("TOFCompBandMaxK", "[0]", 2, 24); TOFCompBandMaxK->SetParameter(0, 2.9);
  SetCompBand(AliPID::kKaon, AliPIDResponse::kTOF, TOFCompBandMinK, TOFCompBandMaxK);

  TF1 *TOFIdBandMinK = new TF1("TOFIdBandMinK", "[0]", 0, 2); TOFIdBandMinK->SetParameter(0, -3.1);
  TF1 *TOFIdBandMaxK = new TF1("TOFIdBandMaxK", "[0]", 0, 2); TOFIdBandMaxK->SetParameter(0, 2.9);
  SetIdBand(AliPID::kKaon, AliPIDResponse::kTOF, TOFIdBandMinK, TOFIdBandMaxK);

  //TOF pi: shift by -0.15
  TF1 *TOFCompBandMinPi = new TF1("TOFCompBandMinPi", "[0]", 2, 24); TOFCompBandMinPi->SetParameter(0, -3.15);
  TF1 *TOFCompBandMaxPi = new TF1("TOFCompBandMaxPi", "[0]", 2, 24); TOFCompBandMaxPi->SetParameter(0, 2.85);
  SetCompBand(AliPID::kPion, AliPIDResponse::kTOF, TOFCompBandMinPi, TOFCompBandMaxPi);

  TF1 *TOFIdBandMinPi = new TF1("TOFIdBandMinPi", "[0]", 0, 2); TOFIdBandMinPi->SetParameter(0, -3.15);
  TF1 *TOFIdBandMaxPi = new TF1("TOFIdBandMaxPi", "[0]", 0, 2); TOFIdBandMaxPi->SetParameter(0, 2.85);
  SetIdBand(AliPID::kPion, AliPIDResponse::kTOF, TOFIdBandMinPi, TOFIdBandMaxPi);
}

//------------------
void AliAODPidHF::SetIdAsymmetricPID() {
  /// Set identification bands

  SetMatch(10);
  SetTPC(kTRUE);
  SetTOF(kTRUE);

  //TPC K
  Double_t TPCIdBandMinKBins[] = {0, 0.4, 0.5, 0.6, 0.9, 24};
  TH1F *TPCIdBandMinK = new TH1F("TPCIdBandMinK", "TPC Id Band Min K", 5, TPCIdBandMinKBins);
  TPCIdBandMinK->SetBinContent(1, -3); //0   -  0.4
  TPCIdBandMinK->SetBinContent(2, -2); //0.4 -  0.5
  TPCIdBandMinK->SetBinContent(3, -3); //0.5 -  0.6
  TPCIdBandMinK->SetBinContent(4, -2); //0.6 -  0.9
  TPCIdBandMinK->SetBinContent(5, -3); //0.9 - 24

  Double_t TPCIdBandMaxKBins[] = {0, 0.6, 0.7, 24};
  TH1F *TPCIdBandMaxK = new TH1F("TPCIdBandMaxK", "TPC Id Band Max K", 3, TPCIdBandMaxKBins);
  TPCIdBandMaxK->SetBinContent(1, 3); //0   -  0.6
  TPCIdBandMaxK->SetBinContent(2, 2); //0.6 -  0.7
  TPCIdBandMaxK->SetBinContent(3, 3); //0.7 - 24

  SetIdBand(AliPID::kKaon, AliPIDResponse::kTPC, TPCIdBandMinK, TPCIdBandMaxK);
  GetIdBandMin(AliPID::kKaon, AliPIDResponse::kTPC)->SetNpx(5000);
  GetIdBandMax(AliPID::kKaon, AliPIDResponse::kTPC)->SetNpx(5000);

  //TPC pi
  Double_t TPCIdBandMinpiBins[] = {0, 24};
  TH1F *TPCIdBandMinpi = new TH1F("TPCIdBandMinpi", "TPC Id Band Min pi", 1, TPCIdBandMinpiBins);
  TPCIdBandMinpi->SetBinContent(1, -3); //0 - 24

  Double_t TPCIdBandMaxpiBins[] = {0, 0.7, 0.9, 1.3, 1.4, 24};
  TH1F *TPCIdBandMaxpi = new TH1F("TPCIdBandMaxpi", "TPC Id Band Max pi", 5, TPCIdBandMaxpiBins);
  TPCIdBandMaxpi->SetBinContent(1, 3); //0   -  0.7
  TPCIdBandMaxpi->SetBinContent(2, 2); //0.7 -  0.9
  TPCIdBandMaxpi->SetBinContent(3, 3); //0.9 -  1.3
  TPCIdBandMaxpi->SetBinContent(4, 2); //1.3 -  1.4
  TPCIdBandMaxpi->SetBinContent(5, 3); //1.4 - 24

  SetIdBand(AliPID::kPion, AliPIDResponse::kTPC, TPCIdBandMinpi, TPCIdBandMaxpi);
  GetIdBandMin(AliPID::kPion, AliPIDResponse::kTPC)->SetNpx(5000);
  GetIdBandMax(AliPID::kPion, AliPIDResponse::kTPC)->SetNpx(5000);

  //TOF K
  TF1 *TOFIdBandMinK = new TF1("TOFIdBandMinK", "[0]", 0, 24); TOFIdBandMinK->SetParameter(0, -3);
  TF1 *TOFIdBandMaxK = new TF1("TOFIdBandMaxK", "[0]", 0, 24); TOFIdBandMaxK->SetParameter(0, 3);
  
  SetIdBand(AliPID::kKaon, AliPIDResponse::kTOF, TOFIdBandMinK, TOFIdBandMaxK);

  //TOF pi
  TF1 *TOFIdBandMinPi = new TF1("TOFIdBandMinPi", "[0]", 0, 24); TOFIdBandMinPi->SetParameter(0, -3);
  TF1 *TOFIdBandMaxPi = new TF1("TOFIdBandMaxPi", "[0]", 0, 24); TOFIdBandMaxPi->SetParameter(0, 3);
  
  SetIdBand(AliPID::kPion, AliPIDResponse::kTOF, TOFIdBandMinPi, TOFIdBandMaxPi);
}

//------------------
void AliAODPidHF::SetIdCompAsymmetricPID() {
  /// Set compatibility and identification bands

  SetMatch(10);
  SetTPC(kTRUE);
  SetTOF(kTRUE);

  //TPC K
  TF1 *TPCCompBandMinK = new TF1("TPCCompBandMinK", "[0]", 0, 24); TPCCompBandMinK->SetParameter(0, -3);
  TF1 *TPCCompBandMaxK = new TF1("TPCCompBandMaxK", "[0]", 0, 24); TPCCompBandMaxK->SetParameter(0, 3);
  
  SetCompBand(AliPID::kKaon, AliPIDResponse::kTPC, TPCCompBandMinK, TPCCompBandMaxK);

  Double_t TPCIdBandMinKBins[6] = {0, 0.45, 0.55, 0.7, 1.1, 24};
  TH1F *TPCIdBandMinK = new TH1F("TPCIdBandMinK", "TPC Id Band Min K", 5, TPCIdBandMinKBins);
  TPCIdBandMinK->SetBinContent(1, -2); //0-0.45
  TPCIdBandMinK->SetBinContent(2, -1); //0.45-0.55
  TPCIdBandMinK->SetBinContent(3, -2); //0.55-0.7
  TPCIdBandMinK->SetBinContent(4, -1); //0.7-1.1
  TPCIdBandMinK->SetBinContent(5, -2); //1.1-24
  
  Double_t TPCIdBandMaxKBins[4] = {0, 0.5, 0.7, 24};
  TH1F *TPCIdBandMaxK = new TH1F("TPCIdBandMaxK", "TPC Id Band Max K", 3, TPCIdBandMaxKBins);
  TPCIdBandMaxK->SetBinContent(1, 2); //0-0.5
  TPCIdBandMaxK->SetBinContent(2, 1); //0.5-0.7
  TPCIdBandMaxK->SetBinContent(3, 2); //0.7-24

  SetIdBand(AliPID::kKaon, AliPIDResponse::kTPC, TPCIdBandMinK, TPCIdBandMaxK);
  GetIdBandMin(AliPID::kKaon, AliPIDResponse::kTPC)->SetNpx(5000);
  GetIdBandMax(AliPID::kKaon, AliPIDResponse::kTPC)->SetNpx(5000);

  //TPC pi
  TF1 *TPCCompBandMinpi = new TF1("TPCCompBandMinpi", "[0]", 0, 24); TPCCompBandMinpi->SetParameter(0, -3);
  TF1 *TPCCompBandMaxpi = new TF1("TPCCompBandMaxpi", "[0]", 0, 24); TPCCompBandMaxpi->SetParameter(0, 3);
  
  SetCompBand(AliPID::kPion, AliPIDResponse::kTPC, TPCCompBandMinpi, TPCCompBandMaxpi);

  Double_t TPCIdBandMinpiBins[2] = {0, 24};
  TH1F *TPCIdBandMinpi = new TH1F("TPCIdBandMinpi", "TPC Id Band Min pi", 1, TPCIdBandMinpiBins);
  TPCIdBandMinpi->SetBinContent(1, -2); //0-24
  
  Double_t TPCIdBandMaxpiBins[4] = {0, 0.7, 1.7, 24};
  TH1F *TPCIdBandMaxpi = new TH1F("TPCIdBandMaxpi", "TPC Id Band Max pi", 3, TPCIdBandMaxpiBins);
  TPCIdBandMaxpi->SetBinContent(1, 2); //0-0.7
  TPCIdBandMaxpi->SetBinContent(2, 1); //0.7-1.7
  TPCIdBandMaxpi->SetBinContent(3, 2); //1.7-24
  
  SetIdBand(AliPID::kPion, AliPIDResponse::kTPC, TPCIdBandMinpi, TPCIdBandMaxpi);
  GetIdBandMin(AliPID::kPion, AliPIDResponse::kTPC)->SetNpx(5000);
  GetIdBandMax(AliPID::kPion, AliPIDResponse::kTPC)->SetNpx(5000);

  //TOF K
  TF1 *TOFCompBandMinK = new TF1("TOFCompBandMinK", "[0]", 2, 24); TOFCompBandMinK->SetParameter(0, -3);
  TF1 *TOFCompBandMaxK = new TF1("TOFCompBandMaxK", "[0]", 2, 24); TOFCompBandMaxK->SetParameter(0, 3);

  SetCompBand(AliPID::kKaon, AliPIDResponse::kTOF, TOFCompBandMinK, TOFCompBandMaxK);

  TF1 *TOFIdBandMinK = new TF1("TOFIdBandMinK", "[0]", 0, 2); TOFIdBandMinK->SetParameter(0, -3);
  TF1 *TOFIdBandMaxK = new TF1("TOFIdBandMaxK", "[0]", 0, 2); TOFIdBandMaxK->SetParameter(0, 3);

  SetIdBand(AliPID::kKaon, AliPIDResponse::kTOF, TOFIdBandMinK, TOFIdBandMaxK);

  //TOF pi
  TF1 *TOFCompBandMinpi = new TF1("TOFCompBandMinpi", "[0]", 2, 24); TOFCompBandMinpi->SetParameter(0, -3);
  TF1 *TOFCompBandMaxpi = new TF1("TOFCompBandMaxpi", "[0]", 2, 24); TOFCompBandMaxpi->SetParameter(0, 3);

  SetCompBand(AliPID::kPion, AliPIDResponse::kTOF, TOFCompBandMinpi, TOFCompBandMaxpi);

  TF1 *TOFIdBandMinpi = new TF1("TOFIdBandMinpi", "[0]", 0, 2); TOFIdBandMinpi->SetParameter(0, -3);
  TF1 *TOFIdBandMaxpi = new TF1("TOFIdBandMaxpi", "[0]", 0, 2); TOFIdBandMaxpi->SetParameter(0, 3);

  SetIdBand(AliPID::kPion, AliPIDResponse::kTOF, TOFIdBandMinpi, TOFIdBandMaxpi);
}

//------------------
void AliAODPidHF::EnableNsigmaTPCDataCorr(Int_t run, Int_t system, Bool_t isPass1) {

  fApplyNsigmaTPCDataCorr = kTRUE;
  SetNsigmaTPCDataDrivenCorrection(run, system, fNPbinsNsigmaTPCDataCorr, fPlimitsNsigmaTPCDataCorr, fNEtabinsNsigmaTPCDataCorr, fEtalimitsNsigmaTPCDataCorr, fMeanNsigmaTPCPionData, fMeanNsigmaTPCKaonData, fMeanNsigmaTPCProtonData, fSigmaNsigmaTPCPionData, fSigmaNsigmaTPCKaonData, fSigmaNsigmaTPCProtonData, isPass1);
}

//------------------
void AliAODPidHF::GetNsigmaTPCMeanSigmaData(Float_t &mean, Float_t &sigma, AliPID::EParticleType species, Float_t pTPC, Float_t eta) const {
    
  Int_t bin = TMath::BinarySearch(fNPbinsNsigmaTPCDataCorr,fPlimitsNsigmaTPCDataCorr,pTPC);
  if(bin<0) bin=0; //underflow --> equal to min value
  else if(bin>fNPbinsNsigmaTPCDataCorr-1) bin=fNPbinsNsigmaTPCDataCorr-1; //overflow --> equal to max value

  Int_t etabin = TMath::BinarySearch(fNEtabinsNsigmaTPCDataCorr,fEtalimitsNsigmaTPCDataCorr,TMath::Abs(eta));
  if(etabin<0) etabin=0; //underflow --> equal to min value
  else if(etabin>fNEtabinsNsigmaTPCDataCorr-1) etabin=fNEtabinsNsigmaTPCDataCorr-1; //overflow --> equal to max value

  switch(species) {
    case AliPID::kPion: 
    {
      mean = fMeanNsigmaTPCPionData[etabin][bin];
      sigma = fSigmaNsigmaTPCPionData[etabin][bin];
      break;
    }
    case AliPID::kKaon: 
    {
      mean = fMeanNsigmaTPCKaonData[etabin][bin];
      sigma = fSigmaNsigmaTPCKaonData[etabin][bin];
      break;
    }
    case AliPID::kProton: 
    {
      mean = fMeanNsigmaTPCProtonData[etabin][bin];
      sigma = fSigmaNsigmaTPCProtonData[etabin][bin];
      break;
    }
    default: 
    {
      mean = 0.;
      sigma = 1.;
      break;
    }
  }
}

//___________________________________________________________________________________//
void AliAODPidHF::SetNsigmaTPCDataDrivenCorrection(Int_t run, Int_t system, Int_t &nPbins, Float_t Plims[kMaxPBins+1], Int_t &nEtabins, Float_t absEtalims[kMaxEtaBins+1], 
  vector<vector<Float_t> > &meanNsigmaTPCpion, vector<vector<Float_t> > &meanNsigmaTPCkaon, vector<vector<Float_t> > &meanNsigmaTPCproton, 
  vector<vector<Float_t> > &sigmaNsigmaTPCpion, vector<vector<Float_t> > &sigmaNsigmaTPCkaon, vector<vector<Float_t> > &sigmaNsigmaTPCproton, Bool_t isPass1) 
  {

  meanNsigmaTPCpion.resize(kMaxEtaBins,vector<Float_t>(kMaxPBins,0.));
  meanNsigmaTPCkaon.resize(kMaxEtaBins,vector<Float_t>(kMaxPBins,0.));
  meanNsigmaTPCproton.resize(kMaxEtaBins,vector<Float_t>(kMaxPBins,0.));
  sigmaNsigmaTPCpion.resize(kMaxEtaBins,vector<Float_t>(kMaxPBins,1.));
  sigmaNsigmaTPCkaon.resize(kMaxEtaBins,vector<Float_t>(kMaxPBins,1.));
  sigmaNsigmaTPCproton.resize(kMaxEtaBins,vector<Float_t>(kMaxPBins,1.));

  if(run>=295585 && run<=296623 && system==kPbPb010) { //LHC18q 0-10%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.225065, -0.119186, -0.099263, -0.106165, -0.113968, -0.180788, -0.306835, -0.365154},
                                         {-0.157762, -0.072164, -0.066882, -0.079671, -0.087264, -0.167757, -0.299523, -0.287218},
                                         {-0.013474, 0.011867, -0.013237, -0.026186, -0.048963, -0.119485, -0.211451, -0.209846},
                                         {0.329239, 0.175729, 0.115309, 0.066772, 0.002108, -0.052078, -0.121908, -0.133652},
                                         {0.648693, 0.312451, 0.207059, 0.124000, 0.061547, 0.044737, 0.013488, -0.024819}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.168045, -0.184170, -0.134005, -0.153342, -0.187885, -0.247042, -0.345098, -0.451211},
                                         {0.009426, -0.126109, -0.132931, -0.133621, -0.153672, -0.214908, -0.307614, -0.156097},
                                         {0.137680, -0.022762, -0.072650, -0.100962, -0.138850, -0.190713, -0.275153, -0.089141},
                                         {0.337654, 0.201968, 0.067451, -0.021638, -0.119303, -0.166833, -0.236204, -0.080545},
                                         {0.385523, 0.538211, 0.257916, 0.043189, -0.133722, -0.162471, -0.223222, -0.092731}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.190173, -0.195611, -0.116308, -0.342037, -0.131468, -0.225489, -0.141616, 0.077360},
                                           {-0.104356, -0.083717, -0.024616, -0.244703, -0.164418, -0.195518, -0.107934, 0.287976},
                                           {-0.026761, -0.045308, 0.044370, -0.152323, -0.091770, -0.160313, -0.083499, 0.167515},
                                           {0.014289, -0.011640, 0.150733, -0.063745, 0.055933, -0.083332, -0.015881, 0.185935},
                                           {0.062978, -0.018201, 0.205506, 0.060088, 0.156797, -0.022258, 0.002592, 0.211372}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {1.037687, 1.057941, 1.058772, 1.054653, 1.061417, 1.118399, 1.151653, 1.113338},
                                          {1.028331, 1.049659, 1.045960, 1.047649, 1.053873, 1.104584, 1.135918, 1.082807},
                                          {1.028647, 1.044127, 1.045863, 1.039959, 1.048004, 1.084940, 1.102891, 1.066465},
                                          {1.033039, 1.035785, 1.037631, 1.024840, 1.038347, 1.064217, 1.087659, 1.085300},
                                          {1.030865, 1.017514, 1.027900, 1.024781, 1.037124, 1.052824, 1.080056, 1.080240}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.018388, 1.105860, 1.088447, 1.039056, 1.011861, 0.968054, 0.953466, 0.993792},
                                          {1.032788, 1.071545, 1.087832, 1.053688, 1.012775, 0.970389, 0.954898, 0.935175},
                                          {1.030084, 1.045686, 1.073807, 1.039074, 0.999492, 0.958561, 0.939375, 0.944366},
                                          {0.990590, 1.058523, 1.080466, 1.039387, 1.000431, 0.951797, 0.922074, 0.886299},
                                          {1.000247, 1.037122, 1.077713, 1.036642, 0.998754, 0.950503, 0.911743, 0.819221}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.017302, 1.089266, 1.131230, 1.174314, 1.212759, 1.006333, 1.013844, 1.208998},
                                            {1.006405, 1.082345, 1.118237, 1.146892, 1.086678, 1.007054, 1.017808, 1.268190},
                                            {1.013036, 1.071208, 1.111884, 1.136804, 1.107289, 1.003273, 1.000704, 1.190774},
                                            {1.026807, 1.084265, 1.112386, 1.158228, 1.191843, 1.014501, 0.998355, 1.129179},
                                            {1.006287, 1.067167, 1.113327, 1.198577, 1.208199, 1.034627, 1.002670, 1.109491}
                                           };

    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.656082, -0.604754, -0.63195, -0.669819, -0.708323, -0.746162, -0.800557, -0.893548},
                                              {-0.711512, -0.686848, -0.711177, -0.752377, -0.781678, -0.838341, -0.80682, -0.917012},
                                              {-0.650884, -0.706274, -0.752449, -0.798673, -0.835543, -0.859084, -0.823375, -0.854207},
                                              {-0.479705, -0.700093, -0.815494, -0.895986, -0.958185, -0.980633, -0.967819, -0.996364},
                                              {-0.264033, -0.637537, -0.828537, -0.952223, -1.06091, -1.10456, -1.09789, -1.15579}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.48114, -0.672897, -0.625657, -0.776678, -0.786824, -0.708909, -0.822472, -0.491422},
                                              {-0.432004, -0.71966, -0.708949, -0.870034, -0.856239, -0.825942, -0.871391, -1.17962},
                                              {-0.336167, -0.71892, -0.735327, -0.93073, -0.864011, -0.891611, -0.924604, -0.735026},
                                              {-0.391054, -0.716122, -0.796868, -1.0208, -0.984637, -0.998813, -1.01377, -1.06832},
                                              {-0.551251, -0.696801, -0.815058, -1.05691, -1.06688, -1.06648, -1.07023, -1.05183}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {-0.200581, -0.16751, -0.451043, -0.952266, -0.852847, -0.760682, -0.676723, -0.603716},
                                                {-0.123522, -0.128086, -0.444873, -0.846087, -0.988114, -1.05446, -0.761678, -0.785548},
                                                {-0.100534, -0.1431, -0.448783, -0.8385, -0.843197, -1.05925, -0.878891, -0.69573},
                                                {-0.233023, -0.317509, -0.598837, -0.945108, -1.01043, -1.21354, -0.99634, -0.915479},
                                                {-0.391233, -0.516667, -0.768414, -1.00696, -1.03589, -1.26272, -1.02806, -0.994112}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {0.986692, 1.01991, 1.00333, 0.986744, 0.981785, 1.04139, 1.0638, 1.09162},
                                               {0.968236, 0.999018, 0.984673, 0.963658, 0.963348, 0.991749, 1.00931, 1.07714},
                                               {0.948544, 0.971808, 0.957514, 0.938766, 0.936994, 0.987149, 0.957657, 0.994133},
                                               {0.92104, 0.931385, 0.916291, 0.896921, 0.890271, 0.926601, 0.891902, 0.905638},
                                               {0.909424, 0.89589, 0.881729, 0.860767, 0.842961, 0.873783, 0.83704, 0.758586}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {0.891168, 1.00326, 1.04053, 0.952367, 0.919632, 0.951424, 0.902434, 1.07759},
                                               {0.853805, 0.968212, 1.01839, 0.930692, 0.918841, 0.92428, 0.913715, 0.929855},
                                               {0.830364, 0.911894, 0.987625, 0.912264, 0.939659, 0.906548, 0.893721, 1.00634},
                                               {0.718803, 0.882484, 0.959253, 0.870302, 0.907273, 0.881795, 0.867757, 0.906373},
                                               {0.688955, 0.855596, 0.932222, 0.839553, 0.867639, 0.855212, 0.845177, 0.90448}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {0.771648, 0.841043, 0.917283, 1.12449, 1.0023, 0.952976, 0.963016, 1.01111},
                                                 {0.752951, 0.825488, 0.883897, 1.02998, 1.07061, 0.866346, 0.952794, 0.916068},
                                                 {0.72623, 0.799905, 0.860896, 0.996524, 1.02278, 0.866737, 0.834891, 0.988606},
                                                 {0.70374, 0.759504, 0.811721, 0.944681, 0.96305, 0.81128, 0.85648, 0.900749},
                                                 {0.702538, 0.723393, 0.781419, 0.83867, 0.940137, 0.785817, 0.841202, 0.859564}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else if(run>=295585 && run<=296623 && system==kPbPb1030) { //LHC18q 10-30%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.317593, -0.183524, -0.127651, -0.089234, -0.062979, -0.080720, -0.122664, -0.223814},
                                         {-0.229056, -0.110411, -0.070411, -0.041038, -0.017415, -0.041055, -0.095676, -0.090063},
                                         {-0.110306, -0.048978, -0.035642, -0.020479, -0.011696, -0.038913, -0.080831, -0.073858},
                                         {0.145704, 0.046161, 0.017817, -0.000861, -0.026094, -0.049646, -0.085138, -0.045434},
                                         {0.400589, 0.127285, 0.060913, 0.001900, -0.023850, -0.027839, 0.005667, -0.022088}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.067543, -0.097760, -0.115582, -0.195158, -0.221149, -0.251527, -0.255514, -0.076460},
                                         {0.041216, -0.028869, -0.099025, -0.160005, -0.172629, -0.183296, -0.271547, -0.213141},
                                         {0.171465, 0.010551, -0.055422, -0.128486, -0.158187, -0.183289, -0.283075, -0.100889},
                                         {0.261771, 0.139781, 0.005954, -0.111944, -0.189801, -0.207535, -0.303148, -0.114337},
                                         {0.216975, 0.393007, 0.147375, -0.095322, -0.235820, -0.247256, -0.299180, -0.135729}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.156664, -0.110028, 0.014113, -0.118975, -0.094321, -0.179341, -0.131990, -0.137500},
                                           {-0.053650, -0.026483, 0.091400, -0.038219, -0.059285, -0.118064, -0.030585, 0.013138},
                                           {-0.027011, -0.034813, 0.093307, -0.019499, -0.059152, -0.118164, -0.036743, 0.050302},
                                           {-0.028588, -0.098622, 0.096340, -0.019024, -0.055346, -0.115118, -0.067764, 0.013925},
                                           {-0.064779, -0.159609, 0.072633, 0.060540, -0.037721, -0.109379, -0.109764, 0.066457}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.998493, 1.010638, 1.016402, 1.017521, 1.022144, 1.050891, 1.041638, 0.821282},
                                          {0.988449, 1.005340, 1.008359, 1.010262, 1.010872, 1.031073, 1.015050, 0.908931},
                                          {0.982872, 0.995834, 1.001442, 0.994735, 1.003186, 1.020325, 0.991257, 0.907473},
                                          {0.980495, 0.983644, 0.987519, 0.988086, 0.986011, 0.998937, 0.988648, 0.922585},
                                          {0.979858, 0.967431, 0.976917, 0.979059, 0.979981, 0.990846, 0.968238, 0.929777}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {0.995001, 1.038605, 1.041353, 1.011128, 0.981057, 0.950883, 0.954101, 0.813112},
                                          {0.979466, 1.009251, 1.026508, 1.007649, 0.975699, 0.960884, 0.917550, 0.784654},
                                          {0.985526, 0.990439, 1.013789, 1.002628, 0.969205, 0.948116, 0.898339, 0.796312},
                                          {0.973591, 0.993490, 1.017769, 0.994345, 0.967159, 0.942014, 0.883009, 0.820896},
                                          {0.973191, 0.984407, 1.027771, 0.988651, 0.957722, 0.943232, 0.917779, 0.841613}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.025847, 1.064472, 1.079804, 1.080703, 1.128969, 1.006314, 1.024577, 0.970529},
                                            {1.000369, 1.039794, 1.059678, 1.065775, 1.116239, 1.009144, 1.025479, 0.924459},
                                            {0.989583, 1.046193, 1.062586, 1.059619, 1.061026, 0.992778, 0.993071, 1.051145},
                                            {1.007853, 1.053055, 1.065222, 1.069423, 1.106152, 0.995711, 0.987215, 1.027091},
                                            {0.989768, 1.051719, 1.069527, 1.082615, 1.116969, 1.006839, 0.983334, 1.031133}
                                           };

    // using 0-10% for 10-30% in case of pass1
    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.656082, -0.604754, -0.63195, -0.669819, -0.708323, -0.746162, -0.800557, -0.893548},
                                              {-0.711512, -0.686848, -0.711177, -0.752377, -0.781678, -0.838341, -0.80682, -0.917012},
                                              {-0.650884, -0.706274, -0.752449, -0.798673, -0.835543, -0.859084, -0.823375, -0.854207},
                                              {-0.479705, -0.700093, -0.815494, -0.895986, -0.958185, -0.980633, -0.967819, -0.996364},
                                              {-0.264033, -0.637537, -0.828537, -0.952223, -1.06091, -1.10456, -1.09789, -1.15579}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.48114, -0.672897, -0.625657, -0.776678, -0.786824, -0.708909, -0.822472, -0.491422},
                                              {-0.432004, -0.71966, -0.708949, -0.870034, -0.856239, -0.825942, -0.871391, -1.17962},
                                              {-0.336167, -0.71892, -0.735327, -0.93073, -0.864011, -0.891611, -0.924604, -0.735026},
                                              {-0.391054, -0.716122, -0.796868, -1.0208, -0.984637, -0.998813, -1.01377, -1.06832},
                                              {-0.551251, -0.696801, -0.815058, -1.05691, -1.06688, -1.06648, -1.07023, -1.05183}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {-0.200581, -0.16751, -0.451043, -0.952266, -0.852847, -0.760682, -0.676723, -0.603716},
                                                {-0.123522, -0.128086, -0.444873, -0.846087, -0.988114, -1.05446, -0.761678, -0.785548},
                                                {-0.100534, -0.1431, -0.448783, -0.8385, -0.843197, -1.05925, -0.878891, -0.69573},
                                                {-0.233023, -0.317509, -0.598837, -0.945108, -1.01043, -1.21354, -0.99634, -0.915479},
                                                {-0.391233, -0.516667, -0.768414, -1.00696, -1.03589, -1.26272, -1.02806, -0.994112}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {0.986692, 1.01991, 1.00333, 0.986744, 0.981785, 1.04139, 1.0638, 1.09162},
                                               {0.968236, 0.999018, 0.984673, 0.963658, 0.963348, 0.991749, 1.00931, 1.07714},
                                               {0.948544, 0.971808, 0.957514, 0.938766, 0.936994, 0.987149, 0.957657, 0.994133},
                                               {0.92104, 0.931385, 0.916291, 0.896921, 0.890271, 0.926601, 0.891902, 0.905638},
                                               {0.909424, 0.89589, 0.881729, 0.860767, 0.842961, 0.873783, 0.83704, 0.758586}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {0.891168, 1.00326, 1.04053, 0.952367, 0.919632, 0.951424, 0.902434, 1.07759},
                                               {0.853805, 0.968212, 1.01839, 0.930692, 0.918841, 0.92428, 0.913715, 0.929855},
                                               {0.830364, 0.911894, 0.987625, 0.912264, 0.939659, 0.906548, 0.893721, 1.00634},
                                               {0.718803, 0.882484, 0.959253, 0.870302, 0.907273, 0.881795, 0.867757, 0.906373},
                                               {0.688955, 0.855596, 0.932222, 0.839553, 0.867639, 0.855212, 0.845177, 0.90448}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {0.771648, 0.841043, 0.917283, 1.12449, 1.0023, 0.952976, 0.963016, 1.01111},
                                                 {0.752951, 0.825488, 0.883897, 1.02998, 1.07061, 0.866346, 0.952794, 0.916068},
                                                 {0.72623, 0.799905, 0.860896, 0.996524, 1.02278, 0.866737, 0.834891, 0.988606},
                                                 {0.70374, 0.759504, 0.811721, 0.944681, 0.96305, 0.81128, 0.85648, 0.900749},
                                                 {0.702538, 0.723393, 0.781419, 0.83867, 0.940137, 0.785817, 0.841202, 0.859564}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else if(run>=295585 && run<=296623 && system==kPbPb3050) { //LHC18q 30-50%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.324361, -0.160800, -0.078253, -0.019182, 0.037198, 0.063881, 0.039371, 0.085383},
                                         {-0.224308, -0.076201, -0.008182, 0.041642, 0.084393, 0.103260, 0.091632, 0.128449},
                                         {-0.159696, -0.059762, -0.020894, 0.020020, 0.044053, 0.032295, 0.033186, 0.094967},
                                         {0.001655, -0.074936, -0.067459, -0.066945, -0.062178, -0.076584, -0.077504, -0.070499},
                                         {0.134629, -0.086372, -0.119054, -0.148353, -0.149948, -0.139768, -0.100795, -0.077748}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {0.046151, 0.016388, -0.050481, -0.149588, -0.196056, -0.180752, -0.190079, 0.143470},
                                         {0.153567, 0.057117, -0.026156, -0.086432, -0.123796, -0.107294, -0.149493, 0.055667},
                                         {0.168326, 0.017215, -0.032108, -0.101793, -0.143544, -0.135115, -0.218321, 0.031318},
                                         {0.004905, 0.052538, -0.064678, -0.178787, -0.250630, -0.244642, -0.310286, -0.074095},
                                         {0.035571, 0.225276, -0.017677, -0.244894, -0.375957, -0.375866, -0.392168, -0.232767}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.097519, -0.048271, 0.106258, 0.019264, -0.005813, -0.074018, -0.061853, -0.144826},
                                           {0.011767, 0.049442, 0.186854, 0.099568, 0.036226, -0.018250, 0.042317, 0.082425},
                                           {-0.005354, -0.007578, 0.134062, 0.057955, -0.004145, -0.052637, -0.015561, 0.019098},
                                           {-0.043469, -0.149406, 0.022551, -0.041238, -0.103359, -0.149070, -0.127941, -0.058903},
                                           {-0.114942, -0.274365, -0.003453, -0.032073, -0.171188, -0.245157, -0.265916, -0.183356}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.952339, 0.966346, 0.973800, 0.983150, 0.996656, 0.993748, 0.978304, 0.946653},
                                          {0.942172, 0.957841, 0.966333, 0.974654, 0.982796, 0.977678, 0.955283, 0.846376},
                                          {0.929873, 0.945356, 0.953387, 0.960463, 0.968373, 0.960203, 0.925668, 0.889382},
                                          {1.000642, 0.933064, 0.944367, 0.944275, 0.948921, 0.939550, 0.905768, 0.845528},
                                          {0.926115, 0.920114, 0.934629, 0.935449, 0.940571, 0.932503, 0.901855, 0.846489}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {0.969557, 0.980560, 0.997688, 0.987933, 0.958851, 0.944345, 0.909614, 0.890382},
                                          {0.981124, 0.960085, 0.968321, 0.979377, 0.955287, 0.947831, 0.910214, 0.819880},
                                          {0.938036, 0.935218, 0.966406, 0.960659, 0.941661, 0.934175, 0.858276, 0.801148},
                                          {0.978230, 0.940891, 0.967743, 0.952380, 0.934358, 0.921585, 0.872763, 0.846987},
                                          {0.854273, 0.936372, 0.972363, 0.953146, 0.924515, 0.902529, 0.907656, 0.802625}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.023079, 1.040636, 1.037471, 1.018790, 1.035993, 0.997710, 1.028290, 1.113908},
                                            {0.985836, 1.024006, 1.024505, 1.009177, 1.023907, 0.990433, 1.011937, 1.188029},
                                            {0.990252, 1.031607, 1.025671, 1.008726, 1.007395, 0.978731, 0.972668, 1.012270},
                                            {1.017223, 1.034908, 1.038513, 1.010218, 0.984319, 0.970254, 0.958761, 1.014806},
                                            {1.006925, 1.044808, 0.715285, 1.013759, 0.986755, 0.972059, 0.944017, 0.966145}
                                           };

    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.537046, -0.427744, -0.411915, -0.42242, -0.445157, -0.423209, -0.403354, -0.39011},
                                              {-0.451747, -0.358803, -0.347827, -0.360332, -0.379419, -0.373067, -0.309076, -0.201842},
                                              {-0.37701, -0.342516, -0.352131, -0.375476, -0.397523, -0.380363, -0.348125, -0.334354},
                                              {-0.340843, -0.438716, -0.504516, -0.554322, -0.614602, -0.624125, -0.612949, -0.616095},
                                              {-0.273705, -0.510522, -0.643307, -0.739041, -0.830935, -0.860098, -0.885069, -0.956967}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.226216, -0.427422, -0.414774, -0.562412, -0.521969, -0.467143, -0.414713, -0.330372},
                                              {-0.16309, -0.36387, -0.436445, -0.395585, -0.452207, -0.408419, -0.30696, -0.323571},
                                              {-0.106973, -0.382832, -0.361131, -0.402223, -0.452784, -0.400874, -0.342613, -0.185365},
                                              {-0.252347, -0.480454, -0.500724, -0.573321, -0.642967, -0.616602, -0.546648, -0.482116},
                                              {-0.5916, -0.593699, -0.6563, -0.683526, -0.807421, -0.815922, -0.785543, -0.842433}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {0.00677222, 0.0347718, -0.211127, -0.466866, -0.323172, -0.53392, -0.504211, -0.334974},
                                                {0.0935506, 0.0970568, -0.138627, -0.392521, -0.399267, -0.43474, -0.200821, -0.23501},
                                                {0.075394, 0.0609517, -0.170246, -0.409987, -0.420188, -0.448851, -0.267424, -0.313302},
                                                {-0.133011, -0.210294, -0.42431, -0.562203, -0.459603, -0.673718, -0.649959, -0.520375},
                                                {-0.324865, -0.495658, -0.69697, -0.814164, -0.710279, -0.778491, -0.80033, -0.76221}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {0.915632, 0.93365, 0.932587, 0.931425, 0.922551, 0.92571, 0.881836, 0.796746},
                                               {0.906096, 0.925435, 0.924713, 0.919556, 0.906064, 0.911215, 0.866192, 0.773724},
                                               {0.881485, 0.902513, 0.901312, 0.893701, 0.89239, 0.875522, 0.825261, 0.764695},
                                               {0.835562, 0.850935, 0.853431, 0.846889, 0.834667, 0.819543, 0.798447, 0.774576},
                                               {0.809872, 0.807042, 0.80727, 0.799673, 0.785673, 0.757604, 0.74583, 0.726052}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {0.814089, 0.877054, 0.944152, 0.90886, 0.92537, 0.931201, 0.926482, 0.722974},
                                               {0.798753, 0.841944, 0.952407, 0.952185, 0.929863, 0.921988, 0.932327, 0.828574},
                                               {0.733393, 0.821445, 0.903447, 0.928461, 0.917579, 0.921567, 0.923788, 0.727546},
                                               {0.694873, 0.772357, 0.864553, 0.886463, 0.866449, 0.86353, 0.872874, 0.990303},
                                               {0.687564, 0.747428, 0.792902, 0.87729, 0.824179, 0.814195, 0.813793, 0.901867}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {0.758072, 0.796573, 0.838565, 0.942299, 0.990804, 0.914681, 0.900297, 0.872938},
                                                 {0.733353, 0.776217, 0.823004, 0.922272, 1.00522, 0.91665, 0.986521, 1.05648},
                                                 {0.715624, 0.75819, 0.79931, 0.897033, 0.972885, 0.920213, 0.942288, 0.978577},
                                                 {0.691934, 0.724153, 0.7582, 0.793117, 0.879233, 0.832657, 0.8289, 0.932063},
                                                 {0.695097, 0.694242, 0.719957, 0.790272, 0.855222, 0.820038, 0.789541, 0.834146}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else if(run>=295585 && run<=296623 && system==kPbPb5080) { //LHC18q 50-80%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.292480, -0.109690, -0.006546, 0.074381, 0.138770, 0.182105, -0.010704, 0.042879},
                                         {-0.194435, -0.019495, 0.065133, 0.129627, 0.174390, 0.200264, 0.120990, 0.169733},
                                         {-0.183586, -0.051761, 0.009305, 0.057555, 0.103324, 0.095669, -0.010169, 0.124667},
                                         {-0.139560, -0.154218, -0.115706, -0.097739, -0.088791, -0.082103, -0.084627, 0.002918},
                                         {-0.053947, -0.239237, -0.242256, -0.256584, -0.246452, -0.236423, -0.185405, -0.104908}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {0.101589, 0.061254, 0.005293, -0.072268, -0.102444, -0.044265, 0.000680, -0.005517},
                                         {0.217641, 0.100689, 0.041973, -0.013892, -0.025429, 0.043983, 0.052268, 0.169520},
                                         {0.002684, 0.007286, -0.025621, -0.071628, -0.088362, -0.038641, -0.063845, 0.140703},
                                         {0.084500, -0.027425, -0.129124, -0.209249, -0.266129, -0.237653, -0.277646, -0.065473},
                                         {-0.038251, 0.080861, -0.141516, -0.356501, -0.460316, -0.454084, -0.479202, -0.230774}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.082498, -0.009393, 0.173573, 0.114287, 0.095439, 0.026670, 0.024243, 0.117264},
                                           {0.022474, 0.087644, 0.240853, 0.188877, 0.138869, 0.074114, 0.087164, -0.173301},
                                           {0.015226, 0.032181, 0.155128, 0.112671, 0.052385, -0.005417, 0.020613, 0.016407},
                                           {-0.038886, -0.176086, -0.012949, -0.063391, -0.125628, -0.183581, -0.155633, -0.263095},
                                           {-0.177599, -0.352694, -0.201399, -0.120931, -0.255962, -0.357897, -0.371231, -0.385416}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.907559, 0.925039, 0.937333, 0.944514, 0.955643, 0.931166, 0.924014, 0.696171},
                                          {0.897926, 0.915429, 0.934191, 0.943513, 0.940648, 0.920107, 0.789153, 0.766745},
                                          {0.889649, 0.907328, 0.919384, 0.923388, 0.930808, 0.906477, 0.987064, 0.770154},
                                          {0.889654, 0.900043, 0.910578, 0.918840, 0.911155, 0.887032, 0.842598, 0.832887},
                                          {0.896559, 0.890240, 0.904660, 0.913208, 0.904487, 0.890650, 0.836360, 0.780928}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {0.941586, 0.930276, 0.958003, 0.958733, 0.945082, 0.935744, 0.912533, 1.010677},
                                          {0.906601, 0.904342, 0.937911, 0.949574, 0.945744, 0.938570, 0.911436, 0.753891},
                                          {0.988125, 0.884694, 0.924811, 0.940161, 0.935666, 0.918271, 0.884518, 0.811863},
                                          {0.894352, 0.892879, 0.926814, 0.929934, 0.919698, 0.892208, 0.851234, 0.840733},
                                          {0.920926, 0.894999, 0.932580, 0.927593, 0.901108, 0.879334, 0.859266, 0.848501}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.003336, 1.019500, 1.006336, 0.961634, 0.987398, 0.952386, 0.955390, 0.676238},
                                            {0.979906, 1.008984, 0.985871, 0.963956, 0.976405, 0.961994, 0.940296, 1.255870},
                                            {0.977035, 1.025557, 0.996277, 0.959229, 0.978112, 0.946975, 0.937365, 0.712121},
                                            {0.988391, 1.013419, 0.998945, 0.966861, 0.958531, 0.929531, 0.908587, 0.693262},
                                            {0.968214, 1.037859, 1.013852, 0.966745, 0.949422, 0.925684, 0.899639, 0.778341}
                                           };

    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.290461, -0.157362, -0.132967, -0.150766, -0.177387, -0.154474, -0.159506, -0.0819701},
                                              {-0.133763, -0.00498091, 0.0067985, 0.000219383, -0.0387796, -0.0120673, -0.0189252, -0.104004},
                                              {-0.0546859, 0.0146941, 0.012695, -0.0140586, -0.0365883, -0.0455404, -0.0531616, -0.037389},
                                              {-0.0808276, -0.145407, -0.191546, -0.251847, -0.302195, -0.351974, -0.354307, -0.354091},
                                              {-0.093217, -0.288268, -0.403743, -0.500826, -0.589854, -0.641595, -0.693617, -0.700852}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.105712, -0.226045, -0.197616, -0.208206, -0.165878, -0.117968, -0.110323, -0.0885635},
                                              {0.0171305, -0.106295, -0.0665182, 0.00432982, -0.0121926, 0.068422, 0.0699157, 0.0536083},
                                              {0.0668358, -0.115267, -0.049729, 0.00432982, -0.0170053, 0.0328713, 0.0474682, 0.0357199},
                                              {-0.120693, -0.271152, -0.2492, -0.299562, -0.278587, -0.239348, -0.292299, -0.254613},
                                              {-0.485594, -0.44179, -0.427148, -0.531635, -0.522951, -0.552046, -0.679839, -0.666041}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {0.0665573, 0.139883, -0.081605, -0.203813, -0.211302, -0.171301, -0.123612, -0.22991},
                                                {0.178373, 0.233109, 0.0432268, -0.059794, -0.0699837, -0.0425294, 0.00809618, 0.176972},
                                                {0.173038, 0.194914, 0.017688, -0.0766737, -0.081089, -0.0359465, 0.0191725, -0.022404},
                                                {-0.0596051, -0.105099, -0.263882, -0.342295, -0.335415, -0.312026, -0.271335, -0.17117},
                                                {-0.281277, -0.420024, -0.588705, -0.614508, -0.618193, -0.559688, -0.506958, -0.587841}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {0.892435, 0.910563, 0.910262, 0.909437, 0.907557, 0.893655, 0.852062, 0.841},
                                               {0.887739, 0.90651, 0.912569, 0.896137, 0.897385, 0.857293, 0.799653, 0.747155},
                                               {0.861411, 0.890826, 0.889809, 0.891439, 0.87201, 0.847899, 0.830059, 0.813778},
                                               {0.816587, 0.836887, 0.841028, 0.842631, 0.828161, 0.816797, 0.779706, 0.726596},
                                               {0.78294, 0.787287, 0.797276, 0.791521, 0.762976, 0.750954, 0.706446, 0.627991}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {0.820386, 0.862568, 0.913401, 0.941527, 0.945742, 0.949472, 0.891267, 0.85767},
                                               {0.777885, 0.833124, 0.907605, 0.99764, 0.965709, 0.976485, 0.946178, 0.860913},
                                               {0.730192, 0.80115, 0.883422, 0.99764, 0.946244, 0.955063, 0.911931, 0.815714},
                                               {0.689248, 0.772107, 0.846482, 0.87357, 0.888549, 0.89338, 0.837796, 0.810673},
                                               {0.638027, 0.725104, 0.79249, 0.819129, 0.831166, 0.809165, 0.74976, 0.811842}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {0.761836, 0.798676, 0.831566, 0.866324, 0.951735, 0.896439, 0.949077, 0.983275},
                                                 {0.728604, 0.764111, 0.811052, 0.848357, 0.872833, 0.936189, 0.937075, 0.949396},
                                                 {0.729501, 0.747678, 0.796469, 0.839726, 0.916352, 0.925674, 0.894974, 0.85559},
                                                 {0.709612, 0.721265, 0.763333, 0.788841, 0.82415, 0.870806, 0.857862, 0.823096},
                                                 {0.697076, 0.692054, 0.711974, 0.749047, 0.78658, 0.79145, 0.776506, 0.57546}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else if(run>=296690 && run<=297595 && system==kPbPb010) { //LHC18r 0-10%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.191025, -0.069966, -0.043353, -0.110933, -0.153553, -0.151587, -0.206581, -0.333056},
                                         {-0.050894, 0.036801, 0.060649, 0.034148, -0.005891, -0.100589, -0.197833, -0.278303},
                                         {-0.000655, 0.186502, 0.155733, 0.107712, 0.048714, -0.057381, 0.010111, -0.132732},
                                         {0.517092, 0.363080, 0.272265, 0.199289, 0.108668, 0.031334, -0.061000, -0.050126},
                                         {0.870059, 0.519517, 0.369139, 0.001355, 0.200366, 0.137911, 0.050798, -0.001979}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.172679, -0.149360, -0.118343, -0.111048, -0.173474, -0.272699, -0.399108, -0.199539},
                                         {0.016574, -0.057865, -0.047544, -0.017252, -0.082247, -0.188225, -0.204308, -0.197846},
                                         {0.195277, 0.026100, 0.028864, 0.045774, -0.015603, -0.114453, -0.138327, -0.038434},
                                         {0.358194, 0.281201, 0.177425, 0.122200, -0.002712, -0.083754, -0.046955, 0.080940},
                                         {0.387873, 0.555762, 0.398504, 0.193518, 0.002600, -0.068307, -0.028461, 0.080321}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.174141, -0.170543, -0.108674, -0.327194, -0.190786, -0.221476, -0.153456, 0.177733},
                                           {-0.035613, -0.051632, 0.011180, -0.204947, 0.001394, -0.099281, -0.046563, 0.225993},
                                           {-0.015511, -0.026819, 0.072799, -0.107882, 0.079739, -0.037861, 0.023054, 0.308817},
                                           {0.020282, -0.009704, 0.162426, 0.031426, 0.226394, 0.037932, 0.101407, 0.377277},
                                           {0.023680, -0.023035, 0.202086, 0.063541, 0.369817, 0.105155, 0.160985, 0.389877}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {1.146078, 1.173286, 1.166366, 1.145124, 1.137077, 1.196764, 1.240996, 1.225114},
                                          {1.141977, 1.169552, 1.161680, 1.147941, 1.138459, 1.183093, 1.236119, 1.169912},
                                          {1.000606, 1.173940, 1.167339, 1.142365, 1.135951, 1.168155, 0.992925, 1.121113},
                                          {1.156416, 1.159789, 1.160268, 1.126542, 1.125573, 1.141609, 1.153135, 1.110857},
                                          {1.147336, 1.136964, 1.135186, 1.002048, 1.111511, 1.104835, 1.113125, 1.099621}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.022541, 1.193654, 1.172432, 1.151522, 1.113161, 1.047741, 1.002831, 0.996038},
                                          {1.021465, 1.154598, 1.177754, 1.172970, 1.118825, 1.042621, 1.046318, 0.853378},
                                          {0.997829, 1.099294, 1.160628, 1.165639, 1.123121, 1.047782, 1.046315, 0.901296},
                                          {1.011427, 1.080544, 1.157861, 1.153631, 1.114461, 1.038414, 1.046198, 0.954057},
                                          {0.953685, 1.040043, 1.132564, 1.148961, 1.108676, 1.043449, 0.986137, 0.845575}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.050893, 1.075706, 1.108096, 1.265948, 1.181131, 1.090155, 1.114574, 1.275272},
                                            {1.008678, 1.046991, 1.086957, 1.236351, 1.302660, 1.100068, 1.116820, 1.250746},
                                            {0.999061, 1.034440, 1.079335, 1.204056, 1.301164, 1.100452, 1.100579, 1.274839},
                                            {0.996989, 1.020930, 1.068802, 1.171100, 1.293613, 1.106057, 1.104466, 1.257161},
                                            {0.956814, 0.994379, 1.082364, 1.220614, 1.294769, 1.118968, 1.113743, 1.221480}
                                           };

    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.744242, -0.715713, -0.69801, -0.696472, -0.717912, -0.767909, -0.822175, -0.883157},
                                              {-0.976407, -0.964248, -0.976588, -0.969265, -1.00251, -1.04185, -1.08507, -1.02488},
                                              {-0.938573, -1.02253, -1.0532, -1.06874, -1.09608, -1.11066, -1.07855, -1.06274},
                                              {-0.462091, -0.766549, -0.875959, -0.918783, -0.979887, -0.984493, -0.945828, -0.954307},
                                              {0.154123, -0.361271, -0.568491, -0.667592, -0.782836, -0.751772, -0.732903, -0.749}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.468947, -0.636701, -0.601858, -0.806051, -0.94714, -0.842379, -0.955165, -0.898824},
                                              {-0.588647, -0.883708, -0.894757, -1.09769, -1.11786, -1.08056, -1.15336, -1.44054},
                                              {-0.462369, -0.900829, -0.959231, -1.19209, -1.17182, -1.21788, -1.26831, -1.46315},
                                              {-0.28288, -0.585668, -0.942842, -1.13897, -1.18188, -1.1556, -1.20724, -1.06756},
                                              {-0.0830475, -0.129884, -0.40388, -0.905485, -1.03586, -0.963208, -0.95807, -0.591766}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {-0.43448, -0.41261, -0.468653, -0.766399, -0.906529, -0.87423, -0.925983, -0.834281},
                                                {-0.439377, -0.510631, -0.648086, -0.99403, -1.09146, -1.14373, -1.19123, -0.993241},
                                                {-0.400341, -0.54514, -0.649413, -0.979681, -1.22253, -1.27323, -1.27736, -1.12044},
                                                {-0.295184, -0.441872, -0.470702, -0.910766, -1.04581, -1.17824, -1.17277, -0.978326},
                                                {-0.169806, -0.33929, -0.309714, -0.680191, -0.862101, -0.972894, -0.951602, -0.676351}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {1.19971, 1.20244, 1.19018, 1.18674, 1.19735, 1.27442, 1.31886, 1.35234},
                                               {1.17433, 1.18271, 1.16982, 1.17218, 1.17712, 1.26327, 1.33523, 1.30084},
                                               {1.16732, 1.17134, 1.16173, 1.15583, 1.16336, 1.24219, 1.23569, 1.24103},
                                               {1.1773, 1.16117, 1.14767, 1.14647, 1.19338, 1.22358, 1.21504, 1.20365},
                                               {1.21022, 1.17586, 1.16182, 1.15354, 1.15374, 1.22138, 1.19871, 1.20838}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {1.2537, 1.29365, 1.27439, 1.1386, 1.06835, 1.10702, 1.03569, 1.12542},
                                               {1.20352, 1.25335, 1.24745, 1.12229, 1.12846, 1.11836, 1.0746, 1.33407},
                                               {1.17982, 1.21396, 1.23344, 1.13316, 1.15827, 1.10607, 1.06816, 1.14628},
                                               {1.10113, 1.23381, 1.30511, 1.11268, 1.11029, 1.1049, 1.02285, 1.26782},
                                               {1.09653, 1.23731, 1.28769, 1.10684, 1.09593, 1.11015, 1.04836, 1.12603}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {1.13405, 1.18163, 1.24085, 1.27282, 1.12543, 1.0912, 1.04366, 1.10697},
                                                 {1.09801, 1.16854, 1.22617, 1.26278, 1.25383, 1.07191, 1.04581, 1.11811},
                                                 {1.11623, 1.17665, 1.22384, 1.24119, 1.23884, 1.04855, 1.05348, 1.11713},
                                                 {1.08417, 1.16621, 1.22578, 1.34172, 1.24733, 1.04892, 1.05745, 1.14111},
                                                 {1.07421, 1.15563, 1.23173, 1.33796, 1.18478, 1.07793, 1.05178, 1.18215}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else if(run>=296690 && run<=297595 && system==kPbPb1030) { //LHC18r 10-30%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.362832, -0.214416, -0.144585, -0.167250, -0.171063, -0.105322, -0.107191, 0.009823},
                                         {-0.196601, -0.080849, -0.012932, 0.005030, -0.002019, -0.044924, -0.030651, -0.099999},
                                         {-0.023139, 0.050641, 0.061031, 0.051948, 0.031972, -0.035222, -0.020574, -0.079543},
                                         {0.245698, 0.001871, 0.102860, 0.059846, 0.019072, -0.033840, -0.061866, -0.053360},
                                         {0.530210, 0.244812, 0.139886, 0.088187, 0.057971, 0.013332, -0.019454, 0.023011}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.101626, -0.114690, -0.177914, -0.235403, -0.279152, -0.316038, -0.384300, -0.439604},
                                         {0.052146, -0.020212, -0.068919, -0.104394, -0.171097, -0.203185, -0.169125, 0.064661},
                                         {0.162230, 0.022837, -0.023557, -0.052323, -0.104206, -0.151000, -0.189393, -0.070167},
                                         {0.255472, 0.169922, 0.050023, 0.000945, -0.144242, -0.180917, -0.218982, -0.009103},
                                         {0.212816, 0.379694, 0.210836, -0.030621, -0.177924, -0.219872, -0.278034, -0.081553}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.127842, -0.105951, -0.020682, -0.154455, -0.190553, -0.237237, -0.156459, -0.313918},
                                           {0.013450, 0.030051, 0.097971, -0.026394, -0.045325, -0.098226, -0.057618, 0.220656},
                                           {0.003307, 0.002244, 0.114030, -0.011303, -0.021321, -0.061602, 0.016560, -0.007575},
                                           {-0.036168, -0.084552, 0.081779, 0.001038, -0.013646, -0.059595, 0.005110, 0.102638},
                                           {-0.044226, -0.157910, 0.045550, 0.024584, 0.023794, -0.049708, -0.026488, 0.119731}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {1.099109, 1.116216, 1.108694, 1.098980, 1.095777, 1.130881, 1.143756, 1.188649},
                                          {1.096770, 1.115406, 1.111869, 1.098626, 1.090374, 1.112428, 1.113315, 0.882276},
                                          {1.091988, 1.114343, 1.112129, 1.091874, 1.082058, 1.088811, 1.075185, 1.010403},
                                          {1.096639, 1.000070, 1.097670, 1.076073, 1.060905, 1.059529, 1.024134, 0.976708},
                                          {1.088599, 1.080930, 1.079435, 1.057515, 1.038408, 1.028635, 0.999594, 0.920836}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.008989, 1.119319, 1.121537, 1.105223, 1.062584, 1.044995, 1.002856, 1.074909},
                                          {0.997934, 1.084516, 1.111073, 1.108880, 1.077811, 1.044060, 1.041179, 1.126490},
                                          {0.964140, 1.041325, 1.092620, 1.106544, 1.075877, 1.034870, 1.012593, 0.907434},
                                          {0.936910, 1.022693, 1.092146, 0.999002, 1.069088, 1.027403, 0.972605, 0.961660},
                                          {0.936253, 0.987831, 1.069582, 1.089376, 1.057684, 1.027563, 0.947565, 0.820051}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.008565, 1.057382, 1.066843, 1.161353, 1.119564, 1.084088, 1.146888, 0.902038},
                                            {1.000005, 1.034690, 1.048751, 1.143571, 1.147303, 1.093145, 1.119292, 1.101682},
                                            {0.986915, 1.024597, 1.036498, 1.122593, 1.141231, 1.081989, 1.095717, 1.131775},
                                            {0.980491, 1.022194, 1.040435, 0.999788, 1.142627, 1.080839, 1.102095, 1.186842},
                                            {0.973423, 0.994364, 1.037931, 1.098547, 1.142235, 1.092295, 1.088411, 1.181632}
                                           };

    // using 0-10% for 10-30% in case of pass1
    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.744242, -0.715713, -0.69801, -0.696472, -0.717912, -0.767909, -0.822175, -0.883157},
                                              {-0.976407, -0.964248, -0.976588, -0.969265, -1.00251, -1.04185, -1.08507, -1.02488},
                                              {-0.938573, -1.02253, -1.0532, -1.06874, -1.09608, -1.11066, -1.07855, -1.06274},
                                              {-0.462091, -0.766549, -0.875959, -0.918783, -0.979887, -0.984493, -0.945828, -0.954307},
                                              {0.154123, -0.361271, -0.568491, -0.667592, -0.782836, -0.751772, -0.732903, -0.749}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.468947, -0.636701, -0.601858, -0.806051, -0.94714, -0.842379, -0.955165, -0.898824},
                                              {-0.588647, -0.883708, -0.894757, -1.09769, -1.11786, -1.08056, -1.15336, -1.44054},
                                              {-0.462369, -0.900829, -0.959231, -1.19209, -1.17182, -1.21788, -1.26831, -1.46315},
                                              {-0.28288, -0.585668, -0.942842, -1.13897, -1.18188, -1.1556, -1.20724, -1.06756},
                                              {-0.0830475, -0.129884, -0.40388, -0.905485, -1.03586, -0.963208, -0.95807, -0.591766}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {-0.43448, -0.41261, -0.468653, -0.766399, -0.906529, -0.87423, -0.925983, -0.834281},
                                                {-0.439377, -0.510631, -0.648086, -0.99403, -1.09146, -1.14373, -1.19123, -0.993241},
                                                {-0.400341, -0.54514, -0.649413, -0.979681, -1.22253, -1.27323, -1.27736, -1.12044},
                                                {-0.295184, -0.441872, -0.470702, -0.910766, -1.04581, -1.17824, -1.17277, -0.978326},
                                                {-0.169806, -0.33929, -0.309714, -0.680191, -0.862101, -0.972894, -0.951602, -0.676351}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {1.19971, 1.20244, 1.19018, 1.18674, 1.19735, 1.27442, 1.31886, 1.35234},
                                               {1.17433, 1.18271, 1.16982, 1.17218, 1.17712, 1.26327, 1.33523, 1.30084},
                                               {1.16732, 1.17134, 1.16173, 1.15583, 1.16336, 1.24219, 1.23569, 1.24103},
                                               {1.1773, 1.16117, 1.14767, 1.14647, 1.19338, 1.22358, 1.21504, 1.20365},
                                               {1.21022, 1.17586, 1.16182, 1.15354, 1.15374, 1.22138, 1.19871, 1.20838}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {1.2537, 1.29365, 1.27439, 1.1386, 1.06835, 1.10702, 1.03569, 1.12542},
                                               {1.20352, 1.25335, 1.24745, 1.12229, 1.12846, 1.11836, 1.0746, 1.33407},
                                               {1.17982, 1.21396, 1.23344, 1.13316, 1.15827, 1.10607, 1.06816, 1.14628},
                                               {1.10113, 1.23381, 1.30511, 1.11268, 1.11029, 1.1049, 1.02285, 1.26782},
                                               {1.09653, 1.23731, 1.28769, 1.10684, 1.09593, 1.11015, 1.04836, 1.12603}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {1.13405, 1.18163, 1.24085, 1.27282, 1.12543, 1.0912, 1.04366, 1.10697},
                                                 {1.09801, 1.16854, 1.22617, 1.26278, 1.25383, 1.07191, 1.04581, 1.11811},
                                                 {1.11623, 1.17665, 1.22384, 1.24119, 1.23884, 1.04855, 1.05348, 1.11713},
                                                 {1.08417, 1.16621, 1.22578, 1.34172, 1.24733, 1.04892, 1.05745, 1.14111},
                                                 {1.07421, 1.15563, 1.23173, 1.33796, 1.18478, 1.07793, 1.05178, 1.18215}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else if(run>=296690 && run<=297595 && system==kPbPb3050) { //LHC18r 30-50%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.339214, -0.160684, -0.060256, -0.060043, -0.039870, 0.053763, 0.115945, 0.068671},
                                         {-0.153487, -0.003281, 0.000635, 0.129705, 0.151439, 0.141144, 0.162775, 0.122433},
                                         {-0.030722, 0.082842, 0.120633, 0.002293, 0.129770, 0.095171, 0.124722, 0.097544},
                                         {0.108108, 0.064697, 0.055145, 0.039384, 0.017393, -0.013635, -0.027335, -0.001984},
                                         {0.289296, 0.065340, -0.001456, -0.027409, -0.036075, -0.058948, -0.071766, -0.090492}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.007066, 0.007203, -0.072174, -0.147748, -0.219257, -0.207882, -0.145396, -0.226397},
                                         {0.166493, 0.110886, 0.051342, -0.000800, -0.072489, -0.056804, -0.002682, 0.158805},
                                         {0.198191, 0.065111, 0.039822, 0.016201, -0.031021, -0.035406, -0.052571, 0.152084},
                                         {0.192903, 0.119434, 0.017370, -0.077958, -0.156639, -0.153978, -0.210036, 0.021977},
                                         {0.034760, 0.237741, 0.082392, -0.142398, -0.277304, -0.301306, -0.332630, -0.102745}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.075352, -0.037625, 0.076878, 0.025434, -0.038192, -0.091696, -0.090785, -0.130279},
                                           {0.069977, 0.110108, 0.214011, 0.160622, 0.111492, 0.072403, 0.085344, 0.043529},
                                           {0.034632, 0.035115, 0.165312, 0.104921, 0.078413, 0.063851, 0.117921, 0.120917},
                                           {-0.039762, -0.122970, 0.040104, 0.036571, -0.012770, -0.037092, -0.016714, 0.037006},
                                           {-0.105142, -0.249199, -0.073606, -0.011592, -0.053712, -0.140094, -0.150927, -0.067526}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {1.054076, 1.072062, 1.075482, 1.067662, 1.061529, 1.071348, 1.069763, 0.977920},
                                          {1.048322, 1.072738, 1.000194, 1.066192, 1.057138, 1.058840, 1.048689, 1.073029},
                                          {1.038189, 1.064392, 1.063722, 1.288682, 1.039299, 1.029222, 1.015398, 0.970869},
                                          {1.029756, 1.043109, 1.046094, 1.034296, 1.018694, 1.002364, 0.962781, 0.920044},
                                          {1.024971, 1.024221, 1.024885, 1.012428, 0.991032, 0.965466, 0.914833, 0.869332}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {0.973110, 1.064312, 1.078427, 1.074650, 1.044897, 1.033420, 1.026433, 0.949119},
                                          {0.955018, 1.026558, 1.060793, 1.073535, 1.053483, 1.038486, 1.008743, 0.876533},
                                          {0.942043, 0.987616, 1.038474, 1.072654, 1.047256, 1.028571, 0.991543, 0.879501},
                                          {0.914769, 0.964490, 1.034807, 1.056429, 1.030408, 1.013987, 0.945930, 0.879947},
                                          {0.911786, 0.941459, 1.015171, 1.040208, 1.021548, 0.989027, 0.928063, 0.810423}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.019573, 1.041849, 1.038164, 1.099966, 1.110934, 1.087189, 1.107509, 1.163927},
                                            {0.994682, 1.020569, 1.025007, 1.086028, 1.120928, 1.092608, 1.115262, 1.101274},
                                            {0.982203, 1.021021, 1.011548, 1.057053, 1.088326, 1.074358, 1.086468, 1.142056},
                                            {0.977964, 1.017725, 1.016399, 1.031597, 1.075654, 1.061813, 1.076272, 1.106065},
                                            {0.969742, 0.999319, 1.018558, 1.011554, 1.054023, 1.057564, 1.054873, 1.070879}
                                           };

    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.643992, -0.547379, -0.477054, -0.440121, -0.426498, -0.406638, -0.361916, -0.311107},
                                              {-0.714799, -0.624286, -0.588865, -0.554541, -0.539189, -0.50596, -0.442875, -0.504426},
                                              {-0.641763, -0.623418, -0.607051, -0.580924, -0.573627, -0.534728, -0.520039, -0.457556},
                                              {-0.292997, -0.446284, -0.484493, -0.499929, -0.516507, -0.485561, -0.428434, -0.387974},
                                              {0.151661, -0.179137, -0.303383, -0.371113, -0.418849, -0.413089, -0.380399, -0.442395}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.180292, -0.327542, -0.451762, -0.598542, -0.655667, -0.629728, -0.562158, -0.399934},
                                              {-0.247604, -0.448958, -0.522777, -0.71607, -0.737935, -0.683963, -0.577149, -0.59199},
                                              {-0.189726, -0.48024, -0.575707, -0.742612, -0.765291, -0.709371, -0.607647, -0.28355},
                                              {0.00208301, -0.258596, -0.444814, -0.667385, -0.713196, -0.654455, -0.540685, -0.339524},
                                              {-0.000405731, 0.0385342, -0.162143, -0.523791, -0.612858, -0.553858, -0.45908, -0.401148}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {-0.149563, -0.156293, -0.17301, -0.377551, -0.531683, -0.617193, -0.576007, -0.621884},
                                                {-0.14944, -0.204765, -0.255216, -0.438923, -0.675954, -0.637264, -0.640382, -0.585865},
                                                {-0.149585, -0.263429, -0.274904, -0.428442, -0.50586, -0.734513, -0.659482, -0.473013},
                                                {-0.13168, -0.278353, -0.222243, -0.405964, -0.603998, -0.612839, -0.631597, -0.651136},
                                                {-0.108823, -0.280568, -0.193149, -0.250911, -0.53358, -0.55575, -0.545469, -0.421828}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {1.08474, 1.09912, 1.10605, 1.11466, 1.12312, 1.13388, 1.11915, 0.98612},
                                               {1.0785, 1.09545, 1.10056, 1.10451, 1.11234, 1.12323, 1.10217, 1.14173},
                                               {1.0727, 1.09036, 1.09512, 1.10093, 1.10376, 1.10297, 1.1321, 1.06425},
                                               {1.07348, 1.08058, 1.08902, 1.08921, 1.08795, 1.08588, 1.05973, 1.01606},
                                               {1.08087, 1.06969, 1.08016, 1.08215, 1.08229, 1.07935, 1.04663, 1.02234}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {1.10082, 1.14338, 1.17166, 1.08815, 1.07553, 1.05804, 0.98676, 1.09437},
                                               {1.06003, 1.09105, 1.13176, 1.08161, 1.08362, 1.0915, 1.0352, 1.11064},
                                               {1.04461, 1.08936, 1.14854, 1.08121, 1.09336, 1.09623, 1.04665, 0.882711},
                                               {1.01912, 1.09634, 1.15565, 1.08002, 1.09909, 1.09498, 1.05249, 0.999984},
                                               {1.0683, 1.09023, 1.13342, 1.07979, 1.0902, 1.09837, 1.0282, 1.07592}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {1.06529, 1.10828, 1.14221, 1.18478, 1.10225, 1.07407, 1.08529, 1.14558},
                                                 {1.05081, 1.10053, 1.1237, 1.15429, 1.04817, 1.09967, 1.10276, 1.21743},
                                                 {1.06485, 1.11335, 1.13392, 1.12825, 1.12658, 1.08174, 1.10395, 1.2337},
                                                 {1.05451, 1.11328, 1.13998, 1.17462, 1.19486, 1.09464, 1.07388, 1.07995},
                                                 {1.05442, 1.09606, 1.14899, 1.14143, 1.04405, 1.08258, 1.05422, 1.13655}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else if(run>=296690 && run<=297595 && system==kPbPb5080) { //LHC18r 50-80%
    nPbins = 8;
    vector<Float_t> pTPClims = {0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0};
    nEtabins = 5;
    vector<Float_t> absetalims = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8};

    vector<vector<Float_t> > meanPion = {
                                         {-0.242868, -0.034709, 0.086217, 0.098705, 0.120551, -0.003268, 0.317127, 0.208488},
                                         {-0.052033, 0.129527, 0.245216, 0.301874, 0.310191, 0.330253, 0.390647, 0.530561},
                                         {0.028418, 0.178272, 0.232480, 0.250878, 0.256587, 0.221487, 0.246475, 0.239371},
                                         {0.078943, 0.069304, 0.082843, 0.089406, 0.078423, 0.052210, 0.033915, -0.082944},
                                         {0.188599, -0.001766, -0.050489, -0.046980, -0.047724, -0.065067, -0.060780, -0.102816}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {0.099605, 0.124865, 0.055650, -0.015881, -0.056013, 0.010561, 0.059939, 0.286670},
                                         {0.260183, 0.213671, 0.172762, 0.149620, 0.098456, 0.165291, 0.272468, 0.396365},
                                         {0.254415, 0.116044, 0.123464, 0.128821, 0.109343, 0.141212, 0.204949, 0.378121},
                                         {0.191618, 0.108656, 0.027611, -0.032694, -0.077172, -0.070485, -0.010624, 0.221583},
                                         {0.017739, 0.174228, 0.041171, -0.167471, -0.275308, -0.295296, -0.296820, -0.015189}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.048517, 0.038067, 0.164801, 0.188855, 0.123180, 0.051036, 0.028055, 0.230222},
                                           {0.139039, 0.173595, 0.296026, 0.318409, 0.280748, 0.228276, 0.262635, 0.724960},
                                           {0.069845, 0.077207, 0.230373, 0.231992, 0.213941, 0.207665, 0.240780, 0.420941},
                                           {-0.054984, -0.121467, 0.043724, 0.083089, 0.041915, 0.013356, 0.039466, 0.071266},
                                           {-0.132239, -0.292110, -0.129283, -0.016717, -0.054224, -0.159294, -0.191175, -0.064385}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {1.018834, 1.043857, 1.043243, 1.049907, 1.039843, 1.000559, 0.956136, 0.696171},
                                          {1.013415, 1.041146, 1.044810, 1.035901, 1.011686, 1.009264, 0.922757, 0.696171},
                                          {1.001385, 1.036315, 1.037534, 1.034767, 1.013350, 0.991444, 0.949189, 0.746675},
                                          {0.991827, 1.016293, 1.016864, 1.011489, 0.997962, 0.968234, 0.900729, 0.810880},
                                          {0.992146, 0.989005, 0.998556, 0.991347, 0.964742, 0.929473, 0.864147, 0.711892}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {0.966555, 1.016449, 1.049326, 1.055192, 1.046908, 1.031533, 1.011952, 0.979245},
                                          {0.899451, 0.982152, 1.025515, 1.061450, 1.059582, 1.034840, 1.010883, 0.891627},
                                          {0.921892, 0.936196, 1.008438, 1.053296, 1.051750, 1.014570, 1.001926, 0.849365},
                                          {0.901560, 0.920688, 0.999765, 1.034884, 1.028082, 0.988920, 0.981225, 0.896881},
                                          {0.885412, 0.892530, 0.980230, 1.022434, 1.010559, 0.964885, 0.921518, 0.858829}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {0.987262, 1.030739, 1.014339, 1.048503, 1.062421, 1.078677, 1.052815, 1.255870},
                                            {0.976025, 1.000016, 1.006116, 1.032209, 1.072346, 1.064854, 1.084345, 0.676238},
                                            {0.939797, 1.002774, 0.997320, 1.009033, 1.063993, 1.061939, 1.041126, 0.820689},
                                            {0.964952, 1.023224, 1.002267, 0.989680, 1.030757, 1.029380, 1.037391, 1.171224},
                                            {0.952292, 0.999155, 1.012708, 0.967423, 1.023101, 1.015831, 1.004669, 0.830949}
                                           };

    vector<vector<Float_t> > meanPionPass1 = {
                                              {-0.388593, -0.259655, -0.188052, -0.148457, -0.137755, -0.103713, -0.0103896, 0.0183914},
                                              {-0.360916, -0.239814, -0.206334, -0.151416, -0.152129, -0.117832, -0.0512912, -0.0440523},
                                              {-0.275863, -0.225395, -0.194411, -0.162487, -0.166327, -0.137437, -0.0628822, 0.00493976},
                                              {0.0136961, -0.0937519, -0.123901, -0.12913, -0.147936, -0.131946, -0.0580969, -0.140029},
                                              {0.383871, 0.108436, 0.00656927, -0.0673751, -0.108841, -0.136234, -0.0603665, -0.118457}
                                             };
    vector<vector<Float_t> > meanKaonPass1 = {
                                              {-0.0164804, -0.0972051, -0.117233, -0.268106, -0.286532, -0.192977, -0.13506, -0.131241},
                                              {-0.0437434, -0.144134, -0.1716, -0.282154, -0.252251, -0.229581, -0.118788, -0.05944},
                                              {0.0694593, -0.139885, -0.151188, -0.290645, -0.290291, -0.271112, -0.10308, -0.0566368},
                                              {0.183337, 0.0309278, -0.0603532, -0.246455, -0.269791, -0.28694, -0.073309, 0.0979597},
                                              {0.21312, 0.264178, 0.120222, -0.167352, -0.219633, -0.253195, -0.091609, 0.083244}
                                             };
    vector<vector<Float_t> > meanProtonPass1 = {
                                                {-0.00831157, -0.0246502, -0.00931054, -0.059935, -0.186943, -0.197498, -0.307994, -0.0702371},
                                                {0.0332978, -0.0489401, -0.03383, -0.0776816, -0.178209, -0.222742, -0.271989, -0.520199},
                                                {-0.0323603, -0.0684328, -0.0158219, -0.073291, -0.180924, -0.258898, -0.308916, -0.418997},
                                                {-0.0288743, -0.135427, 0.00271318, -0.0413093, -0.181906, -0.220064, -0.259092, -0.174753},
                                                {-0.00748328, -0.151846, -0.0157611, 0.0119542, -0.0992955, -0.196113, -0.216777, -0.149337}
                                               };
    vector<vector<Float_t> > sigmaPionPass1 = {
                                               {1.05503, 1.0735, 1.08331, 1.09315, 1.09289, 1.09399, 1.07299, 1.19263},
                                               {1.05972, 1.08279, 1.09055, 1.08779, 1.095, 1.09622, -1.02437, 1.01151},
                                               {1.05442, 1.08032, 1.08152, 1.08671, 1.08594, 1.0849, 1.09464, 1.00198},
                                               {1.05033, 1.06876, 1.07743, 1.08634, 1.08502, 1.09015, 1.03311, 1.02371},
                                               {1.05402, 1.05733, 1.07034, 1.06867, 1.07451, 1.07605, 1.02273, 1.0096}
                                              };
    vector<vector<Float_t> > sigmaKaonPass1 = {
                                               {1.08457, 1.1115, 1.10762, 1.11613, 1.10712, 1.11228, 1.07115, 1.00753},
                                               {1.0578, 1.08034, 1.09357, 1.11391, 1.12987, 1.11409, 1.08436, 1.10904},
                                               {1.03893, 1.06877, 1.10364, 1.11311, 1.13444, 1.10691, 1.09994, 1.07071},
                                               {1.01472, 1.07273, 1.10807, 1.10768, 1.13076, 1.09433, 1.08405, 0.96925},
                                               {1.09291, 1.07712, 1.1035, 1.09887, 1.12159, 1.07188, 1.05105, 0.956145}
                                              };
    vector<vector<Float_t> > sigmaProtonPass1 = {
                                                 {1.04921, 1.10968, 1.12128, 1.10103, 1.08091, 1.09558, 1.09125, 0.861456},
                                                 {1.08841, 1.0943, 1.12142, 1.10862, 1.08891, 1.09479, 1.05673, 1.04939},
                                                 {1.11483, 1.10683, 1.11962, 1.09238, 1.07412, 1.10004, 1.08878, 1.18803},
                                                 {1.0769, 1.11629, 1.13705, 1.10901, 1.09598, 1.09247, 1.04659, 1.25114},
                                                 {1.06329, 1.10681, 1.15275, 1.11122, 1.0892, 1.06841, 1.06915, 1.09393}
                                                };

    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);
    
    if(!isPass1) {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPion[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaon[iEta];
        meanNsigmaTPCproton[iEta] = meanProton[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPion[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaon[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProton[iEta];
      }
    }
    else {
      for(int iEta=0; iEta<nEtabins; iEta++) {
        meanNsigmaTPCpion[iEta] = meanPionPass1[iEta];
        meanNsigmaTPCkaon[iEta] = meanKaonPass1[iEta];
        meanNsigmaTPCproton[iEta] = meanProtonPass1[iEta];
        sigmaNsigmaTPCpion[iEta] = sigmaPionPass1[iEta];
        sigmaNsigmaTPCkaon[iEta] = sigmaKaonPass1[iEta];
        sigmaNsigmaTPCproton[iEta] = sigmaProtonPass1[iEta];
      }
    }
  }
  else { //default: no correction applied
    nPbins = 1;
    vector<Float_t> pTPClims = {0.,1000.};  
    nEtabins = 1;
    vector<Float_t> absetalims = {0.,2.};  

    vector<Float_t> meanPion = {0.};
    vector<Float_t> meanKaon = {0.};
    vector<Float_t> meanProton = {0.};
    vector<Float_t> sigmaPion = {1.};
    vector<Float_t> sigmaKaon = {1.};
    vector<Float_t> sigmaProton = {1.};
    std::copy(pTPClims.begin(),pTPClims.end(),Plims);
    std::copy(absetalims.begin(),absetalims.end(),absEtalims);

    meanNsigmaTPCpion[0] = meanPion;
    meanNsigmaTPCkaon[0] = meanKaon;
    meanNsigmaTPCproton[0] = meanProton;
    sigmaNsigmaTPCpion[0] = sigmaPion;
    sigmaNsigmaTPCkaon[0] = sigmaKaon;
    sigmaNsigmaTPCproton[0] = sigmaProton;
  }
}