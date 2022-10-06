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
                                         {-0.184653, 0.361796, -0.063490, -0.084251, -0.153846, -0.206109, -0.343694, -0.387847},
                                         {-0.095127, 0.000115, -0.001310, -0.024197, -0.089811, -0.125071, -0.276731, -0.285040},
                                         {0.044578, 0.090402, 0.056943, 0.027620, -0.043285, -0.074859, 0.007401, -0.192704},
                                         {0.355984, 0.235860, 0.170052, 0.102762, 0.026966, -0.049388, -0.122850, -0.174577},
                                         {0.789917, 0.500074, 0.386115, 0.281008, 0.211180, 0.120366, 0.085009, 0.007671}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.147985, -0.342459, -0.251215, -0.095976, -0.104727, -0.083265, -0.099515, 0.034310},
                                         {0.054144, -0.230234, -0.192788, -0.050954, -0.079672, -0.038875, -0.080497, -0.087566},
                                         {0.134583, -0.134194, -0.116771, -0.009438, -0.077332, -0.007457, -0.071729, -0.139253},
                                         {0.345745, 0.057578, -0.003545, 0.057970, -0.032775, 0.014880, -0.050265, -0.077231},
                                         {0.481903, 0.439627, 0.303858, 0.234938, 0.133452, 0.161088, 0.050665, 0.076908}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.482615, -0.278645, -0.173344, -0.273858, -0.214752, -0.126532, 0.010897, 0.173967},
                                           {-0.250260, -0.095922, -0.027161, -0.174421, -0.182634, -0.071802, 0.042606, 0.371931},
                                           {-0.104595, -0.060710, 0.046608, -0.071916, -0.105739, -0.033560, 0.081187, 0.330923},
                                           {-0.214155, -0.047809, 0.095179, 0.006801, 0.021894, 0.037351, 0.139837, 0.310798},
                                           {-0.098777, 0.035917, 0.213733, 0.187997, 0.194436, 0.171149, 0.304858, 0.481044}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.970034, 0.831474, 0.980311, 0.987594, 1.032489, 1.047586, 1.110746, 1.079898},
                                          {0.959897, 0.967210, 0.981411, 0.988660, 1.033790, 1.028895, 1.093639, 1.070336},
                                          {0.943592, 0.957296, 0.969501, 0.979750, 1.008689, 1.013909, 1.011310, 1.016736},
                                          {0.947637, 0.948060, 0.960019, 0.968100, 0.985978, 1.014542, 1.029408, 1.008752},
                                          {0.981082, 0.966207, 0.982094, 0.991465, 1.005105, 1.042266, 1.055854, 1.025440}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.013581, 1.106096, 1.063785, 1.020895, 0.986511, 0.934191, 0.953692, 1.024298},
                                          {0.928656, 1.062147, 1.070935, 1.008039, 1.004113, 0.937559, 0.955645, 0.880607},
                                          {1.038147, 1.047927, 1.053894, 1.002405, 0.987119, 0.937611, 0.932437, 0.800550},
                                          {0.871896, 1.036853, 1.042414, 1.001309, 0.970306, 0.943127, 0.913817, 0.800550},
                                          {0.839236, 1.014803, 1.074403, 1.028232, 1.011712, 0.989421, 0.895366, 0.800550}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.188218, 1.190626, 1.196265, 1.175971, 1.042160, 0.971461, 0.976100, 1.171911},
                                            {1.188218, 1.190626, 1.179793, 1.163684, 1.032194, 0.975340, 0.969536, 1.187465},
                                            {1.188218, 1.190626, 1.151755, 1.138199, 1.046587, 0.982612, 0.957938, 1.187465},
                                            {1.188218, 1.167721, 1.128573, 1.135381, 1.126784, 0.956273, 0.961016, 1.119599},
                                            {1.141203, 1.162393, 1.150685, 1.180325, 1.108841, 1.039609, 1.022764, 1.130694}
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
                                         {-0.283370, -0.121674, -0.075582, -0.057674, -0.069150, -0.091755, -0.134776, -0.168521},
                                         {-0.167996, -0.033431, 0.009880, 0.027641, 0.031169, 0.012482, -0.049401, 0.015069},
                                         {-0.035972, 0.047140, 0.050941, 0.047186, 0.038553, 0.023658, -0.020481, -0.033686},
                                         {0.211666, 0.135312, 0.097353, 0.061209, 0.031138, -0.008842, -0.058618, -0.056704},
                                         {0.603049, 0.359194, 0.270759, 0.196086, 0.145851, 0.115214, -0.004237, 0.039647}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {0.023791, -0.223672, -0.216660, -0.126673, -0.148604, -0.146778, -0.153914, -0.093507},
                                         {0.254925, -0.123895, -0.151308, -0.066765, -0.091561, -0.072130, -0.091512, -0.331942},
                                         {0.317901, -0.077688, -0.096635, -0.024552, -0.054298, -0.027777, -0.111636, -0.089626},
                                         {0.428857, 0.035461, -0.033212, -0.011295, -0.063673, -0.036686, -0.133811, -0.068885},
                                         {0.269039, 0.316576, 0.208976, 0.130007, 0.056365, 0.032858, -0.031329, 0.049233}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.193946, -0.041774, 0.023529, -0.071472, -0.105326, -0.121311, -0.040184, -0.087833},
                                           {0.036593, 0.120100, 0.183316, 0.051970, -0.016761, -0.035116, 0.074296, 0.033127},
                                           {0.035297, 0.087053, 0.181420, 0.077399, 0.007734, -0.018591, 0.085435, 0.144491},
                                           {-0.049428, -0.011938, 0.117597, 0.050936, 0.011910, -0.010469, 0.058240, 0.118302},
                                           {-176.257718, -0.015425, 0.164762, 0.195750, 0.172665, 0.116777, 0.179934, 0.315194}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.985292, 0.987189, 0.994476, 1.007433, 1.028907, 1.043367, 1.057741, 0.917698},
                                          {0.970163, 0.976811, 0.984904, 1.001540, 1.014640, 1.025355, 1.036020, 0.858593},
                                          {0.960742, 0.963977, 0.975611, 0.986411, 0.999715, 1.012931, 0.995695, 0.919031},
                                          {0.955001, 0.951673, 0.961504, 0.971534, 0.986773, 0.993325, 0.980649, 0.903424},
                                          {0.980420, 0.962716, 0.981773, 0.990081, 0.993558, 1.002033, 0.995097, 0.940409}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.050415, 1.122719, 1.095929, 1.023996, 0.996403, 0.959957, 0.930989, 0.819580},
                                          {1.042299, 1.087111, 1.060323, 1.027561, 0.990205, 0.957006, 0.916065, 0.724410},
                                          {1.013682, 1.047452, 1.052988, 1.010398, 0.979929, 0.962425, 0.895738, 0.851427},
                                          {0.913013, 1.024961, 1.043468, 0.996197, 0.980505, 0.965785, 0.877520, 0.787856},
                                          {0.906431, 1.030372, 1.056423, 1.021965, 1.014340, 0.976730, 0.893161, 0.796926}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.299133, 1.256146, 1.208929, 1.166168, 1.067138, 1.011672, 1.008516, 0.973606},
                                            {1.253432, 1.231618, 1.175634, 1.150024, 1.078027, 1.010735, 1.017405, 0.856500},
                                            {1.232899, 1.218378, 1.145956, 1.128872, 1.050302, 0.991076, 0.988346, 0.969494},
                                            {1.247166, 1.183984, 1.128266, 1.108627, 1.032108, 0.989177, 0.976007, 0.959209},
                                            {1.176895, 1.176780, 1.134964, 1.147337, 1.059227, 1.030374, 1.020561, 1.039563}
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
                                         {-0.258026, -0.098174, -0.053314, -0.037314, 0.002447, -0.050709, -0.066735, -0.039035},
                                         {-0.156067, -0.027153, 0.008964, 0.038060, 0.051093, 0.040771, 0.025998, -0.133829},
                                         {-0.030731, 0.049392, 0.051134, 0.063858, 0.071578, 0.039253, 0.035325, -0.001394},
                                         {0.194090, 0.115733, 0.082588, 0.050550, 0.023395, 0.002422, -0.030163, -0.096925},
                                         {0.558235, 0.314788, 0.235057, 0.159226, 0.118855, 0.103936, 0.108531, 0.067412}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {0.062933, -0.175306, -0.194058, -0.135872, -0.122679, -0.144476, -0.117156, -0.034769},
                                         {0.089457, -0.076151, -0.135359, -0.077355, -0.078607, -0.037006, -0.101786, 0.457094},
                                         {0.330708, -0.074933, -0.080090, -0.030227, -0.064564, -0.031211, -0.033114, 0.023632},
                                         {0.350146, 0.059151, -0.042610, -0.006413, -0.082599, -0.080071, -0.106476, -0.013461},
                                         {0.159045, 0.288349, 0.210506, 0.100822, 0.017286, -0.010636, -0.060329, 0.035565}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.197752, -0.062886, 0.051704, 0.001155, -0.047688, -0.066480, -0.009915, -0.024163},
                                           {0.022707, 0.153074, 0.204196, 0.111217, 0.015884, 0.014938, 0.051181, -0.301615},
                                           {0.063096, 0.135891, 0.208465, 0.118118, 0.032716, 0.018858, 0.046326, 0.124939},
                                           {0.004268, 0.025472, 0.125621, 0.083780, 0.016928, -0.012725, 0.041209, 0.138594},
                                           {0.043650, -0.004038, 0.128990, 0.207122, 0.145365, 0.066872, 0.111628, 0.188975}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {1.004241, 1.005783, 1.015425, 1.025364, 1.000384, 1.037614, 1.007406, 1.294102},
                                          {0.986988, 0.981297, 0.993280, 1.013018, 1.034257, 1.019556, 1.042623, 1.294102},
                                          {0.974453, 0.974747, 0.988831, 0.999763, 1.019304, 1.000209, 0.970286, 0.976079},
                                          {0.963036, 0.962815, 0.974899, 0.979256, 0.989649, 0.992986, 0.950496, 0.852717},
                                          {0.992193, 0.976589, 0.986684, 0.998628, 1.001386, 1.004021, 0.962268, 0.890092}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {0.989440, 1.101694, 1.087020, 1.045485, 0.982930, 0.940821, 0.960965, 0.852444},
                                          {0.938504, 1.051016, 1.054329, 1.021777, 0.996675, 0.976912, 0.933042, 1.208222},
                                          {0.991977, 1.045516, 1.057917, 1.025345, 0.999334, 0.972976, 0.927416, 0.792594},
                                          {1.022000, 1.014262, 1.043935, 1.029136, 1.002343, 0.975134, 0.916799, 0.835378},
                                          {1.287251, 1.013271, 1.043850, 1.037127, 1.012268, 0.976532, 0.898373, 0.769351}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.291343, 1.244125, 1.209578, 1.185410, 1.128722, 1.043337, 1.064010, 1.282262},
                                            {1.244684, 1.243539, 1.196232, 1.140684, 1.084718, 1.037788, 1.002907, 1.282442},
                                            {1.282926, 1.216053, 1.151501, 1.127935, 1.061504, 1.019004, 0.983985, 1.190641},
                                            {0.993341, 1.186000, 1.138557, 1.100977, 1.049063, 1.023188, 0.978208, 1.160878},
                                            {1.237456, 1.203790, 1.145391, 1.132460, 1.060131, 1.053141, 1.005814, 1.023906}
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
                                         {-0.144256, 0.084686, 0.177615, 0.238413, 0.294330, 0.311905, 0.349956, 0.207129},
                                         {-0.012614, 0.001473, 0.278050, 0.335107, 0.008866, 0.403851, 0.329365, 0.350922},
                                         {0.018564, 0.174294, 0.231131, 0.265382, 0.305937, 0.318235, 0.309906, 0.329929},
                                         {0.075079, 0.076743, 0.098365, 0.092403, 0.086010, 0.111934, 0.079393, 0.102709},
                                         {0.304347, 0.136284, 0.105478, 0.070984, 0.055943, 0.063587, 0.054957, 0.111635}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {0.338304, 0.121787, 0.047363, 0.113734, 0.078718, 0.182590, 0.204966, 0.201078},
                                         {0.534398, 0.205224, 0.143748, 0.209708, 0.231642, 0.277948, 0.334410, 0.344802},
                                         {0.490480, 0.100758, 0.078271, 0.169979, 0.174474, 0.207179, 0.256516, 0.210144},
                                         {0.378377, 0.012689, -0.039377, 0.018459, -0.007468, 0.034984, 0.089851, 0.107165},
                                         {0.257277, 0.177077, 0.066469, 0.002672, -0.063495, -0.079059, 0.058076, 0.105035}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.001437, 0.288647, 0.363772, 0.331324, 0.265002, 0.212276, 0.249435, 0.146044},
                                           {0.362644, 0.480733, 0.535720, 0.466559, 0.382644, 0.293615, 0.314685, 0.007801},
                                           {0.297508, 0.344544, 0.404930, 0.377014, 0.294349, 0.231484, 0.262850, 0.217701},
                                           {0.111955, 0.094911, 0.150394, 0.172129, 0.098286, 0.032704, 0.073266, -0.069029},
                                           {0.059628, -0.036833, 0.055879, 0.162304, 0.097287, -0.008325, 0.011429, -0.023460}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.960722, 0.976237, 0.995389, 1.011091, 1.019328, 0.993530, 0.923437, 0.717719},
                                          {0.947595, 0.997167, 0.983090, 1.000414, 0.992189, 0.996424, 0.890520, 0.967983},
                                          {0.934720, 0.947736, 0.967144, 0.986093, 0.983782, 0.970657, 0.944202, 0.849360},
                                          {0.916765, 0.927270, 0.947879, 0.959659, 0.957695, 0.939365, 0.880719, 0.879680},
                                          {0.945606, 0.934747, 0.954892, 0.969119, 0.960052, 0.956078, 0.902627, 0.804151}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.071506, 1.086628, 1.067524, 1.050647, 1.026026, 1.015115, 0.949265, 0.887785},
                                          {1.033074, 1.049844, 1.052376, 1.050813, 1.039503, 1.005878, 0.946087, 0.897053},
                                          {1.019604, 1.012808, 1.033489, 1.023531, 1.026400, 0.984494, 0.936115, 0.867207},
                                          {0.937022, 0.986063, 1.020076, 1.002000, 1.004196, 0.967409, 0.917517, 0.834299},
                                          {0.918030, 0.990247, 1.021756, 1.014253, 0.999047, 0.956159, 0.977406, 0.842177}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {0.980272, 1.251209, 1.183687, 1.121365, 1.073480, 1.023050, 0.984423, 0.766073},
                                            {1.274353, 1.227809, 1.144672, 1.113070, 1.040415, 1.034251, 0.988324, 0.682039},
                                            {1.243528, 1.209943, 1.148926, 1.095807, 1.067192, 1.024940, 0.978711, 0.747443},
                                            {1.242904, 1.198735, 1.124272, 1.070151, 1.017776, 0.988903, 0.943155, 0.759892},
                                            {1.274353, 1.174759, 1.143603, 1.073152, 1.035248, 1.009034, 0.942473, 0.763631}
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
                                         {-0.196828, -0.075040, -0.059540, -0.090147, -0.168801, -0.215170, -0.330530, -0.420427},
                                         {-0.128518, -0.027088, -0.019251, -0.041696, -0.134721, -0.157173, -0.285224, -0.325479},
                                         {0.017672, 0.061632, 0.047054, 0.004491, -0.080555, -0.105739, -0.213883, -0.226184},
                                         {0.316975, 0.197993, 0.139924, 0.069684, -0.020624, -0.087247, -0.148538, -0.168119},
                                         {0.771093, 0.453277, 0.349118, 0.237450, 0.160700, 0.045833, 0.018618, -0.022275}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.250354, -0.156928, -0.116049, -0.100731, -0.138547, -0.136585, -0.112824, 0.066743},
                                         {-0.088011, -0.062177, -0.070153, -0.080326, -0.123888, -0.098164, -0.109129, -0.052343},
                                         {0.117415, 0.000152, -0.037164, -0.030129, -0.131634, -0.074912, -0.106370, -0.151782},
                                         {0.224609, 0.211171, 0.077454, 0.010525, -0.064537, -0.069700, -0.106284, -0.175537},
                                         {0.372422, 0.624324, 0.354456, 0.186364, 0.055980, 0.073980, -0.016471, -0.001540}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.531107, -0.273528, -0.201928, -0.284344, -0.209342, -0.127973, -0.026659, 0.160195},
                                           {-0.285947, -0.109632, -0.048825, -0.169340, -0.186199, -0.079545, 0.005589, 0.219579},
                                           {-0.242163, -0.115458, 0.014787, -0.101414, -0.158060, -0.039287, 0.021998, 0.220548},
                                           {-0.174112, -0.106879, 0.078141, -0.034449, -0.016076, 0.026918, 0.069406, 0.266545},
                                           {-0.086245, -0.062128, 0.180774, 0.156378, 0.144139, 0.154574, -0.008910, 0.423626}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.977727, 0.979226, 0.989133, 1.002220, 1.031504, 1.051924, 1.107578, 1.092703},
                                          {0.956762, 0.965964, 0.981122, 0.990063, 1.025744, 1.033574, 1.099267, 1.040435},
                                          {0.945063, 0.957751, 0.971645, 0.978760, 1.011226, 1.012672, 1.062524, 1.005409},
                                          {0.948883, 0.948067, 0.962529, 0.968142, 0.982658, 1.014397, 1.036750, 0.999503},
                                          {0.934816, 0.966958, 0.987133, 0.996842, 1.010535, 1.050921, 1.067253, 1.057234}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.083170, 1.131660, 1.062033, 1.016024, 0.982287, 0.927062, 0.963633, 1.100638},
                                          {0.943351, 1.094854, 1.056409, 1.010605, 0.992640, 0.931446, 0.974996, 1.004806},
                                          {1.044172, 1.037065, 1.053681, 1.002571, 0.995572, 0.933668, 0.954138, 0.862241},
                                          {0.972285, 1.025684, 1.043636, 1.002153, 0.989476, 0.937025, 0.937823, 0.800550},
                                          {1.010903, 1.049958, 1.054553, 1.026049, 1.004937, 0.983670, 0.928826, 0.800550}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.188218, 1.190626, 1.196265, 1.193219, 1.060627, 0.967922, 0.977149, 1.187465},
                                            {1.188218, 1.190626, 1.178498, 1.162875, 1.074818, 0.967511, 0.969763, 1.162103},
                                            {1.188218, 1.190626, 1.156453, 1.145453, 1.044231, 0.973746, 0.953075, 1.174961},
                                            {1.188218, 1.172788, 1.127426, 1.155375, 1.134030, 0.994062, 0.959919, 1.154276},
                                            {1.097285, 1.143206, 1.146974, 1.206302, 1.098323, 1.044445, 1.001289, 1.159706}
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
                                         {-0.289622, -0.131088, -0.077235, -0.070002, -0.078218, -0.084869, -0.164637, -0.229795},
                                         {-0.186043, -0.050958, -0.001672, 0.000918, -0.008057, -0.019380, -0.052167, -0.097639},
                                         {-0.054174, 0.026595, 0.046332, 0.035514, 0.012822, -0.001037, -0.039317, -0.092511},
                                         {0.183774, 0.108748, 0.088214, 0.033535, -0.009612, -0.044275, -0.061562, -0.102588},
                                         {0.575556, 0.329106, 0.254639, 0.163140, 0.095359, 0.057337, 0.048403, 0.076499}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.064384, -0.028766, -0.094976, -0.147923, -0.186167, -0.189076, -0.244179, 0.257015},
                                         {0.088492, 0.064997, -0.021785, -0.073168, -0.127582, -0.129209, -0.196937, -0.290775},
                                         {0.249562, 0.075254, 0.024555, -0.043983, -0.105794, -0.095302, -0.182824, -0.110057},
                                         {0.278555, 0.207490, 0.067068, -0.047485, -0.107397, -0.121857, -0.228939, -0.087735},
                                         {0.258174, 0.529313, 0.305111, 0.099348, -0.021590, -0.046363, -0.115240, -0.009693}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.226646, -0.062501, 0.016117, -0.063109, -0.093844, -0.114600, -0.062076, -0.159304},
                                           {-0.002790, 0.087481, 0.168390, 0.066403, -0.011168, -0.043094, -0.009507, -0.008350},
                                           {-0.018909, 0.063554, 0.171750, 0.087815, 0.021623, -0.022053, 0.025570, 0.007675},
                                           {-0.045649, -0.058214, 0.124238, 0.056287, 0.013588, -0.024546, 0.008559, 0.051425},
                                           {-0.051875, -0.090621, 0.132572, 0.200793, 0.150273, 0.105279, 0.127645, 0.216836}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.990354, 0.989430, 1.001740, 1.015016, 1.033771, 1.049265, 1.073431, 0.984293},
                                          {0.974059, 0.975168, 0.988757, 0.996874, 1.021845, 1.036287, 1.034883, 0.827742},
                                          {0.963158, 0.968667, 0.978734, 0.985739, 1.008324, 1.007514, 1.006969, 0.980196},
                                          {0.959217, 0.954765, 0.962253, 0.973413, 0.982229, 0.990850, 0.975105, 0.894812},
                                          {0.985232, 0.968861, 0.984010, 0.996133, 1.001125, 1.011135, 1.005278, 0.929941}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.110540, 1.138837, 1.077577, 1.008714, 0.986706, 0.962464, 0.925740, 1.186533},
                                          {1.042350, 1.080231, 1.063390, 1.010297, 0.995529, 0.957666, 0.910992, 0.838504},
                                          {1.025853, 1.043577, 1.046277, 1.003180, 0.981662, 0.953283, 0.888012, 0.819415},
                                          {0.959753, 1.023271, 1.041758, 1.003621, 0.984889, 0.957532, 0.869990, 0.864340},
                                          {0.980846, 1.031960, 1.051989, 1.020264, 1.001430, 0.982136, 0.896090, 0.816460}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.299133, 1.276564, 1.226873, 1.182572, 1.083598, 1.010379, 1.025178, 0.773604},
                                            {1.280357, 1.241696, 1.187989, 1.147748, 1.098511, 1.011225, 0.996644, 1.052293},
                                            {1.279895, 1.214740, 1.159997, 1.123810, 1.051801, 1.003962, 0.976248, 0.992783},
                                            {1.253531, 1.186051, 1.129703, 1.107305, 1.023628, 0.987529, 0.988231, 1.084366},
                                            {1.254557, 1.167449, 1.131274, 1.153119, 1.062385, 1.045470, 1.024286, 1.097600}
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
                                         {-0.261147, -0.097552, -0.040036, -0.034878, -0.016643, -0.028297, -0.080007, -0.143438},
                                         {-0.165304, -0.031117, 0.023352, 0.032519, 0.017905, 0.034123, 0.008130, -0.033841},
                                         {-0.043188, 0.038370, 0.064119, 0.058428, 0.049340, 0.059173, 0.028136, -0.022080},
                                         {0.002647, 0.097816, 0.083242, 0.037402, 0.000587, -0.017343, -0.014407, -0.050207},
                                         {0.536993, 0.300864, 0.002182, 0.141599, 0.085799, 0.057342, 0.070835, 0.027963}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {-0.053502, 0.040452, -0.039256, -0.113213, -0.138300, -0.200325, -0.233384, -0.324480},
                                         {0.199417, 0.088825, 0.016996, -0.054013, -0.093688, -0.118856, -0.168093, -0.209238},
                                         {0.231309, 0.136212, 0.055351, -0.031316, -0.071528, -0.072771, -0.128634, -0.029958},
                                         {0.362749, 0.225385, 0.094825, -0.048980, -0.107167, -0.130925, -0.213373, -0.075526},
                                         {0.482420, 0.544269, 0.293814, 0.093279, -0.017018, -0.067444, -0.118263, 0.005495}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.170205, -0.002647, 0.065535, 0.029479, -0.023646, -0.038943, -0.053199, -0.104667},
                                           {0.018332, 0.155277, 0.203296, 0.140974, 0.045770, 0.021113, 0.004024, -0.040219},
                                           {-0.048238, 0.114970, 0.200758, 0.134594, 0.064641, 0.018614, 0.047825, 0.062535},
                                           {-0.016896, -0.034899, 0.128837, 0.103848, 0.019736, -0.012493, -0.003445, -0.017372},
                                           {-0.036684, -0.064504, 0.117689, 0.248536, 0.129954, 0.063029, 0.076023, 0.157235}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {1.008552, 1.010482, 1.019703, 1.021704, 1.035648, 1.050802, 1.051815, 0.941403},
                                          {0.983639, 0.983046, 0.994809, 1.015928, 1.030123, 1.033464, 1.014765, 0.983297},
                                          {0.976606, 0.977820, 0.985901, 0.994227, 1.007506, 0.699220, 0.991977, 0.944983},
                                          {1.000137, 0.964249, 0.974332, 0.981834, 0.985157, 0.996496, 0.952734, 0.878634},
                                          {1.000179, 0.975619, 1.000191, 0.997088, 1.006094, 0.998703, 0.963105, 0.928066}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.055595, 1.123185, 1.088752, 1.037199, 0.999139, 0.969519, 0.936457, 0.978947},
                                          {1.032048, 1.081092, 1.068877, 1.034409, 1.002630, 0.962778, 0.894452, 0.703162},
                                          {1.000203, 1.045034, 1.057223, 1.025865, 1.000571, 0.970833, 0.906630, 0.874453},
                                          {0.892144, 1.030377, 1.038588, 1.010058, 1.011044, 0.965594, 0.877744, 0.803302},
                                          {1.177154, 1.026071, 1.051548, 1.032313, 1.026743, 0.990480, 0.922000, 0.797010}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {1.291343, 1.286387, 1.232884, 1.188643, 1.097036, 1.044273, 1.046515, 1.128127},
                                            {1.289300, 1.239678, 1.197778, 1.156046, 1.077950, 1.038566, 1.034417, 1.017770},
                                            {1.277299, 1.215380, 1.161371, 1.120366, 1.067689, 1.015390, 0.986347, 0.967165},
                                            {1.264481, 1.187517, 1.132495, 1.098047, 1.055154, 1.022212, 0.999335, 1.001760},
                                            {1.256911, 1.195016, 1.140113, 1.132385, 1.067037, 1.054896, 1.034283, 1.064760}
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
                                         {-0.090930, 0.132853, 0.230463, 0.282549, 0.327439, 0.356671, 0.355729, 0.529217},
                                         {0.000722, 0.241683, 0.343128, 0.390189, 0.402642, 0.468007, 0.442862, 0.515279},
                                         {0.091472, 0.241903, 0.306940, 0.334514, 0.358154, 0.364978, 0.371802, 0.313826},
                                         {0.143984, 0.147455, 0.176516, 0.172212, 0.151600, 0.164284, 0.163136, 0.185352},
                                         {0.384527, 0.222017, 0.198837, 0.150235, 0.109312, 0.119336, 0.139504, 0.082313}
                                        };
    vector<vector<Float_t> > meanKaon = {
                                         {0.312939, 0.366516, 0.268415, 0.187053, 0.137318, 0.170548, 0.193186, 0.311248},
                                         {-0.002826, 0.446725, 0.341936, -0.009771, 0.228717, 0.255151, 0.337719, 0.261362},
                                         {0.480352, 0.352318, 0.299527, 0.239067, 0.179967, 0.207959, 0.273301, 0.233576},
                                         {0.427801, -0.008707, 0.182990, 0.080928, 0.037711, 0.039527, 0.103358, 0.144510},
                                         {0.258103, 0.495955, 0.264893, 0.081950, -0.017562, -0.052678, 0.078696, 0.100991}
                                        };
    vector<vector<Float_t> > meanProton = {
                                           {-0.001437, 0.302682, 0.426951, 0.430816, 0.796084, 0.243965, 0.232297, 0.568155},
                                           {-0.001437, 0.468372, 0.561573, 0.560089, 0.458092, 0.361275, 0.367101, 0.642212},
                                           {0.269382, 0.373516, 0.473773, 0.465137, 0.393140, 0.313099, 0.297818, 0.579754},
                                           {0.116385, 0.097116, 0.245755, 0.270444, 0.195629, 0.121001, 0.105098, 0.054838},
                                           {-0.001437, -0.031265, 0.112972, 0.264868, 0.175428, 0.086603, 0.049322, 0.099330}
                                          };
    vector<vector<Float_t> > sigmaPion = {
                                          {0.978110, 0.992582, 1.007336, 1.027124, 1.026068, 1.005216, 0.942417, 0.717719},
                                          {0.999747, 0.977090, 0.993331, 1.012788, 0.994568, 1.007884, 0.899233, 0.749029},
                                          {0.940901, 0.957752, 0.968580, 0.989714, 0.984572, 0.988587, 0.915765, 0.819728},
                                          {0.919462, 0.930581, 0.949072, 0.963085, 0.958612, 0.955708, 0.906270, 0.730735},
                                          {0.949050, 0.928446, 0.954328, 0.968505, 0.970685, 0.960509, 0.903462, 0.834110}
                                         };
    vector<vector<Float_t> > sigmaKaon = {
                                          {1.071575, 1.095412, 1.086576, 1.054975, 1.047588, 1.018820, 0.946009, 0.988692},
                                          {0.894046, 1.056707, 1.060878, 1.000530, 1.040252, 0.999974, 0.955241, 0.905667},
                                          {1.024133, 1.017663, 1.037939, 1.032572, 1.013853, 0.978752, 0.944053, 0.892756},
                                          {0.957289, 0.990632, 1.018180, 1.002266, 1.002538, 0.966059, 0.944323, 0.840634},
                                          {0.920151, 0.986735, 1.020171, 1.012013, 1.006504, 0.961092, 0.969147, 0.834466}
                                         };
    vector<vector<Float_t> > sigmaProton = {
                                            {0.980272, 1.268935, 1.194479, 1.139570, 0.809393, 1.069264, 0.950762, 0.682040},
                                            {0.980272, 1.232782, 1.180821, 1.103780, 1.107370, 1.037780, 0.985345, 0.682039},
                                            {1.274353, 1.207787, 1.149272, 1.083776, 1.059558, 1.027917, 0.957288, 0.682039},
                                            {1.239954, 1.185397, 1.118727, 1.068887, 1.035531, 0.987174, 0.951478, 0.871668},
                                            {0.980272, 1.178614, 1.154873, 1.082578, 1.042558, 1.000969, 0.947090, 0.832644}
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