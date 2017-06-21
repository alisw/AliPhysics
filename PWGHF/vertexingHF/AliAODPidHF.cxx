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
fDefaultPriors(kTRUE)
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
fDefaultPriors(pid.fDefaultPriors)
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
      nSigmaTPC=fPidResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)specie);
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
      return fPidResponse->NumberOfSigmasITS(track, specie);
      break;
    case AliPIDResponse::kTPC:
      return fPidResponse->NumberOfSigmasTPC(track, specie);
      break;
    case AliPIDResponse::kTOF:
      return fPidResponse->NumberOfSigmasTOF(track, specie);
      break;
    default:
      return -999.;
      break;
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
