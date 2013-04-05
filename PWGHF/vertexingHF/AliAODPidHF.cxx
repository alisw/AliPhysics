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

//***********************************************************
// Class AliAODPidHF
// class for PID with AliAODRecoDecayHF
// Authors: D. Caffarri caffarri@pd.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de P. Antonioli pietro.antonioli@bo.infn.it
//***********************************************************
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


ClassImp(AliAODPidHF)

//------------------------------
AliAODPidHF::AliAODPidHF():
  AliAODPid(),
  fnNSigma(5),
  fnSigma(),
  fTOFSigma(160.),
  fnPriors(5),
  fPriors(),
  fnPLimit(2),
  fPLimit(),
  fAsym(kFALSE),
  fTPC(kFALSE),
  fTOF(kFALSE),
  fITS(kFALSE),
  fTRD(kFALSE),
  fMatch(0),
  fCompat(kFALSE),
  fPCompatTOF(1.5),
  fnNSigmaCompat(2),
  fnSigmaCompat(),
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
  fUseCombined(kFALSE)
{
 //
 // Default constructor
 //
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

}
//----------------------
AliAODPidHF::~AliAODPidHF()
{
      // destructor
 //   if(fPLimit) delete fPLimit;
 //   if(fnSigma)  delete fnSigma;
 //   if(fPriors)  delete fPriors;
  delete fTPCResponse;
  for (Int_t ispecies=0;ispecies<AliPID::kSPECIES;++ispecies) {
    delete fPriorsH[ispecies];
  }
}
//------------------------
AliAODPidHF::AliAODPidHF(const AliAODPidHF& pid) :
  AliAODPid(pid),
  fnNSigma(pid.fnNSigma),
  fnSigma(0),
  fTOFSigma(pid.fTOFSigma),
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
  fCompat(pid.fCompat),
  fPCompatTOF(pid.fPCompatTOF),
  fnNSigmaCompat(pid.fnNSigmaCompat),
  fnSigmaCompat(pid.fnSigmaCompat),
  fMC(pid.fMC),
  fOnePad(pid.fOnePad),
  fMCLowEn2011(pid.fMCLowEn2011),
  fppLowEn2011(pid.fppLowEn2011),
  fPbPb(pid.fPbPb),
  fTOFdecide(pid.fTOFdecide),
  fOldPid(pid.fOldPid),
  fPtThresholdTPC(pid.fPtThresholdTPC),
  fMaxTrackMomForCombinedPID(pid.fMaxTrackMomForCombinedPID),
  fPidResponse(pid.fPidResponse),
  fPidCombined(pid.fPidCombined),
  fTPCResponse(0x0),
  fCombDetectors(pid.fCombDetectors),
  fUseCombined(pid.fUseCombined)
{
  
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
  fPriors = new Double_t[fnPriors];
  for(Int_t i=0;i<fnPriors;i++){
    fPriors[i]=pid.fPriors[i];
  }
  for(Int_t i=0;i<AliPID::kSPECIES;i++){
    fPriorsH[i]=pid.fPriorsH[i];
  }

  if(pid.fTPCResponse) fTPCResponse = new AliTPCPIDResponse(*(pid.fTPCResponse));
  //fPidResponse = new AliPIDResponse(*(pid.fPidResponse));
  //fPidCombined = new AliPIDCombined(*(pid.fPidCombined));  
    
}
//----------------------
Int_t AliAODPidHF::RawSignalPID(AliAODTrack *track, TString detector) const{
// raw PID for single detectors, returns the particle type with smaller sigma
   Int_t specie=-1;
   if(detector.Contains("ITS")) return ApplyPidITSRaw(track,specie);
   if(detector.Contains("TPC")) return ApplyPidTPCRaw(track,specie);
   if(detector.Contains("TOF")) return ApplyPidTOFRaw(track,specie);

  return specie;

}
//---------------------------
Bool_t AliAODPidHF::IsKaonRaw(AliAODTrack *track, TString detector) const{
// checks if the track can be a kaon, raw PID applied for single detectors
 Int_t specie=0;

 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,3);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,3);
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,3);

 if(specie==3) return kTRUE;
 return kFALSE;
}
//---------------------------
Bool_t AliAODPidHF::IsPionRaw (AliAODTrack *track, TString detector) const{
// checks if the track can be a pion, raw PID applied for single detectors

 Int_t specie=0;

 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,2);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,2);
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,2);

 if(specie==2) return kTRUE;
 return kFALSE;
}
//---------------------------
Bool_t AliAODPidHF::IsProtonRaw (AliAODTrack *track, TString detector) const{
// checks if the track can be a proton raw PID applied for single detectors

 Int_t specie=0;
 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,4);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,4); 
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,4);

 if(specie==4) return kTRUE;

 return kFALSE;
}
//--------------------------
Bool_t AliAODPidHF::IsElectronRaw(AliAODTrack *track, TString detector) const{
// checks if the track can be an electron raw PID applied for single detectors

 Int_t specie=-1;
 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,0);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,0);
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,0);

 if(specie==0) return kTRUE;

 return kFALSE;
}
//--------------------------
Int_t AliAODPidHF::ApplyPidTPCRaw(AliAODTrack *track,Int_t specie) const{
// n-sigma cut, TPC PID

  Double_t nsigma;
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
// truncated mean, ITS PID

  if(!CheckITSPIDStatus(track)) return 0;
  Int_t pid=-1;

  if(fOldPid){
  Double_t mom=track->P();
  AliAODPid *pidObj = track->GetDetPid();

  Double_t dedx=pidObj->GetITSsignal();
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
   if(clumap&(1<<i)) ++nPointsForPid;
  }

  Bool_t isSA=kTRUE;
  if(track->GetStatus() & AliESDtrack::kTPCin) isSA = kFALSE;

  AliITSPIDResponse itsResponse;
  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fnSigma[4];
   for(Int_t ipart=0;ipart<5;ipart++){
    AliPID::EParticleType type=AliPID::EParticleType(ipart);
    Double_t nsigma = TMath::Abs(itsResponse.GetNumberOfSigmas(mom,dedx,type,nPointsForPid,isSA));
    if((nsigma<nsigmaMax) && (nsigma<fnSigma[4])) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
   AliPID::EParticleType type=AliPID::EParticleType(specie);
    Double_t nsigma = TMath::Abs(itsResponse.GetNumberOfSigmas(mom,dedx,type));
   if (nsigma>fnSigma[4]) {
    pid=-1; 
   }else{
    pid=specie;
   }
  }
 }else{ // old pid

  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fnSigma[4];
   for(Int_t ipart=0;ipart<5;ipart++){
    AliPID::EParticleType type=AliPID::EParticleType(ipart);
    Double_t nsigma = TMath::Abs(fPidResponse->NumberOfSigmasITS(track,type));
    if((nsigma<nsigmaMax) && (nsigma<fnSigma[4])) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
   AliPID::EParticleType type=AliPID::EParticleType(specie);
    Double_t nsigma = TMath::Abs(fPidResponse->NumberOfSigmasITS(track,type));
   if (nsigma>fnSigma[4]) {
    pid=-1;
   }else{
    pid=specie;
   }
  }
 } //new pid

 return pid; 
}
//----------------------------
Int_t AliAODPidHF::ApplyPidTOFRaw(AliAODTrack *track,Int_t specie) const{
// n-sigma cut, TOF PID

  Double_t nsigma;
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
    if(GetnSigmaTOF(track,specie,nsigma)==1){
      nsigma=TMath::Abs(nsigma);
      if (nsigma>fnSigma[3]) pid=-1; 
      else pid=specie;
    }
  }
  return pid; 
  /*
 Double_t time[AliPID::kSPECIESN];
 Double_t sigmaTOFPid[AliPID::kSPECIES];
 AliAODPid *pidObj = track->GetDetPid();
 pidObj->GetIntegratedTimes(time);
 Double_t sigTOF=pidObj->GetTOFsignal();

 AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();
 if (event) {
   AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
   if (tofH && fPidResponse) { // reading new AOD with new aliroot
     AliTOFPIDResponse TOFres = (AliTOFPIDResponse)fPidResponse->GetTOFResponse();
     sigTOF -= TOFres.GetStartTime(track->P());
     if (specie<0) {
       for (Int_t ipart = 0; ipart<5; ipart++) {
	 sigmaTOFPid[ipart]=TOFres.GetExpectedSigma(track->P(),time[ipart],AliPID::ParticleMass(ipart));
       }
     }
     else sigmaTOFPid[specie]=TOFres.GetExpectedSigma(track->P(),time[specie],AliPID::ParticleMass(specie)); //fTOFResponse is set in InitialiseEvent
   } else  pidObj->GetTOFpidResolution(sigmaTOFPid); // reading old AOD with new aliroot
 } else  pidObj->GetTOFpidResolution(sigmaTOFPid);  //reading old AOD with old aliroot

 Int_t pid=-1;

  if(specie<0){  
   Double_t sigmaTOFtrack;
   if (sigmaTOFPid[4]>0) sigmaTOFtrack=sigmaTOFPid[4];
   else sigmaTOFtrack=fTOFSigma;
   Double_t nsigmaMax=sigmaTOFtrack*fnSigma[3];
   for(Int_t ipart=0;ipart<5;ipart++){
    Double_t nsigma=TMath::Abs(sigTOF-time[ipart]);
    if (sigmaTOFPid[ipart]>0) sigmaTOFtrack=sigmaTOFPid[ipart]; 
    else sigmaTOFtrack=fTOFSigma;  // backward compatibility for old AODs
    if((nsigma<nsigmaMax) && (nsigma<fnSigma[3]*sigmaTOFtrack)) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
    Double_t nsigma=TMath::Abs(sigTOF-time[specie]);
    Double_t sigmaTOFtrack;
    if (sigmaTOFPid[specie]>0) sigmaTOFtrack=sigmaTOFPid[specie]; 
    else sigmaTOFtrack=fTOFSigma;  // backward compatibility for old AODs
    if (nsigma>fnSigma[3]*sigmaTOFtrack) {
      pid=-1; 
    }else{
      pid=specie;
    }
  }
 return pid; 
  */
}
//------------------------------
void AliAODPidHF::CombinedProbability(AliAODTrack *track,Bool_t *type) const{
// combined PID stored inside the AOD track

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
  // ITS PID quality cuts
  if ((track->GetStatus()&AliESDtrack::kITSin)==0) return kFALSE;
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nPointsForPid;
  }
  if(nPointsForPid<3) return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckTPCPIDStatus(AliAODTrack *track) const{
  // TPC PID quality cuts
  if ((track->GetStatus()&AliESDtrack::kTPCin )==0) return kFALSE;
  UShort_t nTPCClus=track->GetTPCClusterMap().CountBits();
  if (nTPCClus<70) return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckTOFPIDStatus(AliAODTrack *track) const{
  // TOF PID quality cuts
  if ((track->GetStatus()&AliESDtrack::kTOFout )==0)    return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTIME )==0)     return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTOFpid )==0)   return kFALSE;
  if (!(track->GetStatus()&AliESDtrack::kTOFmismatch)==0)    return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckTRDPIDStatus(AliAODTrack *track) const{
  // TRD PID quality cuts
  if ((track->GetStatus()&AliESDtrack::kTRDout )==0)   return kFALSE;
  return kTRUE;
}
//--------------------------------
Bool_t AliAODPidHF::CheckStatus(AliAODTrack *track,TString detectors) const{

// Quality cuts on the tracks, detector by detector

 if(detectors.Contains("ITS")){
  if ((track->GetStatus()&AliESDtrack::kITSin)==0) return kFALSE;
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
   if(clumap&(1<<i)) ++nPointsForPid;
  }
  if(nPointsForPid<3) return kFALSE;
 }

 if(detectors.Contains("TPC")){
   if ((track->GetStatus()&AliESDtrack::kTPCin )==0) return kFALSE;
   UShort_t nTPCClus=track->GetTPCClusterMap().CountBits();
   if (nTPCClus<70) return kFALSE;
 }

 if(detectors.Contains("TOF")){
   if ((track->GetStatus()&AliESDtrack::kTOFout )==0)    return kFALSE;
   if ((track->GetStatus()&AliESDtrack::kTIME )==0)     return kFALSE;
   if ((track->GetStatus()&AliESDtrack::kTOFpid )==0)   return kFALSE;
   if (!(track->GetStatus()&AliESDtrack::kTOFmismatch)==0)    return kFALSE;
 }


 if(detectors.Contains("TRD")){
  if ((track->GetStatus()&AliESDtrack::kTRDout )==0)   return kFALSE;
  //UChar_t ntracklets = track->GetTRDntrackletsPID();
  //if(ntracklets<4) return kFALSE;
 }

 return kTRUE;
}
//--------------------------------------------
Bool_t AliAODPidHF::TPCRawAsym(AliAODTrack* track,Int_t specie) const{
// TPC nsigma cut PID, different sigmas in different p bins

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
  // combination of the PID info coming from TPC and TOF

  Double_t ptrack=track->P();
  if(ptrack>fMaxTrackMomForCombinedPID) return 1;

  Bool_t okTPC=CheckTPCPIDStatus(track);
  Bool_t okTOF=CheckTOFPIDStatus(track);
  if(ptrack>fPtThresholdTPC) okTPC=kFALSE;

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
	  Double_t sig0tmp=fnSigma[3];
	  SetSigma(3,fnSigmaCompat[1]);
	  if(ApplyPidTOFRaw(track,specie)==specie) tTOFinfo=0;
	  SetSigma(3,sig0tmp);
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

  return -1;

}
//----------------------------------   
Int_t AliAODPidHF::MakeRawPid(AliAODTrack *track, Int_t specie){
// general method to compute PID
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
  // TPC bethe bloch parameters
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
  // Set Bethe Bloch Parameters

 Double_t alephParameters[5];
 GetTPCBetheBlochParams(alephParameters);
 fTPCResponse->SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2],alephParameters[3],alephParameters[4]);

 return;
}
//-----------------------
Bool_t AliAODPidHF::IsTOFPiKexcluded(AliAODTrack *track,Double_t nsigmaK){
  // TOF proton compatibility

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

	//
	// method setting the prior distributions to the AliPIDCombined object of the AliAODPidHF data member
	// all the checks are done directly in the AliPIDCombined object
	//

	GetPidCombined()->SetPriorDistribution(type,prior);
}
//--------------------------------------------------------------------------
void AliAODPidHF::DrawPrior(AliPID::EParticleType type){

	//
	// Drawing prior distribution for type "type"

	new TCanvas();
	GetPidCombined()->GetPriorDistribution(type)->Draw();
}

//--------------------------------------------------------------------------
Int_t AliAODPidHF::GetnSigmaTPC(AliAODTrack *track, Int_t species, Double_t &nsigma) const{
  // get n sigma for TPC 

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
  // get n sigma for TOF

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
  // Exclude a given hypothesis (labelTracks) in detector

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

    if (!(CheckTOFPIDStatus(track))) return kFALSE;
    Double_t nsigma=0.;
    if (GetnSigmaTOF(track,labelTrack,nsigma)==1){
      if(nsigma>nsigmaCut) return kTRUE;
    }
    return kFALSE;

  }
  return kFALSE;

}
//-----------------------------
void AliAODPidHF::SetPriorsHistos(TString priorFileName){
  // Set histograms with priors

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
  // Configuration of combined Bayesian PID

 fPidCombined->SetSelectedSpecies(AliPID::kSPECIES);
  for (Int_t ispecies=0;ispecies<AliPID::kSPECIES;++ispecies) {
    fPidCombined->SetPriorDistribution(static_cast<AliPID::EParticleType>(ispecies),fPriorsH[ispecies]);
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

Int_t AliAODPidHF::GetnSigmaITS(AliAODTrack *track,Int_t species, Double_t &nsigma) const{
  // get n sigma for ITS

  Double_t nsigmaITS=-999;

  if (!CheckITSPIDStatus(track)) return -1;

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

