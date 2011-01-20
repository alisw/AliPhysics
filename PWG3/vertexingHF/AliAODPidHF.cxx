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

//***********************************************************
// Class AliAODPidHF
// class for PID with AliAODRecoDecayHF
// Authors: D. Caffarri caffarri@pd.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de P. Antonioli pietro.antonioli@bo.infn.it
//***********************************************************
#include "AliAODPidHF.h"
#include "AliAODPid.h"
#include "AliPID.h"
#include "AliTPCPIDResponse.h"
#include "AliITSPIDResponse.h"
#include "AliTOFPIDResponse.h"
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
  fMC(kFALSE),
  fOnePad(kFALSE),
  fPbPb(kFALSE)
{
 //
 // Default constructor
 //
 fPLimit=new Double_t[fnPLimit];
 fnSigma=new Double_t[fnNSigma];
 fPriors=new Double_t[fnPriors];

 for(Int_t i=0;i<fnNSigma;i++){
  fnSigma[i]=0.;
 }
 for(Int_t i=0;i<fnPriors;i++){
  fPriors[i]=0.;
 }
 for(Int_t i=0;i<fnPLimit;i++){
  fPLimit[i]=0.;
 }

}
//----------------------
AliAODPidHF::~AliAODPidHF()
{
      // destructor
 //   if(fPLimit) delete fPLimit;
 //   if(fnSigma)  delete fnSigma;
 //   if(fPriors)  delete fPriors;
}
//------------------------
AliAODPidHF::AliAODPidHF(const AliAODPidHF& pid) :
  AliAODPid(pid),
  fnNSigma(pid.fnNSigma),
  fnSigma(pid.fnSigma),
  fTOFSigma(pid.fTOFSigma),
  fnPriors(pid.fnPriors),
  fPriors(pid.fPriors),
  fnPLimit(pid.fnPLimit),
  fPLimit(pid.fPLimit),
  fAsym(pid.fAsym),
  fTPC(pid.fTPC),
  fTOF(pid.fTOF),
  fITS(pid.fITS),
  fTRD(pid.fTRD),
  fMatch(pid.fMatch),
  fCompat(pid.fCompat),
  fMC(pid.fMC),
  fOnePad(pid.fOnePad),
  fPbPb(pid.fPbPb)
  {
  
  for(Int_t i=0;i<5;i++){
    fPriors[i]=pid.fPriors[i];
  }
  
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

  if(!CheckStatus(track,"TPC")) return 0;
  AliAODPid *pidObj = track->GetDetPid();
  
  Double_t dedx=pidObj->GetTPCsignal();
  Double_t mom = pidObj->GetTPCmomentum();
  AliTPCPIDResponse tpcResponse;
  SetBetheBloch(tpcResponse); 

  Int_t pid=-1;
  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fnSigma[0];
   for(Int_t ipart=0;ipart<5;ipart++){
    AliPID::EParticleType type=AliPID::EParticleType(ipart);
    Double_t nsigma = TMath::Abs(tpcResponse.GetNumberOfSigmas(mom,dedx,track->GetTPCNcls(),type));
    if((nsigma<nsigmaMax) && (nsigma<fnSigma[0])) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
   AliPID::EParticleType type=AliPID::EParticleType(specie);
    Double_t nsigma = TMath::Abs(tpcResponse.GetNumberOfSigmas(mom,dedx,track->GetTPCNcls(),type));
   if (nsigma>fnSigma[0]) {
    pid=-1; 
   }else{
    pid=specie;
   }
  }

 return pid;

}
//----------------------------
Int_t AliAODPidHF::ApplyPidITSRaw(AliAODTrack *track,Int_t specie) const{
// truncated mean, ITS PID

  if(!CheckStatus(track,"ITS")) return 0;

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
  Int_t pid=-1;
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
 return pid; 
}
//----------------------------
Int_t AliAODPidHF::ApplyPidTOFRaw(AliAODTrack *track,Int_t specie) const{
// n-sigma cut, TOF PID

 if(!CheckStatus(track,"TOF")) return 0;

 Double_t time[AliPID::kSPECIESN];
 Double_t sigmaTOFPid[AliPID::kSPECIES];
 AliAODPid *pidObj = track->GetDetPid();
 pidObj->GetIntegratedTimes(time);
 Double_t sigTOF=pidObj->GetTOFsignal();
 pidObj->GetTOFpidResolution(sigmaTOFPid);

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
//--------------------
void AliAODPidHF::BayesianProbability(AliAODTrack *track,Double_t *pid) const{
// bayesian PID for single detectors or combined

  if(fITS && !fTPC && !fTOF) {BayesianProbabilityITS(track,pid);return;}
  if(fTPC && !fITS && !fTOF) {BayesianProbabilityTPC(track,pid);return;}
  if(fTOF && !fITS && !fTPC) {BayesianProbabilityTOF(track,pid);return;}

    Double_t probITS[5]={1.,1.,1.,1.,1.};
    Double_t probTPC[5]={1.,1.,1.,1.,1.};
    Double_t probTOF[5]={1.,1.,1.,1.,1.};
    if(fITS) BayesianProbabilityITS(track,probITS);
    if(fTPC) BayesianProbabilityTPC(track,probTPC);
    if(fTOF) BayesianProbabilityTOF(track,probTOF);
    Double_t probTot[5]={0.,0.,0.,0.,0.};
    for(Int_t i=0;i<5;i++){
     probTot[i]=probITS[i]*probTPC[i]*probTOF[i];
    }
    for(Int_t i2=0;i2<5;i2++){
     pid[i2]=probTot[i2]*fPriors[i2]/(probTot[0]*fPriors[0]+probTot[1]*fPriors[1]+probTot[2]*fPriors[2]+probTot[3]*fPriors[3]+probTot[4]*fPriors[4]);
    }

 return;

}
//------------------------------------
void AliAODPidHF::BayesianProbabilityITS(AliAODTrack *track,Double_t *prob) const{

// bayesian PID for ITS
 AliAODpidUtil pid;
 Double_t itspid[AliPID::kSPECIES];
 pid.MakeITSPID(track,itspid);
 for(Int_t ind=0;ind<AliPID::kSPECIES;ind++){
  if(fTOF || fTPC || fTRD){
   prob[ind]=itspid[ind];
  }else{
   prob[ind]=itspid[ind]*fPriors[ind]/(itspid[0]*fPriors[0]+itspid[1]*fPriors[1]+itspid[2]*fPriors[2]+itspid[3]*fPriors[3]+itspid[4]*fPriors[4]);
  }
 }
 return;

}
//------------------------------------
void AliAODPidHF::BayesianProbabilityTPC(AliAODTrack *track,Double_t *prob) const{
// bayesian PID for TPC

 AliAODpidUtil pid;
 Double_t tpcpid[AliPID::kSPECIES];
 pid.MakeTPCPID(track,tpcpid);
 for(Int_t ind=0;ind<AliPID::kSPECIES;ind++){
  if(fTOF || fITS || fTRD){
   prob[ind]=tpcpid[ind];
  }else{
   prob[ind]=tpcpid[ind]*fPriors[ind]/(tpcpid[0]*fPriors[0]+tpcpid[1]*fPriors[1]+tpcpid[2]*fPriors[2]+tpcpid[3]*fPriors[3]+tpcpid[4]*fPriors[4]);
 }
}
 return;

}
//------------------------------------
void AliAODPidHF::BayesianProbabilityTOF(AliAODTrack *track,Double_t *prob) const{
// bayesian PID for TOF

 AliAODpidUtil pid;
 Double_t tofpid[AliPID::kSPECIES];
 pid.MakeTOFPID(track,tofpid);
 for(Int_t ind=0;ind<AliPID::kSPECIES;ind++){
  if(fTPC || fITS || fTRD){
   prob[ind]=tofpid[ind];
  }else{
  prob[ind]=tofpid[ind]*fPriors[ind]/(tofpid[0]*fPriors[0]+tofpid[1]*fPriors[1]+tofpid[2]*fPriors[2]+tofpid[3]*fPriors[3]+tofpid[4]*fPriors[4]);
 }
}
 return;

}
//---------------------------------
void AliAODPidHF::BayesianProbabilityTRD(AliAODTrack *track,Double_t *prob) const{
// bayesian PID for TRD

 AliAODpidUtil pid;
 Double_t trdpid[AliPID::kSPECIES];
 pid.MakeTRDPID(track,trdpid);
 for(Int_t ind=0;ind<AliPID::kSPECIES;ind++){
  if(fTPC || fITS || fTOF){
   prob[ind]=trdpid[ind];
  }else{
   prob[ind]=trdpid[ind]*fPriors[ind]/(trdpid[0]*fPriors[0]+trdpid[1]*fPriors[1]+trdpid[2]*fPriors[2]+trdpid[3]*fPriors[3]+trdpid[4]*fPriors[4]);
 }
}
  return;

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
  AliAODPid *pidObj = track->GetDetPid();
  Float_t *mom=pidObj->GetTRDmomentum();
  Int_t ntracklets=0;
  for(Int_t iPl=0;iPl<6;iPl++){
   if(mom[iPl]>0.) ntracklets++;
  }
   if(ntracklets<4) return kFALSE;
 }

 return kTRUE;
}
//--------------------------------------------
Bool_t AliAODPidHF::TPCRawAsym(AliAODTrack* track,Int_t specie) const{
// TPC nsigma cut PID, different sigmas in different p bins

  if(!CheckStatus(track,"TPC")) return kFALSE;
  AliAODPid *pidObj = track->GetDetPid();
  Double_t mom = pidObj->GetTPCmomentum();
  Double_t dedx=pidObj->GetTPCsignal();
  
  AliTPCPIDResponse tpcResponse;
  SetBetheBloch(tpcResponse); 
  AliPID::EParticleType type=AliPID::EParticleType(specie);
  Double_t nsigma = TMath::Abs(tpcResponse.GetNumberOfSigmas(mom,dedx,track->GetTPCNcls(),type)); 

  if(mom<fPLimit[0] && nsigma<fnSigma[0]) return kTRUE;
  if(mom<fPLimit[1] && mom>fPLimit[0] && nsigma<fnSigma[1]) return kTRUE;
  if(mom>fPLimit[1] && nsigma<fnSigma[2]) return kTRUE;

 return kFALSE;
}
//------------------
Int_t AliAODPidHF::MatchTPCTOF(AliAODTrack *track,Int_t mode,Int_t specie,Bool_t compat){
// combination of the PID info coming from TPC and TOF
 if(mode==1){
  //TOF || TPC (a la' Andrea R.)
 // convention: 
 // for the single detectors: -1 = kFALSE, 1 = kTRUE, 0 = compatible
 // the method returns the sum of the response of the 2 detectors
  if(fTPC && fTOF) {if(!CheckStatus(track,"TPC") && !CheckStatus(track,"TOF")) return 0;}

  
  Int_t tTPCinfo=0;
  if(fTPC){
  if(CheckStatus(track,"TPC")) {
   if(fAsym) {
    if(TPCRawAsym(track,specie)) {
      tTPCinfo=1;
     }else{
      tTPCinfo=-1;
     }
   }else{
    if(specie==2 && IsPionRaw(track,"TPC")) {
     tTPCinfo=1;
    }else{
     tTPCinfo=-1;
    }
    if(specie==3 && IsKaonRaw(track,"TPC")) {
     tTPCinfo=1;
    }else{
     tTPCinfo=-1;
    }
    if(specie==4 && IsProtonRaw(track,"TPC")) {
     tTPCinfo=1;
    }else{
     tTPCinfo=-1;
    }

   }


   if(compat && tTPCinfo<0){
    Double_t sig0tmp=fnSigma[0];
    SetSigma(0,3.);
    if(specie==2 && IsPionRaw(track,"TPC")) tTPCinfo=0;
    if(specie==3 && IsKaonRaw(track,"TPC")) tTPCinfo=0;
    if(specie==4 && IsProtonRaw(track,"TPC")) tTPCinfo=0;
    SetSigma(0,sig0tmp);
   }

  }
 }

 Int_t tTOFinfo=0;
 if(fTOF){
  if(!CheckStatus(track,"TOF") && fTPC) return tTPCinfo;

  tTOFinfo=-1;
 
  if(specie==2 && IsPionRaw(track,"TOF")) tTOFinfo=1;
  if(specie==3 && IsKaonRaw(track,"TOF")) tTOFinfo=1;
  if(specie==4 && IsProtonRaw(track,"TOF")) tTOFinfo=1;

  if(compat && tTOFinfo>0){
   Double_t ptrack=track->P();
   if(ptrack>1.5) tTOFinfo=0;
  }
 }

 if(tTPCinfo+tTOFinfo==0 && fITS){
  if(!CheckStatus(track,"ITS")) return tTPCinfo+tTOFinfo;
  Int_t tITSinfo = -1;
  if(specie==2 && IsPionRaw(track,"TOF")) tITSinfo=1;
  if(specie==3 && IsKaonRaw(track,"TOF")) tITSinfo=1;
  if(specie==4 && IsProtonRaw(track,"TOF")) tITSinfo=1;
  return tITSinfo;
 }

 return tTPCinfo+tTOFinfo;
}
 if(mode==2){
  //TPC & TOF (a la' Yifei)
 // convention: -1 = kFALSE, 1 = kTRUE, 0 = not identified
  Int_t tTPCinfo=0; 
  
  if(fTPC && CheckStatus(track,"TPC")) {
   tTPCinfo=1;
   if(fAsym){
    if(!TPCRawAsym(track,specie)) tTPCinfo=-1;
   }else{
    if(specie==2 && !IsPionRaw(track,"TPC")) tTPCinfo=-1;
    if(specie==3 && !IsKaonRaw(track,"TPC")) tTPCinfo=-1;
    if(specie==4 && !IsProtonRaw(track,"TPC")) tTPCinfo=-1;
   }
  }

  Int_t tTOFinfo=1;
  if(fTOF){
   if(fTPC && !CheckStatus(track,"TOF")) return tTPCinfo;

   if(specie==2 && !IsPionRaw(track,"TOF")) tTOFinfo=-1;
   if(specie==3 && !IsKaonRaw(track,"TOF")) tTOFinfo=-1;
   if(specie==4 && !IsProtonRaw(track,"TOF")) tTOFinfo=-1;
  }

 if(tTOFinfo==1 && tTPCinfo==1) return 1;

 if(tTPCinfo+tTOFinfo==0 && fITS){
  if(!CheckStatus(track,"ITS")) return tTPCinfo+tTOFinfo;
  Int_t tITSinfo = -1;
  if(specie==2 && IsPionRaw(track,"TOF")) tITSinfo=1;
  if(specie==3 && IsKaonRaw(track,"TOF")) tITSinfo=1;
  if(specie==4 && IsProtonRaw(track,"TOF")) tITSinfo=1;
  return tITSinfo;
 }

   return -1;

}

 if(mode==3){
 //TPC for p<fPLimit[0], TOF for p>=fPLimit[0] (a la' Andrea A.)
 // convention (temporary): -1 = kFALSE, 1 = kTRUE, 0 = not identified
  
  if(fTPC && fTOF) if(!CheckStatus(track,"TPC") && !CheckStatus(track,"TOF")) return 0;

  Double_t ptrack=track->P();
  

  Int_t tTPCinfo=-1;

   if(ptrack>=fPLimit[0] && ptrack<fPLimit[1] && fTPC) {  
    if(!CheckStatus(track,"TPC")) return 0;
    if(fAsym) {
     if(TPCRawAsym(track,specie)) tTPCinfo=1;
    }else{
     if(specie==2 && IsPionRaw(track,"TPC")) tTPCinfo=1;
     if(specie==3 && IsKaonRaw(track,"TPC")) tTPCinfo=1;
     if(specie==4 && IsProtonRaw(track,"TPC")) tTPCinfo=1;
    } 
    return tTPCinfo;
   }

   Int_t tTOFinfo=-1;
   if(ptrack>=fPLimit[1] && fTOF){
    if(!CheckStatus(track,"TOF")) return 0;
    if(specie==2 && IsPionRaw(track,"TOF")) tTOFinfo=1;
    if(specie==3 && IsKaonRaw(track,"TOF")) tTOFinfo=1;
    if(specie==4 && IsProtonRaw(track,"TOF")) tTOFinfo=1;
    return tTOFinfo;
   }

   Int_t tITSinfo=-1;
   if(ptrack<fPLimit[0] && fITS){
    if(!CheckStatus(track,"ITS")) return 0;
    if(specie==2 && IsPionRaw(track,"ITS")) tITSinfo=1;
    if(specie==3 && IsKaonRaw(track,"ITS")) tITSinfo=1;
    if(specie==4 && IsProtonRaw(track,"ITS")) tITSinfo=1;
    return tITSinfo;
   }

 }

 return -1;

}
//----------------------------------   
Int_t AliAODPidHF::MakeRawPid(AliAODTrack *track, Int_t specie){
// general method to compute PID
 if(fMatch>0){
  return MatchTPCTOF(track,fMatch,specie,fCompat); 
 }else{
  if(fTPC && !fTOF && !fITS) {
   Int_t tTPCres=ApplyPidTPCRaw(track,specie);
   if(tTPCres==specie){return 1;}else{return tTPCres;};
  }else{
   AliError("You should enable just one detector if you don't want to match");
   return 0;
  }
  if(fTOF && !fTPC && !fITS) {
   Int_t tTOFres=ApplyPidTOFRaw(track,specie); 
   if(tTOFres==specie){return 1;}else{return tTOFres;};
  }else{
   AliError("You should enable just one detector if you don't want to match");
   return 0;
  }

  if(fITS && !fTPC && !fTOF) {
   Int_t tITSres=ApplyPidITSRaw(track,specie);
   if(tITSres==specie){return 1;}else{return tITSres;};
  }else{
   AliError("You should enable just one detector if you don't want to match");
   return 0;
  }
 } 
  
}
//--------------------------------------------
void AliAODPidHF::SetBetheBloch(AliTPCPIDResponse tpcResp) const{

 Double_t alephParameters[5];
 if(fMC){
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;

 }else{

  if(fOnePad) {

   alephParameters[0] =1.34490e+00/50; 
   alephParameters[1] = 2.69455e+01; 
   alephParameters[2] = TMath::Exp(-2.97552+01); 
   alephParameters[3] = 2.35339e+00; 
   alephParameters[4] = 5.98079e+00;

  } else if(fPbPb) {

   alephParameters[0] = 1.25202/50.; 
   alephParameters[1] = 2.74992e+01; 
   alephParameters[2] = TMath::Exp(-3.31517e+01); 
   alephParameters[3] = 2.46246; 
   alephParameters[4] = 6.78938;

  } else {
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;

  }

 }

 tpcResp.SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2],alephParameters[3],alephParameters[4]);


 return;

}
//-----------------------
Bool_t AliAODPidHF::IsTOFPiKexcluded(AliAODTrack *track,Double_t nsigmaK){


 if(!CheckStatus(track,"TOF")) return 0;

  Double_t time[AliPID::kSPECIESN];
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  AliAODPid *pidObj = track->GetDetPid();
  pidObj->GetIntegratedTimes(time);
  Double_t sigTOF=pidObj->GetTOFsignal();
  pidObj->GetTOFpidResolution(sigmaTOFPid);
  Double_t sigmaTOFtrack;
  if (sigmaTOFPid[3]>0) sigmaTOFtrack=sigmaTOFPid[3];
  else sigmaTOFtrack=fTOFSigma;  // backward compatibility for old AODs
  
  if((sigTOF-time[3])>nsigmaK*sigmaTOFtrack)return kTRUE;// K, Pi excluded (->LIKELY A PROTON)
  
  return kFALSE;

}

