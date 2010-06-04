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
// Authors: D. Caffarri caffarri@pd.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de
//***********************************************************
#include "AliAODPidHF.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
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
  fAsym(kFALSE)
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
  fAsym(pid.fAsym)
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

  if(!CheckStatus(track,"TPC")) return -1;
  AliAODPid *pidObj = track->GetDetPid();
  
  Double_t dedx=pidObj->GetTPCsignal();
  Double_t mom = pidObj->GetTPCmomentum();
  AliTPCPIDResponse tpcResponse;
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

  if(!CheckStatus(track,"ITS")) return -1;

  Double_t mom=track->P();
  AliAODPid *pidObj = track->GetDetPid();

  Double_t dedx=pidObj->GetITSsignal();
  AliITSPIDResponse itsResponse;
  Int_t pid=-1;
  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fnSigma[4];
   for(Int_t ipart=0;ipart<5;ipart++){
    AliPID::EParticleType type=AliPID::EParticleType(ipart);
    Double_t nsigma = TMath::Abs(itsResponse.GetNumberOfSigmas(mom,dedx,type));
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

 if(!CheckStatus(track,"TOF")) return -1;

 Double_t time[AliPID::kSPECIESN];
 AliAODPid *pidObj = track->GetDetPid();
 pidObj->GetIntegratedTimes(time);
 Double_t sigTOF=pidObj->GetTOFsignal();
// AliTOFPIDResponse tofResponse;
 Int_t pid=-1;

  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fTOFSigma*fnSigma[3];
   for(Int_t ipart=0;ipart<5;ipart++){
    //AliPID::EParticleType type=AliPID::EParticleType(ipart);
    //Double_t nsigma = tofResponse.GetExpectedSigma(track->P(),time[type],AliPID::ParticleMass(type));
    Double_t nsigma=TMath::Abs(sigTOF-time[ipart]);
    if((nsigma<nsigmaMax) && (nsigma<fnSigma[3]*fTOFSigma)) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
   //AliPID::EParticleType type=AliPID::EParticleType(specie);
   //Double_t nsigma = TMath::Abs(tofResponse.GetExpectedSigma(track->P(),time[type],AliPID::ParticleMass(type)));
    Double_t nsigma=TMath::Abs(sigTOF-time[specie]);
   if (nsigma>fnSigma[3]*fTOFSigma) {
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
void AliAODPidHF::BayesianProbability(AliAODTrack *track,TString detectors,Double_t *pid) const{
// bayesian PID for single detectors or combined

  if(detectors.Contains("ITS")) {BayesianProbabilityITS(track,pid);return;}
  if(detectors.Contains("TPC")) {BayesianProbabilityTPC(track,pid);return;}
  if(detectors.Contains("TOF")) {BayesianProbabilityTOF(track,pid);return;}

  if(detectors.Contains("All")) {
    Double_t probITS[5]={0.,0.,0.,0.,0.};
    Double_t probTPC[5]={0.,0.,0.,0.,0.};
    Double_t probTOF[5]={0.,0.,0.,0.,0.};
    BayesianProbabilityITS(track,probITS);
    BayesianProbabilityTPC(track,probTPC);
    BayesianProbabilityTOF(track,probTOF);
    Double_t probTot[5]={0.,0.,0.,0.,0.};
    for(Int_t i=0;i<5;i++){
     probTot[i]=probITS[i]*probTPC[i]*probTOF[i];
    }
    for(Int_t i2=0;i2<5;i2++){
     pid[i2]=probTot[i2]*fPriors[i2]/(probTot[0]*fPriors[0]+probTot[1]*fPriors[1]+probTot[2]*fPriors[2]+probTot[3]*fPriors[3]+probTot[4]*fPriors[4]);
    }
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
  prob[ind]=itspid[ind]*fPriors[ind]/(itspid[0]*fPriors[0]+itspid[1]*fPriors[1]+itspid[2]*fPriors[2]+itspid[3]*fPriors[3]+itspid[4]*fPriors[4]);
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
  if(tpcpid[ind]>0.) {
   prob[ind]=tpcpid[ind]*fPriors[ind]/(tpcpid[0]*fPriors[0]+tpcpid[1]*fPriors[1]+tpcpid[2]*fPriors[2]+tpcpid[3]*fPriors[3]+tpcpid[4]*fPriors[4]);
  }else{
   prob[ind]=0.;
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
  prob[ind]=tofpid[ind]*fPriors[ind]/(tofpid[0]*fPriors[0]+tofpid[1]*fPriors[1]+tofpid[2]*fPriors[2]+tofpid[3]*fPriors[3]+tofpid[4]*fPriors[4]);
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
  if(trdpid[ind]>0.) {
   prob[ind]=trdpid[ind]*fPriors[ind]/(trdpid[0]*fPriors[0]+trdpid[1]*fPriors[1]+trdpid[2]*fPriors[2]+trdpid[3]*fPriors[3]+trdpid[4]*fPriors[4]);
  }else{
   prob[ind]=0.;
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
  if(nPointsForPid<3) return kFALSE;// track not to be used for PID purposes
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
 // convention (temporary): 
 // for the single detectors: -1 = kFALSE, 1 = kTRUE, 0 = compatible
 // the method returns the sum of the response of the 2 detectors
  if(!CheckStatus(track,"TPC") && !CheckStatus(track,"TOF")) return 0;

  Int_t TPCinfo=0;
  if(CheckStatus(track,"TPC")) {
   if(fAsym) {
    if(TPCRawAsym(track,specie)) TPCinfo=1;
   }else{
    if(specie==2 && IsPionRaw(track,"TPC")) TPCinfo=1;
    if(specie==3 && IsKaonRaw(track,"TPC")) TPCinfo=1;
    if(specie==4 && IsProtonRaw(track,"TPC")) TPCinfo=1;
   }


  if(compat && TPCinfo<0){
   SetSigma(0,3.);
   if(specie==2 && IsPionRaw(track,"TPC")) TPCinfo=0;
   if(specie==3 && IsKaonRaw(track,"TPC")) TPCinfo=0;
   if(specie==4 && IsProtonRaw(track,"TPC")) TPCinfo=0;
  }

 }

 if(!CheckStatus(track,"TOF")) return TPCinfo;

 Int_t TOFinfo=-1;
 
  if(specie==2 && IsPionRaw(track,"TOF")) TOFinfo=1;
  if(specie==3 && IsKaonRaw(track,"TOF")) TOFinfo=1;
  if(specie==4 && IsProtonRaw(track,"TOF")) TOFinfo=1;

 if(compat && TOFinfo>0){
  Double_t ptrack=track->P();
  if(ptrack>1.5) TOFinfo=0;
 }

 return TPCinfo+TOFinfo;
 }
 if(mode==2){
  //TPC & TOF (a la' Yifei)
 // convention (temporary): -1 = kFALSE, 1 = kTRUE, 0 = not identified
  Int_t TPCinfo=1; 
  if(CheckStatus(track,"TPC")) {
   if(fAsym){
    if(!TPCRawAsym(track,specie)) TPCinfo=-1;
   }else{
    if(specie==2 && !IsPionRaw(track,"TPC")) TPCinfo=-1;
    if(specie==3 && !IsKaonRaw(track,"TPC")) TPCinfo=-1;
    if(specie==4 && !IsProtonRaw(track,"TPC")) TPCinfo=-1;
   }
  }

  Int_t TOFinfo=1;
  if(!CheckStatus(track,"TOF")) return TPCinfo;

   if(specie==2 && !IsPionRaw(track,"TOF")) TOFinfo=-1;
   if(specie==3 && !IsKaonRaw(track,"TOF")) TOFinfo=-1;
   if(specie==4 && !IsProtonRaw(track,"TOF")) TOFinfo=-1;

   if(TOFinfo==1 && TPCinfo==1) return 1;
   return -1;

 }

 if(mode==3){
 //TPC for p<fPLimit[0], TOF for p>=fPLimit[0] (a la' Andrea A.)
 // convention (temporary): -1 = kFALSE, 1 = kTRUE, 0 = not identified
  
  if(!CheckStatus(track,"TPC") && !CheckStatus(track,"TOF")) return 0;

  Double_t ptrack=track->P();
  

  Int_t TPCinfo=-1;

   if(ptrack<fPLimit[0]) {  
    if(!CheckStatus(track,"TPC")) return 0;
    if(fAsym) {
     if(TPCRawAsym(track,specie)) TPCinfo=1;
    }else{
     if(specie==2 && IsPionRaw(track,"TPC")) TPCinfo=1;
     if(specie==3 && IsKaonRaw(track,"TPC")) TPCinfo=1;
     if(specie==4 && IsProtonRaw(track,"TPC")) TPCinfo=1;
    } 
    return TPCinfo;
   }

   Int_t TOFinfo=-1;
   if(ptrack>=fPLimit[0]){
    if(!CheckStatus(track,"TOF")) return 0;
    if(specie==2 && IsPionRaw(track,"TOF")) TOFinfo=1;
    if(specie==3 && IsKaonRaw(track,"TOF")) TOFinfo=1;
    if(specie==4 && IsProtonRaw(track,"TOF")) TOFinfo=1;
    return TOFinfo;
   }

 }

 return -1;
}
   
