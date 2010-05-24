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
// Authors: D. Caffarri caffarri@bo.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de
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
  fSigma(3.),
  fPriors()
{
 //
 // Default constructor
 //

}
//----------------------
AliAODPidHF::~AliAODPidHF()
{
      // destructor
}
//------------------------
AliAODPidHF::AliAODPidHF(const AliAODPidHF& pid) :
  AliAODPid(pid),
  fSigma(pid.fSigma),
  fPriors(pid.fPriors)
  {
  
  for(Int_t i=0;i<5;i++){
    fPriors[i]=pid.fPriors[i];
  }
  
  }

//----------------------
Int_t AliAODPidHF::RawSignalPID(AliAODTrack *track, TString detector){

   Int_t specie=-1;
   if(detector.Contains("ITS")) return ApplyPidITSRaw(track,specie);
   if(detector.Contains("TPC")) return ApplyPidTPCRaw(track,specie);
   if(detector.Contains("TOF")) return ApplyPidTOFRaw(track,specie);

  return specie;

}
//---------------------------
Bool_t AliAODPidHF::IsKaonRaw (AliAODTrack *track, TString detector){

 Int_t specie=0;

 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,3);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,3);
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,3);

 if(specie==3) return kTRUE;
 return kFALSE;
}
//---------------------------
Bool_t AliAODPidHF::IsPionRaw (AliAODTrack *track, TString detector){

 Int_t specie=0;

 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,2);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,2);
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,2);

 if(specie==2) return kTRUE;
 return kFALSE;
}
//---------------------------
Bool_t AliAODPidHF::IsProtonRaw (AliAODTrack *track, TString detector){

 Int_t specie=0;
 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,4);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,4); 
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,4);

 if(specie==4) return kTRUE;

 return kFALSE;
}
//--------------------------
Bool_t AliAODPidHF::IsElectronRaw (AliAODTrack *track, TString detector){

 Int_t specie=-1;
 if(detector.Contains("ITS")) specie=ApplyPidITSRaw(track,0);
 if(detector.Contains("TPC")) specie=ApplyPidTPCRaw(track,0);
 if(detector.Contains("TOF")) specie=ApplyPidTOFRaw(track,0);

 if(specie==0) return kTRUE;

 return kFALSE;
}
//--------------------------
Int_t AliAODPidHF::ApplyPidTPCRaw(AliAODTrack *track,Int_t specie){

  if(!CheckStatus(track,"TPC")) return -1;
  AliAODPid *pidObj = track->GetDetPid();
  
  Double_t dedx=pidObj->GetTPCsignal();
  Double_t mom = pidObj->GetTPCmomentum();
  AliTPCPIDResponse tpcResponse;
  Int_t pid=-1;
  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fSigma;
   for(Int_t ipart=0;ipart<5;ipart++){
    AliPID::EParticleType type=AliPID::EParticleType(ipart);
    Double_t nsigma = TMath::Abs(tpcResponse.GetNumberOfSigmas(mom,dedx,track->GetTPCNcls(),type));
    if((nsigma<nsigmaMax) && (nsigma<fSigma)) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
   AliPID::EParticleType type=AliPID::EParticleType(specie);
    Double_t nsigma = TMath::Abs(tpcResponse.GetNumberOfSigmas(mom,dedx,track->GetTPCNcls(),type));
   if (nsigma>fSigma) {
    pid=-1; 
   }else{
    pid=specie;
   }
  }

 return pid;

}
//----------------------------
Int_t AliAODPidHF::ApplyPidITSRaw(AliAODTrack *track,Int_t specie){

  if(!CheckStatus(track,"ITS")) return -1;

  Double_t mom=track->P();
  AliAODPid *pidObj = track->GetDetPid();

  Double_t dedx=pidObj->GetITSsignal();
  AliITSPIDResponse itsResponse;
  Int_t pid=-1;
  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fSigma;
   for(Int_t ipart=0;ipart<5;ipart++){
    AliPID::EParticleType type=AliPID::EParticleType(ipart);
    Double_t nsigma = TMath::Abs(itsResponse.GetNumberOfSigmas(mom,dedx,type));
    if((nsigma<nsigmaMax) && (nsigma<fSigma)) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
   AliPID::EParticleType type=AliPID::EParticleType(specie);
    Double_t nsigma = TMath::Abs(itsResponse.GetNumberOfSigmas(mom,dedx,type));
   if (nsigma>fSigma) {
    pid=-1; 
   }else{
    pid=specie;
   }
  }
 return pid; 
}
//----------------------------
Int_t AliAODPidHF::ApplyPidTOFRaw(AliAODTrack *track,Int_t specie){

 if(!CheckStatus(track,"TOF")) return -1;

 Double_t time[AliPID::kSPECIESN];
 AliAODPid *pidObj = track->GetDetPid();
 pidObj->GetIntegratedTimes(time);
 AliTOFPIDResponse tofResponse;
 Int_t pid=-1;

  if(specie<0){  // from RawSignalPID : should return the particle specie to wich the de/dx is closer to the bethe-block curve -> performance to be checked
   Double_t nsigmaMax=fSigma;
   for(Int_t ipart=0;ipart<5;ipart++){
    AliPID::EParticleType type=AliPID::EParticleType(ipart);
    Double_t nsigma = tofResponse.GetExpectedSigma(track->P(),time[type],AliPID::ParticleMass(type));
    if((nsigma<nsigmaMax) && (nsigma<fSigma)) {
     pid=ipart;
     nsigmaMax=nsigma;
    }
   }
  }else{ // asks only for one particle specie
   AliPID::EParticleType type=AliPID::EParticleType(specie);
   Double_t nsigma = TMath::Abs(tofResponse.GetExpectedSigma(track->P(),time[type],AliPID::ParticleMass(type)));
   if (nsigma>fSigma) {
    pid=-1; 
   }else{
    pid=specie;
   }
  }
 return pid; 

}
//------------------------------
void AliAODPidHF::CombinedProbability(AliAODTrack *track,Bool_t *type){

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
void AliAODPidHF::BayesianProbability(AliAODTrack *track,TString detectors,Double_t *pid){

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
void AliAODPidHF::BayesianProbabilityITS(AliAODTrack *track,Double_t *prob){

 AliAODpidUtil pid;
 Double_t itspid[AliPID::kSPECIES];
 pid.MakeITSPID(track,itspid);
 for(Int_t ind=0;ind<AliPID::kSPECIES;ind++){
  prob[ind]=itspid[ind]*fPriors[ind]/(itspid[0]*fPriors[0]+itspid[1]*fPriors[1]+itspid[2]*fPriors[2]+itspid[3]*fPriors[3]+itspid[4]*fPriors[4]);
 }
 return;

}
//------------------------------------
void AliAODPidHF::BayesianProbabilityTPC(AliAODTrack *track,Double_t *prob){

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
void AliAODPidHF::BayesianProbabilityTOF(AliAODTrack *track,Double_t *prob){

 AliAODpidUtil pid;
 Double_t tofpid[AliPID::kSPECIES];
 //pid.MakeTOFPID(track,tofpid);
 for(Int_t ind=0;ind<AliPID::kSPECIES;ind++){
  prob[ind]=tofpid[ind]*fPriors[ind]/(tofpid[0]*fPriors[0]+tofpid[1]*fPriors[1]+tofpid[2]*fPriors[2]+tofpid[3]*fPriors[3]+tofpid[4]*fPriors[4]);
 }
 return;

}
//---------------------------------
void AliAODPidHF::BayesianProbabilityTRD(AliAODTrack *track,Double_t *prob){

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
Bool_t AliAODPidHF::CheckStatus(AliAODTrack *track,TString detectors){


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
