#include "AliESDEvent.h"
#include "AliITSPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"
#include "AliPIDtools.h"

std::map<Int_t, AliTPCPIDResponse *> AliPIDtools::pidTPC;     /// we should use better hash map
std::map<Int_t, AliPIDResponse *> AliPIDtools::pidAll;        /// we should use better hash map
AliESDtrack  AliPIDtools::dummyTrack;/// dummy value to save CPU - unfortunately PID object use AliVtrack - for the moment create global variable t avoid object constructions

AliPIDResponse* AliPIDtools::GetPID(Int_t hash ) {return pidAll[hash];}
AliTPCPIDResponse& AliPIDtools::GetTPCPID(Int_t hash ) {return pidAll[hash]->GetTPCResponse();}
AliITSPIDResponse& AliPIDtools::GetITSPID(Int_t hash ) {return pidAll[hash]->GetITSResponse();}
AliTOFPIDResponse& AliPIDtools::GetTOFPID(Int_t hash ) {return pidAll[hash]->GetTOFResponse();}

Int_t AliPIDtools::GetHash(Int_t run, Int_t passNumber, TString recoPass,Bool_t isMC){
  recoPass+=run;
  recoPass+=passNumber;
  recoPass+=isMC;
  return recoPass.Hash();
}

Double_t AliPIDtools::BetheBlochAleph(Int_t hash, Double_t bg){
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  if (tpcPID) return tpcPID->Bethe(bg);
  return 0;
}

///  AliPIDtools::BetheBlochITS
/// \param hash   - hash value
/// \param p      - momentum (where?)
/// \param mass   - mass
/// \return
Double_t AliPIDtools::BetheBlochITS(Int_t hash, Double_t p, Double_t mass){
  if (pidAll[hash]== nullptr) return 0;
  AliITSPIDResponse &itsPID=GetITSPID(hash);
  return itsPID.Bethe(p, mass);
}

/// AliPIDtools::GetExpectedITSSignal(
/// \param hash   - hash value
/// \param p      - momentum (where?)
/// \return
Double_t AliPIDtools::GetExpectedITSSignal(Int_t hash, Double_t p, Int_t  particle){
  if (pidAll[hash]== nullptr) return 0;
  AliITSPIDResponse &itsPID=GetITSPID(hash);
  return itsPID.Bethe(p, (AliPID::EParticleType)particle);
}


/// GetExpectedTPCSignal
/// \param hash       - hash value of the PID version
/// \param p          - momenta
/// \param particle   - particle type
/// \return           - mean TPCdedx
Double_t AliPIDtools::GetExpectedTPCSignal(Int_t hash, Double_t p, Int_t  particle) {
  Double_t xyz[3] = {0., 0., 0.};
  Double_t pxyz[3] = {0, 0., 0.};
  Double_t cv[21] = {0.}; // dummy parameters for dummy tracks
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  if (tpcPID==0) return 0;
  pxyz[0]=p;
  dummyTrack.Set(xyz, pxyz, cv, 1);
  Double_t dEdx = tpcPID->GetExpectedSignal(&dummyTrack, (AliPID::EParticleType)particle, AliTPCPIDResponse::kdEdxDefault, kFALSE, kTRUE);
  return dEdx;
}
/// Load and reguster PID objects in hash maps
/// \param run
/// \param passNumber
/// \param recoPass
/// \param isMC
/// \return  - hash value of PID
Int_t AliPIDtools::LoadPID(Int_t run, Int_t passNumber, TString recoPass, Bool_t isMC){
  // Int_t run=246751, passNumber=1; TString recoPass("pass1"); Bool_t isMC=0;
  AliESDEvent ev;
  AliPIDResponse *pid = new AliPIDResponse(isMC);
  pid->SetUseTPCMultiplicityCorrection();
  pid->SetOADBPath("$ALICE_PHYSICS/OADB");
  pid->InitialiseEvent(&ev,passNumber, recoPass, run);
  AliTPCPIDResponse &tpcpid=pid->GetTPCResponse();
  // pid.InitFromOADB(246751,1,"pass1");
  Int_t  hash=GetHash(run,passNumber, recoPass,isMC);
  pidAll[hash]=pid;     /// we should clone them
  pidTPC[hash]=&tpcpid;  ///
  return hash;
}

Double_t AliPIDtools::GetExpectedTOFSigma(Int_t hash, Float_t mom, Int_t  type){
  Double_t dummyTime=0;
  if (pidAll[hash]== nullptr) return 0;
  AliTOFPIDResponse &tofPID=GetTOFPID(hash);
  return tofPID.GetExpectedSigma(mom,dummyTime,(AliPID::EParticleType)type);

}
Double_t AliPIDtools::GetExpectedTOFSignal(Int_t hash, const AliVTrack *track, Int_t type){
  if (pidAll[hash]== nullptr) return 0;
  AliTOFPIDResponse &tofPID=GetTOFPID(hash);
  return tofPID.GetExpectedSignal(track, (AliPID::EParticleType)type);
}

