#include "AliESDEvent.h"
#include "AliTPCPIDResponse.h"
#include "AliPIDResponse.h"
#include "AliESDtrack.h"
#include "AliPIDtools.h"

std::map<Int_t, AliTPCPIDResponse *> AliPIDtools::pidTPC;     /// we should use better hash map
std::map<Int_t, AliPIDResponse *> AliPIDtools::pidAll;        /// we should use better hash map
AliESDtrack  AliPIDtools::dummyTrack;/// dummy value to save CPU - unfortunately PID object use AliVtrack - for the moment create global varaible t avoid object constructions

AliTPCPIDResponse* AliPIDtools::GetTPCPID(Int_t hash ) {return pidTPC[hash];}
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
/// GetExpectedTPCSignal
/// \param hash       - hash value of the PID version
/// \param p          - momenta
/// \param particle   - particle type
/// \return           - mean TPCdedx
Double_t AliPIDtools::GetExpectedTPCSignal(Int_t hash, Double_t p, AliPID::EParticleType particle) {
  Double_t xyz[3] = {0., 0., 0.};
  Double_t pxyz[3] = {0, 0., 0.};
  Double_t cv[21] = {0.}; // dummy parameters for dummy tracks
  AliTPCPIDResponse *tpcPID=pidTPC[hash];
  if (tpcPID==0) return 0;
  pxyz[0]=p;
  dummyTrack.Set(xyz, pxyz, cv, 1);
  Double_t dEdx = tpcPID->GetExpectedSignal(&dummyTrack, particle, AliTPCPIDResponse::kdEdxDefault, kFALSE, kTRUE);
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
