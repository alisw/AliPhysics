#include "AliESEHelpers.h"
#include "AliAODEvent.h"
#include "AliNanoAODTrackMapping.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVParticle.h"
#include "AliPIDResponse.h"
#include "AliNanoAODTrack.h"
#include "AliNanoAODHeader.h"

ClassImp(AliESEEvtCut)
ClassImp(AliESETrkCut)
ClassImp(AliAnalysisESESetter)


void AliAnalysisESESetter::SetNanoAODHeader(const AliAODEvent */* event*/   , AliNanoAODHeader * head  ) {

  // if(fHelperPID->GetPIDType()==kBayes)fHelperPID->GetPIDCombined()->SetDefaultTPCPriors();//FIXME maybe this can go in the UserCreateOutputObject?

  //  if(!isMC){
  Double_t QvecA=fEventCuts->GetqV0A();
  Double_t QvecC=fEventCuts->GetqV0C();
    //  }

  Double_t Cent = fEventCuts->GetCent();
  // TODO: MAPPING FOR HEADERS!
  head->SetVar(0, QvecA);
  head->SetVar(1, QvecC);
  head->SetVar(2, Cent);


}

void AliAnalysisESESetter::SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * spTrack) {
  // Set custom variables in the special track
  // 1. Cache the indexes
  
  static  Int_t kcstNSigmaTPCPi  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTPCPi");
  static  Int_t kcstNSigmaTPCKa  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTPCKa");
  static  Int_t kcstNSigmaTPCPr  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTPCPr");

  static  Int_t kcstNSigmaTOFPi  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTOFPi");
  static  Int_t kcstNSigmaTOFKa  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTOFKa");
  static  Int_t kcstNSigmaTOFPr  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstNSigmaTOFPr");

  // TODO: set Bayes vars in special track

  // static const Int_t kcstBayesTPCPi  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("kcstBayesTPCPi");
  // static const Int_t kcstBayesTPCKa  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("kcstBayesTPCKa");
  // static const Int_t kcstBayesTPCPr  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("kcstBayesTPCPr");

  // static const Int_t kcstBayesTOFPi  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstBayesTOFPi");
  // static const Int_t kcstBayesTOFKa  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstBayesTOFKa");
  // static const Int_t kcstBayesTOFPr  = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstBayesTOFPr");

  // Get the PID info
  static AliPIDResponse * pidResponse = 0;
  if(!pidResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    pidResponse = inputHandler->GetPIDResponse();
  }

  const AliVParticle *inEvHMain = dynamic_cast<const AliVParticle *>(aodTrack);

  Double_t nsigmaTPCkProton = pidResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton);
  Double_t nsigmaTPCkKaon   = pidResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon); 
  Double_t nsigmaTPCkPion   = pidResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion); 

  Double_t nsigmaTOFkProton = pidResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton);
  Double_t nsigmaTOFkKaon   = pidResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon); 
  Double_t nsigmaTOFkPion   = pidResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion); 



  

  spTrack->SetVar(kcstNSigmaTPCPi, nsigmaTPCkPion);
  spTrack->SetVar(kcstNSigmaTPCKa, nsigmaTPCkKaon);
  spTrack->SetVar(kcstNSigmaTPCPr, nsigmaTPCkProton);

  spTrack->SetVar(kcstNSigmaTOFPi, nsigmaTOFkPion);
  spTrack->SetVar(kcstNSigmaTOFKa, nsigmaTOFkKaon);
  spTrack->SetVar(kcstNSigmaTOFPr, nsigmaTOFkProton);
  //TODO: set the bayes vars
  

}
