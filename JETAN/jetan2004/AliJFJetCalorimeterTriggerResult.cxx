// $Id$

#include <Riostream.h>
#include "AliJFJetCalorimeterTriggerResult.h"

ClassImp(AliJFJetCalorimeterTriggerResult)

AliJFJetCalorimeterTriggerResult::AliJFJetCalorimeterTriggerResult() 
                             : AliJFJetTriggerResult(),
                               fCenterPhi(0),fCenterEta(0),fCenterEt(0),fThreshEt(0),
			       fInnerEt(0),fOuterEt(0),fRin(0),fRmid(0),
			       fRout(0),fPhiBins(0),fEtaBins(0)
{
}

AliJFJetCalorimeterTriggerResult::AliJFJetCalorimeterTriggerResult(Float_t PtMin,Float_t PtMax,Float_t PhiMin,Float_t PhiMax,Float_t EtaMin,Float_t EtaMax,Bool_t Neutral,Bool_t EM,Bool_t Charged,Float_t CPhi,Float_t CEta,Float_t CEt,Float_t CThr,Float_t Inner,Float_t Outer,Float_t Rin,Float_t Rmid,Float_t Rout,Int_t PhiBins,Int_t EtaBins) 
                             : AliJFJetTriggerResult(PtMin,PtMax,PhiMin,PhiMax,EtaMin,EtaMax,Neutral,EM,Charged),
			       fCenterPhi(CPhi),fCenterEta(CEta),fCenterEt(CEt),fThreshEt(CThr),
			       fInnerEt(Inner),fOuterEt(Outer),fRin(Rin),fRmid(Rmid),
			       fRout(Rout),fPhiBins(PhiBins),fEtaBins(EtaBins)
{
}


AliJFJetCalorimeterTriggerResult::~AliJFJetCalorimeterTriggerResult()
{
}
