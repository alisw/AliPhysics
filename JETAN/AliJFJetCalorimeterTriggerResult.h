// $Id$

#ifndef ALIJFJETCALORIMETERTRIGGERRESULTH
#define ALIJFJETCALORIMETERTRIGGERRESULTH

#include "AliJFJetTriggerResult.h"

class AliJFJetCalorimeterTriggerResult : public AliJFJetTriggerResult
{
 public:
  AliJFJetCalorimeterTriggerResult();
  AliJFJetCalorimeterTriggerResult(Float_t PtMin,Float_t PtMax,Float_t PhiMin,Float_t PhiMax,Float_t EtaMin,Float_t EtaMax,Bool_t Neutral,Bool_t EM,Bool_t Charged,Float_t CPhi,Float_t CEta,Float_t CEt,Float_t CThr,Float_t Inner,Float_t Outer,Float_t Rin,Float_t Rmid,Float_t Rout,Int_t PhiBins,Int_t EtaBins);

  virtual ~AliJFJetCalorimeterTriggerResult();

  Float_t GetCenterPhi() const {return fCenterPhi;}
  Float_t GetCenterEta() const {return fCenterEta;}
  Float_t GetCenterEt()  const {return fCenterEt;}
  Float_t GetThreshEt()  const {return fThreshEt;}
  Float_t GetInnerEt()   const {return fInnerEt;}
  Float_t GetOuterEt()   const {return fOuterEt;}
  Float_t GetRin()       const {return fRin;}
  Float_t GetRmid()      const {return fRmid;}
  Float_t GetRout()      const {return fRout;}
  Int_t GetPhiBins()     const {return fPhiBins;}
  Int_t GetEtaBins()     const {return fEtaBins;}

 protected:
  Float_t fCenterPhi;
  Float_t fCenterEta;
  Float_t fCenterEt;
  Float_t fThreshEt;
  Float_t fInnerEt;
  Float_t fOuterEt;

  Float_t fRin;
  Float_t fRmid;
  Float_t fRout;

  Int_t fPhiBins;
  Int_t fEtaBins;


  ClassDef(AliJFJetCalorimeterTriggerResult,1) //AliJFJetCalorimeterTriggerResult class
};

#endif /*ALIJFJETCALORIMETERTRIGGERRESULTH*/
