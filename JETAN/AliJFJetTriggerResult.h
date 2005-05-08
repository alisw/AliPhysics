// $Id$

#ifndef ALIJFJETTRIGGERRESULTH
#define ALIJFJETTRIGGERRESULTH

#ifndef ROOT_TObject 
#include <TObject.h>
#endif


class AliJFJetTriggerResult : public TObject
{
 public:
  AliJFJetTriggerResult();
  AliJFJetTriggerResult(Float_t PtMin,Float_t PtMax,Float_t EtaMin,Float_t EtaMax,Float_t PhiMin,Float_t PhiMax,Bool_t Neutral,Bool_t EM,Bool_t Charged);
  virtual ~AliJFJetTriggerResult();

  Float_t GetPtMin()  const {return fPtMin;}
  Float_t GetPtMax()  const {return fPtMax;}
  Float_t GetEtaMin() const {return fEtaMin;}
  Float_t GetEtaMax() const {return fEtaMax;}
  Float_t GetPhiMin() const {return fPhiMin;}
  Float_t GetPhiMax() const {return fPhiMax;}
  Bool_t IsNeutral()  const {return fNeutral;}
  Bool_t IsCharged()  const {return fCharged;}
  Bool_t IsEM()       const {return fEM;}


 protected:

  Float_t fPtMin;
  Float_t fPtMax;
  Float_t fEtaMin;
  Float_t fEtaMax;
  Float_t fPhiMin;
  Float_t fPhiMax;
  Bool_t fNeutral;
  Bool_t fCharged;
  Bool_t fEM;

  ClassDef(AliJFJetTriggerResult,1) //AliJFJetTriggerResult class
};

#endif /*ALIJFJETTRIGGERRESULTH*/
