// $Id$

#ifndef ALIJFJETCALORIMETERTRIGGERH
#define ALIJFJETCALORIMETERTRIGGERH

#include "AliJFJetTrigger.h"
#include <TObject.h>

class TH2F;

class AliJFJetCalorimeterTrigger : public AliJFJetTrigger
{
 public:
  AliJFJetCalorimeterTrigger(Int_t n=50);
  virtual ~AliJFJetCalorimeterTrigger();

  Int_t Init(TClonesArray *particles);
  Int_t Run();

  void Clean();
  void Debug(){};

  void SetRadius(Float_t in=0.2,Float_t mid=0.7,Float_t out=1.0);
  void SetBins(Int_t phis=25,Int_t etas=25);

 protected:
  Float_t fRin;
  Float_t fRmid;
  Float_t fRout;
  Float_t fAin;
  Float_t fAout;

  Int_t fPhiBins;
  Int_t fEtaBins;

  TH2F *fSeedPlane; //!

  ClassDef(AliJFJetCalorimeterTrigger,1) //AliJFJetCalorimeterTrigger class
};

#endif /*ALIJFJETCALORIMETERTRIGGERH*/
