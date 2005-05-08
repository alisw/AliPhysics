// $Id$

#ifndef ALIJFPARTICLESCUTH
#define ALIJFPARTICLESCUTH

#include <TMath.h>

class TParticle;
class TClonesArray;
class TH2F;

class AliJFParticlesCut
{
 public:
  AliJFParticlesCut(TClonesArray *p=0);
  virtual ~AliJFParticlesCut(){};

  void SetPtCut(Float_t ptmin=0, Float_t ptmax=1000);
  void SetPhiCut(Float_t phi=2*TMath::Pi()){SetPhiCut(0,phi);}
  void SetPhiCut(Float_t phimin, Float_t phimax);
  void SetEtaCut(Float_t e=1){SetEtaCut(-e,e);}
  void SetEtaCut(Float_t emin, Float_t emax);
  void SetNeutral(Bool_t b=kTRUE){fNeutral=b;}
  void SetCharged(Bool_t b=kTRUE){fCharged=b;}
  void SetEM(Bool_t b=kTRUE){fEM=b;}

  void SetParticles(TClonesArray *p){fParts=p;}
  TClonesArray* GetParticles(TClonesArray */*p*/){return fParts;}

  Int_t Cut();
  Int_t Cut(TClonesArray *p);

  TH2F* CreateHistogram(Char_t *title,Char_t *text,Int_t phibins=100,Int_t etabins=100);

  Bool_t IsAcceptedParticle(TParticle *p);
 
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

  TClonesArray *fParts;

  ClassDef(AliJFParticlesCut,1) //AliJFParticlesCut class

};

#endif /*ALIJFPARTICLESCUTH*/
