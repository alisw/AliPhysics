// $Id$

#ifndef ALIJFJETFINDERH
#define ALIJFJETFINDERH

#include <TObject.h>

class TCollection;
class TParticle;
class JFTower;

class AliJFJetFinder
{
 public:
  AliJFJetFinder(Int_t n=50);
  virtual ~AliJFJetFinder();

  virtual Int_t Init(TClonesArray */*particles*/){return 0;}
  virtual Int_t Run(){return 0;}

  inline TObjArray*  GetJets()  {return &fJets;}
  inline Int_t const GetNJets() const {return fNJets;}

  virtual Bool_t IsAcceptedParticle(TParticle *p);
  virtual Bool_t IsAcceptedTower(JFTower*){return kFALSE;}

  virtual void Debug(){};
  virtual void Clean(){fJets.Delete();fNJets=0;}

  void SetPtCut(Float_t ptmin=0, Float_t ptmax=1000);
  void SetPhiCut(Float_t phi=6.4){SetPhiCut(0,phi);}
  void SetPhiCut(Float_t phimin, Float_t phimax);
  void SetEtaCut(Float_t e=1){SetEtaCut(-e,e);}
  void SetEtaCut(Float_t emin, Float_t emax);
  void SetNeutral(Bool_t b=kTRUE){fNeutral=b;}
  void SetCharged(Bool_t b=kTRUE){fCharged=b;}
  void SetEM(Bool_t b=kTRUE){fEM=b;}

 protected:
  Int_t fNJets;
  Int_t fNJetsMax;
  TObjArray fJets;

  Float_t fPtMin;
  Float_t fPtMax;
  Float_t fEtaMin;
  Float_t fEtaMax;
  Float_t fPhiMin;
  Float_t fPhiMax;
  Bool_t fNeutral;
  Bool_t fCharged;
  Bool_t fEM;

  ClassDef(AliJFJetFinder,1) //AliJFJetFinder class
};

#endif /*ALIJFJETFINDERH*/
