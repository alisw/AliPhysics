// $Id$

#ifndef ALITKTOWERV2_H
#define ALITKTOWERV2_H

#include <TParticle.h>
#include <TClonesArray.h>

class AliTkTowerV2 : public TObject {
 public:
  AliTkTowerV2();
  AliTkTowerV2(const AliTkTowerV2 &t);
  ~AliTkTowerV2();

  Float_t getEta()       const {return fEta;}
  Float_t getPhi()       const {return fPhi;}
  Float_t getEt()        const {return fEt;}
  Int_t getNParticles()  const {return fNParticles;}
  Bool_t getUpdate()     const {return fUpdate;}
  Float_t getEtCharged() {if(fUpdate) update(); return fEtCharged;}
  Float_t getEtEM()      {if(fUpdate) update(); return fEtEM;};

  TClonesArray *getParticleList() const {return fParticles;}
  TClonesArray *getChargedParticleList() const;
  TClonesArray *getNeutralParticleList() const;

  void setEta(Float_t Eta) {fEta = Eta;}
  void setPhi(Float_t Phi) {fPhi = Phi;}
  void setEt(Float_t Et)   {fEt = Et;}
  void addParticle(const TParticle *particle);
  void setParticleList(TClonesArray *ptr);

  void Print(Option_t *) const {
    cout << "AliTkTower " << fEt << " " << fEta << " " << fPhi << endl;
  }
  ULong_t Hash() const {return 0;}
  Bool_t IsEqual(const TObject */*obj*/) const {return kFALSE;}
  Bool_t IsSortable() const {return kTRUE;}
  Int_t  Compare(const TObject *obj) const;
  void Clear(Option_t *option="");

 private:
  Bool_t isChargedParticle(const TParticle *particle) const;
  Bool_t isEMParticle(const TParticle *particle) const;
  void update();

  // member variables
  Float_t fEta;
  Float_t fPhi;
  Float_t fEt;
  Int_t fNParticles;
  Float_t fEtCharged;
  Float_t fEtEM;
  Bool_t fUpdate;
  TClonesArray *fParticles; //->

  ClassDef(AliTkTowerV2,2) //AliTkTower class
};
#endif
