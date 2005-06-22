// $Id$

#ifndef ALITKCONEJET_H
#define ALITKCONEJET_H

#include "AliTkTowerV2.h"

class AliTkConeJet : public TObject {
 public:
  AliTkConeJet();
  AliTkConeJet(Float_t pt,Float_t eta,Float_t phi,Int_t type=0);
  AliTkConeJet(const AliTkConeJet &j);
  ~AliTkConeJet();

  void addTower(AliTkTowerV2 *tower);

  Float_t getEtEM() const;
  Float_t getEtMarked(Float_t ptCut=0.) const;
  Float_t getEtMarkedFrac(Float_t ptCut=0.) const;
  Float_t getEtCharged() const;
  Float_t getEtChargedMarked(Float_t ptCut=0.) const;
  Float_t getE(Float_t ptCut=0.) const;
  Float_t getECharged(Float_t ptCut=0.) const;

  TClonesArray *getParticles(Float_t ptCut=0.) const;
  TClonesArray *getChargedParticles(Float_t ptCut=0.) const;
  TClonesArray *getNeutralParticles(Float_t ptCut=0.) const;

  Int_t getNParticles() const;
  Int_t getNChargedParticles() const;
  Int_t getNNeutralParticles() const;
  Int_t getNParticles(Float_t ptCut) const;
  Int_t getNChargedParticles(Float_t ptCut) const;
  Int_t getNNeutralParticles(Float_t ptCut) const;

  void getAxis(Float_t &x,Float_t &y,Float_t &z,Float_t ptcut=0.) const;
  void getChAxis(Float_t &x,Float_t &y,Float_t &z,Float_t ptcut=0.) const;
  TParticle *getLeadingPart(Float_t ptcut=0.) const;
  TParticle *getLeadingChPart(Float_t ptcut=0.) const;

  //get Et, eta and phi from particles in towers
  void calculateFromParticles(Float_t &et, Float_t &eta, Float_t &phi, Float_t ptcut=0.); 

  //calculate interesting values
  void calculateValues(Float_t ptcut=0.);

  void Print(Option_t *) const {
    cout << "AliTkConeJet " << fEt << " " << fEta << " " << fPhi << endl;
  }

  ULong_t Hash() const {return 0;}
  Bool_t IsEqual(const TObject */*obj*/) const {return kFALSE;}
  Bool_t IsSortable() const {return kTRUE;}
  Int_t  Compare(const TObject *obj) const;
  void Clear(Option_t *option="");

  void setEta(Float_t eta) {fEta = eta;}
  void setPhi(Float_t phi) {fPhi = phi;}
  void setEt(Float_t et)   {fEt = et;}
  void setType(Int_t t)    {fType = t;}
  Float_t getEta()   const {return fEta;}
  Float_t getPhi()   const {return fPhi;}
  Float_t getEt()    const {return fEt;}
  Int_t getType()    const {return fType;}
  Int_t getNTowers() const {return fNTowers;}
  Float_t getPtCut()       const {return fPtCut;}
  Float_t getCEta()        const {return fCEta;}
  Float_t getCPhi()        const {return fCPhi;}
  Float_t getCEt()         const {return fCEt;}
  Float_t getPLength()     const {return fPLength;}
  Float_t getXAxis()       const {return fXAxis;}
  Float_t getYAxis()       const {return fYAxis;}
  Float_t getZAxis()       const {return fZAxis;}
  Float_t getPtLead()      const {return fPtLead;}
  TParticle* getLeadPart() const {return fLeadPart;}
  Int_t getNParts()        const {return fNParts;}
  TClonesArray *getParts() const {return fParts;}
  TClonesArray *getTowerList() const {return fTowers;}

  static Float_t Diff(const AliTkConeJet *jet1, const AliTkConeJet *jet2, Float_t &etdiff, Float_t &phidiff, Float_t &etadiff);

 private:
  Float_t fEta;
  Float_t fPhi;
  Float_t fEt;
  Int_t fType;    //indicate quark, gluon or other information

  Int_t fNTowers; //number of towers
  TClonesArray *fTowers; //->

  /* properties here are calculated for particles
     regardless of their charge */
  Float_t fPtCut; //pt cut at which properties where calculated
  Float_t fCEta;  //calculated from particles
  Float_t fCPhi;
  Float_t fCEt;
  Float_t fPLength; //length of axis
  Float_t fXAxis;   //axis
  Float_t fYAxis;
  Float_t fZAxis;
  Float_t fPtLead;  //leading particle pt
  TParticle *fLeadPart; //->
  Int_t fNParts;    //number of parts in all towers
  TClonesArray *fParts; //->

  ClassDef(AliTkConeJet,3) //AliTkConeJet class
};

ostream& operator<<(ostream& s,const AliTkConeJet& t);
#endif
