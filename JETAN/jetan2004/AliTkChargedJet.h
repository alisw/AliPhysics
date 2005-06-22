// $Id$

#ifndef ALITKCHARGEDJET_H
#define ALITKCHARGEDJET_H

#include <Riostream.h>
#include <list>
#include <TObject.h>
#include "AliTkChargedJetFinder.h"

class TClonesArray;
class TParticle;

class AliTkChargedJet : public TObject {

 public:
  AliTkChargedJet();
  AliTkChargedJet(jet);
  ~AliTkChargedJet();

  Float_t getEta() const      { return fEta; }
  Float_t getPhi() const      { return fPhi; }
  Float_t getPt()  const      { return fPt; }
  Int_t getNParticles() const { return fNParticles; }

  TClonesArray *getParticles() const { return fParticles; }

  Float_t getPtInRadius(Float_t r)      const;
  Int_t getParticlesInRadius(Float_t r) const;
  Int_t getNChargedInRadius(Float_t r)  const;
  Int_t getNNeutralInRadius(Float_t r)  const;

 private:

  Float_t fEta;
  Float_t fPhi;
  Float_t fPt;
  Int_t fNParticles;
  TClonesArray *fParticles; //->

  ClassDef(AliTkChargedJet,1)
};
#endif
