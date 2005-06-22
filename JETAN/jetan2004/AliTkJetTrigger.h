// $Id$

/****************************************************************************
 * AliTkJetTrigger.h                                                           *
 * Test version of jet-trigger, includes many analysis functions to tune    *
 * paramters                                                                *
 * Thorsten Kollegger <kollegge@ikf.physik.uni-frankfurt.de>                *
 ***************************************************************************/

#ifndef TKJETTRIGGER_H
#define TKJETTRIGGER_H

#include <Riostream.h>
#include <list>
#include <TParticle.h>
#include "AliTkJetTriggerEvent.h"

class AliTkJetTrigger : public TObject {
 public:
  AliTkJetTrigger() : TObject() {  
    mEvOutName=0;
    setDefaults();
  }
  virtual ~AliTkJetTrigger(){if(mEvOutName) delete[] mEvOutName;}

  void init();
  void initEvent();
  void initEvent(TClonesArray *particles, Int_t type);
  void initEvent(TClonesArray *particles);
  void make();
  void run() { make(); }
  void finishEvent();
  void finish();

  void setDefaults();
  void setEvOutFilename(Char_t *filename);
  Char_t *getEvOutFilename() { return mEvOutName; }

  void setParticleType(Int_t type);
  Int_t getParticleType() { return mParticleType; }
  void addConeRadius(Float_t radius);
  void addConeNParticles(UInt_t nParticles);
  void addConePtThreshold(Float_t ptThreshold);

  void mixParticles(TClonesArray *particles, Int_t type);

  void setParticleMinPt(Float_t f){fParticleMinPt=f;}
  void setParticleMinPtSeed(Float_t f){fParticleMinPtSeed=f;}

  void setExpectedParticleNumber(UInt_t n) {mExpectedNParticles = n; }
  UInt_t getExpectedParticleNumber() { return mExpectedNParticles; }

 private:

  void createSeedPoints();
  
  Bool_t isParticleAccepted(TParticle *particle);
  Bool_t isParticleAccepted(TParticle *particle, Int_t type);
  Bool_t isParticleAcceptedPythia(TParticle *particle);
  Bool_t isParticleAcceptedPythiaAliceGeo(TParticle *particle);
  Bool_t isParticleAcceptedHijing(TParticle *particle);
  Bool_t isParticleAcceptedHijingAliceGeo(TParticle *particle);

  void clearParticles();

  list<TParticle*> *getParticlesInCone(AliTkEtaPhiVector center, Float_t radius);
  list<TParticle*> *getParticlesInCone(AliTkEtaPhiVector center,
				      Float_t radius,
				      TClonesArray *particles);
  Bool_t decide(list<TParticle*> *particles, Float_t ptThreshold,
		UInt_t nParticles);

  Int_t mParticleType;
  TClonesArray *mParticles;
  UInt_t mExpectedNParticles;
  
  list<AliTkEtaPhiVector> seedPoints;
  list<Float_t> mConeRadia;
  list<UInt_t> mConeNParticles;
  list<Float_t> mConePtThreshold;

  Char_t *mEvOutName; //!
  TFile *mEvOutFile; //!
  TTree *mEvOutTree; //!
  AliTkJetTriggerEvent *mEvOutEvent; //!
  TObjArray *mDecisions;

  Float_t fParticleMinPt;
  Float_t fParticleMinPtSeed;

  ClassDef(AliTkJetTrigger,0)
};


#endif
