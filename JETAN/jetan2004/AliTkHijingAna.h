// $Id$

#ifndef ALITKHIJINGANA_H
#define ALITKHIJINGANA_H

#include <Riostream.h>
#include <TObject.h>
#include "AliTkEtaPhiVector.h"

class AliTkHijingAna : public TObject {
 public:
  // how to use it: - call setCenter
  //                - call setParticles at beginning (not mandatory)
  //                - call getParticlesInRadius with maximum radius
  //                - call analysis functions
  //                - repeat last 2 steps with decreasing radia

  // center of cone
  void setCenter(Float_t eta, Float_t phi);
  Float_t getCenterEta() const;
  Float_t getCenterPhi() const;

  // helper functions
  Bool_t isParticleInRadius(TParticle *particle,Float_t radius);
  Bool_t isParticleAccepted(TParticle *particle, Float_t pt=0.1, Bool_t neutral=kFALSE);
  Bool_t isParticleAcceptedALICE(TParticle *particle, Float_t pt=0.1, Bool_t neutral=kFALSE);
  void clear();
  void clearCone();

  // functions to create local copies of particles
  void setParticles(TClonesArray *particles);
  void setParticles(TClonesArray *particles,Float_t pt,Bool_t neutral);
  TClonesArray *getParticlesInRadius(Float_t radius);
  TClonesArray *getParticlesInRadius(TClonesArray *particles, Float_t radius);

  // analysis functions
  Int_t getNParticlesInRadius(Float_t ptCut);
  Int_t getNParticlesInRadius() { return getNParticlesInRadius(0); }
  Float_t getEtInRadius(Float_t ptCut);
  Float_t getEtInRadius() { return getEtInRadius(0); }


  //--- old stuff, don't use anymore
  Float_t getEtInRadius(TClonesArray *particles,Float_t radius);
  //--- end of old stuff

 private:
  AliTkEtaPhiVector center;
  TClonesArray *oParticles;
  TClonesArray *mParticles;
  Float_t mRadius;
  Bool_t updated;

  ClassDef(AliTkHijingAna,0)
};

#endif
