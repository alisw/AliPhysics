#ifndef ALIESDCALOTRACK_H
#define ALIESDCALOTRACK_H

#include "TObject.h"

class AliPHOSRecParticle;

class AliESDCaloTrack : public TObject {
 public:
  AliESDCaloTrack(AliPHOSRecParticle * part = 0x0) {
    fPHOSRecParticle = part;
  }
  virtual ~AliESDCaloTrack (){}

 private:
  AliPHOSRecParticle * fPHOSRecParticle; // pointer to PHOS particle
  ClassDef(AliESDCaloTrack,1)  //ESD calorimeter track class 
};

#endif 
