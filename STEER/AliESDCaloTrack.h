#ifndef ALIESDCALOTRACK_H
#define ALIESDCALOTRACK_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//   Class AliESDCaloTrack
//   This is the class to deal with during the physical analysis of data
//   It converts calorimeter (PHOS or EMCAL) reconstructed particles   
//   into event summary data object
//-------------------------------------------------------------------------

#include "TObject.h"
class AliPHOSRecParticle;

class AliESDCaloTrack : public TObject {

public:
  AliESDCaloTrack() {}
  virtual ~AliESDCaloTrack() {}
  AliESDCaloTrack(AliPHOSRecParticle* recpart);
  Float_t Px() { return fPx; }
  Float_t Py() { return fPy; }
  Float_t Pz() { return fPz; }

private:
  Float_t fPx; // x-component of PHOS rec.particle
  Float_t fPy; // y-component of PHOS rec.particle
  Float_t fPz; // z-component of PHOS rec.particle

  ClassDef(AliESDCaloTrack,1)  //ESD calorimeter track class 
};

#endif 
