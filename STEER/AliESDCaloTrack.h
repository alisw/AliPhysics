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
#include "TParticle.h"

class AliESDCaloTrack : public TObject {

public:
  AliESDCaloTrack() {}
  virtual ~AliESDCaloTrack() {}
  AliESDCaloTrack(TParticle* recpart);
  Float_t Px() { return fRecParticle->Px(); }
  Float_t Py() { return fRecParticle->Py(); }
  Float_t Pz() { return fRecParticle->Pz(); }

private:
  TParticle *fRecParticle; // reconstructed particle from PHOS or EMCAL

  ClassDef(AliESDCaloTrack,2)  //ESD calorimeter track class 
};

#endif 
