#ifndef ALIBODY_H
#define ALIBODY_H
/* Copyright(c) 1998-2005, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////
//  Manager class for detector:  AliEMCALWsuCosmicRaySetUp        //
//   This is the envelop for Alice                                //
///////////////////////////////////////////////////////////////////
 
#include "AliModule.h"

class AliEMCALWsuCosmicRaySetUp : public AliModule {
 
public:
  AliEMCALWsuCosmicRaySetUp();
  AliEMCALWsuCosmicRaySetUp(const char *name, const char *title);
  virtual     ~AliEMCALWsuCosmicRaySetUp() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 0;}
  void  DrawWSUC(float cxy=0.025) const; // *MENU*

  ClassDef(AliEMCALWsuCosmicRaySetUp,1)  // Class manager for the Wsu Cosmic Ray SetUp
};

#endif
