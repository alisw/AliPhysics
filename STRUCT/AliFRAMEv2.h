#ifndef ALIFRAMEV2_H
#define ALIFRAMEV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id */

/////////////////////////////////////////////////////////
//  Manager and class for detector: FRAME  version 2    //
/////////////////////////////////////////////////////////
 
#include "AliFRAME.h"

class AliFRAMEv2 : public AliFRAME {
  
public:
  AliFRAMEv2();
  AliFRAMEv2(const char *name, const char *title);
  virtual       ~AliFRAMEv2() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const;
  virtual void   SetHoles(Int_t flag=0) {fHoles = flag;}
  virtual Int_t  Holes() {return fHoles;}
 public:
  Int_t  fHoles;
  
   ClassDef(AliFRAMEv2,2)  //Class for FRAME version 2
};
 
#endif
