#ifndef ALIFRAMEV1_H
#define ALIFRAMEV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and class for detector: FRAME  version 0    //
/////////////////////////////////////////////////////////
 
#include "AliFRAME.h"

class AliFRAMEv1 : public AliFRAME {
  
public:
  AliFRAMEv1();
  AliFRAMEv1(const char *name, const char *title);
  virtual       ~AliFRAMEv1() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 1;}
   
   ClassDef(AliFRAMEv1,1)  //Class for FRAME version 0
};
 
#endif
