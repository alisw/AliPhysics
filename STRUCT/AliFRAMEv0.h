#ifndef FRAMEv0_H
#define FRAMEv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and class for detector: FRAME  version 0    //
/////////////////////////////////////////////////////////
 
#include "AliFRAME.h"

class AliFRAMEv0 : public AliFRAME {
  
public:
  AliFRAMEv0();
  AliFRAMEv0(const char *name, const char *title);
  virtual       ~AliFRAMEv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   Init();
  virtual Int_t  IsVersion() const {return 0;}
   
   ClassDef(AliFRAMEv0,1)  //Class for FRAME version 0
};
 
#endif
