#ifndef ALIPIPEV0_H
#define ALIPIPEV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and class for detector: PIPE  version 0    //
/////////////////////////////////////////////////////////
 
#include "AliPIPE.h"

class AliPIPEv0 : public AliPIPE {
  
public:
  AliPIPEv0();
  AliPIPEv0(const char *name, const char *title);
  virtual       ~AliPIPEv0() {}
  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const {return 0;}
  virtual void  DrawModule();
   
   ClassDef(AliPIPEv0,1)  //Class for PIPE version 0
};
 
#endif
