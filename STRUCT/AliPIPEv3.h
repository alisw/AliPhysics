#ifndef ALIPIPEV3_H
#define ALIPIPEV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: PIPE          //
////////////////////////////////////////////////
 
#include "AliPIPE.h"
 
 
class AliPIPEv3 : public AliPIPE {
 
public:
  AliPIPEv3();
  AliPIPEv3(const char *name, const char *title);
  virtual      ~AliPIPEv3() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 3;}
  virtual void  DrawModule();
  virtual void  Undulation(char *undul, Float_t pitch, Float_t thick, Float_t zundul, Float_t rundul,
                           char (*cone)[5]);
  ClassDef(AliPIPEv3,1)  //Class for PIPE version 3
};

#endif
