#ifndef PIPEv1_H
#define PIPEv1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: PIPE          //
////////////////////////////////////////////////
 
#include "AliPIPE.h"
 
 
class AliPIPEv1 : public AliPIPE {
 
public:
  AliPIPEv1();
  AliPIPEv1(const char *name, const char *title);
  virtual      ~AliPIPEv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  DrawModule();
  
  ClassDef(AliPIPEv1,1)  //Class for PIPE version 1
};

#endif
