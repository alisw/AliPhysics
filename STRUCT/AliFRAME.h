#ifndef FRAME_H
#define FRAME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: FRAME         //
////////////////////////////////////////////////
 
#include "AliModule.h"


class AliFRAME : public AliModule {
  
public:
  AliFRAME();
  AliFRAME(const char *name, const char *title);
  virtual      ~AliFRAME() {}
  virtual void   Init() {}
  virtual Int_t IsVersion() const =0;
 
   ClassDef(AliFRAME,1)  //Class for Space Frame
};

#endif



