#ifndef ALIDIPO_H
#define ALIDIPO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for Module: DIPO          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliDIPO : public AliModule {
 
public:
  AliDIPO();
  AliDIPO(const char *name, const char *title);
  virtual      ~AliDIPO() {}
  virtual void  Init();
  
  ClassDef(AliDIPO,1)  //Class for the dipole magnet
};

#endif
