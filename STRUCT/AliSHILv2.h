#ifndef ALISHILv2_H
#define ALISHILv2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for Module: SHIL          //
////////////////////////////////////////////////
 
#include "AliSHIL.h"
 
 
class AliSHILv2 : public AliSHIL {
  
public:
  AliSHILv2();
  AliSHILv2(const char *name, const char *title);
  virtual      ~AliSHILv2() {}
  virtual void  CreateGeometry();
  virtual void  Init();
  virtual Int_t IsVersion() const {return 0;}
  virtual void  SetPbCone(Bool_t flag = kTRUE)         {fPbCone        = flag;}
  virtual void  SetWriteGeometry(Bool_t flag = kFALSE)  {fWriteGeometry = flag;}	  
 protected:
  Bool_t fPbCone;           // outer Pb cone option flag 
  Bool_t fWriteGeometry;    // flag to write out the fluka geometry
  ClassDef(AliSHILv2,1)     // Muon Shield Class (Open Geometry)
};

#endif
