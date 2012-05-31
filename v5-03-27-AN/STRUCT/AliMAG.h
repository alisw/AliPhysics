#ifndef ALIMAG_H
#define ALIMAG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: MAG           //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliMAG : public AliModule {
 
public:
   AliMAG();
   AliMAG(const char *name, const char *title);
   virtual      ~AliMAG() {}
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual void  Init();
   virtual Int_t IsVersion() const {return 0;}
 
   ClassDef(AliMAG,1)  //Class manager for detector:MAG
};

#endif
