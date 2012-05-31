#ifndef ALIBODY_H
#define ALIBODY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: BODY          //
//   This is the envelop for Alice            //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliBODY : public AliModule {
 
public:
   AliBODY();
   AliBODY(const char *name, const char *title);
   virtual      ~AliBODY() {}
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual Int_t IsVersion() const {return 0;}

   ClassDef(AliBODY,1)  //Class manager for the ALICE body
};

#endif
