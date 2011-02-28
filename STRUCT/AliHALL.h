#ifndef ALIHALL_H
#define ALIHALL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for detector: HALL          //
////////////////////////////////////////////////
 
#include "AliModule.h"
 
 
class AliHALL : public AliModule {
 
public:
   AliHALL();
   AliHALL(const char *name, const char *title);
   virtual      ~AliHALL() {}
   virtual void  CreateGeometry();
   virtual void  CreateMaterials();
   virtual void  Init();
   virtual Int_t IsVersion() const {return 0;}
 
   ClassDef(AliHALL,1)  //Class for ALICE experimental hall
};

#endif
