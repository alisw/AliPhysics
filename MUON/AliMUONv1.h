#ifndef ALIMUONV1_H
#define ALIMUONV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
 
#include "AliMUON.h"

class AliMUONv1 : public AliMUON {
public:
   AliMUONv1();
   AliMUONv1(const char *name, const char *title);
   virtual  ~AliMUONv1() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();
   virtual Int_t  IsVersion() const {return 1;}
   virtual void   StepManager();
private:
   ClassDef(AliMUONv1,1)  // MUON Detector class Version 1
};
#endif







