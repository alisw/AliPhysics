#ifndef MUONV0_H
#define MUONV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include "AliMUONClusterFinder.h"
class AliMUONv0 : public AliMUON {
public:
   AliMUONv0();
   AliMUONv0(const char *name, const char *title);
   virtual  ~AliMUONv0() {}
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   Init();
   virtual Int_t  IsVersion() const {return 0;}
   virtual void   StepManager();
   virtual void   Trigger(Float_t (*)[4], Float_t (*)[4], Int_t& iflag);
private:
   ClassDef(AliMUONv0,1)  //Hits manager for set:MUON version 0
};
#endif







