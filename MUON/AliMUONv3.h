#ifndef ALIMUONV3_H
#define ALIMUONV3_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 3    //
/////////////////////////////////////////////////////////
 
#include "AliMUONv1.h"

class AliMUONv3 : public AliMUONv1 {
public:
   AliMUONv3();
   AliMUONv3(const char *name, const char *title);
   virtual  ~AliMUONv3() {}
   virtual void   StepManager();
protected:
   
private:

   ClassDef(AliMUONv3,1)  // MUON Detector class Version 3
};
#endif







