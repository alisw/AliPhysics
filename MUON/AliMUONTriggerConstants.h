#ifndef ALIMUONTRIGGERCONSTANTS_H
#define ALIMUONTRIGGERCONSTANTS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>

class AliMUONTriggerConstants :
public TNamed {    
 public:
  AliMUONTriggerConstants();
  ~AliMUONTriggerConstants();
  
public:
  static Int_t Nmodule();              
  static Int_t ModuleId(Int_t imodule);
  static Int_t NstripX(Int_t imodule);
  static Int_t NstripY(Int_t imodule);
  static Float_t XcMin(Int_t imodule); 
  static Float_t XcMax(Int_t imodule);
  static Int_t CircuitId(Int_t icircuit);
  
private:
  static const Int_t fgNmodule;        // total number of module
  static const Int_t fgModuleId[126];  // module Id. number
  static const Int_t fgNstripX[126];   // number of X strips 
  static const Int_t fgNstripY[126];   // number of Y strips 
  static const Float_t fgXcMin[126];   // min X pos of module
  static const Float_t fgXcMax[126];   // max X poa of module 
  static const Int_t fgCircuitId[234]; // circuit Id. number   
  
  ClassDef(AliMUONTriggerConstants,1) // Trigger Constants class

};
#endif













