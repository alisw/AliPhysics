#ifndef TRDv1_H
#define TRDv1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 1    //
////////////////////////////////////////////////////////
 
#include "AliTRD.h"
             
class AliTRDv1 : public AliTRD {

public:
  AliTRDv1() {}
  AliTRDv1(const char *name, const char *title);
  virtual      ~AliTRDv1() {}
  virtual void  CreateGeometry();
  virtual void  CreateMaterials();
  virtual Int_t IsVersion() const {return 1;}
  virtual void  SetHits(Int_t ihit = 1) { fHitsOn = ihit; };
  virtual void  StepManager();
  virtual void  Init();

protected:
  Int_t        fIdSens;     // Sensitive volume identifier

  Int_t        fIdSpace1;   // Spaceframe volume identifier
  Int_t        fIdSpace2;   // 
  Int_t        fIdSpace3;   // 

  Int_t        fIdChamber1; // Driftchamber volume identifier
  Int_t        fIdChamber2; // 
  Int_t        fIdChamber3; // 

  Int_t        fHitsOn;     // Used to switch hits on
            
  ClassDef(AliTRDv1,1)      // Transition Radiation Detector version 1

};

#endif
