#ifndef ALIZDC_H
#define ALIZDC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and classes for set ZDC           //
////////////////////////////////////////////////
 
#include "AliDetector.h"

 
class AliZDC : public AliDetector {

public:
  AliZDC();
  AliZDC(const char *name, const char *title);
  virtual       ~AliZDC();
  virtual void  AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual void  BuildGeometry();
  virtual void  CreateGeometry() {}
  virtual void  CreateMaterials() {}
  Int_t         DistancetoPrimitive(Int_t px, Int_t py);
  virtual Int_t IsVersion() const =0;
  virtual void  ResetDigits(); 
  virtual void  StepManager();
  
  // Switching off the shower development in ZDCs
  void  NoShower(){fNoShower=1;}
  void  Shower()  {fNoShower=0;}
  
protected:
  // TClonesArray of stored hits -> not reset et finish event
  // 	 (for digitization at the end of the event)
  TClonesArray *fStHits;
  Int_t fNStHits;
  
  Int_t   fNPrimaryHits;	// Number of primary particles

  Int_t   fNoShower;		// Flag to switch off the shower	

  ClassDef(AliZDC,1)  	// Zero Degree Calorimeter base class
};
 
#endif
