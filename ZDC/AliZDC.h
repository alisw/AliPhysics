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
  virtual void  StepManager();
  virtual void  ResetDigits(); 

protected:

  // Parameters for hadronic calorimeters geometry
  Float_t fDimZN[3];  // Dimensions of neutron detector
  Float_t fDimZP[3];  // Dimensions of proton detector
  Float_t fPosZN[3];  // Position of neutron detector
  Float_t fPosZP[3];  // Position of proton detector
  Float_t fFibZN[3];  // Fibers for neutron detector
  Float_t fFibZP[3];  // Fibers for proton detector
  Float_t fGrvZN[3];  // Grooves for neutron detector
  Float_t fGrvZP[3];  // Grooves for proton detector
  Int_t   fDivZN[3];  // Division for neutron detector
  Int_t   fDivZP[3];  // Division for proton detector
  Int_t   fTowZN[2];  // Tower for neutron detector
  Int_t   fTowZP[2];  // Tower for proton detector

  // Parameters for EM calorimeter geometry
  Float_t fDimZEMPb;  // z-dimension of the Pb slice
  Float_t fDimZEMAir; // scotch
  Float_t fFibRadZEM; // External fiber radius (including cladding)
  Float_t fFibZEM[3]; // Fibers for EM calorimeter
  Float_t fDimZEM[6]; // Dimensions of EM detector
  Float_t fPosZEM[3]; // Position of EM detector
  Int_t   fDivZEM[3]; // Divisions for EM detector
  
  // TClonesArray of stored hits -> not reset et finish event
  // 	(for digitization at the end of the event)
  TClonesArray *fStHits;
  Int_t fNStHits;
  
  Int_t   fNPrimaryHits;

   ClassDef(AliZDC,1)  // Zero Degree Calorimeter base class
};
 
#endif
