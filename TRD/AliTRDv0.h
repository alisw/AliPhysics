#ifndef TRDv0_H
#define TRDv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 0    //
////////////////////////////////////////////////////////
 
#include "AliTRD.h"

//_____________________________________________________________________________ 
class AliTRDv0 : public AliTRD {

 public:

  AliTRDv0() {};
  AliTRDv0(const char *name, const char *title);
  ~AliTRDv0() {};
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual Int_t   IsVersion() const { return 0; };
  virtual void    StepManager();
  virtual void    Init();

  virtual void    SetHits(Int_t ihit = 1) { fHitsOn = ihit; };

          void    SetSensChamber(Int_t ichamber) { };
          void    SetSensPlane(Int_t iplane)     { };
          void    SetSensSector(Int_t isector)   { };

          Int_t   GetSensChamber() { return 0; };
          Int_t   GetSensPlane()   { return 0; };
          Int_t   GetSensSector()  { return 0; };

 protected:

  Int_t        fIdSens;     // Sensitive volume identifier

  Int_t        fIdChamber1; // Driftchamber volume identifier
  Int_t        fIdChamber2; // 
  Int_t        fIdChamber3; // 

  Int_t        fHitsOn;     // Used to switch hits on

  ClassDef(AliTRDv0,1)      // Transition Radiation Detector version 0 (fast simulator)

};

#endif
