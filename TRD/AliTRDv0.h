#ifndef ALITRDV0_H
#define ALITRDV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 0    //
////////////////////////////////////////////////////////
 
#include "AliTRD.h"

class AliTRDsim;

//_____________________________________________________________________________ 
class AliTRDv0 : public AliTRD {

 public:

  AliTRDv0();
  AliTRDv0(const char *name, const char *title);
  virtual ~AliTRDv0();

  virtual void       CreateGeometry();
  virtual void       CreateMaterials();
  virtual Int_t      IsVersion() const           { return 0; };
  virtual void       StepManager();
  virtual void       Init();

  virtual void       SetHits()                   { fHitsOn = 1; };

          void       SetSensChamber(Int_t )      { };
          void       SetSensPlane(Int_t )        { };
          void       SetSensSector(Int_t )       { };
          void       SetSensSector(Int_t ,Int_t) { };

          Int_t      GetSensChamber() const      { return 0; };
          Int_t      GetSensPlane() const        { return 0; };
          Int_t      GetSensSector() const       { return 0; };
          Int_t      GetSensSectorRange() const  { return 0; };

          AliTRDsim *CreateTR()                  { return 0; };
          AliTRDsim *GetTR() const               { return 0; };

 protected:

  Int_t        fIdSens;     // Sensitive volume identifier

  Int_t        fIdChamber1; // Driftchamber volume identifier
  Int_t        fIdChamber2; // Driftchamber volume identifier
  Int_t        fIdChamber3; // Driftchamber volume identifier

  Int_t        fHitsOn;     // Used to switch hits on

  ClassDef(AliTRDv0,1)      // Transition Radiation Detector version 0 (fast simulator)

};

#endif
