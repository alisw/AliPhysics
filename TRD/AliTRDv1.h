#ifndef TRDv1_H
#define TRDv1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 1    //
////////////////////////////////////////////////////////

// Energy spectrum of the delta-rays 
Double_t Ermilova(Double_t *x, Double_t *par);
 
#include <TF1.h> 

#include "AliTRD.h"

//_____________________________________________________________________________
class AliTRDv1 : public AliTRD {

 public:

  AliTRDv1() {};
  AliTRDv1(const char *name, const char *title);
  ~AliTRDv1();
  virtual void    CreateGeometry();
  virtual void    CreateMaterials();
  virtual Int_t   IsVersion() const { return 1; };
  virtual void    StepManager();
  virtual void    Init();

          void    SetSensPlane(Int_t iplane = 0);
          void    SetSensChamber(Int_t ichamber = 0);
          void    SetSensSector(Int_t isector = 0);

          Int_t   GetSensPlane()   { return fSensPlane;   };
          Int_t   GetSensChamber() { return fSensChamber; };
          Int_t   GetSensSector()  { return fSensSector;  };

 protected:

  Int_t        fIdSens;                 // Sensitive volume identifier

  Int_t        fIdChamber1;             // Driftchamber volume identifier
  Int_t        fIdChamber2;             // 
  Int_t        fIdChamber3;             // 

  Int_t        fSensSelect;             // Switch to select only parts of the detector
  Int_t        fSensPlane;              // Sensitive detector plane
  Int_t        fSensChamber;            // Sensitive detector chamber
  Int_t        fSensSector;             // Sensitive detector sector

 private:

  virtual Double_t BetheBloch(Double_t bg);

  TF1         *fDeltaE;                 // Energy distribution of the delta-electrons
   
  ClassDef(AliTRDv1,1)                  // Transition Radiation Detector version 1 (slow simulator)

};

#endif
