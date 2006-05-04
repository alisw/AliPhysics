#ifndef ALITRDV1_H
#define ALITRDV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////
//  Manager and hits classes for set:TRD version 1    //
////////////////////////////////////////////////////////

// Energy spectrum of the delta-rays 
Double_t Ermilova(Double_t *x, Double_t *par);
Double_t IntSpecGeant(Double_t *x, Double_t *par);
 
#include "AliTRD.h"

class TF1;
class TTree;
class TFile;
class AliTRDsim;

//_____________________________________________________________________________
class AliTRDv1 : public AliTRD {

 public:

  AliTRDv1();
  AliTRDv1(const char *name, const char *title);
  AliTRDv1(const AliTRDv1 &trd);
  virtual ~AliTRDv1();
  AliTRDv1 &operator=(const AliTRDv1 &trd);

  virtual void       Copy(TObject &trd) const;
  virtual void       CreateGeometry();
  virtual void       CreateMaterials();
  virtual void       CreateTRhit(Int_t det);
  virtual Int_t      IsVersion() const          { return 1; };
  virtual void       StepManager();
  virtual void       Init();

          void       StepManagerErmilova();
          void       StepManagerGeant();
          void       StepManagerFixedStep();
          void       SelectStepManager(Int_t t);
          void       SetStepSize(Double_t s)    { fStepSize = s; };

          void       SetSensPlane(Int_t iplane = 0);
          void       SetSensChamber(Int_t ichamber = 0);
          void       SetSensSector(Int_t isector);
          void       SetSensSector(Int_t isector, Int_t nsector);

          void       SetTR(Bool_t kTRUE)        { fTRon = kTRUE; };

          Int_t      GetSensPlane() const       { return fSensPlane;       };
          Int_t      GetSensChamber() const     { return fSensChamber;     };
          Int_t      GetSensSector() const      { return fSensSector;      };
          Int_t      GetSensSectorRange() const { return fSensSectorRange; };

          Bool_t     GetTR() const              { return fTRon;            };
          AliTRDsim *GetTRDsim() const          { return fTR;              };

 protected:

  void        *StepManagerEntity();

  Int_t        fSensSelect;             // Switch to select only parts of the detector
  Int_t        fSensPlane;              // Sensitive detector plane
  Int_t        fSensChamber;            // Sensitive detector chamber
  Int_t        fSensSector;             // Sensitive detector sector 
  Int_t        fSensSectorRange;        // Sensitive detector range

  Bool_t       fTRon;                   // Switch for TR simulation
  AliTRDsim   *fTR;                     // TR simulator

  Int_t        fTypeOfStepManager;      // Type of Step Manager.
  Double_t     fStepSize;               // Used for the fixed step size

 private:

  Double_t     BetheBloch(Double_t bg);
  Double_t     BetheBlochGeant(Double_t bg);
  void 	       Stepping();
  
  TF1         *fDeltaE;                 // Energy distribution of the delta-electrons (Ermilova)
  TF1         *fDeltaG;                 // Energy distribution of the
	                                // Delta-electrons (GEANT) for StepManagerGeant
  Float_t      fTrackLength0;         	// Save the track length at chamber entrance  
  Int_t	       fPrimaryTrackPid;        // Save the id of the primary track  

  ClassDef(AliTRDv1,4)                  // Transition Radiation Detector version 1 (slow simulator)

};

#endif
