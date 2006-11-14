#ifndef ALITRD_H
#define ALITRD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Manager and hits classes for set: TRD                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TLorentzVector.h>
#include <TVirtualMC.h>

#include "AliDetector.h"

#include "AliTRDTriggerL1.h"

class TFile;
class TLorentzVector;

class AliRun;
class AliDigit;

class AliTRDhit;
class AliTRDsim;
class AliTRDgeometry;

class AliTRD : public AliDetector {

 public:

  AliTRD();
  AliTRD(const char *name, const char *title);
  virtual         ~AliTRD();

  virtual void     Init();
  virtual Int_t    IsVersion() const          = 0;
  virtual void     ResetDigits();     

  virtual void     BuildGeometry();
  virtual void     CreateGeometry();
  virtual void     CreateMaterials();
  virtual void     DrawModule() const;
  virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
  virtual void     LoadPoints(Int_t track);    

  virtual void     Hits2Digits();
  virtual void     Hits2SDigits();
  virtual void     SDigits2Digits();
  virtual void     Digits2Raw();
  virtual void     MakeBranch(Option_t *option);

  virtual void     AddHit(Int_t, Int_t*, Float_t*)       { }; 
  virtual void     AddHit(Int_t track, Int_t det, Float_t *hits, Int_t q, Bool_t inDrift); 

  virtual void     SetTreeAddress();

  virtual void     StepManager()              = 0; 
  virtual void     StepManagerErmilova()      = 0;
  virtual void     StepManagerGeant()         = 0;
  virtual void     StepManagerFixedStep()     = 0;
  virtual void     SelectStepManager(Int_t t) = 0;

  virtual void     SetStepSize(Double_t s)    = 0;
  virtual void     SetHits()                             { };
  virtual void     SetDrawTR(Int_t idraw = 1)            { fDrawTR      = idraw;         };
  virtual void     SetDisplayType(Int_t type = 0)        { fDisplayType = type;          };
  virtual void     SetTR(Bool_t )             = 0;

  virtual Bool_t   GetTR() const              = 0;
  AliTRDgeometry  *GetGeometry() const                   { return fGeometry;             };

  virtual AliDigitizer    *CreateDigitizer(AliRunDigitizer *manager) const; 
  virtual AliTRDTriggerL1 *CreateTriggerDetector() const { return new AliTRDTriggerL1(); };

 protected:

  AliTRDgeometry  *fGeometry;           //  The TRD geometry

          Float_t  fGasDensity;         //  The density of the drift gas
          Float_t  fFoilDensity;        //  The density of the entrance window foil

          Int_t    fDrawTR;             //  Switches marking the TR photons in the display
          Int_t    fDisplayType;        //  Display type (0: normal, 1: detailed) 

 private:
  AliTRD(const AliTRD &trd);
  AliTRD  &operator=(const AliTRD &trd);

  ClassDef(AliTRD,9)                    //  Transition Radiation Detector base class

};

#endif
