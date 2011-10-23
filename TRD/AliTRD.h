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


#include "AliDetector.h"
#include "AliTRDTrigger.h"

class AliRawReader;

class AliTRDgeometry;
class AliTriggerDetector;

class AliTRD : public AliDetector {

 public:

  AliTRD();
  AliTRD(const char *name, const char *title);
  virtual         ~AliTRD();

  virtual void     Init();
  virtual Int_t    IsVersion() const          = 0;
  virtual void     ResetDigits();     

  virtual void     CreateGeometry();
  virtual void     CreateMaterials();

  virtual void     Hits2Digits();
  virtual void     Hits2SDigits();
  virtual void     SDigits2Digits();
  virtual void     Digits2Raw();
  virtual Bool_t   Raw2SDigits(AliRawReader* rawReader);

  virtual void     AddHit(Int_t, Int_t*, Float_t*)       { }; 
  virtual void     AddHit(Int_t track, Int_t det, Float_t *hits
                        , Int_t q, Float_t time, Bool_t inDrift); 

  virtual void     SetTreeAddress();

  virtual void     StepManager()              = 0; 

  virtual void     SetStepSize(Double_t s)    = 0;
  virtual void     SetHits()                             { };
  virtual void     SetTR(Bool_t )             = 0;

  virtual Bool_t   GetTR() const              = 0;

          AliTRDgeometry     *GetGeometry() const           { return fGeometry; };
  virtual AliDigitizer       *CreateDigitizer(AliDigitizationInput* digInput) const; 
  virtual AliLoader          *MakeLoader(const char* topfoldername);
  virtual AliTriggerDetector *CreateTriggerDetector() const { return new AliTRDTrigger(); }
  void    SetPrimaryIonisation(Bool_t flag = kTRUE) {fPrimaryIonisation = flag;}  
 protected:

  AliTRDgeometry       *fGeometry;             //  The TRD geometry

  Float_t               fGasDensity;           //  The density of the drift gas
  Float_t               fFoilDensity;          //  The density of the entrance window foil
  Float_t               fGasNobleFraction;     //  The fraction of noble gas in the mixture
  Bool_t                fPrimaryIonisation;    //  switch between Fluka(true) and geant3(false)
 private:

  AliTRD(const AliTRD &trd);
  AliTRD  &operator=(const AliTRD &trd);

  ClassDef(AliTRD,12)                          //  Transition Radiation Detector base class

};

#endif
