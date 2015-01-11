#ifndef ALIFIT_H
#define ALIFIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager and hits classes for set: FIT 
// Alla.Maevskaya@cern.ch 
////////////////////////////////////////////////
 
#include <AliDetector.h>
#include <TTree.h>
#include "AliFITDigit.h"
#include "AliFITHits.h"
class TDirectory;
class TFile;
class AliFIT: public AliDetector {
  
 public:
  AliFIT();
  AliFIT(const char *name, const char *title);
  virtual       ~AliFIT();
  virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits);
  //  virtual void   CreateGeometry() {};
  // virtual void   CreateMaterials() {};
   void AddDigit(Int_t npmt, 
  		Int_t timeCFD, Int_t timeLED, Int_t timeQT0, Int_t timeQT1, Int_t *labels) ;
   virtual void AddDigit(Int_t*, Int_t*) {};
 virtual Int_t  IsVersion()const {return 0;}
  virtual void   Init();
  //  virtual void   DefineOpticalProperties() {};
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void SetHitsAddressBranch(TBranch *b1)
  {b1->SetAddress(&fHits);}
  virtual void   StepManager() {};
  virtual void   ResetHits(); 
  virtual void   ResetDigits();
  virtual void   SetTreeAddress();
  virtual AliDigitizer* CreateDigitizer(AliDigitizationInput* digInput) const;
  virtual  void  Digits2Raw ();
  virtual void  Raw2Digits (AliRawReader *reader,TTree* digitsTree);
  
  virtual void  Raw2Digits() {}
 protected:
  Int_t           fIdSens;    // Sensetive Cherenkov photocathode
  //  AliFITDigit     *fDigits;    // pointer to T0digits
 
 private:
  AliFIT(const AliFIT&);
  AliFIT& operator=(const AliFIT&);
  
  ClassDef(AliFIT,1)  //Base class for the FIT detector
    };


//_____________________________________________________________________________
 
#endif































