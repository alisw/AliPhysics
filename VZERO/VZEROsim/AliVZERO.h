#ifndef ALIVZERO_H
#define ALIVZERO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//////////////////////////////////////////////////
//  Manager and hits classes for set : VZERO    //
//////////////////////////////////////////////////

/*
#include "AliRun.h"
#include "AliMC.h"
#include "AliDetector.h"
#include "AliVZEROLoader.h"
 
#include <TNamed.h>
#include <TTree.h>
*/
#include "AliDetector.h"
#include "AliVZEROTrigger.h"

class TNamed;
class TTree;
class TF1;

class AliVZEROLoader;
class AliVZEROhit; 
class AliVZEROdigit;
class AliVZEROCalibData;
class AliVZERORecoParam;
  
class AliVZERO : public AliDetector {
 
public:

  AliVZERO();
  AliVZERO(const char *name, const char *title);
  virtual       ~AliVZERO();

  virtual void   CreateGeometry();
  virtual void   CreateMaterials();
  virtual void   AddAlignableVolumes() const {}
  virtual Int_t  IsVersion() const = 0;
  virtual void   Init();
  virtual AliLoader* MakeLoader(const char* topfoldername);
  virtual void   Hits2Digits();
  virtual void   Hits2SDigits();
  virtual void   Digits2Raw();
  virtual Bool_t Raw2SDigits(AliRawReader*);
  virtual void   SetTreeAddress();  
  virtual void   MakeBranch(Option_t *option) =0;
  virtual void   StepManager() {};
// Trigger 
  virtual AliTriggerDetector* CreateTriggerDetector() const 
  { return new AliVZEROTrigger(); }
  
  virtual void   SetThickness(Float_t thick)  {fThickness = thick;};
  virtual void   SetThickness1(Float_t thick) {fThickness1 = thick;};
// Set Stepping Parameters
  virtual void   SetMaxStepQua(Float_t p1);
  virtual void   SetMaxStepAlu(Float_t p1);
  virtual void   SetMaxDestepQua(Float_t p1);
  virtual void   SetMaxDestepAlu(Float_t p1);

  AliDigitizer*  CreateDigitizer(AliDigitizationInput* digInput) const;

  void           GetCalibData();
  Float_t        CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const;
  double         SignalShape(double *x, double *par);

protected:

   Int_t   fIdSens1;      // Sensitive volume  in VZERO
   Float_t fThickness;    // Total thickness of box holding Right detector V0R i.e. 4.1 cm
   Float_t fThickness1;   // Thickness of elementary cells i.e. 0.7 cm
  
// Stepping Parameters
   Float_t fMaxStepQua;   // Maximum step size inside the quartz volumes
   Float_t fMaxStepAlu;   // Maximum step size inside the  aluminum volumes
   Float_t fMaxDestepQua; // Maximum relative energy loss in quartz
   Float_t fMaxDestepAlu; // Maximum relative energy loss in aluminum

private:
   AliVZERO(const AliVZERO& /*vzero*/); 
   AliVZERO& operator = (const AliVZERO& /*vzero*/); 

   AliVZEROCalibData *fCalibData;      //! Pointer to the calibration object
   Int_t              fNBins[64];      //! Number of bins in each SDigit
   Float_t            fBinSize[64];    //! Bin size in each SDigit
   TF1*               fTimeSlewing;    //! Function for time slewing correction
   TF1*               fSignalShape;    //! Function for signal shape used in Raw->SDigits
   AliVZERORecoParam *fRecoParam;      //! Reco params used in Raw->SDigits

  ClassDef(AliVZERO,2)  //Class for the VZERO detector
};

//____________________________________________________________

#endif
