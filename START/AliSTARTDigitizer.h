#ifndef ALISTARTDIGITIZER_H
#define ALISTARTDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigitizer.h>
#include <AliLoader.h>
#include <AliRunLoader.h>

#include <AliRunDigitizer.h>
class AliSTART;
class AliSTARThit;
class AliSTARThitPhoton;
class AliSTARTdigit;

class AliSTARTDigitizer : public AliDigitizer {
 public:
  
  AliSTARTDigitizer();
  AliSTARTDigitizer(AliRunDigitizer * manager);
  virtual ~AliSTARTDigitizer();
  virtual Bool_t Init();
  TClonesArray *Hits() const {return fHits;}
  TClonesArray *Photons() const {return fPhotons;}
  TArrayI *timeRightTDC() {return ftimeRightTDC;} //for slow simulation
  TArrayI *timeLeftTDC() {return ftimeLeftTDC;}
  TArrayI *RightADC() {return fRightADC;} //for slow simulation
  TArrayI *LeftADC() {return fLeftADC;}
  // Do the main work
  void Exec (Option_t* /*option=0*/) ;
  Bool_t RegisterPhotoE(/*AliSTARThitPhoton *hit*/);//!!!
  Bool_t GetDebug() const {return fManager->GetDebug();}
  enum {kBgTag = -1};

private:

  AliSTART *fSTART;
  TClonesArray *fPhotons   ; //Number of Cherenkov photons		      
  TClonesArray *fHits      ; // List of hits
  AliSTARTdigit *fdigits   ; // digits
  TArrayI *ftimeRightTDC    ; //array of TDC signal from right sida
  TArrayI *ftimeLeftTDC     ; ////array of TDC signal from left side
  TArrayI *fRightADC    ;   //array of ADC signal from right sida 
  TArrayI *fLeftADC     ;//array of ADC signal from left sida
  
    ClassDef(AliSTARTDigitizer,0)
};    
#endif

