#ifndef ALISTARTDIGITIZER_H
#define ALISTARTDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliDigitizer.h"
#include "AliLoader.h"
#include "AliRunLoader.h"

class AliRunDigitizer;
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
  //  TArrayI *timeRightADC() {return ftimeRightADC;}
  // TArrayI *timeLeftADC() {return ftimeLeftADC;}
  // Do the main work
  void Exec(Option_t* option=0) ;
  Bool_t RegisterPhotoE(AliSTARThitPhoton *hit);			//!!!

  enum {kBgTag = -1};

private:

  AliSTART *START;
  TClonesArray *fPhotons   ; 						//!!! 
  TClonesArray *fHits      ; // List of summable digits
  AliSTARTdigit *fdigits   ; // digits
  TArrayI *ftimeRightTDC    ;
  TArrayI *ftimeLeftTDC     ;
  TArrayI *fRightADC    ;
  TArrayI *fLeftADC     ;
     
    ClassDef(AliSTARTDigitizer,0)
};    
#endif

