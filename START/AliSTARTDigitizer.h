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
class AliSTARTdigit;

class AliSTARTDigitizer : public AliDigitizer {
 public:
  
  AliSTARTDigitizer();
  AliSTARTDigitizer(AliRunDigitizer * manager);
  virtual ~AliSTARTDigitizer();
  virtual Bool_t Init();
  TClonesArray *Hits() const {return fHits;}
  TArrayI *timeTDC() {return ftimeTDC;}
  TArrayI * ADC() {return fADC;} //for slow simulation
  // Do the main work
  void Exec (Option_t* /*option=0*/) ;
  Bool_t RegisterPhotoE(Double_t e);
  enum {kBgTag = -1};

private:

  AliSTART *fSTART;          //!
  TClonesArray *fPhotons   ; //! Number of Cherenkov photons		      
  TClonesArray *fHits      ; //! List of hits
  AliSTARTdigit *fdigits   ; //! digits
  TArrayI *ftimeTDC    ; //! array of TDC signal from right side
  TArrayI *fADC     ;//! array of ADC signal from left sida
  TH1*     fEff;    //! efficiency histogram
  
    ClassDef(AliSTARTDigitizer,1)
};    
#endif

