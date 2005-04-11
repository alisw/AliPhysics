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
  TArrayI * ADC() {return fADC;} 
  TArrayI *timeTDCAmp() {return ftimeTDCAmp;}
  TArrayI * ADCAmp() {return fADCAmp;} 
  TArrayI *SumMult() {return fSumMult;}
  // Do the main work
  void Exec (Option_t* /*option=0*/) ;
  Bool_t RegisterPhotoE(Double_t energy);
  enum {kBgTag = -1};

private:

  AliSTART *fSTART;          //!
  TClonesArray *fHits      ; //! List of hits
  AliSTARTdigit *fdigits   ; //! digits
  TArrayI *ftimeTDC    ; //! array of TDC signal from right side
  TArrayI *fADC     ;//! array of ADC signal from left sida
  TArrayI *ftimeTDCAmp    ; //! array of TDC amplified signal from right side
  TArrayI *fADCAmp     ;//! array of ADC amplified signal from left sida
  TArrayI *fSumMult; // multiplicity
  TH1*     fEff;    //! efficiency histogram
  
    ClassDef(AliSTARTDigitizer,1)
};    
#endif

