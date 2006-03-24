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
  TArrayI *timeCFD() {return ftimeCFD;}
  TArrayI *timeLED() {return ftimeLED;}
  TArrayI * ADC() {return fADC;} 
   TArrayI * ADC0() {return fADC0;} 

  // Do the main work
  void Exec (Option_t* /*option=0*/) ;
  Bool_t RegisterPhotoE(Int_t impt, Double_t energy);
  enum {kBgTag = -1};
 
private:

  AliSTART *fSTART;          //!
  TClonesArray *fHits      ; //! List of hits
  AliSTARTdigit *fdigits   ; //! digits
  TArrayI *ftimeCFD    ; //! array of CFD signal 
  TArrayI *ftimeLED    ; //! array of (LED-GFD) time (amplitude)
  TArrayI *fADC     ;//! array of QTC signals (main amplitude)
  TArrayI *fADC0     ;//! array of QTC signals (main amplitude)
  Int_t fSumMult; // multiplicity
  TObjArray fEffPMT; //pmt registration effeicincy

  
    ClassDef(AliSTARTDigitizer,1)
};    
#endif

