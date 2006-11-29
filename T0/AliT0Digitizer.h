#ifndef ALIT0DIGITIZER_H
#define ALIT0DIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigitizer.h>
#include <AliLoader.h>
#include <AliRunLoader.h>

#include <AliRunDigitizer.h>
class AliT0;
class AliT0hit;
class AliT0digit;

class AliT0Digitizer : public AliDigitizer {
 public:
  
  AliT0Digitizer();
  AliT0Digitizer(AliRunDigitizer * manager);
  virtual ~AliT0Digitizer();
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

  AliT0 *fT0;          //!
  TClonesArray *fHits      ; //! List of hits
  AliT0digit *fdigits   ; //! digits
  TArrayI *ftimeCFD    ; //! array of CFD signal 
  TArrayI *ftimeLED    ; //! array of (LED-GFD) time (amplitude)
  TArrayI *fADC     ;//! array of QTC signals (main amplitude)
  TArrayI *fADC0     ;//! array of QTC signals (main amplitude)
  Int_t fSumMult; // multiplicity
  TObjArray fEffPMT; //pmt registration effeicincy

  AliT0Digitizer(const AliT0Digitizer&);
  AliT0Digitizer& operator=(const AliT0Digitizer);

  
    ClassDef(AliT0Digitizer,1)
};    

typedef AliT0Digitizer AliSTARTDigitizer; // for backward compatibility

#endif

