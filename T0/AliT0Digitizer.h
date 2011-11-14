#ifndef ALIT0DIGITIZER_H
#define ALIT0DIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/******************************************************************
 *    Produde digits from hits
 * Alla.Maevskaya@cern.ch 
 ********************************************************************/

#include <AliDigitizer.h>
#include "AliT0Parameters.h"

#include <AliDigitizationInput.h>
class AliT0;
class AliT0hit;
class AliT0digit;

class AliT0Digitizer : public AliDigitizer {
 public:
  
  AliT0Digitizer();
  AliT0Digitizer(AliDigitizationInput * digInput);
  virtual ~AliT0Digitizer();
  virtual Bool_t Init();
  TClonesArray *Hits() const {return fHits;}
  TArrayI *TimeCFD() {return ftimeCFD;}
  TArrayI *TimeLED() {return ftimeLED;}
  TArrayI * ADC() {return fADC;} 
  TArrayI * ADC0() {return fADC0;} 

  // Do the main work
  void Digitize(Option_t* /*option=0*/) ;
  //  Bool_t RegisterPhotoE(Int_t impt, Double_t energy);
  enum {kBgTag = -1};
 
private:

  AliT0 *fT0;            //!
  TClonesArray *fHits;   //! List of hits
  AliT0digit *fdigits;   //! digits
  TArrayI *ftimeCFD;     //! array of CFD signal 
  TArrayI *ftimeLED;     //! array of (LED-GFD) time (amplitude)
  TArrayI *fADC;         //! array of QTC signals (main amplitude)
  TArrayI *fADC0;        //! array of QTC signals (main amplitude)
  Int_t fSumMult;        // multiplicity
  TObjArray fAmpLED;     // amplitude  CFD-LED dependence #channel -> #MIPs
  TObjArray fAmpQTC;     // amplitude  QTC dependence #channel -> #MIPs
 
  AliT0Parameters  *fParam;           //pointer to T0 parameters class     


  AliT0Digitizer(const AliT0Digitizer&);
  AliT0Digitizer& operator=(const AliT0Digitizer&);

    ClassDef(AliT0Digitizer,4)
};    

typedef AliT0Digitizer AliSTARTDigitizer; // for backward compatibility

#endif

