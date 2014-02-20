#ifndef ALIFITDIGITIZER_H
#define ALIFITDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/******************************************************************
 *    Produde digits from hits
 * Alla.Maevskaya@cern.ch 
 ********************************************************************/

#include <AliDigitizer.h>

#include <AliDigitizationInput.h>
class AliFIT;
class AliFITHits;
class AliFITDigit;

class AliFITDigitizer : public AliDigitizer {
 public:
  
  AliFITDigitizer();
  AliFITDigitizer(AliDigitizationInput * digInput);
  virtual ~AliFITDigitizer();
  virtual Bool_t Init();
  TClonesArray *Hits() const {return fHits;}
  TClonesArray *Digits() const {return fDigits;}
 // Do the main work
  void Digitize(Option_t* /*option=0*/) ;
  void AddDigit(Int_t* digits, Int_t*);
  void AddSDigit(Int_t npmt,
		 Int_t timeCFD, Int_t timeLED, Int_t timeQT0, Int_t timeQT1) ;
 
  void ResetDigits();
  void DigitizeHits();
private:

  AliFIT *fFIT;            //!
  TClonesArray *fHits;   //! List of hits
  TClonesArray *fDigits;   //! digits
  Int_t    fNdigits;                //! Number of digits


  AliFITDigitizer(const AliFITDigitizer&);
  AliFITDigitizer& operator=(const AliFITDigitizer&);

    ClassDef(AliFITDigitizer,1)
};    

#endif

