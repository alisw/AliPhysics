#ifndef ALIACORDEDigitizer_H
#define ALIACORDEDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  
///_________________________________________________________________________
///
///  Class for making Digits in ACORDE 
///_________________________________________________________________________   


// --- Standard library ---

// --- AliRoot header files ---

#include "AliDigitizer.h"

class TClonesArray;
class AliDigitizationInput;
class AliCDBManager;
class AliCDBStorage;
class AliACORDECalibData;

class AliACORDEDigitizer: public AliDigitizer {

 public:

   AliACORDEDigitizer() ;                       // constructor
   AliACORDEDigitizer(AliDigitizationInput* digInput);// constructor
   virtual ~AliACORDEDigitizer() ;              // destructor
  
   virtual Bool_t Init();
   virtual void   Digitize(Option_t* option=0);

   void AddDigit(Int_t* track, Int_t module, Float_t time);
   void AddDigit(Int_t* modul, Float_t time);	
   void ResetDigit();
  
   AliACORDECalibData *GetCalibData() const;
  
 protected:
 
   AliACORDECalibData *fCalibData;  //! calibration data
 
 private:
 
   AliACORDEDigitizer(const AliACORDEDigitizer& /*digitizer*/); 
      
   AliACORDEDigitizer& operator = (const AliACORDEDigitizer& /*digitizer*/); 
  

   Int_t   fNdigits;                //! Number of digits
   TClonesArray *fDigits;           //! List of digits

   ClassDef(AliACORDEDigitizer,1)    // digitizer for ACORDE

};

#endif // AliACORDEDigitizer_H
