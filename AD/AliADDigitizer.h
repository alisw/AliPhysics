#ifndef ALIADDigitizer_H
#define ALIADDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
  
///_________________________________________________________________________
///
///  Class for making Digits in AD 
///_________________________________________________________________________   


// --- Standard library ---

// --- AliRoot header files ---

#include "AliDigitizer.h"

// #include "AliADConst.h"

class TClonesArray;
class AliDigitizationInput;
class AliCDBManager;
class AliCDBStorage;

class AliADDigitizer: public AliDigitizer {

public:
                    AliADDigitizer() ;                              // default constructor
                    AliADDigitizer(AliDigitizationInput* digInput); // constructor
  virtual          ~AliADDigitizer() ;              // destructor

  virtual Bool_t    Init();
  virtual void      Digitize(Option_t* option=0);
  void AddDigit(Int_t* track, Int_t module, Float_t cell);
  void AddDigit(Int_t* modul, Float_t cell);	
 

            void    ResetDigit();
            
 
private:
 
                    AliADDigitizer(const AliADDigitizer& /*digitizer*/); 
                    AliADDigitizer& operator = (const AliADDigitizer& /*digitizer*/); 
  

           Int_t    fNdigits;           //! Number of digits
    TClonesArray*   fDigits;            //! List of digits

   ClassDef(AliADDigitizer,1)     // digitizer for AD

};

#endif // AliADDigitizer_H
