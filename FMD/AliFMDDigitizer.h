#ifndef ALIFMDDIGITIZER_H
#define ALIFMDDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliDigitizer.h"

class AliRunDigitizer;

class AliFMDDigitizer : public AliDigitizer {
 public:
  
  AliFMDDigitizer();
  AliFMDDigitizer(AliRunDigitizer * manager);
  virtual ~AliFMDDigitizer();
  virtual Bool_t Init();
    
  
  // Do the main work
  void Exec(Option_t* option=0) ;
  Int_t PutNoise(Int_t charge){return (Int_t)(gRandom->Gaus(charge,500));}
  TClonesArray *SDigits() const {return fSDigits;}
 

  void ReadDigit(Int_t a[][50][300], Int_t);
  
  enum {kBgTag = -1};
      
   
 private:
    TClonesArray *fDigits;               // ! array with digits
    TClonesArray *fSDigits      ; // List of summable digits
     
    ClassDef(AliFMDDigitizer,0)
};    
#endif

