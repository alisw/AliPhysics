#ifndef ALISTARTDIGITIZER_H
#define ALISTARTDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliDigitizer.h"

class AliRunDigitizer;
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
  
  // Do the main work
  void Exec(Option_t* option=0) ;
  
  enum {kBgTag = -1};

private:

  AliSTART *START;
  TClonesArray *fHits      ; // List of summable digits
  AliSTARTdigit *fdigits   ; // digits
     
    ClassDef(AliSTARTDigitizer,0)
};    
#endif

