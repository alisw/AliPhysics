#ifndef AliRICHDigitizer_h
#define AliRICHDigitizer_h
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliDigitizer.h"

class AliRICHDigitizer : public AliDigitizer 
{
public:
           AliRICHDigitizer();
           AliRICHDigitizer(AliRunDigitizer * manager);
  virtual ~AliRICHDigitizer();
        
  void   Exec(Option_t* option=0);                //virtual
protected:
  ClassDef(AliRICHDigitizer,0)
};    
#endif
