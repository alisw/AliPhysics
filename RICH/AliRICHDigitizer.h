#ifndef AliRICHDigitizer_h
#define AliRICHDigitizer_h
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <AliDigitizer.h>
#include <AliRunDigitizer.h>
#include <AliRun.h>

class AliRICH;

class AliRICHDigitizer : public AliDigitizer 
{
public:
           AliRICHDigitizer();
           AliRICHDigitizer(AliRunDigitizer * manager);
  virtual ~AliRICHDigitizer();
        
  void     Exec(Option_t* option=0);                //virtual
  Bool_t   GetDebug() const {return gAlice->GetDebug();}
  AliRICH* Rich()     const {return fRich;}
protected:
  AliRICH* fRich; //pointer to main RICH object
  ClassDef(AliRICHDigitizer,0)
};    
#endif
