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
           AliRICHDigitizer()                                                {;}
           AliRICHDigitizer(AliRunDigitizer * manager):AliDigitizer(manager) {if(GetDebug())Info("main ctor","Start.");}
  virtual ~AliRICHDigitizer()                                                {if(GetDebug())Info("dtor","Start.");}

        
  void     Exec(Option_t* option=0);                //virtual
  Bool_t   Init();                                  //virtual
  Bool_t   GetDebug() const {return ((gAlice) ? gAlice->GetDebug() : kFALSE);}
  AliRICH* Rich()     const {return fRich;}
protected:
  AliRICH* fRich; //pointer to main RICH object
  ClassDef(AliRICHDigitizer,0)
};    
#endif
