#ifndef ALITOFDigitizer_H
#define ALITOFDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//  Task Class for making Digits in TOF
//  Comment:
//
// -- Author: F. Pierella (Bologna University) pierella@bo.infn.it


#include "TTask.h"
#include "TString.h"

class AliTOFDigitizer: public TTask {

public:
  AliTOFDigitizer() ;          // ctor
  AliTOFDigitizer(const char* HeaderFile,const char* digitsTitle = 0) ; 
  virtual ~AliTOFDigitizer() ; // dtor
  virtual void  Exec(Option_t* option); 
  
  Float_t  GetTimeRes() const {return fTimeRes;}
  Float_t  GetChrgRes() const {return fChrgRes;}
  char*    GetDigitsBranch()const{return (char*) fDigitsTitle.Data();}  

  virtual void Print(Option_t* option) const ;

  void     SetTimeRes(Float_t timeRes)  {fTimeRes = timeRes ;}
  void     SetChrgRes(Float_t chrgRes)  {fChrgRes = chrgRes ;}
  void     SetDigitsBranch(const char* title ) ;

  Bool_t   operator == (const AliTOFDigitizer & sd) const ;

private:
  void     Init() ;
  void     PrintDigits(Option_t* option) ;

private:
  Float_t fTimeRes;                // Time Resolution
  Float_t fChrgRes;                // ADC parameter
  Int_t   fNevents ;               // Number of events to digitize
  TString fDigitsTitle ;           // title of Digits branch
  TString fHeadersFile ;           // input file
  Bool_t  fIsInitialized ;         // kTRUE if Digitizer is initialized
  TClonesArray* fDigits ;          // list of Digits
  TClonesArray* fHits ;            // list of Hits


  ClassDef(AliTOFDigitizer,1)  // Task Class for making Digits in TOF

};

#endif // AliTOFDigitizer_H
