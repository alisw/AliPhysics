#ifndef ALIPHOSSDigitizer_H
#define ALIPHOSSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Task Class for making SDigits in PHOS      
// A Summable Digits is the sum of all hits originating 
// from one primary in one active cell
//*--
//*-- Author: Dmitri Peressounko(SUBATECH & KI)


// --- ROOT system ---
#include "TTask.h"
#include "TString.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSSDigitizer: public TTask {

public:
  AliPHOSSDigitizer() ;          // ctor
  AliPHOSSDigitizer(const char* headerFile, const char *sdigitsTitle = "Default", const Bool_t toSplit = kFALSE) ; 
  virtual ~AliPHOSSDigitizer() ; // dtor

  Float_t  Calibrate(Int_t amp)const {return (amp - fA)/fB ; }
  Int_t    Digitize(Float_t Energy)const { return (Int_t ) ( fA + Energy*fB); }
  virtual void   Exec(Option_t *option); 
  const char *   GetSDigitsBranch()const{return GetName();}  
  const Int_t    GetSDigitsInRun() const {return fSDigitsInRun ;}  
  virtual void Print(Option_t* option) const ;
  void SetSDigitsBranch(const char * title ) ;
  void UseHitsFrom(const char * filename) ;      
  Bool_t operator == (const AliPHOSSDigitizer & sd) const ;

private:
  void     Init() ;
  void     InitParameters() ;
  void     PrintSDigits(Option_t * option) ;

private:
  Float_t fA ;              // Pedestal parameter
  Float_t fB ;              // Slope Digitizition parameters
  Float_t fPrimThreshold ;  // To store primari if Elos > threshold
  Bool_t  fDefaultInit;     //! Says if the task was created by defaut ctor (only parameters are initialized)
  Int_t   fSDigitsInRun ;   //! Total number of sdigits in one run
  TFile * fSplitFile ;      //! file in which SDigits will eventually be stored
  Bool_t  fToSplit ;        //! Says that sigits should be written into splip file

  ClassDef(AliPHOSSDigitizer,1)  // description 

};

#endif // AliPHOSSDigitizer_H
