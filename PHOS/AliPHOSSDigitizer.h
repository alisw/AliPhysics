#ifndef ALIPHOSSDigitizer_H
#define ALIPHOSSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Task Class for making SDigits in PHOS      
//                  
//*-- Author: Dmitri Peressounko(SUBATECH & KI)


// --- ROOT system ---
#include "TTask.h"
#include "TString.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSSDigitizer: public TTask {

public:
  AliPHOSSDigitizer() ;          // ctor
  AliPHOSSDigitizer(const char* HeaderFile,const char *SdigitsTitle = 0) ; 
  virtual ~AliPHOSSDigitizer() ; // dtor

  Float_t  Calibrate(Int_t amp)const {return (amp - fA)/fB ; }
  Int_t    Digitize(Float_t Energy)const { return (Int_t ) ( fA + Energy*fB); }

  virtual void  Exec(Option_t *option); 
  
  Float_t  GetPedestalParameter()const {return fA;}
  Float_t  GetCalibrationParameter()const{return fB;}
  char *   GetSDigitsBranch()const{return (char*) fSDigitsTitle.Data();}  

  virtual void Print(Option_t* option) const ;

  void     SetPedestalParameter(Float_t A){fA = A ;}
  void     SetSlopeParameter(Float_t B){fB = B ;}
  void     SetSDigitsBranch(const char * title ) ;

  Bool_t   operator == (const AliPHOSSDigitizer & sd) const ;

private:
  void     Init() ;
  void     PrintSDigits(Option_t * option) ;

private:
  Float_t fA ;              //Pedestal parameter
  Float_t fB ;              //Slope Digitizition parameters
  Int_t   fNevents ;        // Number of events to digitize
  Float_t fPrimThreshold ;  // To store primari if Elos > threshold
  TString fSDigitsTitle ;   // title of SDigits branch
  TString fHeadersFile ;    //input file
  Bool_t         fIsInitialized ; 
  TClonesArray * fSDigits ; //! list of SDigits
  TClonesArray * fHits ;    //!


  ClassDef(AliPHOSSDigitizer,1)  // description 

};

#endif // AliPHOSSDigitizer_H
