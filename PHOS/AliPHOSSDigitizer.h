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
  AliPHOSSDigitizer(char* HeaderFile,char *SdigitsFile = 0) ; 

  virtual ~AliPHOSSDigitizer() ; // dtor
  Float_t  Calibrate(Int_t amp){return (amp - fA)/fB ; }
  Int_t    Digitize(Float_t Energy){ return (Int_t ) ( fA + Energy*fB); }

  Float_t GetPedestalParameter(){return fA;}
  Float_t GetCalibrationParameter(){return fB;}
  char *GetSDigitsFile()const{return (char*) fSDigitsFile.Data();}  
  virtual void  Exec(Option_t *option); 
  void SetNEvents(Int_t Nevents){fNevents = Nevents;}
  void SetPedestalParameter(Float_t A){fA = A ;}
  void SetSlopeParameter(Float_t B){fB = B ;}
  void SetSDigitsFile(char * file ) ;
  virtual void Print(Option_t* option) const ;
  Bool_t operator == (const AliPHOSSDigitizer & sd) const ;

private:
  Float_t fA ;              //Pedestal parameter
  Float_t fB ;              //Slope Digitizition parameters
  Int_t   fNevents ;        // Number of events to digitize
  Float_t fPrimThreshold ;  // To store primari if Elos > threshold
  TString fSDigitsFile ;    //output file 
  TString fHeadersFile ;    //input file


  ClassDef(AliPHOSSDigitizer,1)  // description 

};

#endif // AliPHOSSDigitizer_H
