#ifndef ALIPHOSDigitizer_H
#define ALIPHOSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Task Class for making SDigits in PHOS      
//                  
//*-- Author: Dmitri Peressounko(SUBATECH & KI)


// --- ROOT system ---
#include "TTask.h"
#include "TObjString.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSSDigitizer ;

class AliPHOSDigitizer: public TTask {

public:
  AliPHOSDigitizer() ;          // ctor
  AliPHOSDigitizer(char *HeaderFile,char * SDigitsBrancheFile = 0) ; 
  virtual ~AliPHOSDigitizer() ;       
  Bool_t  Combinator() ;                         // makes all desirable combination sig+Bg
  void    Exec(Option_t *option);                // Does the job
  void    Digitize(Option_t *option);            // Digitizes SDigits stored in fSDigits
  Float_t GetPinNoise()     const { return fPinNoise;}
  Float_t GetEMCThreshold() const { return fEMCDigitThreshold;}
  Float_t GetCPVNoise()     const { return fCPVNoise ;}
  Float_t GetCPVThreshold() const { return fCPVDigitThreshold ;}
  Float_t GetPPSDNoise()    const { return fPPSDNoise ;}
  Float_t GetPPSDThreshold()const { return fPPSDDigitThreshold ;}
  void    MixWith(char* HeaderFile, Int_t isOutFile = 1, char* SDigitsFile =0) ; // Add another one file to mix
  void    Print(Option_t* option) ;
  void    SetPinNoise(Float_t PinNoise )         {fPinNoise = PinNoise;}
  void    SetEMCThreshold(Float_t EMCThreshold)  {fEMCDigitThreshold = EMCThreshold;}
  void    SetCPVNoise(Float_t CPVNoise)          {fCPVNoise = CPVNoise;}
  void    SetCPVThreshold(Float_t CPVThreshold)  {fCPVDigitThreshold= CPVThreshold;}
  void    SetPPSDNoise(Float_t PPSDNoise)        {fPPSDNoise = PPSDNoise;}
  void    SetPPSDThreshold(Float_t PPSDThreshold){fPPSDDigitThreshold = PPSDThreshold;}

private:
  void   Init(Int_t isOutFile);                   
  Bool_t ReadSDigits() ;            // Read sdigits for particular events
  void   WriteDigits() ;            // Writes Digits for particular event

private:
  TClonesArray * fSDigitsFiles ;    // Names of sdigits branches
  TClonesArray * fHeaderFiles ;     // Names of files with headers to merge
  TString        fDigitsFile ;      // Name of the Digits Branch  
  Int_t          fOutFileNumber ;   // Number of the header file into which Digits are written
  TClonesArray * fSDigits ;         // ! Lists of SDigits 
  TClonesArray * fDigits ;          // ! Final list of digits
  AliPHOSSDigitizer * fSDigitizer ; // ! SDigitizer to extarct some digitizing parameters
  Int_t   fNinputs ;                // Number of input files
  Bool_t  fInitialized ;            // if background file already read?
  TArrayI * fIevent ;               // events to read at the next ReadSDigits() call
  TArrayI * fIeventMax ;            // Maximal number of events in each input file

  Float_t fPinNoise ;               // Electronics noise in EMC
  Float_t fEMCDigitThreshold  ;     // Threshold for storing digits in EMC
  Float_t fCPVNoise ;               // Noise in CPV
  Float_t fCPVDigitThreshold  ;     // Threshold for storing digits in CPV
  Float_t fPPSDNoise ;              // Noise in PPSD
  Float_t fPPSDDigitThreshold ;     // Threshold for storing digits in PPSD


  ClassDef(AliPHOSDigitizer,1)  // description 

};


#endif // AliPHOSDigitizer_H
