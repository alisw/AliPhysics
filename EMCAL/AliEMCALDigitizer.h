#ifndef ALIEMCALDigitizer_H
#define ALIEMCALDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Task Class for making Digits in EMCAL      
//                  
//*-- Author: Sahal Yacoob (LBL)
// based on : AliPHOSDigit
//_________________________________________________________________________ 


// --- ROOT system ---
#include "TTask.h"
#include "TObjString.h"
class TArrayI ;
// --- Standard library ---

// --- AliRoot header files ---
class AliEMCALSDigitizer ;


class AliEMCALDigitizer: public TTask {

public:
  AliEMCALDigitizer() ;          // ctor
  AliEMCALDigitizer(const char *headerFile,const char * sDigitsBranchTitle = "Default") ; 
  virtual ~AliEMCALDigitizer() ;       

  void    Digitize(Option_t *option);            // Make Digits from SDigits stored in fSDigits
  void    Exec(Option_t *option);                // Supervising method

  Float_t GetEMCThreshold() const { return fEMCDigitThreshold;}
  Float_t GetPedestal()     const { return fPedestal; }
  Float_t GetPinNoise()     const { return fPinNoise;}
  Float_t GetSlope()        const { return fSlope; }
  char *  GetDigitsBranch ()const { return (char*)fDigitsTitle.Data() ;}
  TClonesArray * GetHeadersFiles(){ return fHeaderFiles ;}
  TArrayI*    GetCurrentEvents()  { return fIevent ;}
  Int_t   DigitizeEnergy(Float_t energy, Int_t absId) ;
  void    MixWith(char* HeaderFile, char* SDigitsTitle =0) ; // Add another one file to mix
  virtual void    Print(Option_t* option)const ;
  void    Reset() ;   //restarts starts event processing from 0 event(s)

  void    SetEMCThreshold(Float_t EMCThreshold)  {fEMCDigitThreshold = EMCThreshold;}
  void    SetPinNoise(Float_t PinNoise )         {fPinNoise = PinNoise;}

  void    SetDigitsBranch (const char* file) ;
  void    SetSDigitsBranch(const char* file) ;

private:
  Bool_t  Combinator() ;                         // makes all desirable combination sig+Bg
  void    Init();                   
  Bool_t  ReadSDigits() ;            // Read sdigits for particular events
  void    WriteDigits() ;            // Writes Digits for particular event
  void    PrintDigits(Option_t * option) ;

private:
  TClonesArray * fSDigitsTitles ;   // Titles of sdigits branches 
  TClonesArray * fHeaderFiles ;     // Names of files with headers to merge 
  TString        fDigitsTitle ;     // Title of the Digits Branch  
  TClonesArray * fSDigits ;         // ! Lists of SDigits 
  TClonesArray * fDigits ;          // ! Final list of digits
  AliEMCALSDigitizer * fSDigitizer ; // ! SDigitizer to extract some sdigitizing parameters
  Int_t   fNinputs ;                // Number of input files
  Bool_t  fInitialized ;            // 
  TArrayI * fIevent ;               // events to read at the next ReadSDigits() call
  TArrayI * fIeventMax ;            // Maximal number of events in each input file

  Float_t fPedestal ;                // Calibration parameters 
  Float_t fSlope ;                   // read from SDigitizer

  Float_t fPinNoise ;               // Electronics noise in EMC
  Float_t fEMCDigitThreshold  ;     // Threshold for storing digits in EMC


  ClassDef(AliEMCALDigitizer,1)  // description 

};


#endif // AliEMCALDigitizer_H
