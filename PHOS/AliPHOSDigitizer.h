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
class TObjString ;
class TArrayI ;
// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSSDigitizer ;


class AliPHOSDigitizer: public TTask {

public:
  AliPHOSDigitizer() ;          // ctor
  AliPHOSDigitizer(const char *headerFile,const char * sDigitsBranchTitle = 0) ; 
  virtual ~AliPHOSDigitizer() ;       

  void    Digitize(Option_t *option);            // Make Digits from SDigits stored in fSDigits
  void    Exec(Option_t *option);                // Supervising method

  Float_t GetCPVNoise()     const { return fCPVNoise ;}
  Float_t GetCPVThreshold() const { return fCPVDigitThreshold ;}
  Float_t GetEMCThreshold() const { return fEMCDigitThreshold;}
  Float_t GetPedestal()     const { return fPedestal; }
  Float_t GetPinNoise()     const { return fPinNoise;}
  Float_t GetPPSDNoise()    const { return fPPSDNoise ;}
  Float_t GetPPSDThreshold()const { return fPPSDDigitThreshold ;}
  Float_t GetSlope()        const { return fSlope; }
  char *  GetDigitsBranch ()const { return (char*)fDigitsTitle.Data() ;}
  TClonesArray * GetHeadersFiles(){ return fHeaderFiles ;}
  TArrayI*    GetCurrentEvents()  { return fIevent ;}

  void    MixWith(char* HeaderFile, char* SDigitsTitle =0) ; // Add another one file to mix
  virtual void    Print(Option_t* option)const ;
  void    Reset() ;   //restarts starts event processing from 0 event(s)

  void    SetCPVNoise(Float_t CPVNoise)          {fCPVNoise = CPVNoise;}
  void    SetCPVThreshold(Float_t CPVThreshold)  {fCPVDigitThreshold= CPVThreshold;}
  void    SetEMCThreshold(Float_t EMCThreshold)  {fEMCDigitThreshold = EMCThreshold;}
  void    SetPinNoise(Float_t PinNoise )         {fPinNoise = PinNoise;}
  void    SetPPSDNoise(Float_t PPSDNoise)        {fPPSDNoise = PPSDNoise;}
  void    SetPPSDThreshold(Float_t PPSDThreshold){fPPSDDigitThreshold = PPSDThreshold;}

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
  AliPHOSSDigitizer * fSDigitizer ; // ! SDigitizer to extarct some sdigitizing parameters
  Int_t   fNinputs ;                // Number of input files
  Bool_t  fInitialized ;            // 
  TArrayI * fIevent ;               // events to read at the next ReadSDigits() call
  TArrayI * fIeventMax ;            // Maximal number of events in each input file

  Float_t fPedestal ;                // Calibration parameters 
  Float_t fSlope ;                   // read from SDigitizer

  Float_t fPinNoise ;               // Electronics noise in EMC
  Float_t fEMCDigitThreshold  ;     // Threshold for storing digits in EMC
  Float_t fCPVNoise ;               // Noise in CPV
  Float_t fCPVDigitThreshold  ;     // Threshold for storing digits in CPV
  Float_t fPPSDNoise ;              // Noise in PPSD
  Float_t fPPSDDigitThreshold ;     // Threshold for storing digits in PPSD


  ClassDef(AliPHOSDigitizer,1)  // description 

};


#endif // AliPHOSDigitizer_H
