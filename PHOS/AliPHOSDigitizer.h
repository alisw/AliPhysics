#ifndef ALIPHOSDigitizer_H
#define ALIPHOSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Task Class for making SDigits in PHOS      
// Class performs digitization of Summable digits (in the PHOS case this is just
// sum of contributions of all primary particles into given cell). 
// In addition it performs mixing of summable digits from different events.
//                  
//*-- Author: Dmitri Peressounko(SUBATECH & KI)


// --- ROOT system ---
#include "TTask.h"
#include "TObjString.h"
class TArrayI ;
// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSSDigitizer ;


class AliPHOSDigitizer: public TTask {

public:
  AliPHOSDigitizer() ;          // ctor
  AliPHOSDigitizer(const char *headerFile, const char * name = "No Name") ; 
  AliPHOSDigitizer(const AliPHOSDigitizer & dtizer) {( (AliPHOSDigitizer &)dtizer ).Copy(*this) ;} // cpy ctor
  virtual ~AliPHOSDigitizer() ;       

  void    Digitize(const Int_t event) ;            // Make Digits from SDigits 
  void    Exec(Option_t *option);                // Supervising method

  const Float_t GetCPVNoise()     const { return fCPVNoise ;}
  const Float_t GetCPVThreshold() const { return fCPVDigitThreshold ;}
  const Float_t GetEMCThreshold() const { return fEMCDigitThreshold;}
  const Float_t GetPedestal()     const { return fPedestal; }
  const Float_t GetPinNoise()     const { return fPinNoise;}
  const Float_t GetPPSDNoise()    const { return fPPSDNoise ;}
  const Float_t GetPPSDThreshold()const { return fPPSDDigitThreshold ;}
  const Float_t GetSlope()        const { return fSlope; }
  //  const TArrayI      * GetCurrentEvents()const { return fIevent ;}

  void    MixWith(const char* HeaderFile) ; // Add another one file to mix
  virtual void    Print(Option_t* option)const ;
  void    Reset() ;   //restarts starts event processing from 0 event(s)

  void    SetCPVNoise(Float_t CPVNoise)          {fCPVNoise = CPVNoise;}
  void    SetCPVThreshold(Float_t CPVThreshold)  {fCPVDigitThreshold= CPVThreshold;}
  void    SetEMCThreshold(Float_t EMCThreshold)  {fEMCDigitThreshold = EMCThreshold;}
  void    SetPinNoise(Float_t PinNoise )         {fPinNoise = PinNoise;}
  void    SetPPSDNoise(Float_t PPSDNoise)        {fPPSDNoise = PPSDNoise;}
  void    SetPPSDThreshold(Float_t PPSDThreshold){fPPSDDigitThreshold = PPSDThreshold;}

  void    SetSDigitsBranch(const char* file) ;

  AliPHOSDigitizer & operator = (const AliPHOSDigitizer & rvalue)  {
    // assignement operator requested by coding convention but not needed
    abort() ;
    return *this ; 
  }

private:
  void    Init() ;                   
  void    PrintDigits(Option_t * option) ;
  Bool_t  ReadSDigits(Int_t evt) ;            // Read sdigits for particular events
  void    WriteDigits(Int_t evt) ;            // Writes Digits for particular event

private:

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
