#ifndef ALIEMCALSDIGITIZER_H
#define ALIEMCALSDIGITIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
//_________________________________________________________________________
// This is a class that makes SDigits out of Hits
// A Summable Digits is the sum of all hits originating 
// from one in one tower of the EMCAL 
// A threshold for assignment of the primary to SDigit is applied 
//
// SDigits need to hold the energy sum of the hits, but AliEMCALDigit
// can (should) only store amplitude.  Therefore, the SDigit energy is
// "digitized" before being stored and must be "calibrated" back to an
// energy before SDigits are summed to form true Digits
//
//
//*-- Author: Sahal Yacoob (LBL)
// based on : AliPHOSSDigitizer
//_________________________________________________________________________
 
// --- ROOT system ---
#include "TNamed.h"
class TBrowser;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h"

class AliEMCALSDigitizer: public TNamed {

public:
  AliEMCALSDigitizer() ;          // ctor
  AliEMCALSDigitizer(const char * alirunFileName, const char * eventFolderName = AliConfig::GetDefaultEventFolderName()) ; 
  AliEMCALSDigitizer(const AliEMCALSDigitizer & sd) ;
  virtual ~AliEMCALSDigitizer(); // dtor

  Float_t       Digitize(Float_t energy)const; //convert energy in GeV to int amplitude
  Float_t       Calibrate(Float_t amp)const;  //opposite of Digitize()

  virtual void  Digitize(Option_t *option=""); 
  Int_t         GetSDigitsInRun() const {return fSDigitsInRun ;}  
  virtual void  Print(Option_t *option="") const;
  void          Print1(Option_t *option="all");  // *MENU*
  void          SetEventFolderName(TString name) { fEventFolderName = name ; }
  void          SetEventRange(Int_t first=0, Int_t last=-1) {fFirstEvent=first; fLastEvent=last; }

  Bool_t operator == (const AliEMCALSDigitizer & sd) const ;
  const AliEMCALSDigitizer & operator = (const AliEMCALSDigitizer & /*sd*/) {return *this ;}

  virtual void Browse(TBrowser* b);

private:
  void     Init() ;
  void     InitParameters() ; 
  void     PrintSDigits(Option_t * option) ;
  void     Unload() const ;

private:
  Float_t fA ;                     // Pedestal parameter
  Float_t fB ;                     // Slope Digitizition parameters
  Float_t fECPrimThreshold ;       // To store primary if EC Shower Elos > threshold
  Bool_t  fDefaultInit;            //! Says if the object was created by defaut ctor (only parameters are initialized)
  TString fEventFolderName;        // event folder name
  Bool_t  fInit ;                  //! tells if initialisation went OK, will revent exec if not
  Int_t   fSDigitsInRun ;          //! Total number of sdigits in one run
  Int_t   fFirstEvent;             // first event to process
  Int_t   fLastEvent;              // last  event to process
  Float_t fSampling;               // See AliEMCALGeometry
  TClonesArray* fHits;             //-> Temporal array with hits
	
  ClassDef(AliEMCALSDigitizer,8)  // description 
};

#endif // AliEMCALSDIGITIZER_H

