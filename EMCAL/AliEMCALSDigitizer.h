#ifndef ALIEMCALSDigitizer_H
#define ALIEMCALSDigitizer_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
//_________________________________________________________________________
//  Task Class for making SDigits in EMCAL      
//                  
//*-- Author: Sahal Yacoob (LBL)
// based on : AliPHOSSDigitizer
//_________________________________________________________________________
//
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//  July   2003 Yves Schutz: new  IO (à la PHOS)
 
// --- ROOT system ---
#include "TTask.h"
class TFile ;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h"

class AliEMCALSDigitizer: public TTask {

public:
  AliEMCALSDigitizer() ;          // ctor
  AliEMCALSDigitizer(const char * alirunFileName, const char * eventFolderName = AliConfig::fgkDefaultEventFolderName) ; 
  AliEMCALSDigitizer(const AliEMCALSDigitizer & sd) ;
  virtual ~AliEMCALSDigitizer() {;} // dtor

  Float_t       Calibrate(Int_t amp)const {return (amp - fA)/fB ; }
  Int_t         Digitize(Float_t Energy)const { return (Int_t ) ( fA + Energy*fB); }
  virtual void  Exec(Option_t *option); 
  Int_t         GetSDigitsInRun() const {return fSDigitsInRun ;}  
  virtual void  Print(Option_t*) const ;
  void          SetEventFolderName(TString name) { fEventFolderName = name ; }

  Bool_t operator == (const AliEMCALSDigitizer & sd) const ;
  AliEMCALSDigitizer & operator = (const AliEMCALSDigitizer & sd) {return *this ;}


private:
  void     Init() ;
  void     InitParameters() ; 
  void     PrintSDigits(Option_t * option) ;
  void     Unload() const ;

private:
  Float_t fA ;                     // Pedestal parameter
  Float_t fB ;                     // Slope Digitizition parameters
  Float_t fPREPrimThreshold ;      // To store primary if Pre Shower Elos > threshold
  Float_t fECPrimThreshold ;       // To store primary if EC Shower Elos > threshold
  Float_t fHCPrimThreshold ;       // To store primary if HC Shower Elos > threshold
  Bool_t  fDefaultInit;            //! Says if the task was created by defaut ctor (only parameters are initialized)
  TString fEventFolderName;        // event folder name
  Bool_t  fInit ;                  //! tells if initialisation wennt OK, will revent exec if not
  Int_t   fSDigitsInRun ;          //! Total number of sdigits in one run

  ClassDef(AliEMCALSDigitizer,5)  // description 

};

#endif // AliEMCALSDigitizer_H

