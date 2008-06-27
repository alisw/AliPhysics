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
class TBrowser;

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h"

class AliEMCALSDigitizer: public TTask {

public:
  AliEMCALSDigitizer() ;          // ctor
  AliEMCALSDigitizer(const char * alirunFileName, const char * eventFolderName = AliConfig::GetDefaultEventFolderName()) ; 
  AliEMCALSDigitizer(const AliEMCALSDigitizer & sd) ;
  virtual ~AliEMCALSDigitizer(); // dtor

  Int_t         Digitize(Float_t energy)const; //convert energy in GeV to int amplitude
  Float_t       Calibrate(Int_t amp)const;  //opposite of Digitize()

  virtual void  Exec(Option_t *option); 
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
  Bool_t  fDefaultInit;            //! Says if the task was created by defaut ctor (only parameters are initialized)
  TString fEventFolderName;        // event folder name
  Bool_t  fInit ;                  //! tells if initialisation went OK, will revent exec if not
  Int_t   fSDigitsInRun ;          //! Total number of sdigits in one run
  Int_t   fFirstEvent;             // first event to process
  Int_t   fLastEvent;              // last  event to process
  Float_t fSampling;               // See AliEMCALGeometry

  ClassDef(AliEMCALSDigitizer,6)  // description 
};

#endif // AliEMCALSDigitizer_H

