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
//                           of new  IO (à la PHOS)
 
// --- ROOT system ---
#include "TTask.h"
#include "TString.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALSDigitizer: public TTask {

public:
  AliEMCALSDigitizer() ;          // ctor
  AliEMCALSDigitizer(const char* headerFile,const char* hdigitsTitle = "Default", const Bool_t toSplit = kFALSE) ; 
  virtual ~AliEMCALSDigitizer() ; // dtor

  Float_t  Calibrate(Int_t amp)const {return (amp - fA)/fB ; }
  Int_t    Digitize(Float_t Energy)const { return (Int_t ) ( fA + Energy*fB); }
  virtual void  Exec(Option_t *option); 
  const char *  GetSDigitsBranch()const{return GetName();}  
  const Int_t   GetSDigitsInRun() const {return fSDigitsInRun ;}  
  const Float_t GetPedestalParameter()const {return fA;}
  const Float_t GetCalibrationParameter()const{return fB;}
  virtual void  Print(Option_t* option) const ;
  void SetSDigitsBranch(const char * title ) ;
  void SetPedestalParameter(Float_t A){fA = A ;}
  void SetSlopeParameter(Float_t B){fB = B ;}
  void UseHitsFrom(const char * filename) ;      
  Bool_t operator == (const AliEMCALSDigitizer & sd) const ;
//   const Int_t Segment2TowerID(Int_t SegmentID){
//     return Layer2TowerID(SegmentID,0);
// }

private:
  void     Init() ;
  void     InitParameters() ; 
  void     PrintSDigits(Option_t * option) ;

private:
  Float_t fA ;                     // Pedestal parameter
  Float_t fB ;                     // Slope Digitizition parameters
  Float_t fPREPrimThreshold ;      // To store primary if Pre Shower Elos > threshold
  Float_t fECPrimThreshold ;       // To store primary if EC Shower Elos > threshold
  Float_t fHCPrimThreshold ;       // To store primary if HC Shower Elos > threshold
  Bool_t  fDefaultInit;            //! Says if the task was created by defaut ctor (only parameters are initialized)
  Int_t   fSDigitsInRun ;          //! Total number of sdigits in one run
  TFile * fSplitFile ;             //! file in which SDigits will eventually be stored
  Bool_t  fToSplit ;               //! Says that sigits should be written into splip file

  ClassDef(AliEMCALSDigitizer,4)  // description 

};

#endif // AliEMCALSDigitizer_H

