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


// --- ROOT system ---
#include "TTask.h"
#include "TString.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALSDigitizer: public TTask {

public:
  AliEMCALSDigitizer() ;          // ctor
  AliEMCALSDigitizer(const char* HeaderFile,const char *SdigitsTitle = "Default") ; 
  virtual ~AliEMCALSDigitizer() ; // dtor

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

  Bool_t   operator == (const AliEMCALSDigitizer & sd) const ;

private:
  void     Init() ;
  void     PrintSDigits(Option_t * option) ;
  Int_t    Layer2TowerID(Int_t,Bool_t) ;
private:
  Float_t fA ;              //Pedestal parameter
  Float_t fB ;              //Slope Digitizition parameters
  Float_t fLayerRatio ;     //Factor that takes into account difference in light collection between 2 first layers and rest of layers

  Int_t   fNevents ;        // Number of events to digitize
  Float_t fPrimThreshold ;  // To store primary if Elos > threshold
  TString fSDigitsTitle ;   // title of SDigits branch
  TString fHeadersFile ;    //input file
  Bool_t         fIsInitialized ; 
  TClonesArray * fSDigits ; //! list of SDigits
  TClonesArray * fHits ;    //! 

  ClassDef(AliEMCALSDigitizer,1)  // description 

};

#endif // AliEMCALSDigitizer_H
