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
  
  // void TestTowerID(void) ;
  virtual void  Exec(Option_t *option); 
  
  Float_t  GetPedestalParameter()const {return fA;}
  Float_t  GetCalibrationParameter()const{return fB;}
  char *   GetSDigitsBranch()const{return (char*) fSDigitsTitle.Data();}  
  void     SetSplitFile(const TString splitFileName = "EMCAL.SDigits.root" ) ;
  virtual void Print(Option_t* option) const ;

  void     SetPedestalParameter(Float_t A){fA = A ;}
  void     SetSlopeParameter(Float_t B){fB = B ;}
  void     SetSDigitsBranch(const char * title ) ;

  Bool_t   operator == (const AliEMCALSDigitizer & sd) const ;
  Int_t    Segment2TowerID(Int_t SegmentID){
  return Layer2TowerID(SegmentID,kFALSE);
}

private:
  void     Init() ;
  void     InitParameters() ; 
  void     PrintSDigits(Option_t * option) ;
  Int_t    Layer2TowerID(Int_t,Bool_t) ;

private:
  Float_t fA ;              //Pedestal parameter
  Float_t fB ;              //Slope Digitizition parameters
  Float_t fPhotonElectronFactor ;  // number of photon electrons per GeV
  // should be calculated independently for each layer as : 
  // LightYield*LightCollectionEfficiency*LightAttenuation*APDPhotoElectronEfficiency*APDGain

  Bool_t fDefaultInit;      //! Says if the task was created by defaut ctor (only parameters are initialized)
  Int_t   fNevents ;        // Number of events to digitize
  Float_t fTowerPrimThreshold ;  // To store primary in Tower if Elos > threshold
  Float_t fPreShowerPrimThreshold ;  // To store primary if Pre Shower Elos > threshold
  TString fSDigitsTitle ;   // title of SDigits branch
  TString fHeadersFile ;    //input file
  Bool_t         fIsInitialized ; 
  TClonesArray * fSDigits ; //! list of SDigits
  TClonesArray * fHits ;    //! 
  TFile * fSplitFile ;      //! file in which SDigits will eventually be stored

  ClassDef(AliEMCALSDigitizer,2)  // description 

};

#endif // AliEMCALSDigitizer_H
