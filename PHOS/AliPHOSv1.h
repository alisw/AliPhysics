#ifndef ALIPHOSV1_H
#define ALIPHOSV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Implementation version v1 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
// Layout EMC + CPV  has name IHEP
//*--                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
#include "TClonesArray.h"

class TFile;

// --- AliRoot header files ---
#include "AliPHOSv0.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSPID.h"
#include "AliPHOSCPVModule.h"
#include "AliPHOSCPVHit.h"
#include "AliPHOSCPVDigit.h"

class AliPHOSv1 : public AliPHOSv0 {

public:

  AliPHOSv1(void) ;
  AliPHOSv1(const char *name, const char *title="") ;
  AliPHOSv1(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title="") ;
  AliPHOSv1(const AliPHOSv1 & phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  }
  virtual ~AliPHOSv1(void) ;

  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits, Int_t pid, TLorentzVector p, Float_t *pos) ; 
  Float_t        Calibrate(Int_t amp){ return (amp - fDigitizeA)/fDigitizeB ; }
  Int_t          Digitize(Float_t Energy){ return (Int_t ) (fDigitizeA + Energy*fDigitizeB); }
  //  virtual void   Hit2Digit(Int_t event) ;
  virtual void   Hits2SDigits() ;
  virtual void   MakeBranch(Option_t* opt, char *file=0 ) ;
  void           Reconstruction(AliPHOSReconstructioner * Reconstructioner) ;
  void           ResetClusters(){} ;
  virtual void   ResetHits() ; 
  virtual void   SDigits2Digits() ;  
  virtual Int_t  IsVersion(void) const {
    // Gives the version number 
    return 1 ; 
  }

  virtual void   ResetReconstruction() ; // Reset reconstructed objects
  void           SetReconstructioner(AliPHOSReconstructioner& Reconstructioner) {
    // sets the reconstructionner object to be used
    fReconstructioner = &Reconstructioner ;
  }  
  void           SetDigitThreshold(Float_t th) { fDigitThreshold = th ; } 
  void           SetPpsdEnergyThreshold(Float_t enth)      { fPpsdEnergyThreshold = enth ; } 
  void           SetCpvEnergyThreshold(Float_t enth)       { fCpvEnergyThreshold = enth ; } 

  virtual void   SetTreeAddress(); 
  virtual void   StepManager(void) ;                              
  virtual TString Version(void){ 
    // returns the version number 
    return TString("v1") ; 
  }

  AliPHOSv1 & operator = (const AliPHOSv1 & rvalue)  {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *this ; 
  }

  // IHEP's CPV specific functions

  //  AliPHOSCPVModule &GetEMCModule(int n) { return *(AliPHOSCPVModule*)fEMCModules->operator[](n); }
  //  AliPHOSCPVModule &GetCPVModule(int n) { return *(AliPHOSCPVModule*)fCPVModules->operator[](n); }

  void       CPVDigitize (TLorentzVector p, Float_t *xy, Int_t moduleNumber, TClonesArray *digits) ;
  Float_t    CPVPadResponseFunction(Float_t qhit, Float_t zg, Float_t xg) ;
  Double_t   CPVCumulPadResponse(Double_t x, Double_t y) ;

protected:

  Float_t fDigitThreshold ;                       // Threshold for the digit registration 
  Float_t fPpsdEnergyThreshold;                   //PPSD
  Float_t fCpvEnergyThreshold;                    //CPV  
  Float_t fPinElectronicNoise  ;                  // Electronic Noise in the PIN
  Float_t fDigitizeA ;                            //Parameters of the 
  Float_t fDigitizeB ;                            //digitization 
  Int_t   fnSdigits ; 
  AliPHOSReconstructioner  * fReconstructioner ;  // Clusterization and subtracking procedures
  AliPHOSTrackSegmentMaker * fTrackSegmentMaker ; // Reconstructioner of the PHOS track segment: 2 x PPSD + 1 x EMC
  //  TClonesArray             * fEMCModules;         // Array of EMC modules
  //  TClonesArray             * fCPVModules;         // Array of CPV modules for the IHEP's version of CPV

  ClassDef(AliPHOSv1,1)  // Implementation of PHOS manager class for layout EMC+PPSD

};

#endif // AliPHOSV1_H
