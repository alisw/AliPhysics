#ifndef ALIPHOSV1_H
#define ALIPHOSV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Implementation version v1 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOSv0.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSPID.h"

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

  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits ) ; 
  Int_t          Digitize(Float_t Energy);
  virtual void   FinishEvent(void) ;                               
  Int_t IsVersion(void) const { return 1 ; }
  virtual void   MakeBranch(Option_t* opt) ;
  virtual  AliPHOSRecPoint::RecPointsList *  PpsdRecPoints() {
    // Getting list of PPSD RecPoints
    return fPpsdRecPoints ;
  }
  void           Reconstruction(AliPHOSReconstructioner * Reconstructioner) ;
  void           ResetClusters(){} ;
  virtual void   ResetDigits() ; 
  virtual void   ResetReconstruction() ; // Reset reconstructed objects
  void           SetReconstructioner(AliPHOSReconstructioner& Reconstructioner) {
    // sets the reconstructionner object to be used
    fReconstructioner = &Reconstructioner ;
  }  
  void           SetDigitThreshold(Float_t th) { fDigitThreshold = th ; } 
  virtual void   SetTreeAddress(); 
  virtual void   StepManager(void) ;                              
  virtual TString Version(void){ 
    // returns the version number 
    return TString("v1") ; 
  }

  AliPHOSv1 & operator = (const AliPHOSv1 & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }

protected:

  Float_t fDigitThreshold ;                       // Threshold for the digit registration 
  Int_t fNTmpHits ;                               //!  Used internally for digitalization
  Float_t fPinElectronicNoise  ;                  // Electronic Noise in the PIN
  AliPHOSRecPoint::RecPointsList * fPpsdRecPoints ; // The RecPoints (clusters) list in PPSD 
  AliPHOSReconstructioner * fReconstructioner ;   // Reconstrutioner of the PHOS event: Clusterization and subtracking procedures
  TClonesArray * fTmpHits ;                       //!  Used internally for digitalization 
  AliPHOSTrackSegmentMaker * fTrackSegmentMaker ; // Reconstructioner of the PHOS track segment: 2 x PPSD + 1 x EMC

  ClassDef(AliPHOSv1,1)  // Implementation of PHOS manager class for layout EMC+PPSD

};

#endif // AliPHOSV1_H
