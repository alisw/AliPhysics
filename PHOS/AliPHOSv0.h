#ifndef ALIPHOSV0_H
#define ALIPHOSV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
// Implementation version v0 of PHOS Manager class 
// Layout EMC + PPSD has name GPS2  
//                  
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOS.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSPID.h"

class AliPHOSv0 : public AliPHOS {

public:

  AliPHOSv0(void) ;
  AliPHOSv0(const char *name, const char *title="") ;
  AliPHOSv0(AliPHOSReconstructioner * Reconstructioner, const char *name, const char *title="") ;
  virtual ~AliPHOSv0(void) ;

  virtual void   AddHit( Int_t primary, Int_t id, Float_t *hits ) ; // adds a pre-digitilized hit to the hit tree 
  virtual void   BuildGeometry(void) ;                              // creates the geometry for the ROOT display
  void           BuildGeometryforPHOS(void) ;                       // creates the PHOS geometry for the ROOT display
  void           BuildGeometryforPPSD(void) ;                       // creates the PPSD geometry for the ROOT display
  virtual void   CreateGeometry(void) ;                             // creates the geometry for GEANT
  void           CreateGeometryforPHOS(void) ;                      // creates the PHOS geometry for GEANT
  void           CreateGeometryforPPSD(void) ;                      // creates the PPSD geometry for GEANT
  Int_t          Digitize(Float_t Energy);
  void           FinishEvent(void) ;                                // makes the digits from the hits 
  virtual AliPHOSGeometry * GetGeometry() { return fGeom ; }  
  virtual void   Init(void) ;                                       // does nothing
  Int_t IsVersion(void) const { return 0 ; }
  void           MakeBranch(Option_t* opt) ;
  virtual  AliPHOSRecPoint::RecPointsList *  PpsdRecPoints() {
    // Getting list of PPSD RecPoints
    return fPpsdRecPoints ;
  }
  void           Reconstruction(AliPHOSReconstructioner * Reconstructioner) ;
  void           ResetClusters(){} ;
  virtual void   ResetDigits() ; 
  void           SetReconstructioner(AliPHOSReconstructioner& Reconstructioner) {fReconstructioner = &Reconstructioner ;} 
  void           SetDigitThreshold(Float_t th) { fDigitThreshold = th ; } 
  virtual void   SetTreeAddress(); 
  virtual void   StepManager(void) ;                                // does the tracking through PHOS and a preliminary digitalization
  virtual TString Version(void){ return TString("v0"); }
protected:

  Float_t fDigitThreshold ;                       // Threshold for the digit registration 
  AliPHOSGeometry * fGeom ;                       // Geometry definition
  Int_t fNTmpHits ;                               //!  Used internally for digitalization
  Float_t fPinElectronicNoise  ;                  // Electronic Noise in the PIN
  AliPHOSRecPoint::RecPointsList * fPpsdRecPoints ; // The RecPoints (clusters) list in PPSD 
  virtual void               ResetReconstruction() ; // Reset reconstructed objects
  AliPHOSReconstructioner * fReconstructioner ;   // Reconstrutioner of the PHOS event: Clusterization and subtracking procedures
  TClonesArray * fTmpHits ;                       //!  Used internally for digitalization 
  AliPHOSTrackSegmentMaker * fTrackSegmentMaker ; // Reconstructioner of the PHOS track segment: 2 x PPSD + 1 x EMC

  ClassDef(AliPHOSv0,1)  // Implementation of PHOS manager class for layout EMC+PPSD

};

#endif // AliPHOSV0_H
