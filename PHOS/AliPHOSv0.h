#ifndef ALIPHOSV0_H
#define ALIPHOSV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Manager class  for PHOS                   //
//  Version SUBATECH                          //
//  Author  Y. Schutz SUBATECH                //
//       geometry parametrized for any        //  
//       shape of modules                     //
////////////////////////////////////////////////

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
  RecPointsList* EmcClusters() {return fEmcClusters;}               // gets TClonesArray of cluster in the crystals 
  void           FinishEvent(void) ;                                // makes the digits from the hits 
  virtual AliPHOSGeometry * GetGeometry() { return fGeom ; }  
  virtual void   Init(void) ;                                       // does nothing
  Int_t IsVersion(void) const { return 0 ; }
  void           MakeBranch(Option_t* opt) ;
  RecPointsList* PpsdClusters() { return fPpsdClusters ; }          // gets TClonesArray of clusters in the PPSD 
  void           Reconstruction(AliPHOSReconstructioner * Reconstructioner) ;
  RecParticlesList * RecParticles() { return fRecParticles ; }      // gets TClonesArray of reconstructed particles
  void           ResetClusters(){} ;
  void           SetReconstructioner(AliPHOSReconstructioner& Reconstructioner) {fReconstructioner = &Reconstructioner ;} 
  virtual void   StepManager(void) ;                                // does the tracking through PHOS and a preliminary digitalization
  TrackSegmentsList *    TrackSegments(){return fTrackSegments ;}
  
protected:
  Float_t fPINElectronicNoise  ;         // Electronic Noise in the PIN
  RecPointsList * fEmcClusters ;        //!  (!=do not stream)
  AliPHOSGeometry * fGeom ;             // geometry definition
  Int_t fNTmpHits ;                     //!  used internally for digitalization
  RecPointsList * fPpsdClusters ;       //!
  AliPHOSReconstructioner * fReconstructioner ; // Reconstrutioner of the PHOS event: Clusterization and subtracking procedures
  TClonesArray * fTmpHits ;             //!  idem
  AliPHOSTrackSegmentMaker * fTrackSegmentMaker ;
  TrackSegmentsList * fTrackSegments ;  //! idem
  RecParticlesList * fRecParticles ;    //! idem

  ClassDef(AliPHOSv0,1)  // PHOS main class , version subatech

};

#endif // AliPHOSV0_H
