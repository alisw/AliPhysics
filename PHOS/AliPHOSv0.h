#ifndef ALIPHOSXXX_H
#define ALIPHOSXXX_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Short description                         //
//  Author                 SUBATECH           //
//      comment                               //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSxxx {

public:

  virtual ~AliPHOSxxx() ; // dtor

private:

ClassDef(AliPHOSxxx,1)  // description , version 1

};

#endif // AliPHOSXXX_H
//-*-C++-*-
#ifndef ALIPHOSV4_H
#define ALIPHOSV4_H
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

class AliPHOSv0 : public AliPHOS {

public:

  AliPHOSv0(void) ;
  AliPHOSv0(const char *name, const char *title="") ;
  AliPHOSv0(AliPHOSReconstructioner&  Reconstructioner, const char *name, const char *title="") ;
  virtual       ~AliPHOSv0(void) ;

  virtual void   AddHit( Int_t track, Int_t id, Float_t *hits ) ;   // adds a pre-digitilized hit to the hit tree 
  virtual void   BuildGeometry(void) ;                              // creates the geometry for the ROOT display
  void           BuildGeometryforPHOS(void) ;                       // creates the PHOS geometry for the ROOT display
  void           BuildGeometryforPPSD(void) ;                       // creates the PPSD geometry for the ROOT display
  virtual void   CreateGeometry(void) ;                             // creates the geometry for GEANT
  void           CreateGeometryforPHOS(void) ;                      // creates the PHOS geometry for GEANT
  void           CreateGeometryforPPSD(void) ;                      // creates the PPSD geometry for GEANT
  Int_t          Digitize(Float_t Energy);
  RecPointsList* EmcClusters() {return fEmcClusters;}                 // gets TClonesArray of cluster in the crystals 
  void           FinishEvent(void) ;                                // makes the digits from the hits 
  virtual void   Init(void) ;                                       // does nothing
  void           MakeBranch(Option_t* opt) ;
  RecPointsList* PpsdClusters() {return fPpsdClusters;}             // gets TClonesArray of clusters in the PPSD 
  void           Reconstruction(AliPHOSReconstructioner& Reconstructioner) ;
  void           ResetClusters(){} ;
  void           SetReconstructioner(AliPHOSReconstructioner& Reconstructioner) {fReconstructioner = &Reconstructioner;} //
  virtual void   StepManager(void) ;                                // does the tracking through PHOS and a preliminary digitalization
  TObjArray *   TrackSegments(){return fTrackSegments ;}
  // inlines

  virtual AliPHOSGeometry * GetGeometry() { return fGeom ; }  
  Int_t IsVersion(void) const { return 4 ; }

private:

  AliPHOSGeometry  *        fGeom ; // geometry definition
  RecPointsList    *        fEmcClusters;    //!
  Int_t                     fNTmpHits ;     //!  used internally for digitalization (!=do not stream)
  RecPointsList    *        fPpsdClusters;  //!
  TObjArray *               fTrackSegments ;//!
  TClonesArray *            fTmpHits ;      //!  idem
  AliPHOSReconstructioner * fReconstructioner ; // Reconstrutioner of the PHOS event: Clusterization and subtracking procedures
  AliPHOSTrackSegmentMaker       * fTrackSegmentMaker ;
public:

  ClassDef(AliPHOSv0,1)  // PHOS main class , version subatech

};

#endif // AliPHOSV4_H
