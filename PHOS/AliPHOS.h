#ifndef ALIPHOS_H
#define ALIPHOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for PHOS     
//                  
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH)

// --- ROOT system ---
#include "TString.h"

// --- AliRoot header files ---

#include "AliDetector.h"
#include "AliPHOSGeometry.h" 
#include "AliRecPoint.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"

class AliPHOS : public AliDetector {

 public:

  AliPHOS(const char* name, const char* title): AliDetector(name,title) {} 
  AliPHOS() : AliDetector() {
    // default ctor
  } 
  virtual ~AliPHOS() ; 
 
  virtual void CreateMaterials() ;                                     // defines the material of the detector
  virtual AliPHOSGeometry * GetGeometry() = 0 ;  
  AliPHOSRecPoint::RecPointsList* EmcRecPoints(Int_t evt=0) ;          // gets Array of cluster in the crystals 
  AliPHOSRecParticle::RecParticlesList* RecParticles(Int_t evt = 0) ; // gets Array of reconstructed particles
  AliPHOSTrackSegment::TrackSegmentsList* TrackSegments(Int_t evt=0) ;// gets Array of track segments
  virtual AliPHOSRecPoint::RecPointsList* PpsdRecPoints(Int_t evt=0)=0;// gets Array of clusters in the PPSD 
  virtual TString Version() {
    // gives the version number (void for the base class)
    return TString(" ") ; } 

 protected:
  
  AliPHOSRecPoint::RecPointsList * fEmcRecPoints ;         // The RecPoints (clusters) list in EMC 
  AliPHOSTrackSegment::TrackSegmentsList * fTrackSegments ;// The TrackSegment list in PHOS
  AliPHOSRecParticle::RecParticlesList * fRecParticles ;   // The reconstructed particles list in PHOS


  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
