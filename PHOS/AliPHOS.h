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

// --- AliRoot header files ---

#include "AliDetector.h"
#include "AliPHOSGeometry.h" 
#include "AliRecPoint.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"

class AliPHOS : public AliDetector {

 public:

  AliPHOS(const char* name, const char* title): AliDetector(name,title) {} 
  AliPHOS() : AliDetector() {} 
  virtual ~AliPHOS() ; 
 
  virtual void CreateMaterials() ;               // defines the material of the detector
  virtual AliPHOSGeometry * GetGeometry() = 0 ;  
  RecPointsList* EmcRecPoints() {return fEmcClusters;}               // gets Array of cluster in the crystals 
  RecParticlesList * RecParticles() { return fRecParticles ; }      // gets Array of reconstructed particles
  TrackSegmentsList *    TrackSegments(){return fTrackSegments ;} // gets Array of track segments
  virtual RecPointsList* PpsdRecPoints() = 0 ;        // gets Array of clusters in the PPSD 

 protected:
  
  RecPointsList * fEmcClusters ;                  // The RecPoints (clusters) list in EMC 
  TrackSegmentsList * fTrackSegments ;            // The TrackSegment list in PHOS
  RecParticlesList * fRecParticles ;              // The reconstructed particles list in PHOS


  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
