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
  AliPHOS(const AliPHOS & phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  }
  virtual ~AliPHOS() ; 
 
  virtual  void              CreateMaterials() ;                     // defines the material of the detector
  virtual  AliPHOSRecPoint::RecPointsList *  EmcRecPoints() {
    // Getting list of RecPoints
    return fEmcRecPoints ;
  }
  virtual  AliPHOSGeometry * GetGeometry() = 0 ;
  virtual   AliPHOSRecPoint::RecPointsList *    PpsdRecPoints()=0;// gets Array of clusters in the PPSD 
  virtual void               SetTreeAddress();                 // Tree Address for reconstruction lists
  virtual  AliPHOSRecParticle::RecParticlesList *  RecParticles() {
    // Getting list of RecParticles
    return fRecParticles ;
  }
  virtual  AliPHOSTrackSegment::TrackSegmentsList *  TrackSegments() {
    // Getting list of TrackSegments
    return fTrackSegments ;
  }
  virtual TString Version() {return TString(" ") ; } 
 
  AliPHOS & operator = (const AliPHOS & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }
 
 protected:
  
  AliPHOSRecPoint::RecPointsList * fEmcRecPoints ;         // The RecPoints (clusters) list in EMC 
  AliPHOSTrackSegment::TrackSegmentsList * fTrackSegments ;// The TrackSegment list in PHOS
  AliPHOSRecParticle::RecParticlesList * fRecParticles ;   // The reconstructed particles list in PHOS


  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
