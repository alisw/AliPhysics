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

  AliPHOS() ;
  AliPHOS(const char* name, const char* title="");
  AliPHOS(const AliPHOS & phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  }
  virtual ~AliPHOS() ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not used this definition but the one below
    assert(0==1) ; 
  }
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits ) = 0 ;   
  virtual void   CreateMaterials() ;                     
  AliPHOSRecPoint::RecPointsList *  EmcRecPoints() const {
    // Getting list of RecPoints
    return fEmcRecPoints ;
  }
  virtual  AliPHOSGeometry * GetGeometry() = 0 ;

  Int_t   IsVersion(void) const { return -1 ; } 
  AliPHOSRecPoint::RecPointsList * PpsdRecPoints() const {
    // to be redefined when ppsd is present
    return  fPpsdRecPoints ;
  } 
  virtual void  SetTreeAddress();                
  AliPHOSRecParticle::RecParticlesList *  RecParticles() const {
    // Getting list of RecParticles
    return fRecParticles ;
  }
  TClonesArray *SDigits() const {return fSDigits;}

  AliPHOSTrackSegment::TrackSegmentsList *  TrackSegments() const {
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
  TClonesArray                           *fSDigits      ; // List of summable digits
  AliPHOSRecPoint::RecPointsList         *fEmcRecPoints ; // The RecPoints (clusters) list in EMC 
  AliPHOSRecPoint::RecPointsList         *fPpsdRecPoints ;// The RecPoints (clusters) list in PPSD (veto)
  AliPHOSTrackSegment::TrackSegmentsList *fTrackSegments ;// The TrackSegment list in PHOS
  AliPHOSRecParticle::RecParticlesList   *fRecParticles ; // The reconstructed particles list in PHOS

  ClassDef(AliPHOS,2) // Photon Spectrometer Detector (base class)

} ;

#endif // ALIPHOS_H
