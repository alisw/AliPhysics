#ifndef ALIPHOSRECONSTRUCTIONER_H
#define ALIPHOSRECONSTRUCTIONER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Algorithm class for the reconstruction: clusterizer
//                                          track segment maker
//                                          particle identifier   
//*--
//*-- Author: Gines Martinez & Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TObject.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSPID.h"
#include "TClonesArray.h" 

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSReconstructioner : public TObject {

public:

  AliPHOSReconstructioner(){} //ctor            
  AliPHOSReconstructioner(AliPHOSClusterizer * Clusterizer, AliPHOSTrackSegmentMaker * Tracker, 
			  AliPHOSPID * Identifier); //ctor         
  AliPHOSReconstructioner(const AliPHOSReconstructioner & phos) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  }
   
  ~AliPHOSReconstructioner(){} // dtor

  AliPHOSClusterizer * GetClusterizer() { return fClusterizer ; }
  void Init(AliPHOSClusterizer * Clusterizer, AliPHOSTrackSegmentMaker * Tracker, 
			  AliPHOSPID * Identifier) ;  
  void Make(TClonesArray * DL, 
	    AliPHOSRecPoint::RecPointsList * emccl, 
	    AliPHOSRecPoint::RecPointsList * ppsdl, 
	    AliPHOSTrackSegment::TrackSegmentsList * trsl, 
	    AliPHOSRecParticle::RecParticlesList * rpl)    ; // does the job for EMC+PPSD
  void Make(TClonesArray * DL, 
	    AliPHOSRecPoint::RecPointsList * emccl, 
	    AliPHOSRecPoint::RecPointsList * ppsdl) ;        // does the job for EMC+CPV

  void SetDebugReconstruction(Bool_t deb) { fDebugReconstruction = deb; }

  AliPHOSReconstructioner & operator = (const AliPHOSReconstructioner & rvalue)  {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *this ; 
  }
  

private:
  
  Bool_t               fDebugReconstruction;      // For debuging of the Reconstruction procedure
  AliPHOSClusterizer * fClusterizer ;             // Method for clusterization 
  AliPHOSTrackSegmentMaker * fTrackSegmentMaker ; // Method for track segments finding
  AliPHOSPID * fPID ;                             // Method for identifying the type of particle
 
  ClassDef(AliPHOSReconstructioner,1)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIPHOSRECONSTRUCTIONER_H
