#ifndef ALIPHOSRECONSTRUCTIONER_H
#define ALIPHOSRECONSTRUCTIONER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Algorithme class for the reconstruction   //
//                                            //
//  Author Gines MARTINEZ     SUBATECH        //
//                                            //
//  january 2000: added Particle guesser (YS) //
//                                            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSParticleGuesser.h"
#include "TClonesArray.h" 

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSReconstructioner : public TObject {

public:

  AliPHOSReconstructioner(); //ctor            
  AliPHOSReconstructioner(AliPHOSClusterizer * Clusterizer, AliPHOSTrackSegmentMaker * Tracker, 
			  AliPHOSParticleGuesser * Guesser); //ctor            
  ~AliPHOSReconstructioner(); // dtor

  AliPHOSClusterizer * GetClusterizer() { return fClusterizer ; }
  void Init(AliPHOSClusterizer * Clusterizer, AliPHOSTrackSegmentMaker * Tracker, 
			  AliPHOSParticleGuesser * Guesser) ;  
  void Make(TClonesArray * DL, RecPointsList * emccl, RecPointsList * ppsdl, 
	    TrackSegmentsList * trsl, RecParticlesList * rpl) ; // does the job


private:
  
  AliPHOSClusterizer * fClusterizer ;             // Method for clusterization 
  AliPHOSTrackSegmentMaker * fTrackSegmentMaker ; // Method for track segments finding
  AliPHOSParticleGuesser * fParticleGuesser ;     // Method for guessing the type of particle
 
  ClassDef(AliPHOSReconstructioner,1)  // Reconstruction interface , version 1

}; 

#endif // ALIPHOSRECONSTRUCTIONER_H
