#ifndef ALIPHOSRECONSTRUCTIONER_H
#define ALIPHOSRECONSTRUCTIONER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Algorithme class for the reconstruction   //
//                                            //
//  Author Gines MARTINEZ     SUBATECH        //
//                                            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "TClonesArray.h" 

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSReconstructioner : public TObject {

public:

  AliPHOSReconstructioner(); //ctor            
  AliPHOSReconstructioner(AliPHOSClusterizer& Clusterizer, AliPHOSTrackSegmentMaker& Tracker); //ctor            
  ~AliPHOSReconstructioner(); // dtor

  AliPHOSClusterizer * GetClusterizer() { return fClusterizer ; }
  void Make(TClonesArray * DL, RecPointsList * emccl, RecPointsList * ppsdl, TrackSegmentsList * trsl) ; // does the job


private:
  
  AliPHOSClusterizer * fClusterizer ;             // Method for clusterization 
  AliPHOSTrackSegmentMaker * fTrackSegmentMaker ; // Method for track segments finding 
  

public:
  
  ClassDef(AliPHOSReconstructioner,1)  // Reconstruction interface , version 1

}; 

#endif // ALIPHOSRECONSTRUCTIONER_H
