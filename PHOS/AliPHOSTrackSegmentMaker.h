#ifndef ALIPHOSTRACKSEGMENTMAKER_H
#define ALIPHOSTRACKSEGMENTMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////
//  Subtrackin class for PHOS                    //
//  Version SUBATECH                             //
//  Author Dmitri Peressounko RRC Ki             //
//     comment: finds pairs of clusters EMC+PPSD //  
//              performs unfolding.              //
///////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "TObjArray.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"

typedef TObjArray TrackSegmentsList ;

class  AliPHOSTrackSegmentMaker : public TObject {

public:

  AliPHOSTrackSegmentMaker() ;                     
  
  virtual ~ AliPHOSTrackSegmentMaker(){}  // dtor

  virtual void MakeTrackSegments(DigitsList * DL, RecPointsList * emcl, RecPointsList * ppsdl, TrackSegmentsList * trsl ) = 0 ; 
                                         // does the job
  virtual void SetMaxEmcPpsdDistance(Float_t r) = 0 ; 
  virtual void SetUnfoldFlag() = 0 ;
  virtual void UnsetUnfoldFlag() = 0 ;

  ClassDef( AliPHOSTrackSegmentMaker,1)  // subtracking implementation , version 1

};

#endif // ALIPHOSTRACKSEGMENTMAKER_H
