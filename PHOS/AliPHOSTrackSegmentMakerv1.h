#ifndef ALIPHOSTRACKSEGMENTMAKERV1_H
#define ALIPHOSTRACKSEGMENTMAKERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////
//  Track Segment Maker class for PHOS           //
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
#include "AliPHOSTrackSegmentMaker.h"
#include "TMinuit.h" 

class  AliPHOSTrackSegmentMakerv1 : public AliPHOSTrackSegmentMaker {

public:

  AliPHOSTrackSegmentMakerv1() ;                     
  virtual ~ AliPHOSTrackSegmentMakerv1() ; // dtor
  
  Bool_t  FindFit(AliPHOSEmcRecPoint * emcRP, int * MaxAt, Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) ; //Used in UnfoldClusters, calls TMinuit

  void    FillOneModule(DigitsList * Dl, RecPointsList * emcIn, TObjArray * emcOut, RecPointsList * ppsdIn, 
			TObjArray * ppsdOutUp, TObjArray * ppsdOutLow, Int_t &PHOSModule, Int_t & emcStopedAt, 
			Int_t & ppsdStopedAt) ; // Unfolds clusters and fills temporary arrais   

  Float_t GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * EmcClu , AliPHOSPpsdRecPoint * Ppsd , Bool_t & TooFar ) ; // see R0

  void    MakeLinks(TObjArray * EmcRecPoints, TObjArray * PpsdRecPointsUp, TObjArray * PpsdRecPointsLow, 
		    TClonesArray * LinkLowArray, TClonesArray *LinkUpArray) ; //Evaluates distances(links) between EMC and PPSD

  void    MakePairs(TObjArray * EmcRecPoints, TObjArray * PpsdRecPointsUp, TObjArray * PpsdRecPointsLow, 
		    TClonesArray * LinkLowArray, TClonesArray * LinkUpArray, TrackSegmentsList * trsl) ; 
                    //Finds pairs(triplets) with smallest link

  void    MakeTrackSegments(DigitsList * DL, RecPointsList * emcl, RecPointsList * ppsdl, TrackSegmentsList * trsl ) ; // does the job

  void    SetMaxEmcPpsdDistance(Float_t r){ fR0 = r ;} //Radius within which we look for ppsd cluster

  static Double_t ShowerShape(Double_t r) ; // Shape of shower used in unfolding; class member function (not object member function)

  void    UnfoldClusters(DigitsList * DL, RecPointsList * emcIn, AliPHOSEmcRecPoint * iniEmc, Int_t Nmax, 
		         int * maxAt, Float_t * maxAtEnergy, TObjArray * emclist) ; //Unfolds overlaping clusters using TMinuit package

private:

  Float_t fDelta ;    // parameter used for sorting
  Float_t fR0 ;       // Maximal distance between EMC and PPSD clusters of one Track Segment in module plane
  TMinuit * fMinuit ; // Minuit object needed by cluster unfolding

  ClassDef( AliPHOSTrackSegmentMakerv1,1)  // track segment maker implementation , version 1

};

#endif // AliPHOSTRACKSEGMENTMAKERV1_H
