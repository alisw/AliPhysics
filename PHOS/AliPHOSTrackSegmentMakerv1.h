#ifndef ALIPHOSTRACKSEGMENTMAKERV1_H
#define ALIPHOSTRACKSEGMENTMAKERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Implementation version 1 of algorithm class to construct PHOS track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)

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
  AliPHOSTrackSegmentMakerv1(const AliPHOSTrackSegmentMakerv1 & tsm) {
    // cpy ctor: no implementation yet
    // requested by the Coding Convention
    assert(0==1) ; 
  }
   
  virtual ~ AliPHOSTrackSegmentMakerv1() ; // dtor
  
  Bool_t  FindFit(AliPHOSEmcRecPoint * emcRP, int * MaxAt, Float_t * maxAtEnergy, 
		  Int_t NPar, Float_t * FitParametres) ; //Used in UnfoldClusters, calls TMinuit
  void    FillOneModule(AliPHOSRecPoint::RecPointsList * emcIn, 
			TArrayI * emcOut, 
			AliPHOSRecPoint::RecPointsList * ppsdIn, 
			TArrayI * ppsdOutUp, 
			TArrayI * ppsdOutLow, 
			Int_t &PHOSModule, 
			Int_t & emcStopedAt, 
			Int_t & ppsdStopedAt) ; // Fills temporary arrais with clusters from one module EMC+PPSD
  void    FillOneModule(AliPHOSRecPoint::RecPointsList * emcIn, 
			TArrayI * emcOut, 
			AliPHOSRecPoint::RecPointsList * cpvIn, 
			TArrayI * cpvOut, 
			Int_t & PHOSModule, 
			Int_t & emcStopedAt, 
			Int_t & cpvStopedAt) ;  // Fills temporary arrais with clusters from one module EMC+CPV
  Float_t GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * EmcClu , AliPHOSPpsdRecPoint * Ppsd , Bool_t & TooFar ) ; // see R0

  void    MakeLinks(TArrayI * EmcRecPoints, TArrayI * PpsdRecPointsUp, TArrayI * PpsdRecPointsLow, 
		    TClonesArray * LinkLowArray, TClonesArray *LinkUpArray) ; //Evaluates distances(links) between EMC and PPSD
  void    MakePairs(TArrayI * EmcRecPoints, 
		    TArrayI * PpsdRecPointsUp, 
		    TArrayI * PpsdRecPointsLow, 
		    TClonesArray * LinkLowArray, 
		    TClonesArray * LinkUpArray, 
		    AliPHOSTrackSegment::TrackSegmentsList * trsl) ; //Finds pairs(triplets) with smallest link
  void    MakeTrackSegments(DigitsList * DL, 
			    AliPHOSRecPoint::RecPointsList * emcl, 
			    AliPHOSRecPoint::RecPointsList * ppsdl, 
			    AliPHOSTrackSegment::TrackSegmentsList * trsl ) ; // does the job
  void    MakeTrackSegmentsCPV(DigitsList * DL, 
			    AliPHOSRecPoint::RecPointsList * emcl, 
			    AliPHOSRecPoint::RecPointsList * ppsdl ) ; // just unfold EMC and CPV clusters
  virtual void SetMaxEmcPpsdDistance(Float_t r){ fR0 = r ;}
  virtual void    SetUnfoldFlag() { fUnfoldFlag = kTRUE ; } ; 
  static Double_t ShowerShape(Double_t r) ; // Shape of shower used in unfolding; class member function (not object member function)
  void    UnfoldAll(DigitsList * Dl, AliPHOSRecPoint::RecPointsList * emcIn) ; 
                                            // Unfolds and sorts all EMC clusters
  void  UnfoldClusters(DigitsList * DL, 
		       AliPHOSRecPoint::RecPointsList * emcIn, 
		       AliPHOSEmcRecPoint * iniEmc, 
		       Int_t Nmax, 
		       int * maxAt, 
		       Float_t * maxAtEnergy ) ; //Unfolds overlaping clusters using TMinuit package
  virtual void UnsetUnfoldFlag() { fUnfoldFlag = kFALSE ; } 

  AliPHOSTrackSegmentMakerv1 & operator = (const AliPHOSTrackSegmentMakerv1 & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }

private:

  Float_t fDelta ;     // parameter used for sorting
  TMinuit * fMinuit ;  // Minuit object needed by cluster unfolding
  Float_t fR0 ;        // Maximum distance between a EMC RecPoint and a PPSD RecPoint   
  Bool_t fUnfoldFlag ; // Directive to unfold or not the clusters in case of multiple maxima

  ClassDef( AliPHOSTrackSegmentMakerv1,1)  // Implementation version 1 of algorithm class to make PHOS track segments 

};

#endif // AliPHOSTRACKSEGMENTMAKERV1_H
