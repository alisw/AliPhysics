#ifndef ALIMUONSEGMENT_H
#define ALIMUONSEGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

#include <TROOT.h>

class AliMUONHitForRec;
class AliMUONTrackParam;

class AliMUONSegment : public TObject {
 public:
  AliMUONSegment(){
    // Constructor
    ;} // Constructor
  virtual ~AliMUONSegment(){
    // Destructor
    ;} // Destructor
  AliMUONSegment (const AliMUONSegment& AliMUONSegment); // copy constructor
  AliMUONSegment& operator=(const AliMUONSegment& AliMUONSegment); // assignment operator
  AliMUONSegment(AliMUONHitForRec* Hit1, AliMUONHitForRec* Hit2); // Constructor from two HitForRec's

  AliMUONHitForRec* GetHitForRec1(void);
  AliMUONHitForRec* GetHitForRec2(void);
  Double_t GetBendingCoorReso2(void);
  void SetBendingCoorReso2(Double_t BendingCoorReso2);
  Double_t GetNonBendingCoorReso2(void);
  void SetNonBendingCoorReso2(Double_t NonBendingCoorReso2);
  Double_t GetBendingImpact(void);
  Bool_t GetInTrack(void);
  void SetInTrack(Bool_t InTrack);

  AliMUONSegment* CreateSegmentFromLinearExtrapToStation (Int_t Station, Double_t MCSfactor);
  Double_t NormalizedChi2WithSegment(AliMUONSegment* Segment, Double_t Sigma2Cut);
  AliMUONHitForRec* CreateHitForRecFromLinearExtrapToChamber (Int_t Chamber, Double_t MCSfactor);
  void UpdateFromStationTrackParam(AliMUONTrackParam *TrackParam, Double_t MCSfactor, Double_t Dz1, Double_t Dz2);

  // What is necessary for sorting TClonesArray's; sufficient too ????
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(TObject* Segment); // "Compare" function for sorting
 protected:
 private:
  AliMUONHitForRec* fHitForRecPtr1; // pointer to HitForRec in first chamber
  AliMUONHitForRec* fHitForRecPtr2; // pointer to HitForRec in second chamber
  // Bending plane:
  Double_t fBendingCoor; // Coordinate in bending plane
  Double_t fBendingSlope; // Slope in bending plane
  // Covariance in bending plane:
  Double_t fBendingCoorReso2; // Covariance(coordinate C1 in first chamber)
  Double_t fBendingSlopeReso2; // Covariance(slope)
  Double_t fBendingCoorSlopeReso2; // Covariance(C1,slope)
  Double_t fBendingImpact; // Impact parameter in bending plane
  // Non Bending plane:
  Double_t fNonBendingCoor; // Coordinate in non bending plane
  Double_t fNonBendingSlope; // Slope in non bending plane
  // Covariance in non bending plane:
  Double_t fNonBendingCoorReso2; // Covariance(coordinate C1 in first chamber)
  Double_t fNonBendingSlopeReso2; // Covariance(slope)
  Double_t fNonBendingCoorSlopeReso2; // Covariance(C1,slope)
  Double_t fNonBendingImpact; // Impact parameter in non bending plane
  Bool_t fInTrack; // TRUE if segment belongs to one track
  
  ClassDef(AliMUONSegment, 1) // Class definition in ROOT context
    };
	
#endif
