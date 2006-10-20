#ifndef ALIMUONSEGMENT_H
#define ALIMUONSEGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONSegment
/// \brief Segment for reconstruction in ALICE dimuon spectrometer
///
////////////////////////////////////////////////////////////
/// Segment for reconstruction in ALICE dimuon  spectrometer
////////////////////////////////////////////////////////////

#include <TObject.h>

class AliMUONHitForRec;
class AliMUONTrackParam;

class AliMUONSegment : public TObject 
{
 public:
  AliMUONSegment(); // default constructor
  virtual ~AliMUONSegment(){} // Destructor
  AliMUONSegment(AliMUONHitForRec* Hit1, AliMUONHitForRec* Hit2); // Constructor from two HitForRec's

  // Inline functions for Get and Set
  AliMUONHitForRec* GetHitForRec1(void) const {return fHitForRecPtr1;}
  AliMUONHitForRec* GetHitForRec2(void) const {return fHitForRecPtr2;}
  Double_t GetBendingCoor(void) const {return fBendingCoor;}
  void SetBendingCoor(Double_t BendingCoor) {fBendingCoor = BendingCoor;}
  Double_t GetBendingSlope(void) const {return fBendingSlope;}
  void SetBendingSlope(Double_t BendingSlope) {fBendingSlope = BendingSlope;}
  Double_t GetNonBendingCoor(void) const {return fNonBendingCoor;}
  void SetNonBendingCoor(Double_t NonBendingCoor) {fNonBendingCoor = NonBendingCoor;}
  Double_t GetNonBendingSlope(void) const {return fNonBendingSlope;}
  void SetNonBendingSlope(Double_t NonBendingSlope) {fNonBendingSlope = NonBendingSlope;}
  Double_t GetBendingCoorReso2(void) const {return fBendingCoorReso2;}
  void SetBendingCoorReso2(Double_t BendingCoorReso2) {fBendingCoorReso2 = BendingCoorReso2;}
  Double_t GetNonBendingCoorReso2(void) const {return fNonBendingCoorReso2;}
  void SetNonBendingCoorReso2(Double_t NonBendingCoorReso2) {fNonBendingCoorReso2 = NonBendingCoorReso2;}
  Double_t GetZ(void) const {return fZ;}
  
  Double_t GetBendingImpact(void) const {return fBendingImpact;}
  Bool_t GetInTrack(void) const {return fInTrack;}
  void SetInTrack(Bool_t InTrack) {fInTrack = InTrack;}

  AliMUONSegment* CreateSegmentFromLinearExtrapToStation (Double_t z, Double_t MCSfactor) const;
  Double_t NormalizedChi2WithSegment(AliMUONSegment* Segment, Double_t Sigma2Cut) const;
  AliMUONHitForRec* CreateHitForRecFromLinearExtrapToChamber (Double_t z, Double_t MCSfactor) const;
  void UpdateFromStationTrackParam(AliMUONTrackParam *TrackParam, Double_t MCSfactor, Double_t Dz1, Double_t Dz2, Double_t Dz3, Int_t Station, Double_t InverseMomentum);

  // What is necessary for sorting TClonesArray's; sufficient too ????
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(const TObject* Segment) const; // "Compare" function for sorting

  void Print(Option_t* opt="") const;
  
 private:
  AliMUONHitForRec* fHitForRecPtr1; ///< pointer to HitForRec in first chamber
  AliMUONHitForRec* fHitForRecPtr2; ///< pointer to HitForRec in second chamber
  // Bending plane:
  Double_t fBendingCoor; ///< Coordinate in bending plane
  Double_t fBendingSlope; ///< Slope in bending plane
  // Covariance in bending plane:
  Double_t fBendingCoorReso2; ///< Covariance(coordinate C1 in first chamber)
  Double_t fBendingSlopeReso2; ///< Covariance(slope)
  Double_t fBendingCoorSlopeReso2; ///< Covariance(C1,slope)
  Double_t fBendingImpact; ///< Impact parameter in bending plane
  // Non Bending plane:
  Double_t fNonBendingCoor; ///< Coordinate in non bending plane
  Double_t fNonBendingSlope; ///< Slope in non bending plane
  // Covariance in non bending plane:
  Double_t fNonBendingCoorReso2; ///< Covariance(coordinate C1 in first chamber)
  Double_t fNonBendingSlopeReso2; ///< Covariance(slope)
  Double_t fNonBendingCoorSlopeReso2; ///< Covariance(C1,slope)
  Double_t fNonBendingImpact; ///< Impact parameter in non bending plane
  Double_t fZ;                ///< Z of the segment
  Bool_t fInTrack; ///< TRUE if segment belongs to one track
  
  AliMUONSegment (const AliMUONSegment& AliMUONSegment); // copy constructor
  AliMUONSegment& operator=(const AliMUONSegment& AliMUONSegment); // assignment operator

  ClassDef(AliMUONSegment, 1) // Segment for reconstruction in ALICE dimuon spectrometer
};
	
#endif
