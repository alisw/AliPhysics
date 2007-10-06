#ifndef ALIMUONHITFORREC_H
#define ALIMUONHITFORREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONHitForRec
/// \brief Hit for reconstruction in ALICE dimuon spectrometer
///
/// \author J. Gosset

#include <TObject.h>

class AliMUONVCluster;
class AliMUONTrackHit;
class AliMUONTrackParam;

class AliMUONHitForRec : public TObject {
 public:
  AliMUONHitForRec(); // Constructor
  virtual ~AliMUONHitForRec(); // Destructor
  AliMUONHitForRec (const AliMUONHitForRec& AliMUONHitForRec); // copy constructor
  AliMUONHitForRec& operator=(const AliMUONHitForRec& AliMUONHitForRec); // assignment operator
  AliMUONHitForRec(AliMUONVCluster* theRawCluster); // Constructor from raw cluster

  // Inline functions for Get and Set
           /// Return coordinate (cm) in bending plane
  Double_t GetBendingCoor(void) const { return fBendingCoor;}
           /// Set coordinate (cm) in bending plane
  void SetBendingCoor(Double_t BendingCoor) { fBendingCoor = BendingCoor;}
           /// Return coordinate (cm) in non bending plane
  Double_t GetNonBendingCoor(void) const { return fNonBendingCoor;}
           /// Set coordinate (cm) in non bending plane
  void SetNonBendingCoor(Double_t NonBendingCoor) { fNonBendingCoor = NonBendingCoor;}
           /// Return Z coordinate (cm)
  Double_t GetZ(void) const { return fZ;}
           /// Set Z coordinate (cm)
  void SetZ(Double_t Z) { fZ = Z;}
           /// Return resolution**2 (cm**2) on coordinate in bending plane
  Double_t GetBendingReso2(void) const { return fBendingReso2;}
           /// Set resolution**2 (cm**2) on coordinate in bending plane
  void SetBendingReso2(Double_t BendingReso2) { fBendingReso2 = BendingReso2;}
           /// Return resolution**2 (cm**2) on coordinate in non bending plane
  Double_t GetNonBendingReso2(void) const { return fNonBendingReso2;}
           /// Set resolution**2 (cm**2) on coordinate in non bending plane
  void SetNonBendingReso2(Double_t NonBendingReso2) { fNonBendingReso2 = NonBendingReso2;}
           /// Return chamber number (0...)
  Int_t GetChamberNumber(void) const { return fChamberNumber;}
           /// Set chamber number (0...)
  void SetChamberNumber(Int_t ChamberNumber) { fChamberNumber = ChamberNumber;}
           /// Return detection element Id
  Int_t GetDetElemId(void) const {return fDetElemId;}
           /// Set detection element Id
  void SetDetElemId(Int_t id) { fDetElemId = id;}
           /// Return hit number (0...)
  Int_t GetHitNumber(void) const { return fHitNumber;}
           /// Set hit number (0...)
  void SetHitNumber(Int_t HitNumber) { fHitNumber = HitNumber;}
           /// Return track number (0...) in TTR
  Int_t GetTTRTrack(void) const { return fTTRTrack;}
           /// Set track number (0...) in TTR
  void SetTTRTrack(Int_t TTRTrack) { fTTRTrack = TTRTrack;}
           /// Return Track ref. signal (1) or background (0)
  Int_t GetTrackRefSignal(void) const { return fTrackRefSignal;}
           /// Set Track ref. signal (1) or background (0)
  void SetTrackRefSignal(Int_t TrackRefSignal) { fTrackRefSignal = TrackRefSignal;}
           /// Return number of TrackHit's made with HitForRec
  Int_t GetNTrackHits(void) const { return fNTrackHits;}
           /// Set number of TrackHit's made with HitForRec
  void SetNTrackHits(Int_t NTrackHits) { fNTrackHits = NTrackHits;}

  Double_t NormalizedChi2WithHitForRec(AliMUONHitForRec* Hit, Double_t Sigma2Cut) const;

  /// What is necessary for sorting TClonesArray's; sufficient too ????
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(const TObject* HitForRec) const; // "Compare" function for sorting

  virtual void Print(Option_t* opt="") const;
  
 private:
  Double_t fBendingCoor; ///< coordinate (cm) in bending plane
  Double_t fNonBendingCoor; ///< coordinate (cm) in non bending plane
  Double_t fZ; ///< Z coordinate (cm)
  Double_t fBendingReso2; ///< resolution**2 (cm**2) on coordinate in bending plane
  Double_t fNonBendingReso2; ///< resolution**2 (cm**2) on coordinate in non bending plane

  // links back to original hit for various checks
  // ideal would be real link to "hit" or "reconstructed hit"
  // if everything would be in memory ????
  Int_t fChamberNumber; ///< chamber number (0...)
  Int_t fDetElemId; ///< detection element Id   
  Int_t fHitNumber; ///< hit number (0...): RawCluster in "chamber" event of TR or track ref. hit in "track" event of TTR
  Int_t fTTRTrack; ///< track number (0...) in TTR
  Int_t fTrackRefSignal; ///< Track ref. signal (1) or background (0)

  Int_t fNTrackHits; //!<  number of TrackHit's made with HitForRec
  
  ClassDef(AliMUONHitForRec, 2) // Hit for reconstruction in ALICE dimuon spectrometer
    };
	
#endif
