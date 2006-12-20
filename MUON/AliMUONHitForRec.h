#ifndef ALIMUONHITFORREC_H
#define ALIMUONHITFORREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONHitForRec
/// \brief Hit for reconstruction in ALICE dimuon spectrometer

#include <TObject.h>

class AliTrackReference;
class AliMUONRawCluster;
class AliMUONTrackHit;
class AliMUONTrackParam;

class AliMUONHitForRec : public TObject {
 public:
  AliMUONHitForRec(); // Constructor
  virtual ~AliMUONHitForRec(){} // Destructor
  AliMUONHitForRec (const AliMUONHitForRec& AliMUONHitForRec); // copy constructor
  AliMUONHitForRec& operator=(const AliMUONHitForRec& AliMUONHitForRec); // assignment operator
  AliMUONHitForRec(AliTrackReference* mHit); // Constructor from track ref. hit
  AliMUONHitForRec(AliMUONRawCluster* theRawCluster); // Constructor from raw cluster

  // Inline functions for Get and Set
  Double_t GetBendingCoor(void) const { return fBendingCoor;}
  void SetBendingCoor(Double_t BendingCoor) { fBendingCoor = BendingCoor;}
  Double_t GetNonBendingCoor(void) const { return fNonBendingCoor;}
  void SetNonBendingCoor(Double_t NonBendingCoor) { fNonBendingCoor = NonBendingCoor;}
  Double_t GetZ(void) const { return fZ;}
  void SetZ(Double_t Z) { fZ = Z;}
  Double_t GetBendingReso2(void) const { return fBendingReso2;}
  void SetBendingReso2(Double_t BendingReso2) { fBendingReso2 = BendingReso2;}
  Double_t GetNonBendingReso2(void) const { return fNonBendingReso2;}
  void SetNonBendingReso2(Double_t NonBendingReso2) { fNonBendingReso2 = NonBendingReso2;}
  Int_t GetChamberNumber(void) const { return fChamberNumber;}
  void SetChamberNumber(Int_t ChamberNumber) { fChamberNumber = ChamberNumber;}
  Int_t GetDetElemId(void) const {return fDetElemId;}
  void SetDetElemId(Int_t id) { fDetElemId = id;}
  Int_t GetHitNumber(void) const { return fHitNumber;}
  void SetHitNumber(Int_t HitNumber) { fHitNumber = HitNumber;}
  Int_t GetTTRTrack(void) const { return fTTRTrack;}
  void SetTTRTrack(Int_t TTRTrack) { fTTRTrack = TTRTrack;}
  Int_t GetTrackRefSignal(void) const { return fTrackRefSignal;}
  void SetTrackRefSignal(Int_t TrackRefSignal) { fTrackRefSignal = TrackRefSignal;}
  Int_t GetNTrackHits(void) const { return fNTrackHits;}
  void SetNTrackHits(Int_t NTrackHits) { fNTrackHits = NTrackHits;}

  Double_t NormalizedChi2WithHitForRec(AliMUONHitForRec* Hit, Double_t Sigma2Cut) const;

  // What is necessary for sorting TClonesArray's; sufficient too ????
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
