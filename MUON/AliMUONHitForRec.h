#ifndef ALIMUONHITFORREC_H
#define ALIMUONHITFORREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

#include <TROOT.h>

class AliMUONHit;
class AliMUONRawCluster;
class AliMUONTrackHit;
class AliMUONTrackParam;

class AliMUONHitForRec : public TObject {
 public:
  AliMUONHitForRec(){
    // Constructor
    ;} // Constructor
  virtual ~AliMUONHitForRec(){
    // Destructor
    ;} // Destructor
  AliMUONHitForRec (const AliMUONHitForRec& AliMUONHitForRec); // copy constructor
  AliMUONHitForRec& operator=(const AliMUONHitForRec& AliMUONHitForRec); // assignment operator
  AliMUONHitForRec(AliMUONHit* mHit); // Constructor from GEANT hit
  AliMUONHitForRec(AliMUONRawCluster* RawCluster); // Constructor from raw cluster

  Double_t GetBendingCoor(void);
  void SetBendingCoor(Double_t BendingCoor);
  Double_t GetNonBendingCoor(void);
  void SetNonBendingCoor(Double_t NonBendingCoor);
  Double_t GetZ(void);
  void SetZ(Double_t Z);
  Double_t GetBendingReso2(void);
  void SetBendingReso2(Double_t BendingReso2);
  Double_t GetNonBendingReso2(void);
  void SetNonBendingReso2(Double_t NonBendingReso2);
  Int_t GetChamberNumber(void);
  void SetChamberNumber(Int_t ChamberNumber);
  Int_t GetHitNumber(void);
  void SetHitNumber(Int_t HitNumber);
  Int_t GetTHTrack(void);
  void SetTHTrack(Int_t THTrack);
  Int_t GetGeantSignal(void);
  void SetGeantSignal(Int_t GeantSignal);
  Int_t GetIndexOfFirstSegment(void);
  void SetIndexOfFirstSegment(Int_t IndexOfFirstSegment);
  Int_t GetNSegments(void);
  void SetNSegments(Int_t NSegments);
  AliMUONTrackHit* GetFirstTrackHitPtr(void);
  void SetFirstTrackHitPtr(AliMUONTrackHit* FirstTrackHitPtr);
  AliMUONTrackHit* GetLastTrackHitPtr(void);
  void SetLastTrackHitPtr(AliMUONTrackHit* LastTrackHitPtr);
  Int_t GetNTrackHits(void);
  void SetNTrackHits(Int_t NTrackHits);

  Double_t NormalizedChi2WithHitForRec(AliMUONHitForRec* Hit, Double_t Sigma2Cut);
/*   void UpdateFromChamberTrackParam(AliMUONTrackParam *TrackParam, Double_t MCSfactor); */

  // What is necessary for sorting TClonesArray's; sufficient too ????
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(TObject* HitForRec); // "Compare" function for sorting
 protected:
 private:
  Double_t fBendingCoor; // coordinate (cm) in bending plane
  Double_t fNonBendingCoor; // coordinate (cm) in non bending plane
  Double_t fZ; // Z coordinate (cm)
  Double_t fBendingReso2; // resolution**2 (cm**2) on coordinate in bending plane
  Double_t fNonBendingReso2; // resolution**2 (cm**2) on coordinate in non bending plane

  // links back to original hit for various checks
  // ideal would be real link to "hit" or "reconstructed hit"
  // if everything would be in memory ????
  Int_t fChamberNumber; // chamber number (0...)
  Int_t fHitNumber; // hit number (0...): RawCluster in "chamber" event of TR or GEANT hit in "track" event of TH
  Int_t fTHTrack; // track number (0...) in TH
  Int_t fGeantSignal; // Geant signal (1) or background (0)

  // links forward to the segment(s) if HitForRec in first chamber of a station
  Int_t fIndexOfFirstSegment; // index of first Segment
  Int_t fNSegments; // number of Segments

  // links forward to reconstructed track hits
  AliMUONTrackHit *fFirstTrackHitPtr ; // pointer to first TrackHit made with HitForRec
  AliMUONTrackHit *fLastTrackHitPtr ; // pointer to last TrackHit made with HitForRec
  Int_t fNTrackHits; // number of TrackHit's made with HitForRec
  
  ClassDef(AliMUONHitForRec, 1) // Class definition in ROOT context
    };
	
#endif
