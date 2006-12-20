#ifndef ALIMUONVTRACKRECONSTRUCTOR_H
#define ALIMUONVTRACKRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONVTrackReconstructor
/// \brief Virtual class for the MUON track reconstruction
///
/////////////////////////////////////////////
/// Virtual MUON track reconstructor in ALICE
/////////////////////////////////////////////

#include <TObject.h>

class TClonesArray;
class AliMUONData;
class AliMUONTriggerTrack;

class AliMUONVTrackReconstructor : public TObject {

 public:
  AliMUONVTrackReconstructor(AliMUONData* data); // default Constructor
  virtual ~AliMUONVTrackReconstructor(); // Destructor

  // Parameters for track reconstruction: public methods
  // Get and Set, Set to defaults
           /// Return minimum value (GeV/c) of momentum in bending plane
  Double_t GetMinBendingMomentum(void) const {return fMinBendingMomentum;}
           /// Return chamber resolution (cm) in bending plane
  Double_t GetBendingResolution(void) const {return fBendingResolution;}
           /// Return chamber resolution (cm) in non-bending plane
  Double_t GetNonBendingResolution(void) const {return fNonBendingResolution;}

  // Reconstructed tracks
           /// Return number of reconstructed tracks
  Int_t GetNRecTracks() const {return fNRecTracks;} // Number
           /// Set number of reconstructed tracks
  void SetNRecTracks(Int_t NRecTracks) {fNRecTracks = NRecTracks;}
           /// Return array of reconstructed tracks
  TClonesArray* GetRecTracksPtr(void) const {return fRecTracksPtr;} // Array
 
  // Functions
  void EventReconstruct(void);
  void EventReconstructTrigger(void);
  virtual void EventDump(void) = 0;  // dump reconstructed event
  void EventDumpTrigger(void);  // dump reconstructed trigger event
  
          /// Return MUON data
  AliMUONData*  GetMUONData() {return fMUONData;}

          /// Set trigger circuit
  void SetTriggerCircuit(TClonesArray* circuit) {fTriggerCircuit = circuit;}


 protected:

  // Defaults parameters for reconstruction
  static const Double_t fgkDefaultMinBendingMomentum; ///< default min. bending momentum for reconstruction
  static const Double_t fgkDefaultMaxBendingMomentum; ///< default max. bending momentum for reconstruction
  static const Double_t fgkDefaultBendingResolution; ///< default bending coordinate resolution for reconstruction 
  static const Double_t fgkDefaultNonBendingResolution; ///< default non bending coordinate resolution for reconstruction
  static const Double_t fgkDefaultBendingVertexDispersion; ///< default vertex dispersion in bending plane for reconstruction
  static const Double_t fgkDefaultNonBendingVertexDispersion; ///< default vertex dispersion in non bending plane for reconstruction
  static const Double_t fgkDefaultMaxNormChi2MatchTrigger; ///< default maximum normalized chi2 of tracking/trigger track matching
  
  // Parameters for track reconstruction
  Double_t fMinBendingMomentum; ///< minimum value (GeV/c) of momentum in bending plane
  Double_t fMaxBendingMomentum; ///< maximum value (GeV/c) of momentum in bending plane
  Double_t fBendingResolution; ///< chamber resolution (cm) in bending plane
  Double_t fNonBendingResolution; ///< chamber resolution (cm) in non bending plane
  Double_t fBendingVertexDispersion; ///< vextex dispersion (cm) in bending plane
  Double_t fNonBendingVertexDispersion; ///< vextex dispersion (cm) in non bending plane
  Double_t fMaxNormChi2MatchTrigger; ///< maximum normalized chi2 of tracking/trigger track matching
  
  Double_t* fSegmentMaxDistBending; ///< maximum distance (cm) for segments in bending plane
  Double_t* fSegmentMaxDistNonBending; ///< maximum distance (cm) for segments in non bending plane
  
  TClonesArray* fHitsForRecPtr; ///< pointer to the array of hits for reconstruction
  Int_t fNHitsForRec; ///< number of hits for reconstruction
  Int_t* fNHitsForRecPerChamber; ///< number of HitsForRec
  Int_t* fIndexOfFirstHitForRecPerChamber; ///< index (0...) of first HitForRec

  // Reconstructed tracks
  TClonesArray *fRecTracksPtr; ///< pointer to array of reconstructed tracks
  Int_t fNRecTracks; ///< number of reconstructed tracks

  // data container
  AliMUONData* fMUONData; ///< Data container for MUON subsystem 

  // Functions
  AliMUONVTrackReconstructor (const AliMUONVTrackReconstructor& rhs); ///< copy constructor
  AliMUONVTrackReconstructor& operator=(const AliMUONVTrackReconstructor& rhs); ///< assignment operator
  
  void SortHitsForRecWithIncreasingChamber();
  TClonesArray *MakeSegmentsInStation(Int_t station);

  virtual void AddHitsForRecFromRawClusters() = 0;
  virtual void MakeTracks(void) = 0;
  virtual void MakeTrackCandidates(void) = 0;
  virtual void FollowTracks(void) = 0;
  virtual void RemoveDoubleTracks(void) = 0;
  virtual void ExtrapTracksToVertex(void) = 0;
  virtual void FillMUONTrack(void) = 0;

 private:
  
  AliMUONTriggerTrack* fTriggerTrack; ///< Trigger track structure

  TClonesArray* fTriggerCircuit;      //!< trigger circuit array
  
  // Functions
  void SetReconstructionParametersToDefaults(void);
  
  void ResetTracks(void);
  void ResetHitsForRec(void);
  
  void ValidateTracksWithTrigger(void);
  
  Bool_t MakeTriggerTracks(void);


  ClassDef(AliMUONVTrackReconstructor, 0) // MUON track reconstructor in ALICE
};
	
#endif
