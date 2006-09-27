#ifndef ALIMUONTRACKRECONSTRUCTOR_H
#define ALIMUONTRACKRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONTrackReconstructor
/// \brief Class for the MUON track reconstruction
///
/////////////////////////////////////
/// MUON track reconstructor in ALICE
/////////////////////////////////////

#include <TObject.h>
#include "AliMUONConstants.h"

class AliMUONHitForRec;
class AliMUONSegment;
class TClonesArray;
class TFile;
class TTree;
class AliMUONData;
class AliRunLoader;
class AliLoader;
class AliTrackReference;
class AliMUONTriggerTrack;

class AliMUONTrackReconstructor : public TObject {

 public:
  AliMUONTrackReconstructor(AliLoader* loader, AliMUONData* data); // default Constructor
  virtual ~AliMUONTrackReconstructor(void); // Destructor

  // Parameters for track reconstruction: public methods
  // Get and Set, Set to defaults
  Double_t GetMinBendingMomentum(void) const {return fMinBendingMomentum;}
  void SetMinBendingMomentum(Double_t MinBendingMomentum) {fMinBendingMomentum = MinBendingMomentum;}
  Double_t GetMaxBendingMomentum(void) const {return fMaxBendingMomentum;}
  void SetMaxBendingMomentum(Double_t MaxBendingMomentum) {fMaxBendingMomentum = MaxBendingMomentum;}
  Double_t GetMaxChi2(void) const {return fMaxChi2;}
  void SetMaxChi2(Double_t MaxChi2) {fMaxChi2 = MaxChi2;}
  Double_t GetMaxSigma2Distance(void) const {return fMaxSigma2Distance;}
  void SetMaxSigma2Distance(Double_t MaxSigma2Distance) {fMaxSigma2Distance = MaxSigma2Distance;}
  Double_t GetBendingResolution(void) const {return fBendingResolution;}
  void SetBendingResolution(Double_t BendingResolution) {fBendingResolution = BendingResolution;}
  Double_t GetNonBendingResolution(void) const {return fNonBendingResolution;}
  void SetNonBendingResolution(Double_t NonBendingResolution) {fNonBendingResolution = NonBendingResolution;}
  Double_t GetChamberThicknessInX0(void) const {return fChamberThicknessInX0;}
  void SetChamberThicknessInX0(Double_t ChamberThicknessInX0) {fChamberThicknessInX0 = ChamberThicknessInX0;}
  Double_t GetSimpleBValue(void) const {return fSimpleBValue;}
  void SetSimpleBValue(Double_t SimpleBValue) {fSimpleBValue = SimpleBValue;}
  Double_t GetSimpleBLength(void) const {return fSimpleBLength;}
  void SetSimpleBLength(Double_t SimpleBLength) {fSimpleBLength = SimpleBLength;}
  Double_t GetSimpleBPosition(void) const {return fSimpleBPosition;}
  void SetSimpleBPosition(Double_t SimpleBPosition) {fSimpleBPosition = SimpleBPosition;}
  Int_t GetRecTrackRefHits(void) const {return fRecTrackRefHits;}
  void SetRecTrackRefHits(Int_t RecTrackRefHits) {fRecTrackRefHits = RecTrackRefHits;}
  Double_t GetEfficiency(void) const {return fEfficiency;}
  void SetEfficiency(Double_t Efficiency) {fEfficiency = Efficiency;}
  void SetReconstructionParametersToDefaults(void);

  // Parameters for Track Ref. background events
  TFile* GetBkgTrackRefFile(void) const {return fBkgTrackRefFile;}
  void SetBkgTrackRefFile(Text_t *BkgTrackRefFileName); // set background file for track ref. hits
  void NextBkgTrackRefEvent(void); // next event in background file for track ref. hits

  // Hits for reconstruction
  Int_t GetNHitsForRec(void) const {return fNHitsForRec;} // Number

  // Reconstructed tracks
  Int_t GetNRecTracks() const {return fNRecTracks;} // Number
  void SetNRecTracks(Int_t NRecTracks) {fNRecTracks = NRecTracks;}
  TClonesArray* GetRecTracksPtr(void) const {return fRecTracksPtr;} // Array
 
  // Hits on reconstructed tracks
  Int_t GetNRecTrackHits() const {return fNRecTrackHits;} // Number
  void SetNRecTrackHits(Int_t NRecTrackHits) {fNRecTrackHits = NRecTrackHits;}
  TClonesArray* GetRecTrackHitsPtr(void) const {return fRecTrackHitsPtr;} // Array

  // Functions
  Double_t GetImpactParamFromBendingMomentum(Double_t BendingMomentum) const;
  Double_t GetBendingMomentumFromImpactParam(Double_t ImpactParam) const;
  void EventReconstruct(void);
  void EventReconstructTrigger(void);
  void EventDump(void);  // dump reconstructed event
  void EventDumpTrigger(void);  // dump reconstructed trigger event
  //PH  void FillEvent();      // fill and write tree of reconstructed events
  void SetTrackMethod(Int_t iTrackMethod); //AZ
  Int_t GetTrackMethod(void) const {return fTrackMethod;} 
  void FillMUONTrack(void); // set track parameters at hits for Kalman track
  //Int_t fMuons; // AZ - number of muons within acceptance - just for tests

  AliMUONData*  GetMUONData() {return fMUONData;}

 private:

  // Constants which should be elsewhere ????
  static const Int_t fgkMaxMuonTrackingChambers = 10; ///< Max number of Muon tracking chambers
  static const Int_t fgkMaxMuonTrackingStations = 5; ///< Max number of Muon tracking stations

  // Defaults parameters for reconstruction
  static const Double_t fgkDefaultMinBendingMomentum; ///< default min. bending momentum for reconstruction
  static const Double_t fgkDefaultMaxBendingMomentum; ///< default max. bending momentum for reconstruction
  static const Double_t fgkDefaultMaxChi2; ///< default max. track chi2 for reconstruction
  static const Double_t fgkDefaultMaxSigma2Distance; ///< default square of max. distance for window size 
  static const Double_t fgkDefaultBendingResolution; ///< default bending coordinate resolution for reconstruction 
  static const Double_t fgkDefaultNonBendingResolution; ///< default non bending coordinate resolution for reconstruction
  static const Double_t fgkDefaultChamberThicknessInX0; ///< default chamber thickness in X0 for reconstruction
  // Simple magnetic field:
  // Value taken from macro MUONtracking.C: 0.7 T, hence 7 kG
  // Length and Position from reco_muon.F, with opposite sign:
  // Length = ZMAGEND-ZCOIL
  // Position = (ZMAGEND+ZCOIL)/2
  // to be ajusted differently from real magnetic field ????
  static const Double_t fgkDefaultSimpleBValue; ///< default value of magnetic field (dipole)
  static const Double_t fgkDefaultSimpleBLength; ///< default length of magnetic field (dipole)
  static const Double_t fgkDefaultSimpleBPosition; ///< default position of magnetic field (dipole)
  static const Int_t fgkDefaultRecTrackRefHits; ///< default flag for reconstrution track ref. hits or Clusters
  static const Double_t fgkDefaultEfficiency; ///< default chamber efficiency for track ref. hits recontruction


  Int_t fTrackMethod; ///< AZ - tracking method

  // Parameters for track reconstruction
  Double_t fMinBendingMomentum; ///< minimum value (GeV/c) of momentum in bending plane
  // Parameters for track reconstruction
  Double_t fMaxBendingMomentum; ///< maximum value (GeV/c) of momentum in bending plane
  Double_t fMaxChi2; ///< maximum Chi2 per degree of Freedom
  Double_t fMaxSigma2Distance; ///< maximum square distance in units of the variance (maximum chi2)
  Double_t fRMin[fgkMaxMuonTrackingChambers]; ///< minimum radius (cm)
  Double_t fRMax[fgkMaxMuonTrackingChambers]; ///< maximum radius (cm)
  Double_t fSegmentMaxDistBending[fgkMaxMuonTrackingStations]; ///< maximum distance (cm) for segments in bending plane
  Double_t fSegmentMaxDistNonBending[fgkMaxMuonTrackingStations]; ///< maximum distance (cm) for segments in non bending plane
  Double_t fBendingResolution; ///< chamber resolution (cm) in bending plane
  Double_t fNonBendingResolution; ///< chamber resolution (cm) in non bending plane
  Double_t fChamberThicknessInX0; ///< chamber thickness in number of radiation lengths
                                  // how to take it from simulation ????
  Double_t fSimpleBValue; ///< simple magnetic field: value (kG)
  Double_t fSimpleBLength; ///< simple magnetic field: length (cm)
  Double_t fSimpleBPosition; ///< simple magnetic field: Z central position (cm)
  Int_t fRecTrackRefHits; ///< reconstruction from raw clusters (0) or from track ref. hits (1)
  Double_t fEfficiency; ///< chamber efficiency (used for track ref. hits only)

  // Parameters for track ref. background events
  // should be in AliMUON class ????
  TFile *fBkgTrackRefFile; ///< pointer to file
  TTree *fBkgTrackRefTK; ///< pointer to tree TK
  TClonesArray *fBkgTrackRefParticles;   ///< pointer to list of particles in tree TK
  TTree *fBkgTrackRefTTR; ///< pointer to tree TTR
  Int_t fBkgTrackRefEventNumber; ///< event number
  
  // Hits for reconstruction (should be in AliMUON ????)
  TClonesArray *fHitsForRecPtr; ///< pointer to the array of hits for reconstruction
  Int_t fNHitsForRec; ///< number of hits for reconstruction
  // Information per chamber (should be in AliMUONChamber ????)
  Int_t fNHitsForRecPerChamber[fgkMaxMuonTrackingChambers]; ///< number of HitsForRec
  Int_t fIndexOfFirstHitForRecPerChamber[fgkMaxMuonTrackingChambers]; ///< index (0...) of first HitForRec

  // Segments inside a station
  TClonesArray *fSegmentsPtr[fgkMaxMuonTrackingStations]; ///< array of pointers to the segments for each station
  Int_t fNSegments[fgkMaxMuonTrackingStations]; ///< number of segments for each station

  // Reconstructed tracks
  TClonesArray *fRecTracksPtr; ///< pointer to array of reconstructed tracks
  Int_t fNRecTracks; ///< number of reconstructed tracks

  // Track hits on reconstructed tracks
  TClonesArray *fRecTrackHitsPtr; ///< pointer to array of hits on reconstructed tracks
  Int_t fNRecTrackHits; ///< number of hits on reconstructed tracks

  // data container
  AliMUONData* fMUONData; ///< Data container for MUON subsystem 

  // alice loader
  AliLoader* fLoader; ///< MUON loader to get data

  Int_t fMuons; ///< AZ - number of muons within acceptance - just for tests

  AliMUONTriggerTrack* fTriggerTrack; ///< Trigger track structure

  // Functions
  AliMUONTrackReconstructor (const AliMUONTrackReconstructor& rhs); // copy constructor
  AliMUONTrackReconstructor& operator=(const AliMUONTrackReconstructor& rhs); // assignment operator
  void ResetHitsForRec(void);
  void MakeEventToBeReconstructed(void);
  void AddHitsForRecFromTrackRef(TTree *TTR, Int_t Signal);
  AliMUONHitForRec* NewHitForRecFromTrackRef(AliTrackReference* Hit, Int_t TrackNumber, Int_t Signal);
  TClonesArray *CleanTrackRefs(TTree *treeTR);
/*   void AddHitsForRecFromCathodeCorrelations(TTree* TC); */
  void AddHitsForRecFromRawClusters(TTree* TR);
  void SortHitsForRecWithIncreasingChamber();
  void MakeSegments(void);
  void ResetSegments(void);
  void MakeSegmentsPerStation(Int_t Station);
  void MakeTracks(void);
  Bool_t MakeTriggerTracks(void);
  void ResetTrackHits(void);
  void ResetTracks(void);
  Int_t MakeTrackCandidatesWithTwoSegments(AliMUONSegment *BegSegment);
  Int_t MakeTrackCandidatesWithOneSegmentAndOnePoint(AliMUONSegment *BegSegment);
  void MakeTrackCandidates(void);
  void FollowTracks(void);
  void RemoveDoubleTracks(void);
  void UpdateTrackParamAtHit(void);
  void UpdateHitForRecAtHit(void);
  void ValidateTracksWithTrigger(void);


  //AZ - for Kalman Filter
  void MakeTrackCandidatesK(void);
  void FollowTracksK(void);
  void RemoveDoubleTracksK(void);
  void GoToVertex(void);
  Bool_t CheckCandidateK(Int_t icand, Int_t nSeeds) const;


  ClassDef(AliMUONTrackReconstructor, 0) // MUON track reconstructor in ALICE
    };
	
#endif
