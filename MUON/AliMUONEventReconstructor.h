#ifndef ALIMUONEVENTRECONSTRUCTOR_H
#define ALIMUONEVENTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

////////////////////////////////////
// MUON event reconstructor in ALICE
////////////////////////////////////

#include <TObject.h>

class AliMUONHit;
class AliMUONHitForRec;
class AliMUONSegment;
class TClonesArray;
class TFile;
class TTree;
class AliMUONRecoEvent;
class AliMUONData;
class AliRunLoader;
class AliLoader;

class AliMUONEventReconstructor : public TObject {

 public:
  AliMUONEventReconstructor(AliLoader* loader); // default Constructor
  virtual ~AliMUONEventReconstructor(void); // Destructor

  // Parameters for event reconstruction: public methods
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
  Int_t GetRecGeantHits(void) const {return fRecGeantHits;}
  void SetRecGeantHits(Int_t RecGeantHits) {fRecGeantHits = RecGeantHits;}
  Double_t GetEfficiency(void) const {return fEfficiency;}
  void SetEfficiency(Double_t Efficiency) {fEfficiency = Efficiency;}
  Int_t GetPrintLevel(void) const {return fPrintLevel;}
  void SetPrintLevel(Int_t PrintLevel) {fPrintLevel = PrintLevel;}
  void SetReconstructionParametersToDefaults(void);

  // Parameters for GEANT background events
  TFile* GetBkgGeantFile(void) const {return fBkgGeantFile;}
  void SetBkgGeantFile(Text_t *BkgGeantFileName); // set background file for GEANT hits
  void NextBkgGeantEvent(void); // next event in background file for GEANT hits

  // Hits for reconstruction
  Int_t GetNHitsForRec(void) const {return fNHitsForRec;} // Number

  // Reconstructed tracks
  Int_t GetNRecTracks() const {return fNRecTracks;} // Number
  void SetNRecTracks(Int_t NRecTracks) {fNRecTracks = NRecTracks;}
  TClonesArray* GetRecTracksPtr(void) const {return fRecTracksPtr;} // Array
 
 // Reconstructed trigger tracks
   Int_t GetNRecTriggerTracks() const {return fNRecTriggerTracks;} // Number
   void SetNRecTriggerTracks(Int_t NRecTriggerTracks) {fNRecTriggerTracks = NRecTriggerTracks;}
   TClonesArray* GetRecTriggerTracksPtr(void) const {return fRecTriggerTracksPtr;} // Array

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
  void FillEvent();      // fill and write tree of reconstructed events
  void SetTrackMethod(Int_t TrackMethod) {fTrackMethod = TrackMethod;} //AZ
  Int_t GetTrackMethod(void) const {return fTrackMethod;} //AZ
  //Int_t fMuons; // AZ - number of muons within acceptance - just for tests

  AliMUONData*  GetMUONData() {return fMUONData;}

 protected:
  AliMUONEventReconstructor (const AliMUONEventReconstructor& rhs); // copy constructor
  AliMUONEventReconstructor& operator=(const AliMUONEventReconstructor& rhs); // assignment operator

 private:

  // Constants which should be elsewhere ????
  static const Int_t fgkMaxMuonTrackingChambers = 10; // Max number of Muon tracking chambers
  static const Int_t fgkMaxMuonTrackingStations = 5; // Max number of Muon tracking stations

  // Defaults parameters for reconstruction
  static const Double_t fgkDefaultMinBendingMomentum; // default min. bending momentum for reconstruction
  static const Double_t fgkDefaultMaxBendingMomentum; // default max. bending momentum for reconstruction
  static const Double_t fgkDefaultMaxChi2; // default max. track chi2 for reconstruction
  static const Double_t fgkDefaultMaxSigma2Distance; // default square of max. distance for window size 
  static const Double_t fgkDefaultBendingResolution; // default bending coordinate resolution for reconstruction 
  static const Double_t fgkDefaultNonBendingResolution; // default non bending coordinate resolution for reconstruction
  static const Double_t fgkDefaultChamberThicknessInX0; // default chamber thickness in X0 for reconstruction
  // Simple magnetic field:
  // Value taken from macro MUONtracking.C: 0.7 T, hence 7 kG
  // Length and Position from reco_muon.F, with opposite sign:
  // Length = ZMAGEND-ZCOIL
  // Position = (ZMAGEND+ZCOIL)/2
  // to be ajusted differently from real magnetic field ????
  static const Double_t fgkDefaultSimpleBValue; // default value of magnetic field (dipole)
  static const Double_t fgkDefaultSimpleBLength; // default length of magnetic field (dipole)
  static const Double_t fgkDefaultSimpleBPosition; // default position of magnetic field (dipole)
  static const Int_t fgkDefaultRecGeantHits; // default flag for reconstrution GEANT hits or Clusters
  static const Double_t fgkDefaultEfficiency; // default chamber efficiency for GEANT hits recontruction

  static const Int_t fgkDefaultPrintLevel; // default print level

  Int_t fTrackMethod; // AZ - tracking method

  // Parameters for event reconstruction
  Double_t fMinBendingMomentum; // minimum value (GeV/c) of momentum in bending plane
  // Parameters for event reconstruction
  Double_t fMaxBendingMomentum; // maximum value (GeV/c) of momentum in bending plane
  Double_t fMaxChi2; // maximum Chi2 per degree of Freedom
  Double_t fMaxSigma2Distance; // maximum square distance in units of the variance (maximum chi2)
  Double_t fRMin[fgkMaxMuonTrackingChambers]; // minimum radius (cm)
  Double_t fRMax[fgkMaxMuonTrackingChambers]; // maximum radius (cm)
  Double_t fSegmentMaxDistBending[fgkMaxMuonTrackingStations]; // maximum distance (cm) for segments in bending plane
  Double_t fSegmentMaxDistNonBending[fgkMaxMuonTrackingStations]; // maximum distance (cm) for segments in non bending plane
  Double_t fBendingResolution; // chamber resolution (cm) in bending plane
  Double_t fNonBendingResolution; // chamber resolution (cm) in non bending plane
  Double_t fChamberThicknessInX0; // chamber thickness in number of radiation lengths
                                  // how to take it from simulation ????
  Double_t fSimpleBValue; // simple magnetic field: value (kG)
  Double_t fSimpleBLength; // simple magnetic field: length (cm)
  Double_t fSimpleBPosition; // simple magnetic field: Z central position (cm)
  Int_t fRecGeantHits; // reconstruction from raw clusters (0) or from GEANT hits (1)
  Double_t fEfficiency; // chamber efficiency (used for GEANT hits only)
  Int_t fPrintLevel; // print level

  // Parameters for GEANT background events
  // should be in AliMUON class ????
  TFile *fBkgGeantFile; // pointer to file
  TTree *fBkgGeantTK; // pointer to tree TK
  TClonesArray *fBkgGeantParticles;   // pointer to list of particles in tree TK
  TTree *fBkgGeantTH; // pointer to tree TH
  TClonesArray *fBkgGeantHits; // pointer to list of hits in tree TH
  Int_t fBkgGeantEventNumber; // event number
  
  // Hits for reconstruction (should be in AliMUON ????)
  TClonesArray *fHitsForRecPtr; // pointer to the array of hits for reconstruction
  Int_t fNHitsForRec; // number of hits for reconstruction
  // Information per chamber (should be in AliMUONChamber ????)
  Int_t fNHitsForRecPerChamber[fgkMaxMuonTrackingChambers]; // number of HitsForRec
  Int_t fIndexOfFirstHitForRecPerChamber[fgkMaxMuonTrackingChambers]; // index (0...) of first HitForRec

  // Segments inside a station
  TClonesArray *fSegmentsPtr[fgkMaxMuonTrackingStations]; // array of pointers to the segments for each station
  Int_t fNSegments[fgkMaxMuonTrackingStations]; // number of segments for each station

  // Reconstructed tracks
  TClonesArray *fRecTracksPtr; // pointer to array of reconstructed tracks
  Int_t fNRecTracks; // number of reconstructed tracks

  // Reconstructed trigger tracks
  TClonesArray *fRecTriggerTracksPtr; // pointer to array of reconstructed trigger tracks
  Int_t fNRecTriggerTracks; // number of reconstructed trigger tracks

  // Track hits on reconstructed tracks
  TClonesArray *fRecTrackHitsPtr; // pointer to array of hits on reconstructed tracks
  Int_t fNRecTrackHits; // number of hits on reconstructed tracks

  // Objects needed for tree writing
  AliMUONRecoEvent *fRecoEvent; // the reconstructed event
  TTree            *fEventTree; // tree of reconstructed events
  TFile            *fTreeFile;  // file where the tree is outputed

  // data container
  AliMUONData* fMUONData; // Data container for MUON subsystem 

  // alice loader
  AliLoader* fLoader; // MUON loader to get data

  Int_t fMuons; // AZ - number of muons within acceptance - just for tests

  // Functions
  void ResetHitsForRec(void);
  void MakeEventToBeReconstructed(void);
  void AddHitsForRecFromGEANT(TTree *TH);
  void AddHitsForRecFromBkgGEANT(TTree *TH, TClonesArray *Hits);
  AliMUONHitForRec* NewHitForRecFromGEANT(AliMUONHit* Hit, Int_t TrackNumber, Int_t HitNumber, Int_t Signal);
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
  void ResetTriggerTracks(void);
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

  ClassDef(AliMUONEventReconstructor, 0) // MUON event reconstructor in ALICE
    };
	
#endif
