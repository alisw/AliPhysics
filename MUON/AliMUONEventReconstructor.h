#ifndef ALIMUONEVENTRECONSTRUCTOR_H
#define ALIMUONEVENTRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

#include <TROOT.h>

class AliMUONHit;
class AliMUONHitForRec;
class AliMUONSegment;
class TClonesArray;
class TFile;
class TTree;

// Constants which should be elsewhere ????
const Int_t kMaxMuonTrackingChambers = 10;
const Int_t kMaxMuonTrackingStations = kMaxMuonTrackingChambers / 2;

class AliMUONEventReconstructor : public TObject {

 public:
  AliMUONEventReconstructor(void); // Constructor
  virtual ~AliMUONEventReconstructor(void); // Destructor
  AliMUONEventReconstructor (const AliMUONEventReconstructor& Reconstructor); // copy constructor
  AliMUONEventReconstructor& operator=(const AliMUONEventReconstructor& Reconstructor); // assignment operator

  // Parameters for event reconstruction: public methods
  // Get and Set, Set to defaults
  Double_t GetMinBendingMomentum(void) {return fMinBendingMomentum;}
  void SetMinBendingMomentum(Double_t MinBendingMomentum) {fMinBendingMomentum = MinBendingMomentum;}
  Double_t GetMaxSigma2Distance(void) {return fMaxSigma2Distance;}
  void SetMaxSigma2Distance(Double_t MaxSigma2Distance) {fMaxSigma2Distance = MaxSigma2Distance;}
  Double_t GetBendingResolution(void) {return fBendingResolution;}
  void SetBendingResolution(Double_t BendingResolution) {fBendingResolution = BendingResolution;}
  Double_t GetNonBendingResolution(void) {return fNonBendingResolution;}
  void SetNonBendingResolution(Double_t NonBendingResolution) {fNonBendingResolution = NonBendingResolution;}
  Double_t GetChamberThicknessInX0(void) {return fChamberThicknessInX0;}
  void SetChamberThicknessInX0(Double_t ChamberThicknessInX0) {fChamberThicknessInX0 = ChamberThicknessInX0;}
  Double_t GetSimpleBValue(void) {return fSimpleBValue;}
  void SetSimpleBValue(Double_t SimpleBValue) {fSimpleBValue = SimpleBValue;}
  Double_t GetSimpleBLength(void) {return fSimpleBLength;}
  void SetSimpleBLength(Double_t SimpleBLength) {fSimpleBLength = SimpleBLength;}
  Double_t GetSimpleBPosition(void) {return fSimpleBPosition;}
  void SetSimpleBPosition(Double_t SimpleBPosition) {fSimpleBPosition = SimpleBPosition;}
  Int_t GetRecGeantHits(void) {return fRecGeantHits;}
  void SetRecGeantHits(Int_t RecGeantHits) {fRecGeantHits = RecGeantHits;}
  Double_t GetEfficiency(void) {return fEfficiency;}
  void SetEfficiency(Double_t Efficiency) {fEfficiency = Efficiency;}
  Int_t GetPrintLevel(void) {return fPrintLevel;}
  void SetPrintLevel(Int_t PrintLevel) {fPrintLevel = PrintLevel;}
  void SetReconstructionParametersToDefaults(void);

  // Parameters for GEANT background events
  TFile* GetBkgGeantFile(void) {return fBkgGeantFile;}
  void SetBkgGeantFile(Text_t *BkgGeantFileName); // set background file for GEANT hits
  void NextBkgGeantEvent(void); // next event in background file for GEANT hits

  // Hits for reconstruction
  Int_t GetNHitsForRec() {return fNHitsForRec;} // Number

  // Functions
  Double_t GetImpactParamFromBendingMomentum(Double_t BendingMomentum);
  Double_t GetBendingMomentumFromImpactParam(Double_t ImpactParam);
  void EventReconstruct(void);

 protected:

 private:

  // Parameters for event reconstruction
  Double_t fMinBendingMomentum; // minimum value (GeV/c) of momentum in bending plane
  Double_t fMaxSigma2Distance; // maximum square distance in units of the variance (maximum chi2)
  Double_t fRMin[kMaxMuonTrackingChambers]; // minimum radius (cm)
  Double_t fRMax[kMaxMuonTrackingChambers]; // maximum radius (cm)
  Double_t fSegmentMaxDistBending[kMaxMuonTrackingStations]; // maximum distance (cm) for segments in bending plane
  Double_t fSegmentMaxDistNonBending[kMaxMuonTrackingStations]; // maximum distance (cm) for segments in bending plane
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
  TClonesArray *fBkgGeantHits;   // pointer to list of hits in tree TH
  Int_t fBkgGeantEventNumber;   // event number
  
  // Hits for reconstruction (should be in AliMUON ????)
  TClonesArray *fHitsForRecPtr; // pointer to the array of hits for reconstruction
  Int_t fNHitsForRec; // number of hits for reconstruction
  // Information per chamber (should be in AliMUONChamber ????)
  Int_t fNHitsForRecPerChamber[kMaxMuonTrackingChambers]; // number of HitsForRec
  Int_t fIndexOfFirstHitForRecPerChamber[kMaxMuonTrackingChambers]; // index (0...) of first HitForRec

  // Segments inside a station
  TClonesArray *fSegmentsPtr[kMaxMuonTrackingStations]; // array of pointers to the segments for each station
  Int_t fNSegments[kMaxMuonTrackingStations]; // number of segments for each station

  // Tracks
  TClonesArray *fRecTracksPtr; // pointer to array of reconstructed tracks
  Int_t fNRecTracks; // number of reconstructed tracks

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
  void ResetTracks(void);
  Int_t MakeTrackCandidatesWithTwoSegments(AliMUONSegment *BegSegment);
  Int_t MakeTrackCandidatesWithOneSegmentAndOnePoint(AliMUONSegment *BegSegment);
  void MakeTrackCandidates(void);
  void FollowTracks(void);

  ClassDef(AliMUONEventReconstructor, 1) // Class definition in ROOT context
    };
	
#endif
