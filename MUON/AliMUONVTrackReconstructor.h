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
           /// Return chamber thickness in number of radiation lengths
  Double_t GetChamberThicknessInX0(void) const {return fChamberThicknessInX0;}

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
  
  // needed in MakeSegmentsPerStation !!!
  Double_t GetImpactParamFromBendingMomentum(Double_t BendingMomentum) const;
  Double_t GetBendingMomentumFromImpactParam(Double_t ImpactParam) const;
  
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
  static const Double_t fgkDefaultMaxSigma2Distance; ///< default square of max. distance for window size 
  // Simple magnetic field:
  // Value taken from macro MUONtracking.C: 0.7 T, hence 7 kG
  // Length and Position from reco_muon.F, with opposite sign:
  // Length = ZMAGEND-ZCOIL
  // Position = (ZMAGEND+ZCOIL)/2
  // to be ajusted differently from real magnetic field ????
  static const Double_t fgkDefaultSimpleBValue; ///< default value of magnetic field (dipole)
  static const Double_t fgkDefaultSimpleBLength; ///< default length of magnetic field (dipole)
  static const Double_t fgkDefaultSimpleBPosition; ///< default position of magnetic field (dipole)
  
  // Parameters for track reconstruction
  Double_t fMinBendingMomentum; ///< minimum value (GeV/c) of momentum in bending plane
  Double_t fMaxBendingMomentum; ///< maximum value (GeV/c) of momentum in bending plane
  Double_t fBendingResolution; ///< chamber resolution (cm) in bending plane
  Double_t fNonBendingResolution; ///< chamber resolution (cm) in non bending plane
  Double_t fMaxSigma2Distance; ///< maximum square distance in units of the variance (maximum chi2)
  Double_t fChamberThicknessInX0; ///< chamber thickness in number of radiation lengths
  Double_t fSimpleBValue; ///< simple magnetic field: value (kG)
  Double_t fSimpleBLength; ///< simple magnetic field: length (cm)
  Double_t fSimpleBPosition; ///< simple magnetic field: Z central position (cm)
  
  Double_t* fSegmentMaxDistBending; ///< maximum distance (cm) for segments in bending plane
  Double_t* fSegmentMaxDistNonBending; ///< maximum distance (cm) for segments in non bending plane
  
  // Hits for reconstruction (should be in AliMUON ????)
  TClonesArray* fHitsForRecPtr; ///< pointer to the array of hits for reconstruction
  Int_t fNHitsForRec; ///< number of hits for reconstruction
  // Information per chamber (should be in AliMUONChamber ????)
  Int_t* fNHitsForRecPerChamber; ///< number of HitsForRec
  Int_t* fIndexOfFirstHitForRecPerChamber; ///< index (0...) of first HitForRec

  // Segments inside a station
  TClonesArray** fSegmentsPtr; ///< array of pointers to the segments for each station
  Int_t* fNSegments; ///< number of segments for each station

  // Reconstructed tracks
  TClonesArray *fRecTracksPtr; ///< pointer to array of reconstructed tracks
  Int_t fNRecTracks; ///< number of reconstructed tracks

  // data container
  AliMUONData* fMUONData; ///< Data container for MUON subsystem 

  // Functions
  AliMUONVTrackReconstructor (const AliMUONVTrackReconstructor& rhs); ///< copy constructor
  AliMUONVTrackReconstructor& operator=(const AliMUONVTrackReconstructor& rhs); ///< assignment operator
  
  void SortHitsForRecWithIncreasingChamber();
  void MakeSegmentsPerStation(Int_t Station);

  virtual void AddHitsForRecFromRawClusters() = 0;
  virtual void MakeSegments(void) = 0;
  virtual void MakeTracks(void) = 0;
  virtual void MakeTrackCandidates(void) = 0;
  virtual void FollowTracks(void) = 0;
  virtual void RemoveDoubleTracks(void) = 0;

 private:
  
  AliMUONTriggerTrack* fTriggerTrack; ///< Trigger track structure

  TClonesArray* fTriggerCircuit;      //!< trigger circuit array
  
  // Functions
  void SetReconstructionParametersToDefaults(void);
  
  void ResetTracks(void);
  void ResetSegments(void);
  void ResetHitsForRec(void);
  
  void ValidateTracksWithTrigger(void);
  
  Bool_t MakeTriggerTracks(void);


  ClassDef(AliMUONVTrackReconstructor, 0) // MUON track reconstructor in ALICE
    };
	
#endif
