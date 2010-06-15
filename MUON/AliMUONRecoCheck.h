#ifndef ALIMUONRECOCHECK_H
#define ALIMUONRECOCHECK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup evaluation
/// \class AliMUONRecoCheck
/// \brief Utility class to check reconstruction

#include <TObject.h>

class TClonesArray;
class TFile;
class TTree;
class AliESDEvent;
class AliMCEventHandler;
class AliMUONVTrackStore;
class AliMUONVTriggerTrackStore;
class AliMUONTrack;
class AliMUONTriggerTrack;
class AliMUONGeometryTransformer;
class AliMUONTriggerCircuit;
class AliMUONLocalTrigger;

class AliMUONRecoCheck : public TObject 
{
public:
  AliMUONRecoCheck(const Char_t *chLoader, const Char_t *pathSim = "./");
  AliMUONRecoCheck(AliESDEvent *esdEvent, AliMCEventHandler *mcEventHandler);
  virtual ~AliMUONRecoCheck();

  /// Return the list of reconstructed tracks
  AliMUONVTrackStore* ReconstructedTracks(Int_t event, Bool_t refit = kTRUE);

  /// Return the list of reconstructed trigger tracks
  AliMUONVTriggerTrackStore* TriggeredTracks(Int_t event);

  void TriggerToTrack(const AliMUONLocalTrigger& locTrg, AliMUONTriggerTrack& triggerTrack);
	
  /// Return reference muon tracks
  AliMUONVTrackStore* TrackRefs(Int_t event);

  /// Return triggerable reference tracks
  AliMUONVTriggerTrackStore* TriggerableTracks(Int_t event);
	
  /// Return reconstructible reference tracks
  AliMUONVTrackStore* ReconstructibleTracks(Int_t event, UInt_t requestedStationMask = 0x1F, Bool_t request2ChInSameSt45 = kTRUE);

	
  /// Return the run number of the current ESD event
  Int_t GetRunNumber();
  
  /// Return the total number of events.
  Int_t NumberOfEvents() const;
  
  /// Return the reconstructed data of current event
  const AliESDEvent* GetESDEvent() { return fESDEvent; }
  
  /// Return the interface to the Monte Carlo data of current event
  const AliMCEventHandler* GetMCEventHandler() { return fMCEventHandler; }
  
  /// Return the track from the store matched with the given track (or 0x0) and the fraction of matched clusters.
  static AliMUONTrack* FindCompatibleTrack(AliMUONTrack &track, AliMUONVTrackStore &trackStore,
					   Int_t &nMatchClusters, Bool_t useLabel = kFALSE, Double_t sigmaCut = 10.);
  
  /// Return the trigger track from the store matched with the given track (or 0x0)
  static AliMUONTriggerTrack* FindCompatibleTrack(AliMUONTriggerTrack &track,
                                                  const AliMUONVTriggerTrackStore &triggerTrackStore,
                                                  Double_t sigmaCut = 10.);
  
private:
  /// Not implemented
  AliMUONRecoCheck(const AliMUONRecoCheck& rhs);
  /// Not implemented
  AliMUONRecoCheck& operator = (const AliMUONRecoCheck& rhs);

  void ResetStores();
  
  void MakeReconstructedTracks(Bool_t refit);
	
  void MakeTriggeredTracks();
  
  void MakeTrackRefs();
  
  void CleanMuonTrackRef(const AliMUONVTrackStore *tmpTrackRefStore);
  
  void MakeReconstructibleTracks(UInt_t requestedStationMask, Bool_t request2ChInSameSt45 = kTRUE);
	
  void MakeTriggerableTracks();
	
  Bool_t InitCircuit();

private:
  AliMCEventHandler* fMCEventHandler; ///< to access MC truth information
  AliESDEvent* fESDEvent; ///< ESD event to access MUON data
  TTree* fESDTree;        ///< ESD tree to access MUON data
  TFile* fESDFile;        ///< ESD file to access MUON data
  
  Int_t fCurrentEvent; ///< current event number
  
  AliMUONVTrackStore* fTrackRefStore;     ///< current simulated tracks (owner)
  AliMUONVTrackStore* fRecoTrackRefStore; ///< current reconstructible tracks (owner)
  AliMUONVTriggerTrackStore* fRecoTriggerRefStore; ///< current triggerable tracks (owner)
  AliMUONVTrackStore* fRecoTrackStore;    ///< current reconstructed tracks (owner)
  AliMUONVTriggerTrackStore* fRecoTriggerTrackStore;    ///< current reconstructed trigger tracks (owner)
	
  AliMUONGeometryTransformer* fGeometryTransformer; ///< geometry transformer
  AliMUONTriggerCircuit* fTriggerCircuit; ///< trigger circuit
  
  Bool_t fESDEventOwner;         ///< using constructor from the analysis task

  ClassDef(AliMUONRecoCheck, 0)   //Utility class to check reconstruction
};

#endif

