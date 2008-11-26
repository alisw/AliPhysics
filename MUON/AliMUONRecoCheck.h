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

class AliMUONRecoCheck : public TObject 
{
public:
  AliMUONRecoCheck(const Char_t *chLoader, const Char_t *pathSim = "./");
  AliMUONRecoCheck(AliESDEvent *esdEvent, AliMCEventHandler *mcEventHandler);
  virtual ~AliMUONRecoCheck();

  /// Return the list of reconstructed tracks
  AliMUONVTrackStore* ReconstructedTracks(Int_t event);
  
  /// Return reference muon tracks
  AliMUONVTrackStore* TrackRefs(Int_t event);

  /// Return reconstructible reference tracks
  AliMUONVTrackStore* ReconstructibleTracks(Int_t event);
  
  /// Return the total number of events.
  Int_t NumberOfEvents() const;
  
  /// Return the reconstructed data of current event
  const AliESDEvent* GetESDEvent() { return fESDEvent; }
  
  /// Return the interface to the Monte Carlo data of current event
  const AliMCEventHandler* GetMCEventHandler() { return fMCEventHandler; }
  
private:
  /// Not implemented
  AliMUONRecoCheck(const AliMUONRecoCheck& rhs);
  /// Not implemented
  AliMUONRecoCheck& operator = (const AliMUONRecoCheck& rhs);

  void ResetStores();
  
  void MakeReconstructedTracks();
  
  void MakeTrackRefs();
  
  void CleanMuonTrackRef(const AliMUONVTrackStore *tmpTrackRefStore);
  
  void MakeReconstructibleTracks();

private:
  AliMCEventHandler* fMCEventHandler; ///< to access MC truth information
  AliESDEvent* fESDEvent; ///< ESD event to access MUON data
  TTree* fESDTree;        ///< ESD tree to access MUON data
  TFile* fESDFile;        ///< ESD file to access MUON data
  
  Int_t fCurrentEvent; ///< current event number
  
  AliMUONVTrackStore* fTrackRefStore;     ///< current simulated tracks (owner)
  AliMUONVTrackStore* fRecoTrackRefStore; ///< current reconstructible tracks (owner)
  AliMUONVTrackStore* fRecoTrackStore;    ///< current reconstructed tracks (owner)
  
  Bool_t fESDEventOwner;         ///< using constructor from the analysis task

  ClassDef(AliMUONRecoCheck, 0)   //Utility class to check reconstruction
};

#endif

