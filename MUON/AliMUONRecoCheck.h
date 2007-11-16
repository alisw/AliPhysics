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
class AliMCEventHandler;
class AliMUONDataInterface;
class AliMUONVTrackStore;
class AliESDEvent;

class AliMUONRecoCheck : public TObject 
{
public:
  AliMUONRecoCheck(Char_t *chLoader, Char_t *pathSim = "./");
  virtual ~AliMUONRecoCheck();

  /// Return the list of reconstructed tracks
  AliMUONVTrackStore* ReconstructedTracks(Int_t event);
  
  /// Create and return a list of reconstructed tracks from ESD data.
  static AliMUONVTrackStore* ReconstructedTracks(AliESDEvent* esd, Bool_t padMissing = kFALSE);
  
  /// Return reference muon tracks
  AliMUONVTrackStore* TrackRefs(Int_t event);

  /// Return reconstructible reference tracks
  AliMUONVTrackStore* ReconstructibleTracks(Int_t event);
  
  /// Return the total number of events.
  Int_t NumberOfEvents() const;
  
  /// Return the interface to the reconstructed data.
  AliMUONDataInterface* GetDataInterface() { return fDataInterface; };
  
  /// Return the interface to the Monte Carlo data.
  AliMCEventHandler* GetMCEventHandler() { return fMCEventHandler; };
  
private:
  /// Not implemented
  AliMUONRecoCheck(const AliMUONRecoCheck& rhs);
  /// Not implemented
  AliMUONRecoCheck& operator = (const AliMUONRecoCheck& rhs);

  void ResetStores();
  
  void MakeTrackRefs();
  
  void CleanMuonTrackRef(const AliMUONVTrackStore *tmpTrackRefStore);
  
  void MakeReconstructibleTracks();

private:
  AliMCEventHandler* fMCEventHandler; ///< to access MC truth information
  AliMUONDataInterface* fDataInterface; ///< to access MUON data
  
  Int_t fCurrentEvent; ///< current event number
  
  AliMUONVTrackStore* fTrackRefStore;     ///< current simulated tracks (owner)
  AliMUONVTrackStore* fRecoTrackRefStore; ///< current reconstructible tracks (owner)
  AliMUONVTrackStore* fRecoTrackStore;    ///< current reconstructed tracks (owner)
  
  ClassDef(AliMUONRecoCheck, 0)   //Utility class to check reconstruction
};

#endif

