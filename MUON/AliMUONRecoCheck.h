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

class AliMUONRecoCheck : public TObject 
{
public:
  AliMUONRecoCheck(Char_t *chLoader, Char_t *pathSim = "./");
  virtual ~AliMUONRecoCheck();

  Int_t NumberOfEvents() const;
  
  AliMUONVTrackStore* ReconstructedTracks(Int_t event);
  
  AliMUONVTrackStore* TrackRefs(Int_t event);

  AliMUONVTrackStore* ReconstructibleTracks(Int_t event);
  
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

