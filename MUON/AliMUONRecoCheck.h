#ifndef ALIMUONRECOCHECK_H
#define ALIMUONRECOCHECK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup evaluation
/// \class AliMUONRecoCheck
/// \brief Utility class to check reconstruction

#include <TObject.h>
#include "AliMUONTrack.h"

class TClonesArray;
class AliMUONMCDataInterface;
class AliMUONDataInterface;
class AliMUONVTrackStore;

class AliMUONRecoCheck : public TObject 
{
public:
  AliMUONRecoCheck(Char_t *chLoader, Char_t *chLoaderSim);
  virtual ~AliMUONRecoCheck();

  /// Return number of reconstructed tracks
  AliMUONVTrackStore* ReconstructedTracks(Int_t event);
  
  /// Return reference muon tracks
  AliMUONVTrackStore* TrackRefs(Int_t event);

  /// Return reconstructible ref tracks
  AliMUONVTrackStore* ReconstructibleTracks(Int_t event);
  
  Int_t NumberOfEvents() const;
  
private:
  /// Not implemented
  AliMUONRecoCheck(const AliMUONRecoCheck& rhs);
  /// Not implemented
  AliMUONRecoCheck& operator = (const AliMUONRecoCheck& rhs);

  AliMUONVTrackStore* MakeReconstructibleTracks(const AliMUONVTrackStore& refTracks);

  AliMUONVTrackStore* MakeTrackRefs(Int_t event);
  
  AliMUONVTrackStore* CleanMuonTrackRef(const AliMUONVTrackStore& refTracks);
  
private:
    
  AliMUONMCDataInterface* fMCDataInterface; ///< to access MC information
  AliMUONDataInterface* fDataInterface; ///< to access MUON data
  
  ClassDef(AliMUONRecoCheck, 0)   //Utility class to check reconstruction
};

#endif








