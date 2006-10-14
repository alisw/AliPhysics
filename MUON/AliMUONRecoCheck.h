#ifndef ALIMUONRECOCHECK_H
#define ALIMUONRECOCHECK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONRecoCheck
/// \brief Utility class to check reconstruction

#include <TObject.h>
#include "AliMUONTrack.h"

class TClonesArray;
class AliMUONData;
class AliRunLoader;


class AliMUONRecoCheck : public TObject 
{
public:
  AliMUONRecoCheck(AliRunLoader *runloader, AliMUONData *muondata);
  virtual          ~AliMUONRecoCheck();

  void MakeTrackRef();
                /// Add track reference
  void AddMuonTrackReference(const AliMUONTrack *muonTrack) 
    {new ((*fMuonTrackRef)[fMuonTrackRef->GetEntriesFast()]) AliMUONTrack(*muonTrack);}

  void PrintEvent() const;
  void ResetTracks() const;
  void CleanMuonTrackRef();
  void ReconstructibleTracks();
                /// Return number of reconstructible tracks
  Int_t GetNumberOfReconstuctibleTracks() {return fReconstructibleTracks;}
                /// Return number of reconstructed tracks
  Int_t GetNumberOfRecoTracks() {return fRecoTracks;}
  TClonesArray *GetTrackReco();
                /// Return reference muon tracks
  TClonesArray *GetMuonTrackRef() {return fMuonTrackRef;}

private:
  
  AliRunLoader* fRunLoader;     ///< alice run loader 
  AliMUONData*  fMUONData;      ///< Data container for MUON subsystem 
  TClonesArray* fMuonTrackRef;  ///< reference muon tracks
  TClonesArray* fTrackReco;     ///< reconstructed muon tracks
  Int_t fReconstructibleTracks; ///< number of reconstructible tracks 
  Int_t fRecoTracks;            ///< number of reconstructed tracks 

  AliMUONRecoCheck(const AliMUONRecoCheck& rhs);
  AliMUONRecoCheck& operator = (const AliMUONRecoCheck& rhs);

  ClassDef(AliMUONRecoCheck, 0)   //Utility class to check reconstruction
};

#endif








