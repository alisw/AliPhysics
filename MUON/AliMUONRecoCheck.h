#ifndef ALIMUONRECOCHECK_H
#define ALIMUONRECOCHECK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONRecoCheck
/// \brief Utility class to check reconstruction

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliMUONRecoCheck                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TObject.h>
#include "AliMUONTrack.h"

class TClonesArray;
class AliMUONData;
class AliRunLoader;


class AliMUONRecoCheck : public TObject 
{
public:
  AliMUONRecoCheck(Char_t *chLoader);
  virtual          ~AliMUONRecoCheck();

  AliMUONData*  GetMUONData() {return fMUONData;}
  void MakeTrackRef();
  void AddMuonTrackReference(const AliMUONTrack *muonTrack) 
    {new ((*fMuonTrackRef)[fMuonTrackRef->GetEntriesFast()]) AliMUONTrack(*muonTrack);}
  void PrintEvent() const;
  void ResetTracks() const;
  AliRunLoader* GetRunLoader() {return fRunLoader;}
  void CleanMuonTrackRef();
  void ReconstructibleTracks();
  Int_t GetNumberOfReconstuctibleTracks() {return fReconstructibleTracks;}
  Int_t GetNumberOfRecoTracks() {return fRecoTracks;}
  TClonesArray *GetTrackReco();
  TClonesArray *GetMuonTrackRef() {return fMuonTrackRef;}

protected:
  AliMUONRecoCheck(const AliMUONRecoCheck& rhs);
  AliMUONRecoCheck& operator = (const AliMUONRecoCheck& rhs);

private:
  
  AliRunLoader* fRunLoader;     // alice run loader 
  AliMUONData*  fMUONData;      // Data container for MUON subsystem 
  TClonesArray* fMuonTrackRef;  // reference muon tracks
  TClonesArray* fTrackReco;     // reconstructed muon tracks
  Int_t fReconstructibleTracks; // number of reconstructible tracks 
  Int_t fRecoTracks; // number of reconstructed tracks 

  ClassDef(AliMUONRecoCheck, 0)   //Utility class to check reconstruction
};

#endif








