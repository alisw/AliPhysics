#ifndef AliMuonForwardTrack_H
#define AliMuonForwardTrack_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup MFTrec
/// \class AliMuonForwardTrack
/// \brief ALICE muon forward track, combining the information of the Muon Spectrometer and the Muon Forward Tracker
///
/// \author Antonio Uras <antonio.uras@cern.ch>, Bogdan Vulpescu <bogdan.vulpescu@clermont.in2p3.fr>
/// \date February 2nd, 2016

#include "AliMUONTrack.h"

class AliMFTCluster;

class AliMuonForwardTrack : public AliMUONTrack {

public:

  AliMuonForwardTrack();
  AliMuonForwardTrack(AliMUONTrack& MUONTrack);

  AliMuonForwardTrack(const AliMuonForwardTrack&);
  AliMuonForwardTrack &operator=(const AliMuonForwardTrack&);

  /// \brief Set the MC label of the attached MFT track
  void SetTrackMCId(Int_t id) { fTrackMCId = id; }
  /// \brief Get the MC label of the attached MFT track
  Int_t GetTrackMCId() { return fTrackMCId; }

  /// \brief overload of the AliMUONTrack function
  void AddTrackParamAtMFTCluster(AliMUONTrackParam &trackParam, AliMUONVCluster &cluster, const Int_t mftid); 

protected:

  Int_t fTrackMCId;   ///< \brief MC label of the attached MFT track

  ClassDef(AliMuonForwardTrack,1)
    
};

//====================================================================================================================================================

#endif



