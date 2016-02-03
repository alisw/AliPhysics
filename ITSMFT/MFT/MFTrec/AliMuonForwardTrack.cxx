/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//====================================================================================================================================================
//
//      Description of an ALICE muon forward track, combining the information of the Muon Spectrometer and the Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliMUONTrackParam.h"
#include "AliMUONRawCluster.h"
#include "AliMuonForwardTrack.h"

ClassImp(AliMuonForwardTrack)

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack():
AliMUONTrack()
{

  // default constructor
  
}

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack(AliMUONTrack& muonTrack):
AliMUONTrack(muonTrack)
{

  fTrackMCId = -1;

}

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack(const AliMuonForwardTrack& track): 
AliMUONTrack(track),
fTrackMCId(track.fTrackMCId)
{

  // copy constructor

}


//====================================================================================================================================================

AliMuonForwardTrack& AliMuonForwardTrack::operator=(const AliMuonForwardTrack& track) 
{

  // assignment operator

  // check assignement to self
  if (this == &track) return *this;

  // base class assignement
  AliMUONTrack::operator=(track);

  fTrackMCId = track.fTrackMCId;
  
}

//====================================================================================================================================================

void AliMuonForwardTrack::AddTrackParamAtMFTCluster(AliMUONTrackParam &trackParam, AliMUONVCluster &muonCluster, const Int_t mftid) {

  trackParam.SetUniqueID(5000000+mftid); 
  AddTrackParamAtCluster(trackParam, muonCluster, kTRUE);

}

