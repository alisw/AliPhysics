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
//      Class for the description of the clusters of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TObject.h"
#include "AliMUONRawCluster.h"
#include "AliMUONVCluster.h"
#include "AliMFTCluster.h"

ClassImp(AliMFTCluster)

//====================================================================================================================================================

AliMFTCluster::AliMFTCluster():
  TObject(),
  fX(0), 
  fY(0), 
  fZ(0),
  fErrX(0), 
  fErrY(0), 
  fErrZ(0),
  fNElectrons(0),
  fNMCTracks(0),
  fPlane(0),
  fSize(0),
  fTrackChi2(0),
  fLocalChi2(0)
{

  // default constructor

  for (Int_t iTrack=0; iTrack<fNMaxMCTracks; iTrack++) fMCLabel[iTrack] = -1;

}

//====================================================================================================================================================

AliMUONRawCluster* AliMFTCluster::CreateMUONCluster() {

  AliMUONRawCluster *cluster = new AliMUONRawCluster();
  
  cluster->SetXYZ(GetX(), GetY(), GetZ());
  cluster->SetErrXY(GetErrX(),GetErrY());
  cluster->SetDetElemId(100);   // to get the cluster compatible with the AliMUONTrack::AddTrackParamAtCluster(...) method

  return cluster;

}

//====================================================================================================================================================
