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

/* $Id$ */

///////////////////////////////////////////////////
//
// Interface
// for
// MUON
// track
// extrapolation
//
///////////////////////////////////////////////////

#include "AliPWG3TrackExtrapInterface.h"
#include "AliMUONTrackExtrap.h"
#include "AliESDMuonTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMagF.h" 

#include <Riostream.h>
#include <TGeoManager.h>

/// \cond CLASSIMP
ClassImp(AliPWG3TrackExtrapInterface) // Class implementation in ROOT context
/// \endcond

  //__________________________________________________________________________
void AliPWG3TrackExtrapInterface::SetMagField(const AliMagF* magField)
{
  /// Set the magnetic field (required for track extrapolation)
  
  AliMUONTrackExtrap::SetField(magField);
  
}

  //__________________________________________________________________________
Bool_t AliPWG3TrackExtrapInterface::SetGeometry(const char* fileName)
{
  /// Set the geometry (required for absorber correction)
  
  if (!gGeoManager) {
    TGeoManager::Import(fileName);
    if (!gGeoManager) {
      cout<<"E-AliPWG3TrackExtrapInterface::ImportGeo: getting geometry from file "<<fileName<<" failed"<<endl;
      return kFALSE;
    }
  }
  return kTRUE;
  
}

  //__________________________________________________________________________
void AliPWG3TrackExtrapInterface::ExtrapToVertexUncorrected(AliESDMuonTrack* muonTrack, Double_t zVtx)
{
  /// Extrapolation to the vertex (at the z position "zVtx") without Branson and energy loss corrections.
  /// Returns the extrapolated track parameters in the current ESDMuonTrack.
  
  AliMUONTrackParam trackParam;
  trackParam.GetParamFromUncorrected(*muonTrack);
  // dummy vertex transverse position (not used anyway)
  Double_t xVtx = 0.;
  Double_t yVtx = 0.;
  AliMUONTrackExtrap::ExtrapToVertex(&trackParam, xVtx, yVtx, zVtx, kFALSE, kFALSE);
  trackParam.SetParamFor(*muonTrack);
  
}

  //__________________________________________________________________________
void AliPWG3TrackExtrapInterface::ExtrapToVertexWithELoss(AliESDMuonTrack* muonTrack, Double_t zVtx)
{
  /// Extrapolation to the vertex (at the z position "zVtx") with energy loss correction only.
  /// Returns the extrapolated track parameters in the current ESDMuonTrack.
  
  AliMUONTrackParam trackParam;
  trackParam.GetParamFromUncorrected(*muonTrack);
  // compute energy loss correction assuming linear propagation
  Double_t xVtx = trackParam.GetNonBendingCoor() + (zVtx - trackParam.GetZ()) * trackParam.GetNonBendingSlope();
  Double_t yVtx = trackParam.GetBendingCoor()    + (zVtx - trackParam.GetZ()) * trackParam.GetBendingSlope();
  AliMUONTrackExtrap::ExtrapToVertex(&trackParam, xVtx, yVtx, zVtx, kFALSE, kTRUE);
  trackParam.SetParamFor(*muonTrack);
  
}

  //__________________________________________________________________________
void AliPWG3TrackExtrapInterface::ExtrapToVertexWithBranson(AliESDMuonTrack* muonTrack, Double_t vtx[3])
{
  /// Extrapolation to the vertex (at the z position "zVtx") with Branson correction only.
  /// Returns the extrapolated track parameters in the current ESDMuonTrack.
  
  AliMUONTrackParam trackParam;
  trackParam.GetParamFromUncorrected(*muonTrack);
  AliMUONTrackExtrap::ExtrapToVertex(&trackParam, vtx[0], vtx[1], vtx[2], kTRUE, kFALSE);
  trackParam.SetParamFor(*muonTrack);
  
}

  //__________________________________________________________________________
void AliPWG3TrackExtrapInterface::ExtrapToVertex(AliESDMuonTrack* muonTrack, Double_t vtx[3])
{
  /// Extrapolation to the vertex (at the z position "zVtx") with Branson and energy loss corrections.
  /// Returns the extrapolated track parameters in the current ESDMuonTrack.
  
  AliMUONTrackParam trackParam;
  trackParam.GetParamFromUncorrected(*muonTrack);
  AliMUONTrackExtrap::ExtrapToVertex(&trackParam, vtx[0], vtx[1], vtx[2], kTRUE, kTRUE);
  trackParam.SetParamFor(*muonTrack);
  
}

  //__________________________________________________________________________
Double_t AliPWG3TrackExtrapInterface::TotalMomentumEnergyLoss(AliESDMuonTrack* muonTrack, Double_t zVtx)
{
  /// Calculate the total momentum energy loss in-between the track position and the vertex assuming a linear propagation
  AliMUONTrackParam trackParam;
  trackParam.GetParamFromUncorrected(*muonTrack);
  Double_t xVtx = trackParam.GetNonBendingCoor() + (zVtx - trackParam.GetZ()) * trackParam.GetNonBendingSlope();
  Double_t yVtx = trackParam.GetBendingCoor()    + (zVtx - trackParam.GetZ()) * trackParam.GetBendingSlope();
  return AliMUONTrackExtrap::TotalMomentumEnergyLoss(&trackParam, xVtx, yVtx, zVtx);
}

  //__________________________________________________________________________
Double_t AliPWG3TrackExtrapInterface::TotalMomentumEnergyLoss(AliESDMuonTrack* muonTrack, Double_t vtx[3])
{
  /// Calculate the total momentum energy loss in-between the track position and the vertex assuming a linear propagation
  AliMUONTrackParam trackParam;
  trackParam.GetParamFromUncorrected(*muonTrack);
  return AliMUONTrackExtrap::TotalMomentumEnergyLoss(&trackParam, vtx[0], vtx[1], vtx[2]);
}

