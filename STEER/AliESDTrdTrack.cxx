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

///////////////////////////////////////////////////////////////////////////////
//
//  Tracks from the TRD Global Tracking Unit (GTU, trigger)
//  The TRD trigger stores the found tracks 
//  as ESDTrdTrack objects in the ESD object
//  Related classes: AliTRDReconstructor, AliESD
//  Author: B.Vulpescu
//
///////////////////////////////////////////////////////////////////////////////

#include "AliESDTrdTrack.h"

ClassImp(AliESDTrdTrack)

//_____________________________________________________________________________
AliESDTrdTrack::AliESDTrdTrack():
  TObject(),
  fYproj(0),
  fZproj(0),
  fSlope(0),
  fDetector(-1),
  fNtracklets(0),
  fNplanes(0),
  fNclusters(0),
  fPt(0),
  fPhi(0),
  fEta(0),
  fLabel(-1),
  fPID(0),
  fIsElectron(kFALSE)
{

  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliESDTrdTrack::AliESDTrdTrack(const AliESDTrdTrack& track):
  TObject(track),
  fYproj(track.fYproj),
  fZproj(track.fZproj),
  fSlope(track.fSlope),
  fDetector(track.fDetector),
  fNtracklets(track.fNtracklets),
  fNplanes(track.fNplanes),
  fNclusters(track.fNclusters),
  fPt(track.fPt),
  fPhi(track.fPhi),
  fEta(track.fEta),
  fLabel(track.fLabel),
  fPID(track.fPID),
  fIsElectron(track.fIsElectron)
{

  //
  // Copy contructor
  //

}

//_____________________________________________________________________________
AliESDTrdTrack& AliESDTrdTrack::operator=(const AliESDTrdTrack& track)
{
  // 
  // Equal operator
  //

  if (this == &track)
    return *this;

  fYproj      = track.fYproj;
  fZproj      = track.fZproj;
  fSlope      = track.fSlope;
  fDetector   = track.fDetector;
  fNtracklets = track.fNtracklets;
  fNplanes    = track.fNplanes;
  fNclusters  = track.fNclusters;
  fPt         = track.fPt;
  fPhi        = track.fPhi;
  fEta        = track.fEta;
  fLabel      = track.fLabel;
  fPID        = track.fPID;
  fIsElectron = track.fIsElectron;

  return *this;

}

