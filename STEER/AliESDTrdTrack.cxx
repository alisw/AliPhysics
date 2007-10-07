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
  fPt(0),
  fPhi(0),
  fEta(0),
  fPID(0),
  fLabel(0),
  fNtracklets(0),
  fNclusters(0),
  fNplanes(0),
  fDetector(0)
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
  fPt(track.fPt),
  fPhi(track.fPhi),
  fEta(track.fEta),
  fPID(track.fPID),
  fLabel(track.fLabel),
  fNtracklets(track.fNtracklets),
  fNclusters(track.fNclusters),
  fNplanes(track.fNplanes),
  fDetector(track.fDetector)
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
  TObject::operator=(track);
  fYproj      = track.fYproj;
  fZproj      = track.fZproj;
  fSlope      = track.fSlope;
  fPt         = track.fPt;
  fPhi        = track.fPhi;
  fEta        = track.fEta;
  fPID        = track.fPID;
  fLabel      = track.fLabel;
  fNtracklets = track.fNtracklets;
  fNclusters  = track.fNclusters;
  fDetector   = track.fDetector;
  fNplanes    = track.fNplanes;

  return *this;

}

