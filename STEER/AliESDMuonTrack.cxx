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

///////////////////////////////////////////////////////////////////////////////
//
//  Class to describe the MUON tracks
//  in the Event Summary Data class
//  This is where the results of reconstruction
//  are stored for the muons
//  Author: G.Martinez
//
///////////////////////////////////////////////////////////////////////////////


#include "AliESDMuonTrack.h"

ClassImp(AliESDMuonTrack)

//_____________________________________________________________________________
AliESDMuonTrack::AliESDMuonTrack ():
  TObject(),
  fInverseBendingMomentum(0),
  fThetaX(0),
  fThetaY(0),
  fZ(0),
  fBendingCoor(0),
  fNonBendingCoor(0),
  fChi2(0),
  fNHit(0),
  fMatchTrigger(0),
  fChi2MatchTrigger(0)
{
  // Default constructor
}


//_____________________________________________________________________________
AliESDMuonTrack::AliESDMuonTrack (const AliESDMuonTrack& MUONTrack):
  TObject(MUONTrack),
  fInverseBendingMomentum(MUONTrack.fInverseBendingMomentum),
  fThetaX(MUONTrack.fThetaX),
  fThetaY(MUONTrack.fThetaY),
  fZ(MUONTrack.fZ),
  fBendingCoor(MUONTrack.fBendingCoor),
  fNonBendingCoor(MUONTrack.fNonBendingCoor),
  fChi2(MUONTrack.fChi2),
  fNHit(MUONTrack.fNHit),
  fMatchTrigger(MUONTrack.fMatchTrigger),
  fChi2MatchTrigger(MUONTrack.fChi2MatchTrigger)
{
  //
  // Copy constructor
  // Deep copy implemented
  //
}

//_____________________________________________________________________________
AliESDMuonTrack& AliESDMuonTrack::operator=(const AliESDMuonTrack& MUONTrack)
{
  // 
  // Equal operator for a deep copy
  //
  if (this == &MUONTrack)
    return *this;

  fInverseBendingMomentum = MUONTrack.fInverseBendingMomentum; 
  fThetaX                 = MUONTrack.fThetaX;           
  fThetaY                 = MUONTrack.fThetaY ;           
  fZ                      = MUONTrack.fZ;                
  fBendingCoor            = MUONTrack.fBendingCoor;      
  fNonBendingCoor         = MUONTrack.fNonBendingCoor;   
  fChi2                   = MUONTrack.fChi2;             
  fNHit                   = MUONTrack.fNHit ; 

  fMatchTrigger           = MUONTrack.fMatchTrigger;  
  fChi2MatchTrigger       = MUONTrack.fChi2MatchTrigger; 
 
  return *this;
}


