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

////////////////////////////////////////////////////////////////////////
//
// AliTOFTrack class
//
// Authors: Bologna-ITEP-Salerno Group
//
// Description: class for TPC data about extracted tracks for TPC-TOF
// matching (cfr. TDR Section 4.4).
// TRD tracking capabilities have been foreseen by including the member
// variables fRxTRD and fRyTRD.
// It is the candidate class to be written in TreeR for each event.
//
// Member variable summary description: 
// - track momentum and position in the last TPC padrow 
// - track length as given by the geometrical propagation
// - reconstructed mass from time of flight
//
////////////////////////////////////////////////////////////////////////

#include "AliTOFTrack.h"

ClassImp(AliTOFTrack)

AliTOFTrack::AliTOFTrack() 
{
  //
  // Default ctor
  //
  fTrack=0;
  fP=0.;
  fPdgCode=0;
  flTPC=0.;
  fRxTPC=0.;
  fRyTPC=0.;
  fRzTPC=0.;
  fPxTPC=0.;
  fPyTPC=0.;
  fPzTPC=0.;
  fRxTRD=0.;
  fRyTRD=0.;
  fPad=0;
  fMatching=0;
  fLength=-1;
  fTof=-1;
  fMassTOF=-1;
}

//_____________________________________________________________________________
AliTOFTrack::AliTOFTrack(Int_t track, Float_t vtxMom, Int_t pdgcode, Float_t tpcTrackLen, Float_t* tpcpos, Float_t* tpcmom, Float_t* trdpos, Int_t pad, Int_t matchflag, Float_t length, Float_t tof, Float_t massTof)
{
  //
  // par ctor
  //
  fTrack  = track;
  fP      = vtxMom;
  fPdgCode= pdgcode;
  flTPC   = tpcTrackLen;
  fRxTPC  = tpcpos[0];
  fRyTPC  = tpcpos[1];
  fRzTPC  = tpcpos[2];
  fPxTPC  = tpcmom[0];
  fPyTPC  = tpcmom[1];
  fPzTPC  = tpcmom[2];
  fRxTRD  = trdpos[0];
  fRyTRD  = trdpos[1];
  fPad    = pad;
  fMatching= matchflag;
  fLength = length;
  fTof    = tof;
  fMassTOF= massTof;
}

//_____________________________________________________________________________
void AliTOFTrack::SetTrack(Int_t track, Float_t vtxMom, Int_t pdgcode, Float_t tpcTrackLen, Float_t* tpcpos, Float_t* tpcmom, Float_t* trdpos)
{
  //
  // setter
  //
  fTrack  = track;
  fP      = vtxMom;
  fPdgCode= pdgcode;
  flTPC   = tpcTrackLen;
  fRxTPC  = tpcpos[0];
  fRyTPC  = tpcpos[1];
  fRzTPC  = tpcpos[2];
  fPxTPC  = tpcmom[0];
  fPyTPC  = tpcmom[1];
  fPzTPC  = tpcmom[2];
  fRxTRD  = trdpos[0];
  fRyTRD  = trdpos[1];
}
