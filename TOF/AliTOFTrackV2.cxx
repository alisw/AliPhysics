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
// AliTOFTrackV2 class
//
// Author: F. Pierella | pierella@bo.infn.it
//
// Description: output of AliTOFReconstructionerV2
// TRD tracking capabilities have been foreseen by including the member
// variables fxTRD, fyTRD, fzTRD and momentum components fPxTRD, fPyTRD, fPzTRD.
// Class to be written in TreeR for each event.
//
// Member variable summary description: 
// - track momentum and position in the last TPC padrow 
// - track length as given by the geometrical propagation
// - reconstructed mass from time of flight and time of flight itself
// - fMatchingStatus
//             -2 backpropagation goes out of the z acceptance of the TOF
//             -1 failed backpropagation on TOF inner radius
//              0 for fake tracks
//              1 for tracks matched with no signal on TOF (failed DigitFinder)
//              3 for tracks matched with the actual digit
//              4 for tracks matched with a wrong (not its own) TOF digit
////////////////////////////////////////////////////////////////////////

#include "AliTOFTrackV2.h"

ClassImp(AliTOFTrackV2)

AliTOFTrackV2::AliTOFTrackV2() 
{
  //
  // Default ctor
  //
  fTrackLabel=-1;
  fTOFDigitTrackLabel=-1;
  fPTPC=-1;
  fPdgCode=-1;
  fdEdX=-1;
  fxTPC=-1;
  fyTPC=-1;
  fzTPC=-1;
  fPtTPC=-1;
  fPzTPC=-1;
  fxTRD=-1;
  fyTRD=-1;
  fzTRD=-1;
  fPxTRD=-1;
  fPyTRD=-1;
  fPzTRD=-1;
  fMatchingStatus=-1;
  fLength=-1;
  fTof=-1;
  fMassTOF=-1;
  // vertex variables from reconstruction
  fXRecVtx=-1;
  fYRecVtx=-1;
  fZRecVtx=-1;
  fPxRecVtx=-1;
  fPyRecVtx=-1;
  fPzRecVtx=-1;
  fRecTrackLength=-1;
}

//_____________________________________________________________________________
AliTOFTrackV2::AliTOFTrackV2(Int_t trackLabel, Int_t matchingStatus, 
			     Float_t tpcMom, Float_t dEdX, Float_t* tpcXYZ, 
			     Float_t* tpcPtPz, Float_t* /* trdXYZ */, 
			     Float_t* /* trdPxPyPz */)
{
  //
  // par ctor
  //

  fTrackLabel=trackLabel;
  fTOFDigitTrackLabel=-1;
  fPTPC=tpcMom;
  fPdgCode=-1;

  fdEdX=dEdX;
  fxTPC=tpcXYZ[0];
  fyTPC=tpcXYZ[1];
  fzTPC=tpcXYZ[2];
  fPtTPC=tpcPtPz[0];
  fPzTPC=tpcPtPz[1];

  fxTRD=-1;
  fyTRD=-1;
  fzTRD=-1;
  fPxTRD=-1;
  fPyTRD=-1;
  fPzTRD=-1;

  fMatchingStatus=matchingStatus;
  fLength=-1;
  fTof=-1;
  fMassTOF=-1;

  // vertex variables from reconstruction
  fXRecVtx=-1;
  fYRecVtx=-1;
  fZRecVtx=-1;
  fPxRecVtx=-1;
  fPyRecVtx=-1;
  fPzRecVtx=-1;
  fRecTrackLength=-1;
}

//_____________________________________________________________________________
void AliTOFTrackV2::UpdateTrack(Int_t tofDigitTrackLabel, Int_t matching, Float_t tof)
{
  //
  // update the track after the TOF digit assignment
  //
  fTOFDigitTrackLabel=tofDigitTrackLabel;
  fMatchingStatus=matching;
  fTof=tof;
}

//_____________________________________________________________________________
void AliTOFTrackV2::UpdateTrack(Int_t pdgCode, Float_t trackLength)
{
  //
  // update the track 
  //
  fPdgCode=pdgCode;
  fLength=trackLength;
}
