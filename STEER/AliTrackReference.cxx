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


#include "AliTrackReference.h"
#include "TParticle.h"
#include "AliRun.h"
#include "TLorentzVector.h"

// 
// Track Reference object is created every time particle is 
// crossing detector bounds. The object is created by Step Manager
//
// The class stores the following informations:
// track label, 
// track position: X,Y,X
// track momentum px, py, pz
// track length and time of fligth: both in cm
// status bits from Monte Carlo
//


ClassImp(AliTrackReference)

//_______________________________________________________________________
 AliTrackReference::AliTrackReference():
   fTrack(0),
   fX(0),
   fY(0),
   fZ(0),
   fPx(0),
   fPy(0),
   fPz(0),
   fLength(0)
{
  //
  // Default constructor
  // Creates empty object

  for(Int_t i=0; i<16; i++) ResetBit(BIT(i));
}

//_______________________________________________________________________
AliTrackReference::AliTrackReference(Int_t label, TVirtualMC *vMC) {
  //
  // Create Reference object out of label and
  // data in TVirtualMC object
  //
  // Creates an object and fill all parameters 
  // from data in VirtualMC
  //
  // Sylwester Radomski, (S.Radomski@gsi.de)
  // GSI, Jan 31, 2003
  //

  TLorentzVector vec;
  
  fTrack = label;
  fLength = vMC->TrackLength();
  fTime = vMC->TrackTime();

  vMC->TrackPosition(vec);

  fX = vec[0];
  fY = vec[1];
  fZ = vec[2];
  
  vMC->TrackMomentum(vec);
  
  fPx = vec[0];
  fPy = vec[1];
  fPy = vec[2];

  // Set Up status code 
  // Copy Bits from virtual MC

  for(Int_t i=0; i<16; i++) ResetBit(BIT(i));

  SetBit(BIT(0), vMC->IsNewTrack());
  SetBit(BIT(1), vMC->IsTrackAlive());
  SetBit(BIT(2), vMC->IsTrackDisappeared());
  SetBit(BIT(3), vMC->IsTrackEntering());
  SetBit(BIT(4), vMC->IsTrackExiting());
  SetBit(BIT(5), vMC->IsTrackInside());
  SetBit(BIT(6), vMC->IsTrackOut());
  SetBit(BIT(7), vMC->IsTrackStop()); 
}
//_______________________________________________________________________
