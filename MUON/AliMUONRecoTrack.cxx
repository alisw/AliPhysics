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

// ----------------------                                                               
// Class AliMUONRecoTrack 
// ----------------------                                                               
// This class represents a reconstructed muon track
// in the tree of reconstructed events.            
// Authors: Mihaela Gheata, Andrei Gheata 09/10/00

#include <Riostream.h>

#include "AliMUONRecoTrack.h"

/// \cond CLASSIMP
ClassImp(AliMUONRecoTrack)
/// \endcond

//-------------------------------------------------------------------
AliMUONRecoTrack::AliMUONRecoTrack() 
  : TObject(), 
    fSign(0), 
    fFlag(0), 
    fZvr(0.), 
    fChi2r(0.) 
{ 
/// Default constructor
}

//-------------------------------------------------------------------
AliMUONRecoTrack::AliMUONRecoTrack(Bool_t active)
  : TObject(),
    fSign(0),
    fFlag(0),
    fZvr(0.),
    fChi2r(0.)
{
/// Constructor of AliMUONRecoTrack
   if (active) {
   	for (Int_t axis=0; axis<3; axis++) {
      	fPr[axis] = 0.0;
   	}
   	for (Int_t chamber=0; chamber<10; chamber++) {
      	fPosX[chamber] = 0.0;
      	fPosY[chamber] = 0.0;
      	fPosZ[chamber] = 0.0;
   	}
   }
}

//-------------------------------------------------------------------
Double_t AliMUONRecoTrack::Phi() const
{
/// Return trach phi angle
	return TMath::ATan2(fPr[2], fPr[1]);
}

//-------------------------------------------------------------------
Double_t AliMUONRecoTrack::Theta() const
{
/// Return trach theta angle
   return TMath::ACos(fPr[2] / P());
}

//-------------------------------------------------------------------
void AliMUONRecoTrack::SetMomReconstr(Double_t px, Double_t py, Double_t pz)
{
/// Set the track reconstructed momentum 
   fPr[0] = px;
   fPr[1] = py;
   fPr[2] = pz;            
} 
   
//-------------------------------------------------------------------		
void AliMUONRecoTrack::SetHitPosition(Int_t chamber, Double_t x, Double_t y, Double_t z)
{
/// Set hit coordinates in chamber[0..9]
   fPosX[chamber] = x;
   fPosY[chamber] = y;
   fPosZ[chamber] = z;
}

//-------------------------------------------------------------------		
void AliMUONRecoTrack::TrackInfo() const
{
/// Prints momentum info for this track
   cout << "Px=" << GetMomReconstr(0) << " Py=" << GetMomReconstr(1) <<
          " Pz=" << GetMomReconstr(2) << " P=" << P() << endl;
}
