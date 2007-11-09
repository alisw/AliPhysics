/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliMeanVertex.h"

/* $Id$ */

/*****************************************************************************
 *                                                                           *
 * This class contains the coordinates of the mean primary vertex position   *
 * computed by AliITSMeanVertex                                              *
 *                                                                           *
*****************************************************************************/

ClassImp(AliMeanVertex)

//--------------------------------------------------------------------------
AliMeanVertex::AliMeanVertex() :
  AliESDVertex(),
  fTotTracklets(0.),
  fAverTracklets(0.),
  fSigmaOnAverTrack(0.) {
  //
  // Default Constructor, set everything to 0
  //
  for(Int_t i=0;i<3;i++){
    fErrW[i]=0.;
  }
}

//--------------------------------------------------------------------------
AliMeanVertex::AliMeanVertex(Double_t pos[3],Double_t err[3],Double_t cov[6],Int_t nevents, Float_t notracklets, Float_t avertracklets, Float_t signotrackl) :
  AliESDVertex(pos,cov,0.,nevents,"MeanVertex"),
  fTotTracklets(notracklets),
  fAverTracklets(avertracklets),
  fSigmaOnAverTrack(signotrackl) {
  //
  // Standard Constructor
  //
  for(Int_t i=0;i<3;i++){
    fErrW[i]=err[i];
  }
}

//--------------------------------------------------------------------------
void AliMeanVertex::GetErrorsOnPosition(Double_t err[3]) const {
  //
  // Return errors on vertex average position (weighted mean on several events)
  //
  err[0] = fErrW[0];
  err[1] = fErrW[1];
  err[2] = fErrW[2];

  return;
}
