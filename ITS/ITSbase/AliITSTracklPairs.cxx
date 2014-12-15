/**************************************************************************
 * Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////////////
// Helper class for 3D primary vertexing                           //
// Used by AliITSSortTrkl                                          //
// In each object the labels of the 2 paired tracklets             //
// are stored. DCA and crossing point coordinates are also stored. //
// Origin M.Masera (masera@to.infn.it)                             //
/////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include "AliITSTracklPairs.h"

ClassImp(AliITSTracklPairs)

//______________________________________________________________________
AliITSTracklPairs::AliITSTracklPairs():TObject(),
fTrack1(0),
fTrack2(0),
fDCA(0.) {
  // Default constructor

  for(Int_t i=0;i<3;i++)fCross[i]=0.;
}

//______________________________________________________________________
AliITSTracklPairs::AliITSTracklPairs(Int_t t1, Int_t t2, Double_t dca, Double_t *coo):TObject(),
fTrack1(t1),
fTrack2(t2),
fDCA(dca) {
  // Standard constructor

  for(Int_t i=0;i<3;i++)fCross[i]=coo[i];
}

//______________________________________________________________________
AliITSTracklPairs::~AliITSTracklPairs(){
  // destructor
}

//______________________________________________________________________
Double_t AliITSTracklPairs::GetDistance(const AliITSTracklPairs& pair) const {
  // Get distance between the crossing point of pair and this
  Double_t point1[3];
  pair.GetCrossCoord(point1);
  Double_t dist = 0.;
  for(Int_t i=0; i<3; i++)dist+=(point1[i]-fCross[i])*(point1[i]-fCross[i]);
  dist=TMath::Sqrt(dist);
  return dist;
}
