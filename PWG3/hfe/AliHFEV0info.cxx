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
//
// Stores aditional information about the V0 candidates
// author: M.Fasel@gsi.de
//
//
#include "AliHFEV0info.h"

ClassImp(AliHFEV0info)

AliHFEV0info::AliHFEV0info():
  TObject(),
  fTrack(NULL),
  fIDpartnerTrack(-1),
  fIDV0(-1)
{
  //
  // Dummy constructor
  //
}

AliHFEV0info::AliHFEV0info(AliVParticle *track, Int_t idPartnerTrack, Int_t v0id):
  TObject(),
  fTrack(track),
  fIDpartnerTrack(idPartnerTrack),
  fIDV0(v0id)
{
  //
  // Default constructor
  //
}
