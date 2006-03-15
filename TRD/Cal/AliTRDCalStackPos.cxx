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
//                                                                           //
//  TRD calibration class for position parameters of the stacks //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDCalStackPos.h"

ClassImp(AliTRDCalStackPos)

//_____________________________________________________________________________
AliTRDCalStackPos::AliTRDCalStackPos():TNamed()
{
  //
  // AliTRDCalStackPos default constructor
  //

  for (Int_t idet = 0; idet < kNstacks; idet++) {
    for (Int_t i = 0; i < 3; ++i) {
      fStackPos[idet][i] = 0;
      fStackRot[idet][i] = 0;
    }
  }
}

//_____________________________________________________________________________
AliTRDCalStackPos::AliTRDCalStackPos(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDCalStackPos constructor
  //

  for (Int_t idet = 0; idet < kNstacks; idet++) {
    for (Int_t i = 0; i < 3; ++i) {
      fStackPos[idet][i] = 0;
      fStackRot[idet][i] = 0;
    }
  }
}

