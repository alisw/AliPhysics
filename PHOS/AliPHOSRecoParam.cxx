/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

// Base class for the PHOS reconstruction parameters.
// Do not use in the reconstruction; use derivative classes instead.
// Author: Boris Polichtchouk.

// --- AliRoot header files ---
#include "AliPHOSRecoParam.h"

ClassImp(AliPHOSRecoParam)

//-----------------------------------------------------------------------------
AliPHOSRecoParam::AliPHOSRecoParam() : TNamed(),
  fClusteringThreshold(9999),fLocMaxCut(9999),fMinE(9999),fW0(9999),
  fSubtractPedestals(kTRUE),fDecoderVersion("")
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRecoParam::AliPHOSRecoParam(const AliPHOSRecoParam& recoParam):
  TNamed(recoParam),fClusteringThreshold(recoParam.fClusteringThreshold),
  fLocMaxCut(recoParam.fLocMaxCut),fMinE(recoParam.fMinE),fW0(recoParam.fW0),
  fSubtractPedestals(recoParam.fSubtractPedestals),fDecoderVersion(recoParam.fDecoderVersion)
{
  //Copy constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRecoParam& AliPHOSRecoParam::operator = (const AliPHOSRecoParam& recoParam)
{
  //Assignment operator.

  if(this != &recoParam) {
    fClusteringThreshold = recoParam.fClusteringThreshold;
    fLocMaxCut = recoParam.fLocMaxCut;
    fMinE = recoParam.fMinE;
    fW0 = recoParam.fW0;
    fSubtractPedestals = recoParam.fSubtractPedestals;
    fDecoderVersion=recoParam.fDecoderVersion ;
  }

  return *this;
}

