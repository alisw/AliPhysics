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

// This class contains the PHOS EMC reconstruction parameters.
// To get default parameters from any place of your code, use:
//   AliPHOSRecoParam* defPar = AliPHOSRecoParamEmc::GetEmcDefaultParameters();
//   Float_t cluth = defPar->GetClusteringThreshold();
//   ...

// --- AliRoot header files ---
#include "AliPHOSRecoParamEmc.h"

ClassImp(AliPHOSRecoParamEmc)

//-----------------------------------------------------------------------------
AliPHOSRecoParamEmc::AliPHOSRecoParamEmc() : AliPHOSRecoParam()
{
  //Default constructor.

  SetClusteringThreshold(0.2);
  SetLocalMaxCut(0.03);
  SetMinE(0.01);
  SetLogWeight(4.5);
}

//-----------------------------------------------------------------------------
AliPHOSRecoParam* AliPHOSRecoParamEmc::GetEmcDefaultParameters()
{
  //Default parameters for the reconstruction in EMC.

  AliPHOSRecoParam* params = new AliPHOSRecoParamEmc();
  return params;
}
