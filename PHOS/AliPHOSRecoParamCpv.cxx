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

// This class contains the PHOS CPV reconstruction parameters.
// To get default parameters from any place of your code, use:
//   AliPHOSRecoParam* defPar = AliPHOSRecoParamCpv::GetCpvDefaultParameters();
//   Float_t cluth = defPar->GetClusteringThreshold();
//   ...

// --- AliRoot header files ---
#include "AliPHOSRecoParamCpv.h"

ClassImp(AliPHOSRecoParamCpv)

//-----------------------------------------------------------------------------
AliPHOSRecoParamCpv::AliPHOSRecoParamCpv() : AliPHOSRecoParam()
{
  //Default constructor.

  SetClusteringThreshold(0.0);
  SetLocalMaxCut(0.03);
  SetMinE(0.0);
  SetLogWeight(4.0);
}

//-----------------------------------------------------------------------------
AliPHOSRecoParam* AliPHOSRecoParamCpv::GetCpvDefaultParameters()
{
  //Default parameters for the reconstruction in CPV.

  AliPHOSRecoParam* params = new AliPHOSRecoParamCpv();
  return params;
}
