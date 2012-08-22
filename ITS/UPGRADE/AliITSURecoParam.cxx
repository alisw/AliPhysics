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

#include "AliITSURecoParam.h"
#include "AliLog.h"


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with ITS reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



ClassImp(AliITSURecoParam)


//_____________________________________________________________________________
AliITSURecoParam::AliITSURecoParam()
{
  // def c-tor
  SetName("ITS");
  SetTitle("ITS");
}

//_____________________________________________________________________________
AliITSURecoParam::~AliITSURecoParam() 
{
  // destructor
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetHighFluxParam() 
{
  // make default reconstruction  parameters for hig  flux env.
  AliITSURecoParam *param = new AliITSURecoParam(); 
  //
  // put here params
  return param;
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetLowFluxParam() 
{
  // make default reconstruction  parameters for low  flux env.
  AliITSURecoParam *param = new AliITSURecoParam();
  // put here params
  return param;
}

//_____________________________________________________________________________
AliITSURecoParam *AliITSURecoParam::GetCosmicTestParam() 
{
  // make default reconstruction  parameters for cosmics
  AliITSURecoParam *param = new AliITSURecoParam();
  // put here params
  return param;
}

