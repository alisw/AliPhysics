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
:  fNLayers(0)
  ,fTanLorentzAngle(0)
{
  // def c-tor
  SetName("ITS");
  SetTitle("ITS");
}

//_____________________________________________________________________________
AliITSURecoParam::AliITSURecoParam(Int_t nLr)
:  fNLayers(0)
  ,fTanLorentzAngle(0)
{
  // def c-tor
  SetName("ITS");
  SetTitle("ITS");
  SetNLayers(nLr);
}

//_____________________________________________________________________________
AliITSURecoParam::~AliITSURecoParam() 
{
  // destructor
  delete[] fTanLorentzAngle;
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

//_____________________________________________________________________________
void  AliITSURecoParam::SetNLayers(Int_t n)
{
  // set n layers and init all layer dependent arrays
  if (fNLayers>0) AliFatal(Form("Number of layers was already set to %d",fNLayers));
  if (n<1) n = 1; // in case we want to have dummy params
  fNLayers = n;
  fTanLorentzAngle = new Double_t[n];
  //
  for (int i=n;i--;) {
    fTanLorentzAngle[i] = 0;
  }
  //
}

//_____________________________________________________________________________
void  AliITSURecoParam::SetTanLorentzAngle(Int_t lr, Double_t v)
{
  // set Lorentz angle value
  if (lr>=fNLayers) AliFatal(Form("Number of defined layers is %d",fNLayers));
  fTanLorentzAngle[lr] = v;
}
