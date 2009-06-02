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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with PMD reconstruction parameters                                  //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliLog.h"

#include "AliPMDRecoParam.h"

ClassImp(AliPMDRecoParam)

//_____________________________________________________________________________
AliPMDRecoParam::AliPMDRecoParam():
  AliDetectorRecoParam()
{
  //
  // constructor
  //

  for (Int_t i = 0; i < 48; i++)
    {
      fNoiseCut[i] = 0.;
    }

  SetNameTitle("PMD","PMD");
}
//_____________________________________________________________________________
AliPMDRecoParam::AliPMDRecoParam(const AliPMDRecoParam &source):
  AliDetectorRecoParam(source)
{ 
  //copy Ctor

  for (Int_t i = 0; i < 48; i++)
    {
      fNoiseCut[i] = source.fNoiseCut[i];
    }

}
//_____________________________________________________________________________
AliPMDRecoParam& AliPMDRecoParam::operator=(const AliPMDRecoParam &source)
{
  //
  // assign. operator
  //

  if (this != &source)
    {
      for (Int_t i = 0; i < 48; i++)
	{
	  fNoiseCut[i] = source.fNoiseCut[i];
	}
    }

  return *this;
}
//_____________________________________________________________________________
AliPMDRecoParam::~AliPMDRecoParam() 
{
  //
  // destructor
  //  
}

//_____________________________________________________________________________
AliPMDRecoParam *AliPMDRecoParam::GetPbPbParam(){
  //
  // set default reconstruction parameters for PbPb.
  //
  AliPMDRecoParam *param = new AliPMDRecoParam();
    
  for (Int_t i = 0; i < 48; i++)
    {
      param->fNoiseCut[i] = 10.;    // dummy cuts
    }

  return param;
}

//_____________________________________________________________________________
AliPMDRecoParam *AliPMDRecoParam::GetPPParam(){
  //
  // set default reconstruction parameters for PP.
  //
  AliPMDRecoParam *param = new AliPMDRecoParam();
  for (Int_t i = 0; i < 48; i++)
    {
      param->fNoiseCut[i] = 10.;    // dummy cuts
    }

  return param;
}

//_____________________________________________________________________________
AliPMDRecoParam *AliPMDRecoParam::GetCosmicParam(){
  //
  // set default reconstruction parameters for cosmic muon run
  //
  AliPMDRecoParam *param = new AliPMDRecoParam();
  for (Int_t i = 0; i < 48; i++)
    {
      param->fNoiseCut[i] = 15.;    // dummy cuts
    }
  
  return param;
}

//_____________________________________________________________________________
void AliPMDRecoParam::PrintParameters() const
{
  //
  // Printing of the used PMD reconstruction parameters
  //
  for (Int_t i = 0; i < 48; i++)
    {
      AliInfo(Form(" Noise cut in every detector : %f", fNoiseCut[i]));
    }

}
