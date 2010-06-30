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

/////////////////////////////////////////////////////////////////////
//                                                                 //
//      Class with ZDC reconstruction parameters                   //
//      Origin: Chiara.Oppedisano@to.infn.it                       //
//                                                                 //
/////////////////////////////////////////////////////////////////////


#include "AliZDCRecoParamPbPb.h"
#include "AliZDCRecoParampp.h"
#include "AliZDCRecoParam.h"

ClassImp(AliZDCRecoParam)

//_____________________________________________________________________________
AliZDCRecoParam::AliZDCRecoParam() :
  AliDetectorRecoParam(),
  fBeamEnergy(0)
{
  //
  //Default constructor
}

//_____________________________________________________________________________
AliZDCRecoParam::~AliZDCRecoParam()
{
  // destructor

}

//_____________________________________________________________________________
void AliZDCRecoParam::SetGlauberMCDist(Float_t beamEnergy)
{
  // Implemented in AliZDCRecoParamPbPb
  fBeamEnergy = beamEnergy;
  printf("Setting beam energy = %1.0f GeV\n", fBeamEnergy);
}
