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
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store recpoints for ALICE-PMD                      //
//                                                     //
//-----------------------------------------------------//

#include "AliPMDrecpoint1.h"
#include <stdio.h>

ClassImp(AliPMDrecpoint1)

AliPMDrecpoint1::AliPMDrecpoint1()
{
  for (Int_t i = 0; i < 7; i++)
    {
      fClusData[i] = 0.;
    }
}

AliPMDrecpoint1::AliPMDrecpoint1(Float_t *clusdata)
{
  for (Int_t i = 0; i < 7; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
AliPMDrecpoint1::~AliPMDrecpoint1()
{

}

Float_t AliPMDrecpoint1::GetDetector() const
{
  return fClusData[0];
}
Float_t AliPMDrecpoint1::GetSMNumber() const
{
  return fClusData[1];
}
Float_t AliPMDrecpoint1::GetClusX() const
{
  return fClusData[2];
}

Float_t AliPMDrecpoint1::GetClusY() const
{
  return fClusData[3];
}

Float_t AliPMDrecpoint1::GetClusADC() const
{
  return fClusData[4];
}
Float_t AliPMDrecpoint1::GetClusCells() const
{
  return fClusData[5];
}
Float_t AliPMDrecpoint1::GetClusRadius() const
{
  return fClusData[6];
}


