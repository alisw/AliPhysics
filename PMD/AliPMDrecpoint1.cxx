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

ClassImp(AliPMDrecpoint)

AliPMDrecpoint::AliPMDrecpoint()
{
  for (Int_t i = 0; i < 7; i++)
    {
      fClusData[i] = 0.;
    }
}

AliPMDrecpoint::AliPMDrecpoint(Float_t *clusdata)
{
  for (Int_t i = 0; i < 7; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
AliPMDrecpoint::~AliPMDrecpoint()
{

}

Float_t AliPMDrecpoint::GetDetector() const
{
  return fClusData[0];
}
Float_t AliPMDrecpoint::GetSMNumber() const
{
  return fClusData[1];
}
Float_t AliPMDrecpoint::GetClusX() const
{
  return fClusData[2];
}

Float_t AliPMDrecpoint::GetClusY() const
{
  return fClusData[3];
}

Float_t AliPMDrecpoint::GetClusADC() const
{
  return fClusData[4];
}
Float_t AliPMDrecpoint::GetClusCells() const
{
  return fClusData[5];
}
Float_t AliPMDrecpoint::GetClusRadius() const
{
  return fClusData[6];
}


