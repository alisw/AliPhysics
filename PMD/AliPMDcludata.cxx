/***************************************************************************
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
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDcludata.h"
#include <stdio.h>

ClassImp(AliPMDcludata)

AliPMDcludata::AliPMDcludata()
{
  // Default constructor
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = 0.;
    }
}
// --------------------------------------------------------------------- //
AliPMDcludata::AliPMDcludata(Float_t *clusdata)
{
  // Constructor
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
// --------------------------------------------------------------------- //
AliPMDcludata::AliPMDcludata(AliPMDcludata *pmdcludata)
{
  *this = *pmdcludata;
}
// --------------------------------------------------------------------- //

AliPMDcludata::AliPMDcludata(const AliPMDcludata &pmdcludata):
  TObject(pmdcludata)
{
  //Copy Constructor 
  for(Int_t i=0; i<6; i++)
    {
      this->fClusData[i] = pmdcludata.fClusData[i];
    }
}
// --------------------------------------------------------------------- //

AliPMDcludata & AliPMDcludata::operator=(const AliPMDcludata &pmdcludata)
{
  // Assignment operator 
  if(this != &pmdcludata)
    {
      for(Int_t i=0; i<6; i++)
	{
	  this->fClusData[i] = pmdcludata.fClusData[i];
	}
    }
  return *this;
}
// --------------------------------------------------------------------- //

AliPMDcludata::~AliPMDcludata()
{
  // Destructor
}
// --------------------------------------------------------------------- //
Float_t AliPMDcludata::GetClusX() const
{
  return fClusData[0];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcludata::GetClusY() const
{
  return fClusData[1];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcludata::GetClusADC() const
{
  return fClusData[2];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcludata::GetClusCells() const
{
  return fClusData[3];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcludata::GetClusSigmaX() const
{
  return fClusData[4];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcludata::GetClusSigmaY() const
{
  return fClusData[5];
}
// --------------------------------------------------------------------- //
