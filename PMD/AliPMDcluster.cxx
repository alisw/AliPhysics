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
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cluster information                          //
//                                                     //
//-----------------------------------------------------//
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDcluster.h"
#include <stdio.h>

ClassImp(AliPMDcluster)

AliPMDcluster::AliPMDcluster():
  fDet(0),
  fSMN(0)
{
  // Default constructor
  //  fDet = 0;
  //  fSMN = 0;
  for (Int_t i = 0; i < 5; i++)
    {
      fClusData[i] = 0.;
    }
}
AliPMDcluster::AliPMDcluster(Int_t idet, Int_t ismn, Float_t *clusdata)
{
  // Constructor
  fDet = idet;
  fSMN = ismn;
  for (Int_t i = 0; i < 5; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
AliPMDcluster::AliPMDcluster(const AliPMDcluster &pmdcluster):TObject(pmdcluster)
{
  //Copy Constructor 
  if(&pmdcluster == this) return;
  this->fDet = pmdcluster.fDet;
  this->fSMN = pmdcluster.fSMN;
  for(Int_t i=0; i<5; i++)
    {
      this->fClusData[i] = pmdcluster.fClusData[i];
    }
  return;
}

AliPMDcluster & AliPMDcluster::operator=(const AliPMDcluster &pmdcluster)
{
  // Assignment operator 
  if(&pmdcluster == this) return *this;
  this->fDet = pmdcluster.fDet;
  this->fSMN = pmdcluster.fSMN;
  for(Int_t i=0; i<5; i++)
    {
      this->fClusData[i] = pmdcluster.fClusData[i];
    }
  return *this;
}

AliPMDcluster::~AliPMDcluster()
{
  // Destructor
}

Int_t AliPMDcluster::GetDetector() const
{
  return fDet;
}
Int_t AliPMDcluster::GetSMN() const
{
  return fSMN;
}
Float_t AliPMDcluster::GetClusX() const
{
  return fClusData[0];
}
Float_t AliPMDcluster::GetClusY() const
{
  return fClusData[1];
}
Float_t AliPMDcluster::GetClusADC() const
{
  return fClusData[2];
}
Float_t AliPMDcluster::GetClusCells() const
{
  return fClusData[3];
}
Float_t AliPMDcluster::GetClusRadius() const
{
  return fClusData[4];
}
