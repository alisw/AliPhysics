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
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDrecpoint1.h"
#include <stdio.h>

ClassImp(AliPMDrecpoint1)

AliPMDrecpoint1::AliPMDrecpoint1():
  fDet(0),
  fSMN(0)
{
  // Default constructor
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = 0.;
    }
}
// ------------------------------------------------------------------------- //
AliPMDrecpoint1::AliPMDrecpoint1(Int_t idet, Int_t ismn, Float_t *clusdata)
{
  // Constructor
  fDet = idet;
  fSMN = ismn;
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
// ------------------------------------------------------------------------- //
AliPMDrecpoint1::AliPMDrecpoint1(const AliPMDrecpoint1 &pmdrecpoint):TObject(pmdrecpoint)
{
  //Copy Constructor 
  if(&pmdrecpoint == this) return;
  this->fDet = pmdrecpoint.fDet;
  this->fSMN = pmdrecpoint.fSMN;
  for(Int_t i=0; i<6; i++)
    {
      this->fClusData[i] = pmdrecpoint.fClusData[i];
    }
  return;
}
// ------------------------------------------------------------------------- //
AliPMDrecpoint1 & AliPMDrecpoint1::operator=(const AliPMDrecpoint1 &pmdrecpoint)
{
  // Assignment operator 
  if(&pmdrecpoint == this) return *this;
  this->fDet = pmdrecpoint.fDet;
  this->fSMN = pmdrecpoint.fSMN;
  for(Int_t i=0; i<6; i++)
    {
      this->fClusData[i] = pmdrecpoint.fClusData[i];
    }
  return *this;
}
// ------------------------------------------------------------------------- //
AliPMDrecpoint1::~AliPMDrecpoint1()
{
  // Default destructor
}
// ------------------------------------------------------------------------- //
Int_t AliPMDrecpoint1::GetDetector() const
{
  return fDet;
}
// ------------------------------------------------------------------------- //
Int_t AliPMDrecpoint1::GetSMNumber() const
{
  return fSMN;
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecpoint1::GetClusX() const
{
  return fClusData[0];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecpoint1::GetClusY() const
{
  return fClusData[1];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecpoint1::GetClusADC() const
{
  return fClusData[2];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecpoint1::GetClusCells() const
{
  return fClusData[3];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecpoint1::GetClusSigmaX() const
{
  return fClusData[4];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecpoint1::GetClusSigmaY() const
{
  return fClusData[5];
}
// ------------------------------------------------------------------------- //


