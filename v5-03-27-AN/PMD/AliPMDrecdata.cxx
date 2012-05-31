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
//  Date   : May 27 2009                               //
//                                                     //
//                                                     //
//-----------------------------------------------------//
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDrecdata.h"
#include <stdio.h>

ClassImp(AliPMDrecdata)

AliPMDrecdata::AliPMDrecdata():
  fDet(0),
  fSMN(0),
  fTrackNo(0),
  fTrackPid(0)
{
  // Default constructor
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = 0.;
    }
}
// ------------------------------------------------------------------------- //
AliPMDrecdata::AliPMDrecdata(Int_t idet, Int_t ismn, Int_t trno,
			     Int_t trpid, Float_t *clusdata):
  fDet(idet),
  fSMN(ismn),
  fTrackNo(trno),
  fTrackPid(trpid)
{
  // Constructor
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
// ------------------------------------------------------------------------- //
AliPMDrecdata::AliPMDrecdata(AliPMDrecdata *pmdrecpoint):
  fDet(0),
  fSMN(0),
  fTrackNo(0),
  fTrackPid(0)
{
  *this = *pmdrecpoint;
}

// ------------------------------------------------------------------------- //
AliPMDrecdata::AliPMDrecdata(const AliPMDrecdata &pmdrecpoint):
  TObject(pmdrecpoint),
  fDet(pmdrecpoint.fDet),
  fSMN(pmdrecpoint.fSMN),
  fTrackNo(pmdrecpoint.fTrackNo),
  fTrackPid(pmdrecpoint.fTrackPid)
{
  //Copy Constructor 
  for(Int_t i=0; i<6; i++)
    {
      fClusData[i] = pmdrecpoint.fClusData[i];
    }
}
// ------------------------------------------------------------------------- //
AliPMDrecdata & AliPMDrecdata::operator=(const AliPMDrecdata &pmdrecpoint)
{
  // Assignment operator 
  if(this != &pmdrecpoint)
    {
      fDet      = pmdrecpoint.fDet;
      fSMN      = pmdrecpoint.fSMN;
      fTrackNo  = pmdrecpoint.fTrackNo;
      fTrackPid = pmdrecpoint.fTrackPid;
      for(Int_t i=0; i<6; i++)
	{
	  fClusData[i] = pmdrecpoint.fClusData[i];
	}
    }
  return *this;
}
// ------------------------------------------------------------------------- //
AliPMDrecdata::~AliPMDrecdata()
{
  // Default destructor
}
// ------------------------------------------------------------------------- //
Int_t AliPMDrecdata::GetDetector() const
{
  return fDet;
}
// ------------------------------------------------------------------------- //
Int_t AliPMDrecdata::GetSMNumber() const
{
  return fSMN;
}
// ------------------------------------------------------------------------- //
Int_t AliPMDrecdata::GetClusTrackNo() const
{
  return fTrackNo;
}
// ------------------------------------------------------------------------- //
Int_t AliPMDrecdata::GetClusTrackPid() const
{
  return fTrackPid;
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecdata::GetClusX() const
{
  return fClusData[0];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecdata::GetClusY() const
{
  return fClusData[1];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecdata::GetClusADC() const
{
  return fClusData[2];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecdata::GetClusCells() const
{
  return fClusData[3];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecdata::GetClusSigmaX() const
{
  return fClusData[4];
}
// ------------------------------------------------------------------------- //
Float_t AliPMDrecdata::GetClusSigmaY() const
{
  return fClusData[5];
}
// ------------------------------------------------------------------------- //


