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
//  Date   : March 22 2004                             //
//                                                     //
//  Store cluster information                          //
//  after discrimination                               //
//                                                     //
//-----------------------------------------------------//
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDclupid.h"
#include <stdio.h>

ClassImp(AliPMDclupid)

AliPMDclupid::AliPMDclupid():
  fDet(0),
  fSMN(0),
  fTrNo(0),
  fTrPid(0),
  fMstatus(0)
{
  // Default constructor
  for (Int_t i = 0; i < 7; i++)
    {
      fClusData[i] = 0.;
    }
}
// ------------------------------------------------------------------ //
AliPMDclupid::AliPMDclupid(Int_t idet, Int_t ismn, Int_t trno, Int_t trpid,
			   Int_t mstat, Float_t *clusdata):
  fDet(idet),
  fSMN(ismn),
  fTrNo(trno),
  fTrPid(trpid),
  fMstatus(mstat)
{
  // Constructor
  for (Int_t i = 0; i < 7; i++)
    {
      fClusData[i] = clusdata[i];
    }
}
// ------------------------------------------------------------------ //
AliPMDclupid::AliPMDclupid(AliPMDclupid *pmdclupid):
  fDet(0),
  fSMN(0),
  fTrNo(0),
  fTrPid(0),
  fMstatus(0)
{
  *this = *pmdclupid;
}

// ------------------------------------------------------------------ //
AliPMDclupid::AliPMDclupid(const AliPMDclupid &pmdclupid):
  TObject(pmdclupid),
  fDet(pmdclupid.fDet),
  fSMN(pmdclupid.fSMN),
  fTrNo(pmdclupid.fTrNo),
  fTrPid(pmdclupid.fTrPid),
  fMstatus(pmdclupid.fMstatus)
{
  //Copy Constructor 
  for(Int_t i=0; i<7; i++)
    {
      fClusData[i] = pmdclupid.fClusData[i];
    }
}
// ------------------------------------------------------------------ //
AliPMDclupid & AliPMDclupid::operator=(const AliPMDclupid &pmdclupid)
{
  // Assignment operator 
  if(this != &pmdclupid)
    {
      fDet     = pmdclupid.fDet;
      fSMN     = pmdclupid.fSMN;
      fTrNo    = pmdclupid.fTrNo;
      fTrPid   = pmdclupid.fTrPid;
      fMstatus = pmdclupid.fMstatus;
      for(Int_t i=0; i<7; i++)
	{
	  fClusData[i] = pmdclupid.fClusData[i];
	}
    }
  return *this;
}
// ------------------------------------------------------------------ //
AliPMDclupid::~AliPMDclupid()
{
  // Destructor
}
// ------------------------------------------------------------------ //
Int_t AliPMDclupid::GetDetector() const
{
  return fDet;
}
// ------------------------------------------------------------------ //
Int_t AliPMDclupid::GetSMN() const
{
  return fSMN;
}
// ------------------------------------------------------------------ //
Int_t AliPMDclupid::GetClusTrackNo() const
{
  return fTrNo;
}
// ------------------------------------------------------------------ //
Int_t AliPMDclupid::GetClusTrackPid() const
{
  return fTrPid;
}
// ------------------------------------------------------------------ //
Int_t AliPMDclupid::GetClusMatching() const
{
  return fMstatus;
}
// ------------------------------------------------------------------ //
Float_t AliPMDclupid::GetClusX() const
{
  return fClusData[0];
}
// ------------------------------------------------------------------ //
Float_t AliPMDclupid::GetClusY() const
{
  return fClusData[1];
}
// ------------------------------------------------------------------ //
Float_t AliPMDclupid::GetClusADC() const
{
  return fClusData[2];
}
// ------------------------------------------------------------------ //
Float_t AliPMDclupid::GetClusCells() const
{
  return fClusData[3];
}
// ------------------------------------------------------------------ //
Float_t AliPMDclupid::GetClusSigmaX() const
{
  return fClusData[4];
}
// ------------------------------------------------------------------ //
Float_t AliPMDclupid::GetClusSigmaY() const
{
  return fClusData[5];
}
// ------------------------------------------------------------------ //
Float_t AliPMDclupid::GetClusPID() const
{
  return fClusData[6];
}


