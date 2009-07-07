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
#include "AliLog.h"
#include <stdio.h>

ClassImp(AliPMDcluster)

AliPMDcluster::AliPMDcluster():
  fDet(0),
  fSMN(0)
{
  // Default constructor
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = 0.;
    }
  for (Int_t i = 0; i < 19; i++)
    {
      fClusCellDataX[i] = 0;
      fClusCellDataY[i] = 0;
      fClusCellTrack[i] = -1;
      fClusCellPid[i]   = -1;
      fClusCellAdc[i]   = 0;
    }

}
// --------------------------------------------------------------------- //
AliPMDcluster::AliPMDcluster(Int_t idet, Int_t ismn, Float_t *clusdata,
			     Int_t *celldataX, Int_t *celldataY,
			     Int_t *celltrack, Int_t *cellpid,
			     Float_t *celladc):
  fDet(idet),
  fSMN(ismn)
{
  // Constructor
  for (Int_t i = 0; i < 6; i++)
    {
      fClusData[i] = clusdata[i];
    }
  if(fClusData[3]>19){
    AliError(Form("Too many clusters %d >= 19",fClusData[3]));
    fClusData[3] = 19;
  }
	     
  for (Int_t i = 0; i < 19; i++)
    {
      fClusCellDataX[i] = celldataX[i];
      fClusCellDataY[i] = celldataY[i];
      fClusCellTrack[i] = celltrack[i];
      fClusCellPid[i]   = cellpid[i];
      fClusCellAdc[i]   = celladc[i];
    }

}
// --------------------------------------------------------------------- //
AliPMDcluster::AliPMDcluster(AliPMDcluster *pmdcluster):
  fDet(0),
  fSMN(0)
{
  *this = *pmdcluster;
}
// --------------------------------------------------------------------- //

AliPMDcluster::AliPMDcluster(const AliPMDcluster &pmdcluster):
  TObject(pmdcluster),
  fDet(pmdcluster.fDet),
  fSMN(pmdcluster.fSMN)
{
  //Copy Constructor 
  for(Int_t i=0; i<6; i++)
    {
      this->fClusData[i] = pmdcluster.fClusData[i];
    }
  if(fClusData[3]>19){
    AliError(Form("Too many clusters %d >= 19",fClusData[3]));
    fClusData[3] = 19;
  }
  for(Int_t i=0; i<19; i++)
    {
      this->fClusCellDataX[i] = pmdcluster.fClusCellDataX[i];
      this->fClusCellDataY[i] = pmdcluster.fClusCellDataY[i];
      this->fClusCellTrack[i] = pmdcluster.fClusCellTrack[i];
      this->fClusCellPid[i]   = pmdcluster.fClusCellPid[i];
      this->fClusCellAdc[i]   = pmdcluster.fClusCellAdc[i];
    }
}
// --------------------------------------------------------------------- //

AliPMDcluster & AliPMDcluster::operator=(const AliPMDcluster &pmdcluster)
{
  // Assignment operator 
  if(this != &pmdcluster)
    {
      this->fDet = pmdcluster.fDet;
      this->fSMN = pmdcluster.fSMN;
      for(Int_t i=0; i<6; i++)
	{
	  this->fClusData[i] = pmdcluster.fClusData[i];
	}
      if(fClusData[3]>19){
	AliError(Form("Too many clusters %d >= 19",fClusData[3]));
	fClusData[3] = 19;
      }
      for(Int_t i=0; i<19; i++)
	{
	  this->fClusCellDataX[i] = pmdcluster.fClusCellDataX[i];
	  this->fClusCellDataY[i] = pmdcluster.fClusCellDataY[i];
	  this->fClusCellTrack[i] = pmdcluster.fClusCellTrack[i];
	  this->fClusCellPid[i]   = pmdcluster.fClusCellPid[i];
	  this->fClusCellAdc[i]   = pmdcluster.fClusCellAdc[i];
	}
    }
  return *this;
}
// --------------------------------------------------------------------- //

AliPMDcluster::~AliPMDcluster()
{
  // Destructor
}
// --------------------------------------------------------------------- //

Int_t AliPMDcluster::GetDetector() const
{
  return fDet;
}
// --------------------------------------------------------------------- //
Int_t AliPMDcluster::GetSMN() const
{
  return fSMN;
}
// --------------------------------------------------------------------- //
Float_t AliPMDcluster::GetClusX() const
{
  return fClusData[0];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcluster::GetClusY() const
{
  return fClusData[1];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcluster::GetClusADC() const
{
  return fClusData[2];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcluster::GetClusCells() const
{
  return fClusData[3];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcluster::GetClusSigmaX() const
{
  return fClusData[4];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcluster::GetClusSigmaY() const
{
  return fClusData[5];
}
// --------------------------------------------------------------------- //
Int_t AliPMDcluster::GetClusCellX(Int_t i) const
{
  return fClusCellDataX[i];
}
// --------------------------------------------------------------------- //
Int_t AliPMDcluster::GetClusCellY(Int_t i) const
{
  return fClusCellDataY[i];
}
// --------------------------------------------------------------------- //
Int_t AliPMDcluster::GetClusCellTrack(Int_t i) const
{
  return fClusCellTrack[i];
}
// --------------------------------------------------------------------- //
Int_t AliPMDcluster::GetClusCellPid(Int_t i) const
{
  return fClusCellPid[i];
}
// --------------------------------------------------------------------- //
Float_t AliPMDcluster::GetClusCellAdc(Int_t i) const
{
  return fClusCellAdc[i];
}
// --------------------------------------------------------------------- //
