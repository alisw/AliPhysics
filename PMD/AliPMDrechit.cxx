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
//  Date   : February 26 2006                          //
//                                                     //
//  Store cellhits associated to a cluster             //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDrechit.h"

ClassImp(AliPMDrechit)

AliPMDrechit::AliPMDrechit():
  fXcell(0),
  fYcell(0),
  fTrcell(0),
  fPidcell(0),
  fAdccell(0)
{
  // Standard constructor
}

AliPMDrechit::AliPMDrechit(Int_t cellx, Int_t celly, Int_t celltr,
			   Int_t cellpid, Float_t celladc):
  fXcell(cellx),
  fYcell(celly),
  fTrcell(celltr),
  fPidcell(cellpid),
  fAdccell(celladc)

{
  // Constructor
}
AliPMDrechit::AliPMDrechit(AliPMDrechit *pmdrechit):
  fXcell(0),
  fYcell(0),
  fTrcell(0),
  fPidcell(0),
  fAdccell(0)
{
  *this = *pmdrechit;
}

AliPMDrechit::AliPMDrechit(const AliPMDrechit& source):
  TObject(source),
  fXcell(source.fXcell),
  fYcell(source.fYcell),
  fTrcell(source.fTrcell),
  fPidcell(source.fPidcell),
  fAdccell(source.fAdccell)

{
  //Copy Constructor 

}

AliPMDrechit& AliPMDrechit::operator=(const AliPMDrechit& source)
{
  //Copy Constructor 
  if(this != &source)
    {
      fXcell = source.fXcell;
      fYcell = source.fYcell;
      fTrcell = source.fTrcell;
      fPidcell = source.fPidcell;
      fAdccell = source.fPidcell;
    }
  return *this;
}

AliPMDrechit::~AliPMDrechit()
{
  // Default destructor
}
Int_t AliPMDrechit::GetCellX() const
{
  return fXcell;
}
Int_t AliPMDrechit::GetCellY() const
{
  return fYcell;
}
Int_t AliPMDrechit::GetCellTrack() const
{
  return fTrcell;
}
Int_t AliPMDrechit::GetCellPid() const
{
  return fPidcell;
}
Float_t AliPMDrechit::GetCellAdc() const
{
  return fAdccell;
}
