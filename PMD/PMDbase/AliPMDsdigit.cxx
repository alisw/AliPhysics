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
//  used to store the info into TreeS                  //
//                                                     //
//-----------------------------------------------------//
#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDsdigit.h"
#include <stdio.h>

ClassImp(AliPMDsdigit)

AliPMDsdigit::AliPMDsdigit():
  fTrNumber(0),
  fTrPid(0),
  fDet(0),
  fSMN(0),
  fRow(0),
  fColumn(0),
  fEdep(0.)
{
  // Default Constructor
}

AliPMDsdigit::AliPMDsdigit(Int_t trnumber, Int_t trpid, Int_t det, Int_t smn,
			   Int_t irow, Int_t icol, Float_t edep):
  fTrNumber(trnumber),
  fTrPid(trpid),
  fDet(det),
  fSMN(smn),
  fRow(irow),
  fColumn(icol),
  fEdep(edep)
{
  // Constructor
}

AliPMDsdigit::AliPMDsdigit(AliPMDsdigit *pmdsdigit):
  fTrNumber(0),
  fTrPid(0),
  fDet(0),
  fSMN(0),
  fRow(0),
  fColumn(0),
  fEdep(0.)
{
  *this = *pmdsdigit;
}

AliPMDsdigit::AliPMDsdigit(const AliPMDsdigit& pmdsdigit):
  TObject(pmdsdigit),
  fTrNumber(pmdsdigit.fTrNumber),
  fTrPid(pmdsdigit.fTrPid),
  fDet(pmdsdigit.fDet),
  fSMN(pmdsdigit.fSMN),
  fRow(pmdsdigit.fRow),
  fColumn(pmdsdigit.fColumn),
  fEdep(pmdsdigit.fEdep)
{
  //Copy Constructor 
}
AliPMDsdigit & AliPMDsdigit::operator=(const AliPMDsdigit& pmdsdigit)
{
  //Assignment operator 
  if(this != &pmdsdigit)
    {
      fTrNumber   = pmdsdigit.fTrNumber;
      fTrPid      = pmdsdigit.fTrPid;
      fDet        = pmdsdigit.fDet;
      fSMN        = pmdsdigit.fSMN;
      fRow        = pmdsdigit.fRow;
      fColumn     = pmdsdigit.fColumn;
      fEdep       = pmdsdigit.fEdep;
    }
  return *this;
}


AliPMDsdigit::~AliPMDsdigit()
{
  // Default Destructor
}
Int_t AliPMDsdigit::GetTrackNumber() const
{
  return fTrNumber;
}
Int_t AliPMDsdigit::GetTrackPid() const
{
  return fTrPid;
}
Int_t AliPMDsdigit::GetDetector() const
{
  return fDet;
}
Int_t AliPMDsdigit::GetSMNumber() const
{
  return fSMN;
}
Int_t AliPMDsdigit::GetRow() const
{
  return fRow;
}
Int_t AliPMDsdigit::GetColumn() const
{
  return fColumn;
}
Float_t AliPMDsdigit::GetCellEdep() const
{
  return fEdep;
}
