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
  fYcell(0)
{
  // Standard constructor
}

AliPMDrechit::AliPMDrechit(Int_t cellx, Int_t celly):
  fXcell(cellx),
  fYcell(celly)

{
  // Constructor
}
AliPMDrechit::AliPMDrechit(AliPMDrechit *pmdrechit):
  fXcell(0),
  fYcell(0)
{
  *this = *pmdrechit;
}

AliPMDrechit::AliPMDrechit(const AliPMDrechit& source):
  TObject(source),
  fXcell(source.fXcell),
  fYcell(source.fYcell)
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
