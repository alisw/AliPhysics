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
#include "AliPMDisocell.h"

ClassImp(AliPMDisocell)

AliPMDisocell::AliPMDisocell():
  fDet(-1),
  fSmn(-1),
  fRow(-1),
  fCol(-1),
  fAdc(0.)
{
  // Default constructor
}
// --------------------------------------------------------------------- //
AliPMDisocell::AliPMDisocell(Int_t idet, Int_t ismn, Int_t irow, Int_t icol, Float_t cadc):
  fDet(idet),
  fSmn(ismn),
  fRow(irow),
  fCol(icol),
  fAdc(cadc)
{
  // Constructor
}
AliPMDisocell::AliPMDisocell(AliPMDisocell *pmdisocell):
  fDet(-1),
  fSmn(-1),
  fRow(-1),
  fCol(-1),
  fAdc(0.)
{
  *this = *pmdisocell;
}

// --------------------------------------------------------------------- //
// --------------------------------------------------------------------- //
AliPMDisocell::AliPMDisocell(const AliPMDisocell &pmdisocell):
  TObject(pmdisocell),
  fDet(pmdisocell.fDet),
  fSmn(pmdisocell.fSmn),
  fRow(pmdisocell.fRow),
  fCol(pmdisocell.fCol),
  fAdc(pmdisocell.fAdc)
{
  //Copy Constructor 
}
// --------------------------------------------------------------------- //

AliPMDisocell & AliPMDisocell::operator=(const AliPMDisocell &pmdisocell)
{
  // Assignment operator 
  if(this != &pmdisocell)
    {
      fDet = pmdisocell.fDet;
      fSmn = pmdisocell.fSmn;
      fRow = pmdisocell.fRow;
      fCol = pmdisocell.fCol;
      fAdc = pmdisocell.fAdc;

    }
  return *this;
}
// --------------------------------------------------------------------- //

AliPMDisocell::~AliPMDisocell()
{
  // Destructor
}
// --------------------------------------------------------------------- //
Int_t AliPMDisocell::GetDetector() const
{
  return fDet;
}
// --------------------------------------------------------------------- //
Int_t AliPMDisocell::GetSmn() const
{
  return fSmn;
}
// --------------------------------------------------------------------- //
Int_t AliPMDisocell::GetRow() const
{
  return fRow;
}
// --------------------------------------------------------------------- //
Int_t AliPMDisocell::GetCol() const
{
  return fCol;
}
// --------------------------------------------------------------------- //
Float_t AliPMDisocell::GetADC() const
{
  return fAdc;
}
// --------------------------------------------------------------------- //
