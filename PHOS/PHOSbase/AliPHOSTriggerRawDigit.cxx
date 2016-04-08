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
/* 
   4x4 TRU represented by position and amplitude.
   Author: Boris Polishchuk (Boris.Polishchuk@cern.ch)
*/

// --- AliRoot header files ---
#include "AliPHOSTriggerRawDigit.h"

ClassImp(AliPHOSTriggerRawDigit)

AliPHOSTriggerRawDigit::AliPHOSTriggerRawDigit() : 
AliDigitNew(),fType(0),fMod(-1),fXloc(-1),fZloc(-1),fXIdx(-1),fZIdx(-1),fTRURow(-1),fBranch(-1),fL1Threshold(-1)
{
  //Default constructor.
}

AliPHOSTriggerRawDigit::AliPHOSTriggerRawDigit(Int_t module, Int_t xIdx, Int_t zIdx, Int_t TRURow, Int_t branch, Int_t amp) :
  AliDigitNew(),fType(0),fMod(module),fXloc(-1),fZloc(-1),fXIdx(xIdx),fZIdx(zIdx),fTRURow(TRURow),fBranch(branch),fL1Threshold(-1)
{
  fAmp = amp;
}

AliPHOSTriggerRawDigit::AliPHOSTriggerRawDigit(Int_t module, Int_t x, Int_t z, Int_t L1_Threshold, Int_t amp) :
  AliDigitNew(),fType(1),fMod(module),fXloc(x),fZloc(z),fXIdx(-1),fZIdx(-1),fTRURow(-1),fBranch(-1),fL1Threshold(L1_Threshold)
{
  fAmp = amp;
}

AliPHOSTriggerRawDigit::AliPHOSTriggerRawDigit(const AliPHOSTriggerRawDigit & tdigit) :
  AliDigitNew(tdigit),fType(tdigit.fType),fMod(tdigit.fMod),fXloc(tdigit.fXloc),fZloc(tdigit.fZloc),fXIdx(tdigit.fXIdx),fZIdx(tdigit.fZIdx),
  fTRURow(tdigit.fTRURow),fBranch(tdigit.fBranch),fL1Threshold(tdigit.fL1Threshold)
{
  fAmp = tdigit.fAmp;
}

AliPHOSTriggerRawDigit& AliPHOSTriggerRawDigit::operator=(const AliPHOSTriggerRawDigit & tdigit)
{
  if (&tdigit != this) {
    AliDigitNew::operator=(tdigit);
    fType = tdigit.fType;
    fMod = tdigit.fMod;
    fXloc = tdigit.fXloc;
    fZloc = tdigit.fZloc;
    fXIdx = tdigit.fXIdx;
    fZIdx = tdigit.fZIdx;
    fTRURow = tdigit.fTRURow;
    fBranch = tdigit.fBranch;
    fL1Threshold = tdigit.fL1Threshold;
  }

  return *this;
}

void AliPHOSTriggerRawDigit::GetModXZ(Int_t& mod, Int_t& modX, Int_t& modZ)
{
  //Return 4x4 region bottom left cell position.

  if(fType==1) {
    mod = fMod;
    modX = fXloc;
    modZ = fZloc;
    return;
  }
  
  Int_t kN2x2XPrTRURow = 8;
  Int_t kN2x2ZPrBranch = 14;
  
  modX = (fXIdx + fTRURow * kN2x2XPrTRURow)*2;
  modZ = (fZIdx + fBranch * kN2x2ZPrBranch)*2;
  mod = fMod;
  
}
