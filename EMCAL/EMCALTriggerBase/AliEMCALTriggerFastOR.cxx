/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliEMCALTriggerFastOR.h"
#include "AliEMCALGeometry.h"

AliEMCALTriggerFastOR::AliEMCALTriggerFastOR() :
  fAbsId(0),
  fGlobalCol(0),
  fGlobalRow(0),
  fSM(0),
  fCol(0),
  fRow(0),
  fL0Amp(0),
  fL1Amp(0)
{

}

AliEMCALTriggerFastOR::AliEMCALTriggerFastOR(UInt_t L0amp, UInt_t L1amp, Int_t absId, const AliEMCALGeometry* geom) :
  fAbsId(0),
  fGlobalCol(0),
  fGlobalRow(0),
  fSM(0),
  fCol(0),
  fRow(0),
  fL0Amp(L0amp),
  fL1Amp(L1amp)
{
  Initialize(absId, geom);
}

AliEMCALTriggerFastOR::AliEMCALTriggerFastOR(UInt_t L0amp, UInt_t L1amp, Int_t globalRow, Int_t glocalCol, const AliEMCALGeometry* geom) :
  fAbsId(0),
  fGlobalCol(0),
  fGlobalRow(0),
  fSM(0),
  fCol(0),
  fRow(0),
  fL0Amp(L0amp),
  fL1Amp(L1amp)
{
  Initialize(globalRow, glocalCol, geom);
}


void AliEMCALTriggerFastOR::Initialize(UInt_t L0amp, UInt_t L1amp, Int_t absId, const AliEMCALGeometry* geom)
{
  fL0Amp = L0amp;
  fL1Amp = L1amp;

  Initialize(absId, geom);
}

void AliEMCALTriggerFastOR::Initialize(Int_t absId, const AliEMCALGeometry* geom)
{
  Int_t SM = 0;
  Int_t col = 0;
  Int_t row = 0;
  Int_t gcol = 0;
  Int_t grow = 0;

  geom->GetPositionInEMCALFromAbsFastORIndex(absId, gcol, grow);
  geom->GetPositionInSMFromAbsFastORIndex(absId, SM, col, row);

  fCol = col;
  fRow = row;
  fGlobalCol = gcol;
  fGlobalRow = grow;
  fSM = SM;
  fAbsId = absId;
}

void AliEMCALTriggerFastOR::Initialize(UInt_t L0amp, UInt_t L1amp, Int_t globalRow, Int_t glocalCol, const AliEMCALGeometry* geom)
{
  fL0Amp = L0amp;
  fL1Amp = L1amp;

  Initialize(globalRow, glocalCol, geom);
}

void AliEMCALTriggerFastOR::Initialize(Int_t globalRow, Int_t glocalCol, const AliEMCALGeometry* geom)
{
  Int_t absId = 0;
  Int_t SM = 0;
  Int_t col = 0;
  Int_t row = 0;

  geom->GetAbsFastORIndexFromPositionInEMCAL(glocalCol, globalRow, absId);
  geom->GetPositionInSMFromAbsFastORIndex(absId, SM, col, row);

  fCol = col;
  fRow = row;
  fGlobalCol = glocalCol;
  fGlobalRow = globalRow;
  fSM = SM;
  fAbsId = absId;
}
