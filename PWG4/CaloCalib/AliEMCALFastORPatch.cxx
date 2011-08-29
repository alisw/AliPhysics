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

/* $Id$ */

#include "AliEMCALFastORPatch.h"

#include <TObject.h>
#include <TObjArray.h>

#include "AliESDCaloTrigger.h"

ClassImp(AliEMCALFastORPatch)

AliEMCALFastORPatch::AliEMCALFastORPatch()
:TObject(),
fFastORamplitudes(0),
fFastORcolumns(0),
fFastORrows(0),
fSize(0),
fAbsId(-1)
{
  // Dummy constructor
}

AliEMCALFastORPatch::AliEMCALFastORPatch(Int_t aid, Int_t size)
:TObject(),
fFastORamplitudes(new Float_t[size]),
fFastORcolumns(new Int_t[size]),
fFastORrows(new Int_t[size]),
fSize(4),
fAbsId(aid)
{
  // Constructor
  for (Int_t i = 0; i < fSize; i++)
  {
    fFastORamplitudes[i] = -1;
    fFastORrows[i] = -1;
    fFastORcolumns[i] = -1;
  }
}

AliEMCALFastORPatch::~AliEMCALFastORPatch()
{
  // Destructor
  delete[] fFastORamplitudes;
  delete[] fFastORcolumns;
  delete[] fFastORrows;
}

Float_t AliEMCALFastORPatch::GetTotalAmplitude()
{
  Float_t total = 0;
  for (Int_t i = 0; i < GetNumberOfFastOR(); i++)
  {
    Float_t amp = GetFastORamplitude(i);
    if (amp > 0)
      total += amp;
  }
  return total;
}

Float_t AliEMCALFastORPatch::GetFastORamplitude(Int_t i)
{
  if (i < 0 || i > fSize) return -1;
  return fFastORamplitudes[i];
}

Int_t AliEMCALFastORPatch::GetFastORrow(Int_t i)
{
  if (i < 0 || i > fSize) return -1;
  return fFastORrows[i];
}

Int_t AliEMCALFastORPatch::GetFastORcolumn(Int_t i)
{
  if (i < 0 || i > fSize) return -1;
  return fFastORcolumns[i];
}

Int_t AliEMCALFastORPatch::GetNumberOfFastOR()
{
  return fSize;
}

Bool_t AliEMCALFastORPatch::AddFastORat(AliESDCaloTrigger* f, Int_t i)
{
  Float_t amp = 0;
  Int_t gCol = 0, gRow = 0;
  f->GetAmplitude(amp);
  f->GetPosition(gCol, gRow);
  
  return AddFastORat(amp, gCol, gRow, i);
}

Bool_t AliEMCALFastORPatch::AddFastORat(Float_t amp, Int_t gCol, Int_t gRow, Int_t i)
{
  if (i < 0) return kFALSE;
  
  if (i > fSize) Expand(i + 1);
  
  fFastORamplitudes[i] = amp;
  fFastORrows[i] = gRow;
  fFastORcolumns[i] = gCol;
  
  return kTRUE;
}

void AliEMCALFastORPatch::RemoveFastORat(Int_t i)
{
  fFastORamplitudes[i] = -1;
  fFastORrows[i] = -1;
  fFastORcolumns[i] = -1;
}

Int_t AliEMCALFastORPatch::GetAbsId()
{
  return fAbsId;
}

void AliEMCALFastORPatch::SetAbsId(Int_t aid)
{
  fAbsId = aid;
}

void AliEMCALFastORPatch::Expand(Int_t size)
{
  Float_t *newFastORamplitudes = new Float_t[size];
  Int_t *newFastORcolumns = new Int_t[size];
  Int_t *newFastORrows = new Int_t[size];
  
  for (Int_t i = 0; i < fSize; i++)
  {
    newFastORamplitudes[i] = fFastORamplitudes[i];
    newFastORcolumns[i] = fFastORcolumns[i];
    newFastORrows[i] = fFastORrows[i];
  }
  
  delete[] fFastORamplitudes;
  delete[] fFastORcolumns;
  delete[] fFastORrows;
  
  fFastORamplitudes = newFastORamplitudes;
  fFastORcolumns = newFastORcolumns;
  fFastORrows = newFastORrows;
  
  fSize = size;
}

Bool_t AliEMCALFastORPatch::Contains(Int_t row, Int_t col)
{
  for (Int_t i = 0; i < fSize; i++)
  {
    if (col == fFastORcolumns[i] && row == fFastORrows[i])
      return kTRUE;
  }
  return kFALSE;
}