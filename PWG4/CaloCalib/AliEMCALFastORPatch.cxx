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

// --- Root ---
#include <TObject.h>
#include <TArrayI.h>
#include <TArrayF.h>

// --- AliRoot ---
#include "AliESDCaloTrigger.h"

#include "AliEMCALFastORPatch.h"

ClassImp(AliEMCALFastORPatch)

//________________________________________________________________________
AliEMCALFastORPatch::AliEMCALFastORPatch()
  : TObject(),
    fAbsId(-1),
    fFastORamplitudes(new TArrayF(0)),
    fFastORcolumns(new TArrayI(0)),
    fFastORrows(new TArrayI(0))
{
  // Constructor
}

//________________________________________________________________________
AliEMCALFastORPatch::AliEMCALFastORPatch(Int_t aid, Int_t size)
  : TObject(),
    fAbsId(aid),
    fFastORamplitudes(new TArrayF(size)),
    fFastORcolumns(new TArrayI(size)),
    fFastORrows(new TArrayI(size))
{
  // Constructor

  fFastORamplitudes->Reset(-1);
  fFastORrows->Reset(-1);
  fFastORcolumns->Reset(-1);
}

//________________________________________________________________________
AliEMCALFastORPatch::~AliEMCALFastORPatch()
{
  // Destructor
  delete fFastORamplitudes;
  delete fFastORcolumns;
  delete fFastORrows;
}

//________________________________________________________________________
Float_t AliEMCALFastORPatch::GetTotalAmplitude() const
{
  // Return total amplitude of the patch
  
  Float_t total = 0;
  for (Int_t i = 0; i < GetNumberOfFastOR(); i++){
    Float_t amp = GetFastORamplitude(i);
    if (amp > 0)
      total += amp;
  }
  return total;
}

//________________________________________________________________________
Float_t AliEMCALFastORPatch::GetFastORamplitude(Int_t i) const
{
  // Return amplitude of FastOR at position i
  
  if (i < 0 || i > GetNumberOfFastOR()) return -1;
  return fFastORamplitudes->At(i);
}

//________________________________________________________________________
Int_t AliEMCALFastORPatch::GetFastORrow(Int_t i) const
{
  // Return row number of FastOR at position i
  
  if (i < 0 || i > GetNumberOfFastOR()) return -1;
  return fFastORrows->At(i);
}

//________________________________________________________________________
Int_t AliEMCALFastORPatch::GetFastORcolumn(Int_t i) const
{
  // Return column number of FastOR at position i
  
  if (i < 0 || i > GetNumberOfFastOR()) return -1;
  return fFastORcolumns->At(i);
}

//________________________________________________________________________
Int_t AliEMCALFastORPatch::GetNumberOfFastOR() const
{
  // Return total number of FastORs
  
  return fFastORamplitudes->GetSize();
}

//________________________________________________________________________
Bool_t AliEMCALFastORPatch::AddFastORat(AliESDCaloTrigger* f, Int_t i)
{
  // Add a FastOR to the patch and return true if done
  
  Float_t amp = 0;
  Int_t gCol = 0, gRow = 0;
  f->GetAmplitude(amp);
  f->GetPosition(gCol, gRow);
  
  return AddFastORat(amp, gCol, gRow, i);
}

//________________________________________________________________________
Bool_t AliEMCALFastORPatch::AddFastORat(Float_t amp, Int_t gCol, Int_t gRow, Int_t i)
{
  // Add a FastOR to the patch and return true if done
  
  if (i < 0) return kFALSE;
  
  if (i > GetNumberOfFastOR()) Expand(i + 1);
  
  fFastORamplitudes->AddAt(amp, i);
  fFastORrows->AddAt(gRow, i);
  fFastORcolumns->AddAt(gCol, i);
  
  return kTRUE;
}

//________________________________________________________________________
void AliEMCALFastORPatch::RemoveFastORat(Int_t i)
{
  // Remove a FastOR from the patch
  
  fFastORamplitudes->AddAt(-1, i);
  fFastORrows->AddAt(-1, i);
  fFastORcolumns->AddAt(-1, i);
}

//________________________________________________________________________
void AliEMCALFastORPatch::Expand(Int_t size)
{  
  // Expand the arrays
  
  fFastORamplitudes->Set(size);
  fFastORcolumns->Set(size);
  fFastORrows->Set(size);
}

//________________________________________________________________________
Bool_t AliEMCALFastORPatch::Contains(Int_t row, Int_t col) const
{
  // Check if a col, row FastOR is present in the arrays
  
  for (Int_t i = 0; i < GetNumberOfFastOR(); i++){
    if (col == fFastORcolumns->At(i) && row == fFastORrows->At(i))
      return kTRUE;
  }
  return kFALSE;
}
