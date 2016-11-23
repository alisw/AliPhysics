/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliEMCALTriggerPatchADCInfo.h"

ClassImp(AliEMCALTriggerPatchADCInfo)

AliEMCALTriggerPatchADCInfo::AliEMCALTriggerPatchADCInfo():
  TObject(),
  fPatchSize(0),
  fADCValues()
{
}

AliEMCALTriggerPatchADCInfo::AliEMCALTriggerPatchADCInfo(UChar_t patchsize):
  TObject(),
  fPatchSize(patchsize),
  fADCValues()
{
  fADCValues.Allocate(fPatchSize, fPatchSize);
}

AliEMCALTriggerPatchADCInfo::AliEMCALTriggerPatchADCInfo(const AliEMCALTriggerPatchADCInfo &ref):
  TObject(ref),
  fPatchSize(ref.fPatchSize),
  fADCValues()
{
  if(fPatchSize){
    fADCValues.Allocate(fPatchSize, fPatchSize);
    for(UChar_t icol = 0; icol < fPatchSize; icol++){
      for(UChar_t irow = 0; irow < fPatchSize; irow++){
        fADCValues.SetADC(icol, irow, ref.fADCValues(icol, irow));
      }
    }
  }
}

AliEMCALTriggerPatchADCInfo &AliEMCALTriggerPatchADCInfo::operator =(const AliEMCALTriggerPatchADCInfo &ref){
  TObject::operator=(ref);
  if(&ref != this){
    if(ref.fPatchSize){
      if(!fADCValues.IsAllocated()) fADCValues.Allocate(ref.fPatchSize, ref.fPatchSize);
      for(UChar_t icol = 0; icol < ref.fPatchSize; icol++){
        for(UChar_t irow = 0; irow < ref.fPatchSize; irow++){
          fADCValues.SetADC(icol, irow, ref.fADCValues(icol, irow));
        }
      }
    }
    fPatchSize = ref.fPatchSize;
  }
  return *this;
}

void AliEMCALTriggerPatchADCInfo::SetPatchSize(UChar_t patchsize){
  fPatchSize = patchsize;
  fADCValues.Allocate(fPatchSize, fPatchSize);
}

void AliEMCALTriggerPatchADCInfo::SetADC(Int_t adc, UChar_t col, UChar_t row){
  try{
    fADCValues.SetADC(col, row, adc);
  } catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e) {
    // Don't do anything if we are out-of-bounds
  }
}

Int_t AliEMCALTriggerPatchADCInfo::GetADC(UChar_t col, UChar_t row) const {
  try{
    return fADCValues(col, row);
  } catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e){
    return 0;
  }
}

Int_t AliEMCALTriggerPatchADCInfo::GetSumADC() const {
  Int_t adcsum(0);
  for(unsigned int icol = 0; icol < fPatchSize; icol++){
    for(unsigned int irow = 0; irow < fPatchSize; irow++){
      try{
        adcsum += fADCValues(icol, irow);
      } catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e){
      }
    }
  }
  return adcsum;
}

Int_t AliEMCALTriggerPatchADCInfo::GetMaxADC() const {
  Int_t maxadc(0), tmpadc(0);
  for(unsigned int icol = 0; icol < fPatchSize; icol++){
    for(unsigned int irow = 0; irow < fPatchSize; irow++){
      try{
        tmpadc = fADCValues(icol, irow);
      } catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e){
        tmpadc = 0;
      }
      if(tmpadc > maxadc) maxadc = tmpadc;
    }
  }
  return maxadc;
}

Int_t AliEMCALTriggerPatchADCInfo::GetNFastorsContrib() const {
  Int_t ncontrib(0), tmpadc(0);
  for(unsigned int icol = 0; icol < fPatchSize; icol++){
    for(unsigned int irow = 0; irow < fPatchSize; irow++){
      try{
        tmpadc = fADCValues(icol, irow);
      } catch (AliEMCALTriggerDataGrid<Int_t>::OutOfBoundsException &e){
        tmpadc = 0;
      }
      if(tmpadc) ncontrib++;
    }
  }
  return ncontrib;
}
