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

#include "AliQAThresholds.h"

ClassImp(AliQAThresholds)

  AliQAThresholds::AliQAThresholds(Int_t detId): TObject(), fThresholds(), fDetectorId(detId)
{
  // constructor

  fThresholds.SetOwner(kTRUE);
}

AliQAThresholds::~AliQAThresholds()
{
  // destructor
}

Int_t AliQAThresholds::GetDetectorId()
{
  return fDetectorId;
}

void AliQAThresholds::SetDetectorId(Int_t i)
{
  fDetectorId = i;
}

void AliQAThresholds::AddThreshold(TParameter<long>* item)
{
  // Add a threshold at the end of the array of thresholds.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.Add(item);
}
void AliQAThresholds::AddThreshold(TParameter<int>* item)
{
  // Add a threshold at the end of the array of thresholds.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.Add(item);
}
void AliQAThresholds::AddThreshold(TParameter<double>* item)
{
  // Add a threshold at the end of the array of thresholds.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.Add(item);
}
void AliQAThresholds::AddThreshold(TParameter<float>* item)
{
  // Add a threshold at the end of the array of thresholds.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.Add(item);
}

void AliQAThresholds::AddThresholdAt(TParameter<int>* item, Int_t index)
{
  // Add a threshold at index 'index' in the array of thresholds.
  // If index is larger than the current size of the array, expand the array.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.AddAtAndExpand(item, index);
}
void AliQAThresholds::AddThresholdAt(TParameter<long>* item, Int_t index)
{
  // Add a threshold at index 'index' in the array of thresholds.
  // If index is larger than the current size of the array, expand the array.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.AddAtAndExpand(item, index);
}
void AliQAThresholds::AddThresholdAt(TParameter<double>* item, Int_t index)
{
  // Add a threshold at index 'index' in the array of thresholds.
  // If index is larger than the current size of the array, expand the array.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.AddAtAndExpand(item, index);
}
void AliQAThresholds::AddThresholdAt(TParameter<float>* item, Int_t index)
{
  // Add a threshold at index 'index' in the array of thresholds.
  // If index is larger than the current size of the array, expand the array.
  // Ownership of the object is transfered to AliQAThresholds.

  fThresholds.AddAtAndExpand(item, index);
}

TObject* AliQAThresholds::GetThreshold(Int_t i)
{
  // Return the object at position i. Returns 0 if i is out of bounds.

  return fThresholds.At(i);
}

Int_t AliQAThresholds::GetSize()
{
  // Return the number of elements in the thresholds array. 
  // Beware that it is not the number of thresholds, as some elements of the array can be null. 

  return fThresholds.GetSize();
}
