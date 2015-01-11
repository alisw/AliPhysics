/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
 * Global binning definition for the charged particle pt analysis
 *
 *   Author: Markus Fasel
 */
#include <iostream>

#include <TObjArray.h>
#include "AliEMCalTriggerBinningComponent.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerBinningDimension)
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerBinningComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningComponent() :
  TObject(),
  fDimensions(NULL)
{
  /*
   * Main constructor
   */
  fDimensions = new TObjArray;
  fDimensions->SetOwner();
}

//______________________________________________________________________________
AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningComponent(const AliEMCalTriggerBinningComponent& ref) :
  TObject(ref),
  fDimensions(NULL)
{
  /*
   * Copy constructor, creating a deep copy
   */
  fDimensions = new TObjArray;
  fDimensions->SetOwner();
  TIter dimIter(ref.fDimensions);
  AliEMCalTriggerBinningDimension *dim(NULL);
  while((dim = dynamic_cast<AliEMCalTriggerBinningDimension *>(dimIter())))
    fDimensions->Add(new AliEMCalTriggerBinningDimension(*dim));
}

//______________________________________________________________________________
AliEMCalTriggerBinningComponent& AliEMCalTriggerBinningComponent::operator=(const AliEMCalTriggerBinningComponent& ref) {
  /*
   * Assignment operator, doing a deep copy
   */
  TObject::operator=(ref);
  if(&ref != this){
    this->~AliEMCalTriggerBinningComponent();

    fDimensions = new TObjArray;
    fDimensions->SetOwner();
    TIter dimIter(ref.fDimensions);
    AliEMCalTriggerBinningDimension *dim(NULL);
    while((dim = dynamic_cast<AliEMCalTriggerBinningDimension *>(dimIter())))
      fDimensions->Add(new AliEMCalTriggerBinningDimension(*dim));
  }
  return *this;
}

//______________________________________________________________________________
AliEMCalTriggerBinningComponent::~AliEMCalTriggerBinningComponent() {
  /*
   * Destructor
   */
  delete fDimensions;
}

//______________________________________________________________________________
AliEMCalTriggerBinningDimension* AliEMCalTriggerBinningComponent::GetBinning(const char* name) const {
  /*
   * Get binning information for a given axis. Return nullpointer if axis is not yet defined
   *
   * @param name: axis name
   * @return: the axis information
   */
  return dynamic_cast<AliEMCalTriggerBinningDimension *>(fDimensions->FindObject(name));
}

//______________________________________________________________________________
void AliEMCalTriggerBinningComponent::SetBinning(const char* dimname, int nbins, double* binning) {
  /*
   * Set binning for dimension. If not yet existing, create it
   *
   * @param dimname: axis name
   * @param nbins: Number of bins
   * @param binning: array of bin limits (size nbins+1)
   */
  AliEMCalTriggerBinningDimension *dim = GetBinning(dimname);
  if(dim) dim->Set(nbins, binning);
  else {
    dim = new AliEMCalTriggerBinningDimension(dimname, nbins, binning);
    fDimensions->Add(dim);
  }
}

//______________________________________________________________________________
void AliEMCalTriggerBinningComponent::SetBinning(const char* dimname, const TArrayD& binning) {
  /*
   * Set binning for dimension. If not yet existing, create it
   *
   * @param dimname: axis name
   * @param nbins: Number of bins
   * @param binning: array of bin limits (size nbins+1)
   */
  AliEMCalTriggerBinningDimension *dim = GetBinning(dimname);
  if(dim) dim->Set(binning);
  else {
    dim = new AliEMCalTriggerBinningDimension(dimname, binning);
    fDimensions->Add(dim);
  }
}

//______________________________________________________________________________
void AliEMCalTriggerBinningDimension::Print(Option_t * /*option*/) const {
  /*
   * Print the bin limits for a given dimension
   */
  std::cout << "Binning for variable " << GetName() << ":" << std::endl;
  std::cout << "================================================" << std::endl;
  for(int ilim = 0; ilim < fBinning.GetSize(); ilim++){
    std::cout << fBinning[ilim];
    if(ilim < fBinning.GetSize() -1) std::cout << ", ";
  }
  std::cout << std::endl;
}


} /* namespace EMCalTriggerPtAnalysis */

