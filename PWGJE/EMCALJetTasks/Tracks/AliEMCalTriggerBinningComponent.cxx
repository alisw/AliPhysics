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
#include <iostream>

#include <TObjArray.h>
#include "AliEMCalTriggerBinningComponent.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerBinningDimension)
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerBinningComponent)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Main constructor
 */
AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningComponent() :
  TObject(),
  fDimensions(NULL)
{
  fDimensions = new TObjArray;
  fDimensions->SetOwner();
}

/**
 * Copy constructor, creating a deep copy
 * \param ref Reference for the copy
 */
AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningComponent(const AliEMCalTriggerBinningComponent& ref) :
  TObject(ref),
  fDimensions(NULL)
{
  fDimensions = new TObjArray;
  fDimensions->SetOwner();
  TIter dimIter(ref.fDimensions);
  AliEMCalTriggerBinningDimension *dim(NULL);
  while((dim = dynamic_cast<AliEMCalTriggerBinningDimension *>(dimIter())))
    fDimensions->Add(new AliEMCalTriggerBinningDimension(*dim));
}

/**
 * Assignment operator, doing a deep copy
 *
 * \param ref Reference for the assignment
 */
AliEMCalTriggerBinningComponent& AliEMCalTriggerBinningComponent::operator=(const AliEMCalTriggerBinningComponent& ref) {
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

/**
 * Destructor
 */
AliEMCalTriggerBinningComponent::~AliEMCalTriggerBinningComponent() {
  delete fDimensions;
}

/**
 * Get binning information for a given axis. Return nullpointer if axis is not yet defined
 *
 * \param name axis name
 * \return the axis information
 */
AliEMCalTriggerBinningDimension* AliEMCalTriggerBinningComponent::GetBinning(const char* name) const {
  return dynamic_cast<AliEMCalTriggerBinningDimension *>(fDimensions->FindObject(name));
}

/**
 * Set binning for dimension. If not yet existing, create it
 *
 * \param dimname: axis name
 * \param nbins: Number of bins
 * \param binning: array of bin limits (size nbins+1)
 */
void AliEMCalTriggerBinningComponent::SetBinning(const char* dimname, int nbins, double* binning) {
  AliEMCalTriggerBinningDimension *dim = GetBinning(dimname);
  if(dim) dim->Set(nbins, binning);
  else {
    dim = new AliEMCalTriggerBinningDimension(dimname, nbins, binning);
    fDimensions->Add(dim);
  }
}

/**
 * Set binning for dimension. If not yet existing, create it
 *
 * \param dimname axis name
 * \param nbins Number of bins
 * \param binning array of bin limits (size nbins+1)
 */
void AliEMCalTriggerBinningComponent::SetBinning(const char* dimname, const TArrayD& binning) {
  AliEMCalTriggerBinningDimension *dim = GetBinning(dimname);
  if(dim) dim->Set(binning);
  else {
    dim = new AliEMCalTriggerBinningDimension(dimname, binning);
    fDimensions->Add(dim);
  }
}

/**
 * Print the bin limits for a given dimension
 */
void AliEMCalTriggerBinningDimension::Print(Option_t * /*option*/) const {
  std::cout << "Binning for variable " << GetName() << ":" << std::endl;
  std::cout << "================================================" << std::endl;
  for(int ilim = 0; ilim < fBinning.GetSize(); ilim++){
    std::cout << fBinning[ilim];
    if(ilim < fBinning.GetSize() -1) std::cout << ", ";
  }
  std::cout << std::endl;
}


} /* namespace EMCalTriggerPtAnalysis */

