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
#include "AliEMCalTriggerBinningFactory.h"

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
 * Copy constructor, creating a deep copy.
 * \param[in] ref Reference for the copy
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
 * Assignment operator, doing a deep copy.
 * \param[in] ref Reference for the assignment
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
 * \param[in] name axis name
 * \return the axis information
 */
AliEMCalTriggerBinningDimension* AliEMCalTriggerBinningComponent::GetBinning(const char* name) const {
  return dynamic_cast<AliEMCalTriggerBinningDimension *>(fDimensions->FindObject(name));
}

/**
 * Set binning for dimension. If not yet existing, create it
 * \param[in] dimname: axis name
 * \param[in] nbins: Number of bins
 * \param[in] binning: array of bin limits (size nbins+1)
 */
void AliEMCalTriggerBinningComponent::SetBinning(const char* dimname, int nbins, const double* binning) {
  AliEMCalTriggerBinningDimension *dim = GetBinning(dimname);
  if(dim) dim->Set(nbins, binning);
  else {
    dim = new AliEMCalTriggerBinningDimension(dimname, nbins, binning);
    fDimensions->Add(dim);
  }
}

/**
 * Set binning for dimension. If not yet existing, create it.
 * \param[in] dimname axis name
 * \param[in] binning array of bin limits (size nbins+1)
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
 * Set a linear binning for dimension. If not yet existing, create it.
 * \param[in] dimname axis name
 * \param[in] nbins Number of bins
 * \param[in] min Minimum of the range (= lowest bin limit)
 * \param[in] max Maximum of the range (= highest bin limit)
 */
void AliEMCalTriggerBinningComponent::SetLinearBinning(const char *dimname, int nbins, double min, double max){
  TArrayD binning;
  AliEMCalTriggerBinningFactory::CreateLinearBinning(binning, nbins, min, max);
  SetBinning(dimname, binning);
}

/**
 * Print the bin limits for a given dimension. Used in the operator<< of
 * AliEMCalTriggerBinningDimension.
 * \param[in] stream Stream to print the information on
 */
void AliEMCalTriggerBinningDimension::PrintStream(std::ostream &stream) const{
  stream << "Binning for variable " << GetName() << ":\n";
  stream << "================================================\n";
  for(int ilim = 0; ilim < fBinning.GetSize(); ilim++){
    stream << fBinning[ilim];
    if(ilim < fBinning.GetSize() -1) stream << ", ";
  }
}

/**
 * Print the bin limits for a given dimension.
 */
void AliEMCalTriggerBinningDimension::Print(Option_t * /*option*/) const {
  std::cout << *this << std::endl;
}

/**
 * Output stream operator for the binning dimension.
 * \param[in] stream Stream to print the information on
 * \param[in] dim Object to be put on the stream
 * \return Stream after pringing
 */
std::ostream &operator<<(std::ostream &stream, const AliEMCalTriggerBinningDimension &dim){
  dim.PrintStream(stream);
  return stream;
}


} /* namespace EMCalTriggerPtAnalysis */

