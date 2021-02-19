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

#include <TCustomBinning.h>
#include <TLinearBinning.h>
#include <TVariableBinning.h>
#include <TObjArray.h>
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerBinningComponent)
ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData)

using namespace PWGJE::EMCALJetTasks;

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningComponent() :
  TObject(),
  fDimensions(nullptr)
{
  fDimensions = new TObjArray;
  fDimensions->SetOwner();
}

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningComponent(const AliEMCalTriggerBinningComponent& ref) :
  TObject(ref),
  fDimensions(nullptr)
{
  fDimensions = new TObjArray;
  fDimensions->SetOwner();
  for(auto dimIter : *(ref.fDimensions))
    fDimensions->Add(new AliEMCalTriggerBinningData(*(static_cast<AliEMCalTriggerBinningData *>(dimIter))));
}

AliEMCalTriggerBinningComponent& AliEMCalTriggerBinningComponent::operator=(const AliEMCalTriggerBinningComponent& ref) {
  TObject::operator=(ref);
  if(&ref != this){
    this->~AliEMCalTriggerBinningComponent();

    fDimensions = new TObjArray;
    fDimensions->SetOwner();
    for(auto dimIter : *(ref.fDimensions))
      fDimensions->Add(new AliEMCalTriggerBinningData(*(static_cast<AliEMCalTriggerBinningData *>(dimIter))));
  }
  return *this;
}

AliEMCalTriggerBinningComponent::~AliEMCalTriggerBinningComponent() {
  delete fDimensions;
}

TBinning* AliEMCalTriggerBinningComponent::GetBinning(const char* name) const {
  AliEMCalTriggerBinningData *data = dynamic_cast<AliEMCalTriggerBinningData *>(fDimensions->FindObject(name));
  if(data) return data->GetBinning();
  return nullptr;
}

void AliEMCalTriggerBinningComponent::SetBinning(const char* dimname, int nbins, const double* binning) {
  SetBinning(dimname, new TVariableBinning(nbins, binning));
}

void AliEMCalTriggerBinningComponent::SetBinning(const char* dimname, const TArrayD& binning) {
  SetBinning(dimname, new TVariableBinning(binning));
}


void AliEMCalTriggerBinningComponent::SetLinearBinning(const char *dimname, int nbins, double min, double max){
  SetBinning(dimname, new TLinearBinning(nbins, min, max));
}

void AliEMCalTriggerBinningComponent::SetBinning(const char *dimname, TBinning *binning){
  AliEMCalTriggerBinningData *dim = FindBinning(dimname);
  if(dim) dim->SetBinning(binning);
  else fDimensions->Add(new AliEMCalTriggerBinningData(dimname, binning));
}

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData *AliEMCalTriggerBinningComponent::FindBinning(const char *dimname) const {
  return dynamic_cast<AliEMCalTriggerBinningData *>(fDimensions->FindObject(dimname));
}

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData::AliEMCalTriggerBinningData():
  TNamed(),
  fBinning(nullptr)
{

}

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData::AliEMCalTriggerBinningData(const char *name, TBinning *binning):
  TNamed(name, ""),
  fBinning(binning)
{

}

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData::AliEMCalTriggerBinningData(const AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData &data):
  TNamed(data),
  fBinning(nullptr)
{
  TBinning *refbinning = data.GetBinning();
  if(refbinning) fBinning = refbinning->MakeCopy();
}

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData &AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData::operator=(
  const AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData &data ) {
  TNamed::operator=(data);
  if(this != &data){
    if(fBinning) delete fBinning;
    fBinning = data.fBinning;
  }
  return *this;
}

AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData::~AliEMCalTriggerBinningData(){
  if(fBinning) delete fBinning;
}

void AliEMCalTriggerBinningComponent::AliEMCalTriggerBinningData::SetBinning(TBinning *binning){
  if(fBinning) delete fBinning;
  fBinning = binning;
}
