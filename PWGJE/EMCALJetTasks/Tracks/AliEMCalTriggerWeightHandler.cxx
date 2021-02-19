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
#include <TF1.h>
#include <TObjArray.h>

#include "AliEMCalTriggerWeightHandler.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerWeightHandler)
ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerPtHardWeight)

using namespace PWGJE::EMCALJetTasks;

AliEMCalTriggerWeightHandler::AliEMCalTriggerWeightHandler() :
  TObject(),
  fWeightModel(nullptr),
  fBinWeights(nullptr),
  fUseCrossSection(kFALSE)
{
  /*
   * See header file for details
   */
}

AliEMCalTriggerWeightHandler::AliEMCalTriggerWeightHandler(const AliEMCalTriggerWeightHandler &ref):
  TObject(ref),
  fWeightModel(nullptr),
  fBinWeights(nullptr),
  fUseCrossSection(ref.fUseCrossSection)
{
  /*
   * See header file for details
   */
  if(ref.fWeightModel) fWeightModel = new TF1(*ref.fWeightModel);
  if(ref.fBinWeights){
    fBinWeights = new TObjArray;
    for(auto o : *(ref.fBinWeights)){
      fBinWeights->Add(o);
    }
  }
}

AliEMCalTriggerWeightHandler &AliEMCalTriggerWeightHandler::operator =(const AliEMCalTriggerWeightHandler &ref){
  this->~AliEMCalTriggerWeightHandler();
  TObject::operator=(ref);
  if(this != &ref){
    if(ref.fWeightModel) fWeightModel = new TF1(*ref.fWeightModel);
    if(ref.fBinWeights){
      fBinWeights = new TObjArray;
      for(auto o : *(ref.fBinWeights)){
        fBinWeights->Add(o);
      }
    }
    fUseCrossSection = ref.fUseCrossSection;
  }
  return *this;
}

AliEMCalTriggerWeightHandler::~AliEMCalTriggerWeightHandler(){
  /*
   * See header file for details
   */
  if(fWeightModel) delete fWeightModel;
  if(fBinWeights) delete fBinWeights;
}

void AliEMCalTriggerWeightHandler::SetWeightForBin(double ptmin, double ptmax, double weight){
  /*
   * See header file for details
   */
  if(!fBinWeights) fBinWeights = new TObjArray;
  // check whether bin already exists
  AliEMCalTriggerPtHardWeight *exist = NULL, *tmp = NULL;
  for(TIter biniter = TIter(fBinWeights).Begin(); biniter != TIter::End(); ++biniter){
    tmp = static_cast<AliEMCalTriggerPtHardWeight *>(*biniter);
    if(ptmin == tmp->GetPtMin() && ptmax == tmp->GetPtMax()){
      exist = tmp;
      break;
    }
  }
  if(exist) exist->SetWeight(weight);
  else fBinWeights->Add(new AliEMCalTriggerPtHardWeight(ptmin, ptmax, weight));
}

double AliEMCalTriggerWeightHandler::GetEventWeight(const AliMCEvent* const event) const {
  /*
   * See header file for details
   */
  const AliGenPythiaEventHeader *header = dynamic_cast<const AliGenPythiaEventHeader *>(event->GenEventHeader());
  if(!header){
    AliError("Event not a pythia event - returning 1");
    return 1.;
  }
  return GetEventWeight(header);
}

double AliEMCalTriggerWeightHandler::GetEventWeight(const AliGenPythiaEventHeader * const header) const {
  /*
   * See header file for details
   */
  double weight = 1.;
  if(fWeightModel) {
    weight = fWeightModel->Eval(header->GetPtHard());
  } else if(fBinWeights){
    const AliEMCalTriggerPtHardWeight *tmp = FindWeight(header->GetPtHard());
    if(tmp) weight = tmp->GetWeight();
  } else if(fUseCrossSection){
    // Using cross section
    weight = header->GetXsection()/static_cast<double>(header->Trials());
  }
  return weight;
}

const AliEMCalTriggerPtHardWeight *AliEMCalTriggerWeightHandler::FindWeight(Double_t pthard) const{
  /*
   * See header file for details
   */
  const AliEMCalTriggerPtHardWeight *result = NULL, *tmp = NULL;
  for(TIter biniter = TIter(fBinWeights).Begin(); biniter != TIter::End(); ++biniter){
    tmp = static_cast<const AliEMCalTriggerPtHardWeight *>(*biniter);
    if(tmp->IsSelected(pthard)){
      result = tmp;
      break;
    }
  }
  return result;
}
