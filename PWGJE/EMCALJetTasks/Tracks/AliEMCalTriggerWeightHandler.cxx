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

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerWeightHandler)
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerPtHardWeight)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Constructor
 */
AliEMCalTriggerWeightHandler::AliEMCalTriggerWeightHandler() :
  fWeightModel(NULL),
  fBinWeights(NULL),
  fUseCrossSection(kFALSE)
{

}

/**
 * Destructor, cleanup memory assigned
 */
AliEMCalTriggerWeightHandler::~AliEMCalTriggerWeightHandler(){
  if(fWeightModel) delete fWeightModel;
  if(fBinWeights) delete fBinWeights;
}

/**
 * Set weight for a given pt-hard bin to the list of weights. Creates the
 * container if not yet existing.
 * \param ptmin Min. \f$ p_{t} \f$ of the \f$ p_{t} \f$-hard bin
 * \param ptmax Max. \f$ p_{t} \f$ of the \f$ p_{t} \f$-hard bin
 * \param weight Bin weight
 */
void AliEMCalTriggerWeightHandler::SetWeightForBin(double ptmin, double ptmax, double weight){
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

/**
 * Get weight for event
 * \param event Input event
 * \return the weight calculated for the event
 */
double AliEMCalTriggerWeightHandler::GetEventWeight(const AliMCEvent* const event) const {
  const AliGenPythiaEventHeader *header = dynamic_cast<const AliGenPythiaEventHeader *>(event->GenEventHeader());
  if(!header){
    AliError("Event not a pythia event - returning 1");
    return 1.;
  }
  return GetEventWeight(header);
}

/**
 * Get weight for event using a given pythia event header
 * \param header Pythia Event Header
 * \return the weight calculated for the event
 */
double AliEMCalTriggerWeightHandler::GetEventWeight(const AliGenPythiaEventHeader * const header) const {
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

/**
 * Find weihgt for pt-hard value in the list of weights
 * \param pthard Pt-hard value to find a bin for
 * \return weight for the pthard bin (if found), NULL otherwise
 */
const AliEMCalTriggerPtHardWeight *AliEMCalTriggerWeightHandler::FindWeight(Double_t pthard) const{
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


} /* namespace EMCalTriggerPtAnalysis */
