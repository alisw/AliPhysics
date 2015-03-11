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
 * Weight handler for event-dependent reweighting
 *
 *   Author: Markus Fasel
 */
#include <TF1.h>

#include "AliEMCalTriggerWeightHandler.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerWeightHandler)

namespace EMCalTriggerPtAnalysis {

///
/// \brief Constructor
///
AliEMCalTriggerWeightHandler::AliEMCalTriggerWeightHandler() :
  fWeightModel(NULL),
  fUsePtHard(kTRUE)
{

}

///
/// \brief Get weight for event
/// \param event: Input event
/// \return the weight calculated for the event
///
double AliEMCalTriggerWeightHandler::GetEventWeight(const AliMCEvent* const event) const {
  if(!fWeightModel) {
    AliError("Weight model not set - returning 1");
    return 1.;
  }
  double weight = 1.;
  if(fUsePtHard){
    const AliGenPythiaEventHeader *header = dynamic_cast<const AliGenPythiaEventHeader *>(event->GenEventHeader());
    if(header)
      weight = fWeightModel->Eval(header->GetPtHard());
    else
      AliError("Event not a pythia event - returning 1");
  } else {
    // Using cross section
    const AliGenPythiaEventHeader *header = dynamic_cast<const AliGenPythiaEventHeader *>(event->GenEventHeader());
    if(header)
      weight = header->GetXsection()/static_cast<double>(header->Trials());
    else
      AliError("Event not a pythia event - returning 1");
  }
  return weight;
}

} /* namespace EMCalTriggerPtAnalysis */
