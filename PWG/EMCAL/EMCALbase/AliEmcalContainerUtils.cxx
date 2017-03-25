//
// Containers EMCal container utility functions
//

#include "AliEmcalContainerUtils.h"

#include <iostream>

#include <AliVEvent.h>
#include <AliLog.h>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

/**
 * Given a container type, it returns the proper default branch name based on the "usedefault" pattern.
 * This is useful to properly handle creating input objects such as AliEmcalContainer derived objects.
 * If returnObjectType is true, it returns the "default" (unlikely to change) object type instead of the
 * branch name. This is useful to properly determine the type of an object for a TClonesArray.
 *
 * This function can also be very useful in places such as an AddTask(). Using it can significantly reduce
 * code duplication!
 *
 * @param[in] objType Type of the input object
 * @param[in] esdMode True if running with an ESD
 * @param[in] returnObjectType Returns the "default" type of the object rather than the branch name
 *
 * @return The name corresponding to the request branch name or object type.
 */
std::string AliEmcalContainerUtils::DetermineUseDefaultName(InputObject_t objType, bool esdMode, bool returnObjectType)
{
  std::string returnValue = "";
  if (objType == kCluster) {
    if (esdMode == true) {
      if (returnObjectType == true) {
        returnValue = "AliESDCaloCluster";
      }
      else {
        returnValue = "CaloClusters";
      }
    }
    else {
      if (returnObjectType == true) {
        returnValue = "AliAODCaloCluster";
      }
      else {
        returnValue = "caloClusters";
      }
    }
  }
  else if (objType == kTrack) {
    if (esdMode == true) {
      if (returnObjectType == true) {
        returnValue = "AliESDtrack";
      }
      else {
        returnValue = "Tracks";
      }
    }
    else {
      if (returnObjectType == true) {
        returnValue = "AliAODTrack";
      }
      else {
        returnValue = "tracks";
      }
    }
  }
  else if (objType == kCaloCells) {
    if (esdMode == true) {
      if (returnObjectType == true) {
        returnValue = "AliESDCaloCells";
      }
      else {
        returnValue = "EMCALCells";
      }
    }
    else {
      if (returnObjectType == true) {
        returnValue = "AliAODCaloCells";
      }
      else {
        returnValue = "emcalCells";
      }
    }
  }
  else {
    // Default to empty if we are given an unrecognized type with "usedefault"
    returnValue = "";
    AliWarningGeneralStream("AliEmcalContainerUtils") << "Unrecognized combination of inputs! Passed values of input object: " <<  objType << ", esdMode: " << std::boolalpha << esdMode << ", returnObjectType " << returnObjectType << std::endl;
  }

  return returnValue;
}

/**
 * Get the proper event based on whether embedding is enabled or not. Useful when determining from which event
 * an input object should be retrieved. It could either be the current input event or an embedded event.
 *
 * @param[in] inputEvent The input event of the analysis. Will be returned if nothing else is requested. Usually just InputEvent().
 * @param[in] isEmbedding True if the event from embedding should be used.
 *
 * @return The input event to be used
 */
AliVEvent * AliEmcalContainerUtils::GetEvent(AliVEvent * inputEvent, bool isEmbedding)
{
  AliVEvent * event = 0;
  if (isEmbedding) {
    const AliAnalysisTaskEmcalEmbeddingHelper* embedding = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (!embedding) return 0;

    event = embedding->GetExternalEvent();
  }
  else {
    event = inputEvent;
  }

  return event;
}

