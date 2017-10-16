//
// Containers EMCal container utility functions
//

#include "AliEmcalContainerUtils.h"

#include <iostream>

#include <AliVEvent.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

/**
 * Determines the "usedefault" pattern using the Analysis Manager to determine the datatype automatically.
 * This will often work fine, but it may not always. Note that the use cause for this function is assumed to
 * be solely about returning the "usedefault" collection name and not the object type.
 *
 * @param[in] objType Type of the input object
 *
 * @return The name corresponding to the request branch name.
 */
std::string AliEmcalContainerUtils::DetermineUseDefaultName(InputObject_t objType)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalSample", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalSample", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  return DetermineUseDefaultName(objType, dataType == kESD, false);
}

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
 * an input object should be retrieved. It could either be the current input event or an embedded event. This is
 * the const version.
 *
 * @param[in] inputEvent The input event of the analysis. Will be returned if nothing else is requested. Usually just InputEvent().
 * @param[in] isEmbedding True if the event from embedding should be used.
 *
 * @return The input event to be used
 */
const AliVEvent * AliEmcalContainerUtils::GetEvent(const AliVEvent * inputEvent, bool isEmbedding)
{
  const AliVEvent * event = nullptr;
  if (isEmbedding) {
    const AliAnalysisTaskEmcalEmbeddingHelper* embedding = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
    if (!embedding) return nullptr;

    // Need the const cast as GetExternalEvent() returns a non-const event
    event = const_cast<const AliVEvent *>(embedding->GetExternalEvent());
  }
  else {
    event = inputEvent;
  }

  return event;
}

/**
 * Get the proper event based on whether embedding is enabled or not. Useful when determining from which event
 * an input object should be retrieved. It could either be the current input event or an embedded event. This is
 * the non-const version.
 *
 * @param[in] inputEvent The input event of the analysis. Will be returned if nothing else is requested. Usually just InputEvent().
 * @param[in] isEmbedding True if the event from embedding should be used.
 *
 * @return The input event to be used
 */
AliVEvent * AliEmcalContainerUtils::GetEvent(AliVEvent * inputEvent, bool isEmbedding)
{
  return const_cast<AliVEvent *>(GetEvent(const_cast<const AliVEvent *>(inputEvent), isEmbedding));
}

