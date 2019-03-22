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

// string to enum map for use with the %YAML config
const std::map <std::string, AliVEvent::EOfflineTriggerTypes> AliEmcalContainerUtils::fgkPhysicsSelectionMap = {
  {"kMB", AliVEvent::kMB},
  {"kINT1", AliVEvent::kINT1},
  {"kINT7", AliVEvent::kINT7},
  {"kMUON", AliVEvent::kMUON},
  {"kHighMult", AliVEvent::kHighMult},
  {"kHighMultSPD", AliVEvent::kHighMultSPD},
  {"kEMC1", AliVEvent::kEMC1},
  {"kCINT5", AliVEvent::kCINT5},
  {"kINT5", AliVEvent::kINT5},
  {"kCMUS5", AliVEvent::kCMUS5},
  {"kMUSPB", AliVEvent::kMUSPB},
  {"kINT7inMUON", AliVEvent::kINT7inMUON},
  {"kMuonSingleHighPt7", AliVEvent::kMuonSingleHighPt7},
  {"kMUSH7", AliVEvent::kMUSH7},
  {"kMUSHPB", AliVEvent::kMUSHPB},
  {"kMuonLikeLowPt7", AliVEvent::kMuonLikeLowPt7},
  {"kMUL7", AliVEvent::kMUL7},
  {"kMuonLikePB", AliVEvent::kMuonLikePB},
  {"kMuonUnlikeLowPt7", AliVEvent::kMuonUnlikeLowPt7},
  {"kMUU7", AliVEvent::kMUU7},
  {"kMuonUnlikePB", AliVEvent::kMuonUnlikePB},
  {"kEMC7", AliVEvent::kEMC7},
  {"kEMC8", AliVEvent::kEMC8},
  {"kMUS7", AliVEvent::kMUS7},
  {"kMuonSingleLowPt7", AliVEvent::kMuonSingleLowPt7},
  {"kPHI1", AliVEvent::kPHI1},
  {"kPHI7", AliVEvent::kPHI7},
  {"kPHI8", AliVEvent::kPHI8},
  {"kPHOSPb", AliVEvent::kPHOSPb},
  {"kEMCEJE", AliVEvent::kEMCEJE},
  {"kEMCEGA", AliVEvent::kEMCEGA},
  {"kHighMultV0", AliVEvent::kHighMultV0},
  {"kCentral", AliVEvent::kCentral},
  {"kSemiCentral", AliVEvent::kSemiCentral},
  {"kDG", AliVEvent::kDG},
  {"kDG5", AliVEvent::kDG5},
  {"kZED", AliVEvent::kZED},
  {"kSPI7", AliVEvent::kSPI7},
  {"kSPI", AliVEvent::kSPI},
  {"kINT8", AliVEvent::kINT8},
  {"kMuonSingleLowPt8", AliVEvent::kMuonSingleLowPt8},
  {"kMuonSingleHighPt8", AliVEvent::kMuonSingleHighPt8},
  {"kMuonLikeLowPt8", AliVEvent::kMuonLikeLowPt8},
  {"kMuonUnlikeLowPt8", AliVEvent::kMuonUnlikeLowPt8},
  {"kMuonUnlikeLowPt0", AliVEvent::kMuonUnlikeLowPt0},
  {"kUserDefined", AliVEvent::kUserDefined},
  {"kTRD", AliVEvent::kTRD},
  {"kMuonCalo", AliVEvent::kMuonCalo},
  {"kFastOnly", AliVEvent::kFastOnly},
  {"kAny", AliVEvent::kAny},
  {"kAnyINT", AliVEvent::kAnyINT}
};

/**
 * Determines the physics selection that is retrieved from a YAML configuration. Note that the result is an OR of
 * all of the individual selections in the input.
 *
 * @return The desired trigger selection. Note that a `UInt_t` is what is used for fOfflineTriggerMask, so it's fine to return it here.
 */
UInt_t AliEmcalContainerUtils::DeterminePhysicsSelectionFromYAML(const std::vector<std::string> & selections)
{
  UInt_t physicsSelection = 0;
  for (auto selection : selections) {
    auto sel = fgkPhysicsSelectionMap.find(selection);
    AliDebugGeneralStream("AliEmcalContainerUtils", 3) << "Adding physics selection: " << selection << "\n";
    if (sel != fgkPhysicsSelectionMap.end()) {
      physicsSelection |= sel->second;
    }
    else {
      AliFatalGeneralF("AliEmcalContainerUtils", "Could not find physics selection with key \"%s\"", selection.c_str());
    }
  }
  return physicsSelection;
}

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

