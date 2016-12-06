
#include "AliEmcalCopyCollection.h"

#include <string>

#include <TString.h>
#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliVCaloCells.h"
#include "AliESDCaloCells.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCopyCollection);
/// \endcond

AliEmcalCopyCollection::AliEmcalCopyCollection():
  AliAnalysisTaskSE("AliEmcalCopyCollection"),
  fInputObjectType(AliEmcalContainerUtils::kNoDefinedInputObject),
  fCollectionToCopyName(""),
  fNewCollectionName(""),
  fIsEmbedding(false),
  fEventInitialized(false),
  fIsEsd(false),
  fEvent(0)
{
}

AliEmcalCopyCollection::AliEmcalCopyCollection(std::string name, AliEmcalContainerUtils::InputObject_t inputObjectType, std::string collectionToCopyName, std::string newCollectionName, bool isEmbedding):
  AliAnalysisTaskSE(name.c_str()),
  fInputObjectType(inputObjectType),
  fCollectionToCopyName(collectionToCopyName),
  fNewCollectionName(newCollectionName),
  fIsEmbedding(isEmbedding),
  fEventInitialized(false),
  fIsEsd(false),
  fEvent(0)
{
}

/**
 * Used to determine the input type (ESD or AOD).
 */
void AliEmcalCopyCollection::UserCreateOutputObjects()
{
  // Determine event type from AnalysisManager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) {
    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (evhand) {
      if (evhand->InheritsFrom("AliESDInputHandler")) {
        fIsEsd = true;
      }
      else {
        fIsEsd = false;
      }
    }
    else {
      AliError("Event handler not found!");
    }
  }
  else {
    AliError("Analysis manager not found!");
  }
}

/**
 * Creates a new collection during the first event and 
 * then copies the collection into it each event.
 */
void AliEmcalCopyCollection::UserExec(Option_t * option)
{
  // Initialize the event if not intialized
  if (!fEventInitialized)
    CreateNewObjectBranch();

  // Only continue if we are initialized successfully
  if (!fEventInitialized)
    return;

  // Copy cells and cluster to new branches if desired
  CopyBranchToNewObject();
}

/**
 * Steers creation of a new collection in the event. It supports the "usedefault"
 * pattern to determine the input collection branch name.
 */
void AliEmcalCopyCollection::CreateNewObjectBranch()
{
  fEvent = AliEmcalContainerUtils::GetEvent(InputEvent(), fIsEmbedding);
  if (!fEvent) {
    AliError("Could not retrieve event! Returning!");
    return;
  }

  // Check that the input or output objects are not empty
  if (fCollectionToCopyName == "") {
    AliFatal("Collection to copy name is empty! Please set a valid name!");
  }

  if (fNewCollectionName == "") {
    AliFatal("New collection name is empty! Please set a valid name!");
  }

  // Determine use default name if applicable
  if (fCollectionToCopyName == "usedefault") {
    fCollectionToCopyName = AliEmcalContainerUtils::DetermineUseDefaultName(fInputObjectType, fIsEsd);
  }
  if (fNewCollectionName == "usedefault") {
    // This name is almost certainly going to conflict, but checked to be certain
    fNewCollectionName = AliEmcalContainerUtils::DetermineUseDefaultName(fInputObjectType, fIsEsd);
  }

  if (fInputObjectType != AliEmcalContainerUtils::kCaloCells &&
    fInputObjectType != AliEmcalContainerUtils::kCluster &&
    fInputObjectType != AliEmcalContainerUtils::kTrack) {
    AliFatal(TString::Format("Unrecognized input object type %d", fInputObjectType));
  }

  // Create new branch
  NewBranch();

  fEventInitialized = true;
}

/**
 * Create a new branch for the new collection. It checks to ensure that the object does not
 * already exist in the input event before creating it.
 *
 * Note that cells and clusters/tracks are handled differently due to the different types
 * used to store each branch (Ali*CaloCells and TClonesArray respectively).
 */
void AliEmcalCopyCollection::NewBranch()
{
  // Check to ensure that we are not trying to create a new branch on top of the old branch
  TObject * existingObject = fEvent->FindListObject(fNewCollectionName.c_str());
  if (existingObject) {
    AliFatal(TString::Format("Attempted to create a new %d branch, \"%s\", with the same name as an existing branch! Check your configuration! Perhaps \"usedefault\" was used incorrectly?", fInputObjectType, fNewCollectionName.c_str()));
  }
  TClonesArray * existingArray = dynamic_cast<TClonesArray *>(existingObject);

  if (fInputObjectType == AliEmcalContainerUtils::kCaloCells)
  {
    // Cells
    // Create new cells branch and add it to the event
    AliVCaloCells * newCells = 0;
    if (fIsEsd) {
      newCells = new AliESDCaloCells(fNewCollectionName.c_str(), fNewCollectionName.c_str(), AliVCaloCells::kEMCALCell);
    }
    else {
      newCells = new AliAODCaloCells(fNewCollectionName.c_str(), fNewCollectionName.c_str(), AliVCaloCells::kEMCALCell);
    }

    // Add to event
    fEvent->AddObject(newCells);
  }
  else if (fInputObjectType == AliEmcalContainerUtils::kCluster || fInputObjectType == AliEmcalContainerUtils::kTrack) {
    // Clusters and tracks

    // Determine relevant values to able to properly create the branch
    // Branch object type name (the type which is going to be created in the TClonesArray)
    std::string newObjectTypeName = AliEmcalContainerUtils::DetermineUseDefaultName(fInputObjectType, fIsEsd, true);

    // Create new branch and add it to the event
    TClonesArray * newArray = new TClonesArray(newObjectTypeName.c_str());
    newArray->SetName(fNewCollectionName.c_str());
    fEvent->AddObject(newArray);
  }
  else {
    // Shouldn't ever get here, but included for completeness
    AliFatal(TString::Format("Unrecognized input object type %d", fInputObjectType));
  }

}

/**
 * Copies from the collection branch to the new collection branch.
 *
 * Cells, clusters, and tracks are all handled individually.
 */
void AliEmcalCopyCollection::CopyBranchToNewObject()
{
  // Cells
  if (fInputObjectType == AliEmcalContainerUtils::kCaloCells)
  {
    AliVCaloCells * cellsToCopy = dynamic_cast<AliVCaloCells *>(fEvent->FindListObject(fCollectionToCopyName.c_str()));
    AliVCaloCells * newCells = dynamic_cast<AliVCaloCells *>(fEvent->FindListObject(fNewCollectionName.c_str()));
    AliDebug(3, Form("Number of old cells: %d", cellsToCopy->GetNumberOfCells()));
    // Copies from cellsToCopy to newCells!
    cellsToCopy->Copy(*newCells);
    // The name was changed to the currentCells name, but we want it to be newCells, so we restore it
    newCells->SetName(fNewCollectionName.c_str());
    newCells->SetTitle(fNewCollectionName.c_str());

    AliDebug(3, Form("Number of old cells: %d \tNumber of new cells: %d", cellsToCopy->GetNumberOfCells(), newCells->GetNumberOfCells() ));
    if (cellsToCopy->GetNumberOfCells() != newCells->GetNumberOfCells()) {
      AliFatal(Form("Number of old cell: %d != number of new cells %d", cellsToCopy->GetNumberOfCells(), newCells->GetNumberOfCells()));
    }
  }

  // Clusters
  if (fInputObjectType == AliEmcalContainerUtils::kCluster)
  {
    TClonesArray * newClusters = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fNewCollectionName.c_str()));
    TClonesArray * currentClusters = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fCollectionToCopyName.c_str()));
    AliDebug(3, Form("before copy:\t currentClusters->GetEntries(): %d \t newClusters->GetEntries(): %d", currentClusters->GetEntries(), newClusters->GetEntries()));
    CopyClusters(currentClusters, newClusters);
    AliDebug(3, Form("after  copy:\t currentClusters->GetEntries(): %d \t newClusters->GetEntries(): %d", currentClusters->GetEntries(), newClusters->GetEntries()));
  }

  // Tracks
  if (fInputObjectType == AliEmcalContainerUtils::kTrack)
  {
    TClonesArray * newTracks = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fNewCollectionName.c_str()));
    TClonesArray * currentTracks = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject(fCollectionToCopyName.c_str()));
    AliDebug(3, Form("before copy:\t currentTracks->GetEntries(): %d \t newTracks->GetEntries(): %d", currentTracks->GetEntries(), newTracks->GetEntries()));
    for (Int_t i = 0; i < currentTracks->GetEntriesFast(); i++)
    {
      if (fIsEsd)
      {
        AliESDtrack *currentTrack = dynamic_cast<AliESDtrack *>(currentTracks->At(i));
        // Calls copy constructor to create a new track at position i in newTracks
        AliESDtrack *newTrack = dynamic_cast<AliESDtrack *>(new ((*newTracks)[i]) AliESDtrack(*currentTrack));

        // Reset properties of the track to fix TRefArray errors later.
        // Particular combination is from AODTrackFilterTask
        // Resetting the TProcessID of the track is not sufficient!
        newTrack->SetUniqueID(0);
        newTrack->ResetBit(TObject::kHasUUID);
        newTrack->ResetBit(TObject::kIsReferenced);
      }
      else
      {
        AliAODTrack *currentTrack = dynamic_cast<AliAODTrack *>(currentTracks->At(i));
        // Calls copy constructor to create a new track at position i in newTracks
        AliAODTrack *newTrack = dynamic_cast<AliAODTrack *>(new ((*newTracks)[i]) AliAODTrack(*currentTrack));

        // Reset properties of the track to fix TRefArray errors later.
        // Particular combination is from AODTrackFilterTask
        // Resetting the TProcessID of the track is not sufficient!
        newTrack->SetUniqueID(0);
        newTrack->ResetBit(TObject::kHasUUID);
        newTrack->ResetBit(TObject::kIsReferenced);
      }
    }
    AliDebug(3, Form("after  copy:\t currentTracks->GetEntries(): %d \t newTracks->GetEntries(): %d", currentTracks->GetEntries(), newTracks->GetEntries()));
  }
}

/**
 * Copies clusters from an one TClonesArray to another.
 *
 * This function is from AliAnalysisTaskEMCALClusterizeFast, but modified to copy ALL
 * clusters so that PHOS clusters can be used in the EMCal framework.
 */
void AliEmcalCopyCollection::CopyClusters(TClonesArray *orig, TClonesArray *dest)
{
  const Int_t Ncls = orig->GetEntries();

  for(Int_t i=0; i < Ncls; ++i) {
    AliVCluster *oc = static_cast<AliVCluster*>(orig->At(i));

    if (!oc)
      continue;

    // We likely want to include PHOS cells as well, so the copy we make should avoid this check
    //if (!oc->IsEMCAL())
    //    continue;

    AliVCluster *dc = static_cast<AliVCluster*>(dest->New(i));
    //dc->SetType(AliVCluster::kEMCALClusterv1);
    dc->SetType(oc->GetType());
    dc->SetE(oc->E());
    Float_t pos[3] = {0};
    oc->GetPosition(pos);
    dc->SetPosition(pos);
    dc->SetNCells(oc->GetNCells());
    dc->SetCellsAbsId(oc->GetCellsAbsId());
    dc->SetCellsAmplitudeFraction(oc->GetCellsAmplitudeFraction());
    dc->SetID(oc->GetID());
    dc->SetDispersion(oc->GetDispersion());
    dc->SetEmcCpvDistance(-1);
    dc->SetChi2(-1);
    dc->SetTOF(oc->GetTOF());     //time-of-flight
    dc->SetNExMax(oc->GetNExMax()); //number of local maxima
    dc->SetM02(oc->GetM02());
    dc->SetM20(oc->GetM20());
    dc->SetDistanceToBadChannel(oc->GetDistanceToBadChannel()); 
    dc->SetMCEnergyFraction(oc->GetMCEnergyFraction());

    //MC
    UInt_t nlabels = oc->GetNLabels();
    Int_t *labels = oc->GetLabels();

    if (nlabels == 0 || !labels)
      continue;

    AliESDCaloCluster *esdClus = dynamic_cast<AliESDCaloCluster*>(dc);
    if (esdClus) {
      TArrayI parents(nlabels, labels);
      esdClus->AddLabels(parents); 
    }
    else {
      AliAODCaloCluster *aodClus = dynamic_cast<AliAODCaloCluster*>(dc);
      if (aodClus) 
        aodClus->SetLabel(labels, nlabels);
    }
  }
}

AliEmcalCopyCollection * AliEmcalCopyCollection::AddTaskEmcalCopyCollection(AliEmcalContainerUtils::InputObject_t inputObjectType,
                          const char * collectionToCopyName,
                          const char * newCollectionName,
                          bool isEmbedding)

{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalCopyCollection", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalCopyCollection", "This task requires an input event handler");
    return 0;
  }

  TString name = "AliEmcalCopyCollection";
  name += TString::Format("_%s_%s", collectionToCopyName, newCollectionName);
  if (isEmbedding) {
    name += "_Embedded";
  }

  AliEmcalCopyCollection* mgrTask = static_cast<AliEmcalCopyCollection *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  // Create the task that manages copying collections
  AliEmcalCopyCollection* copyCollection = new AliEmcalCopyCollection(name.Data(),
                                    inputObjectType,
                                    collectionToCopyName,
                                    newCollectionName,
                                    isEmbedding);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(copyCollection);

  // Create containers for input
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(copyCollection, 0, cInput);

  return copyCollection;
}
