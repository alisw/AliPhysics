// AliEmcalCorrectionCellCloneContainer
//

#include <sstream>

#include <AliAODEvent.h>
#include <AliVCaloCells.h>
#include <AliAODCaloCells.h>
#include <AliESDCaloCells.h>

#include "AliEmcalContainerUtils.h"

#include "AliEmcalCorrectionCellCloneContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellCloneContainer);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellCloneContainer> AliEmcalCorrectionCellCloneContainer::reg("AliEmcalCorrectionCellCloneContainer");

/**
 * Standard constructor for the correction components.
 */
AliEmcalCorrectionCellCloneContainer::AliEmcalCorrectionCellCloneContainer() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellCloneContainer"),
  fClonedCellsBranchName("emcalCellsCloned"),
  fInitializedClonedCells(false),
  fClonedCells(0)
{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
}

/**
 * Destructor.
 */
AliEmcalCorrectionCellCloneContainer::~AliEmcalCorrectionCellCloneContainer()
{
  // Destructor
}

/**
 * Initialize all needed variables from the %YAML configuration.
 *
 * Note that "usedefault" is only applied to the inputCellsBranchName 
 *
 */
Bool_t AliEmcalCorrectionCellCloneContainer::Initialize()
{
  // Initialization
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Initialize();

  GetProperty("clonedCellsBranchName" , fClonedCellsBranchName);
  
  // Check for "usedefault"
  if ( fClonedCellsBranchName == AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCaloCells, fEsdMode) || 
       fClonedCellsBranchName == "usedefault" ) { 
    AliFatal("Cloned cell branch name cannot be the same as default one!\n"); 
  }
  
  return kTRUE;
}

/**
 * Create user output objects.
 */
void AliEmcalCorrectionCellCloneContainer::UserCreateOutputObjects()
{
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
}

/**
 * The cloned cells are setup during ExecOnce() to atempt to make
 * the cell object creation as efficient as possible.
 */
void AliEmcalCorrectionCellCloneContainer::ExecOnce()
{
  SetupClonedCells();
}

/**
 * Setup the cloned cells object and add it to the input event.
 *
 */
void AliEmcalCorrectionCellCloneContainer::SetupClonedCells()
{
  if (fEsdMode) {
    fClonedCells = new AliESDCaloCells(fClonedCellsBranchName.c_str(), fClonedCellsBranchName.c_str(), AliVCaloCells::kEMCALCell);
  }
  else {
    fClonedCells = new AliAODCaloCells(fClonedCellsBranchName.c_str(), fClonedCellsBranchName.c_str(), AliVCaloCells::kEMCALCell);
  }

  // Add it to the input event
  AddObjectToEvent(fClonedCells, fEventManager.InputEvent());

  AliDebugStream(2) << "Added combined calo cells \"" << fClonedCells->GetName() << "\" with " << fClonedCells->GetNumberOfCells() << " cells to the input event" << std::endl;

  fInitializedClonedCells = true;
}

/**
 * Run each event to fill the combined cells from input and external cells.
 * Note that the combined cells object should have already been created.
 */
Bool_t AliEmcalCorrectionCellCloneContainer::Run()
{
  // Run
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Run();

  CloneCells();

  return kTRUE;
}

/**
 * Copy cells from the input event.
 */
void AliEmcalCorrectionCellCloneContainer::CloneCells()
{
  // Create cells.
  // In general, the memory allocation process for cells is not very efficient...
  if (fInitializedClonedCells == false) {
    SetupClonedCells();
  }

  // Delete any previous container
  fClonedCells->DeleteContainer();
  // Create a new container
  fClonedCells->CreateContainer(fCaloCells->GetNumberOfCells());
  // Add internal and external cells to the combined cells
  AddCellsToClonedObject();
}

/**
 * Takes the input cells and copies them to the cloned cells. 
 */
void AliEmcalCorrectionCellCloneContainer::AddCellsToClonedObject()
{
  // Cell properties
  Short_t cellNumber;
  Double_t amplitude, time, eFrac;
  Int_t mcLabel;
  Bool_t cellHighGain;
  Bool_t getCellResult = kFALSE;

  AliDebugStream(3) << "Copy caloCells \"" << fCaloCells->GetName() << "\" of type \"" << fCaloCells->GetType() << "\" with " << fCaloCells->GetNumberOfCells() << std::endl;

  // Loop over the input cells and add them to the combined cells
  for (unsigned int i = 0; i < static_cast<unsigned int>(fCaloCells->GetNumberOfCells()); i++)
  {
    getCellResult = fCaloCells->GetCell(i, cellNumber, amplitude, time, mcLabel, eFrac);
    if (!getCellResult) {
      AliWarning(TString::Format("Could not get cell %i from cell collection %s", i, fCaloCells->GetName()));
    }
    // Get high gain attribute in addition to cell
    // NOTE: GetCellHighGain() uses the cell position, not cell index, and thus should _NOT_ be used!
    cellHighGain = fCaloCells->GetHighGain(i);

    // Set the properties in the combined cell
    fClonedCells->SetCell(i, cellNumber, amplitude, time, mcLabel, eFrac, cellHighGain);
  }
}

/**
 * Add object to event. Adapted from AliAnalysisTaskEmcal.
 * @param[in] obj Object to be added
 * @param[in] event Event to which the object is added
 * @param[in] attempt If true don't handle error
 */
void AliEmcalCorrectionCellCloneContainer::AddObjectToEvent(TObject *obj, AliVEvent * event, Bool_t attempt)
{
  if (!(event->FindListObject(obj->GetName()))) {
    event->AddObject(obj);
  }
  else {
    if (!attempt) {
      AliFatal(Form("Container with name %s already present. Aborting", obj->GetName()));
    }
  }
}

