// AliEmcalCorrectionCellCombineCollections
//

#include <AliAODEvent.h>
#include <AliVCaloCells.h>
#include <AliAODCaloCells.h>
#include <AliESDCaloCells.h>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliEmcalContainerUtils.h"

#include "AliEmcalCorrectionCellCombineCollections.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellCombineCollections);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellCombineCollections> AliEmcalCorrectionCellCombineCollections::reg("AliEmcalCorrectionCellCombineCollections");

/**
 * Standard constructor for the correction components.
 */
AliEmcalCorrectionCellCombineCollections::AliEmcalCorrectionCellCombineCollections() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellCombineCollections"),
  fExternalCellsBranchName("emcalCells"),
  fCreatedCellsBranchName("emcalCellsCombined"),
  fInitializedCombinedCells(false),
  fCombinedCells(0)
{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
}

/**
 * Destructor.
 */
AliEmcalCorrectionCellCombineCollections::~AliEmcalCorrectionCellCombineCollections()
{
  // Destructor
}

/**
 * Initialize all needed variables from the YAML configuration.
 *
 * Note that "usedefault" is only applied to the externalCellsBranchName because
 * the "usedefault" name for the combinedCellsBranchName is not well defined.
 * However, just using the default in the default YAML configuration should
 * serve more or less the same purpose.
 *
 */
Bool_t AliEmcalCorrectionCellCombineCollections::Initialize()
{
  // Initialization
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Initialize();

  GetProperty("externalCellsBranchName", fExternalCellsBranchName);
  GetProperty("combinedCellsBranchName", fCreatedCellsBranchName);

  // Check for "usedefault"
  if (fExternalCellsBranchName == "usedefault") {
    fExternalCellsBranchName = AliEmcalContainerUtils::DetermineUseDefaultName(AliEmcalContainerUtils::kCaloCells, fEsdMode);
  }

  // NOTE: "usedefault" doesn't make sense for the combined cell branch name. The name needs to be distinct!
  //       Thus, we only do it for fExternalCellsBranchName.
  
  return kTRUE;
}

/**
 * Create user output objects.
 */
void AliEmcalCorrectionCellCombineCollections::UserCreateOutputObjects()
{
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
}

/**
 * The combined cells are setup during ExecOnce() to atempt to make
 * the cell object creation as efficient as possible (although it is
 * difficult to do compared to AliEmcalContainer derived classes).
 */
void AliEmcalCorrectionCellCombineCollections::ExecOnce()
{
  SetupCombinedCells();
}

/**
 * Setup the combined cells object and add it to the input event.
 *
 * NOTE: It could have been added to the embedded event, but then classes outside
 * of AliEmcalCorrection{Task,Component} have some difficulty in using them. It is
 * better to put them in the input event, although it is important to ensure that
 * the requested name is not already in use.
 */
void AliEmcalCorrectionCellCombineCollections::SetupCombinedCells()
{
  if (fEsdMode) {
    fCombinedCells = new AliESDCaloCells(fCreatedCellsBranchName.c_str(), fCreatedCellsBranchName.c_str(), AliVCaloCells::kEMCALCell);
  }
  else {
    fCombinedCells = new AliAODCaloCells(fCreatedCellsBranchName.c_str(), fCreatedCellsBranchName.c_str(), AliVCaloCells::kEMCALCell);
  }

  // Add it to the input event
  // While the CorrectionTask can handle cells in the external event, it is not well handled by
  // other classes, so it should stay in the input event to ensure it is easily available.
  AddObjectToEvent(fCombinedCells, fEventManager.InputEvent());

  AliDebugStream(2) << "Added combined calo cells \"" << fCombinedCells->GetName() << "\" with " << fCombinedCells->GetNumberOfCells() << " cells to the input event" << std::endl;

  fInitializedCombinedCells = true;
}

/**
 * Run each event to fill the combined cells from input and external cells.
 * Note that the combined cells object should have already been created.
 */
Bool_t AliEmcalCorrectionCellCombineCollections::Run()
{
  // Run
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Run();

  CreateCombinedCells();
  
  return kTRUE;
}

/**
 * Create cells combined from the input event and the external (embedded) event.
 * This is only necessary for cells because they basic object is a cell "container" rather than the actual
 * objects, which can then be accessed through EMCal containers (like clusters or tracks).
 */
void AliEmcalCorrectionCellCombineCollections::CreateCombinedCells()
{
  // Create cells.
  // In general, the memory allocation process for cells is not very efficient...
  if (fInitializedCombinedCells == false) {
    SetupCombinedCells();
  }

  // Get the external event
  const AliAnalysisTaskEmcalEmbeddingHelper* embedding = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (!embedding) {
    AliFatal("Could not retrieve external event!");
  }

  // Get the cells to copy
  AliVCaloCells * externalCells = dynamic_cast<AliVCaloCells *>(embedding->GetExternalEvent()->FindListObject(fExternalCellsBranchName.c_str()));
  if (!externalCells) {
    AliFatal(TString::Format("Could not retrieve cells \"%s\" from external event!", fExternalCellsBranchName.c_str()));
  }

  // Delete any previous container
  fCombinedCells->DeleteContainer();
  // Create a new container
  fCombinedCells->CreateContainer(fCaloCells->GetNumberOfCells() + externalCells->GetNumberOfCells());
  // Add internal and external cells to the combined cells
  AddCellsToCombinedCellObject(fCaloCells);
  AddCellsToCombinedCellObject(externalCells);
}

/**
 * Takes the input cells and adds them to the combine cells. To do so,
 * the each individual cell must be duplicated.
 * 
 * @param[in] inputCells Cells to be added to the combined cells
 */
void AliEmcalCorrectionCellCombineCollections::AddCellsToCombinedCellObject(AliVCaloCells * inputCells)
{
  // Cell properties
  Short_t cellNumber;
  Double_t ampltidue, time, eFrac;
  Int_t mcLabel;
  Bool_t cellHighGain;
  Bool_t getCellResult = kFALSE;

  AliDebugStream(3) << "Adding caloCells \"" << inputCells->GetName() << "\" with " << inputCells->GetNumberOfCells() << " cells to the combined cells" << std::endl;

  // Loop over the input cells and add them to the combined cells
  for (unsigned int i = 0; i < inputCells->GetNumberOfCells(); i++)
  {
    getCellResult = inputCells->GetCell(i, cellNumber, ampltidue, time, mcLabel, eFrac);
    if (!getCellResult) {
      AliWarning(TString::Format("Could not get cell %i from cell collection %s", i, inputCells->GetName()));
    }
    // Get high gain attribute in addition to cell
    cellHighGain = inputCells->GetCellHighGain(i);

    // Set the properties in the combined cell
    fCombinedCells->SetCell(i, cellNumber, ampltidue, time, mcLabel, eFrac, cellHighGain);
  }
}

/**
 * Add object to event. Adapted from AliAnalysisTaskEmcal.
 * @param[in] obj Object to be added
 * @param[in] event Event to which the object is added
 * @param[in] attempt If true don't handle error
 */
void AliEmcalCorrectionCellCombineCollections::AddObjectToEvent(TObject *obj, AliVEvent * event, Bool_t attempt)
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

