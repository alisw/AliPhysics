// AliEmcalCorrectionCellCombineCollections
//

#include <sstream>

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
  fVerifyCombinedCells(true),
  fInitializedCombinedCells(false),
  fMergeCells(false),
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
 * Initialize all needed variables from the %YAML configuration.
 *
 * Note that "usedefault" is only applied to the externalCellsBranchName because
 * the "usedefault" name for the combinedCellsBranchName is not well defined.
 * However, just using the default in the default %YAML configuration should
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
  GetProperty("verifyCombinedCells", fVerifyCombinedCells);
  GetProperty("mergeCells", fMergeCells);

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
  AddCellsToCombinedCellObject(fCaloCells, 0);
  AddCellsToCombinedCellObject(externalCells, fCaloCells->GetNumberOfCells());

  
  
  if (fVerifyCombinedCells) {
    VerifyCombinedCells({fCaloCells, externalCells});
  }
}

/**
 * Takes the input cells and adds them to the combine cells. To do so,
 * the each individual cell must be duplicated.
 * 
 * @param[in] inputCells Cells to be added to the combined cells.
 * @param[in] indexOffset Offset into combined cell container where the inputCells should be placed.
 */
void AliEmcalCorrectionCellCombineCollections::AddCellsToCombinedCellObject(AliVCaloCells * inputCells, int indexOffset)
{
  // Cell properties
  Short_t cellNumber;
  Double_t amplitude, time, eFrac;
  Int_t mcLabel;
  Bool_t cellHighGain;
  Bool_t getCellResult = kFALSE;
  Int_t j0 = 0;
  Int_t nMerged = 0;
  
  AliDebugStream(3) << "Adding caloCells \"" << inputCells->GetName() << "\" of type \"" << inputCells->GetType() << "\" with " << inputCells->GetNumberOfCells() << " cells to the combined cells" << std::endl;

  // Loop over the input cells and add them to the combined cells
  for (unsigned int i = 0; i < static_cast<unsigned int>(inputCells->GetNumberOfCells()); i++)
  {
    getCellResult = inputCells->GetCell(i, cellNumber, amplitude, time, mcLabel, eFrac);
    if (!getCellResult) {
      AliWarning(TString::Format("Could not get cell %i from cell collection %s", i, inputCells->GetName()));
    }
    // Get high gain attribute in addition to cell
    // NOTE: GetCellHighGain() uses the cell position, not cell index, and thus should _NOT_ be used!
    cellHighGain = inputCells->GetHighGain(i);
           
    Bool_t newCell = kTRUE;
    // Check if the cell to be added already exists in the list of combined cells
    //
    if ( fMergeCells && indexOffset > 0 )
    {
      // If input data, no signal fraction
      if ( indexOffset == 0 ) eFrac = 0;
      else                    eFrac = amplitude; // not the fraction just signal energy as in other places
      
      // Differenciate MC noise with no label and data cell also with no label.
      if ( mcLabel == -1 )
        mcLabel = -2; 
    
      Short_t cellNumberCo;
      Double_t amplitudeCo, timeCo, eFracCo;
      Int_t mcLabelCo;
      Bool_t cellHighGainCo;
      Bool_t getCellResultCo = kFALSE;
      for (unsigned int j = j0; j < static_cast<unsigned int>(fCombinedCells->GetNumberOfCells()); j++)
      {
        getCellResultCo = fCombinedCells->GetCell(j, cellNumberCo, amplitudeCo, timeCo, mcLabelCo, eFracCo);
        if (!getCellResultCo) {
          AliWarning(TString::Format("Could not get cell %i from cell collection %s", j, fCombinedCells->GetName()));
        }
        cellHighGainCo = fCombinedCells->GetHighGain(j);
        
        // Cell already exists
        if ( cellNumberCo == cellNumber )
        {
          //printf("Merge cell ID %d\n",cellNumberCo);
          
          // Assign time and highGain, the value for the highest energy cell of the two
          if ( amplitude > amplitudeCo )
          {
            timeCo = time;
            cellHighGainCo = cellHighGain;
          }
          
          // Update existing cell
          fCombinedCells->SetCell(j, cellNumber, amplitude+amplitudeCo, timeCo, mcLabel, eFrac, cellHighGainCo);
          
          // Next time, start from index j+1 to avoid too many searches, since lists are ordered by absId
          j0 = j+1;
          
          newCell = kFALSE;
          
          nMerged++;
          
          break;
        } // already existing cell
      } // combined cells loop
    }  // merge cells
    
    // Set the properties in the combined cell
    if ( newCell )
      fCombinedCells->SetCell(i + indexOffset - nMerged, cellNumber, amplitude, time, mcLabel, eFrac, cellHighGain);
  }
  
//  if ( fMergeCells && indexOffset > 0 )
//  {
//    printf("Merged cells %d\n",nMerged);
//    
//    for (unsigned int i = 0; i < static_cast<unsigned int>(fCombinedCells->GetNumberOfCells()); i++)
//    {
//      bool getCellResultCo = fCombinedCells->GetCell(i, cellNumber, amplitude, time, mcLabel, eFrac);
//      if (!getCellResultCo) {
//        printf("Could not get cell %i from cell collection %s", i, fCombinedCells->GetName());
//      }
//      
//      printf("cell %d, id %d, amp %2.3f mc %d, eFrac %2.3f\n",i,cellNumber,amplitude,mcLabel,eFrac);
//    }
//  }
}

/**
 * Explicitly check that the cells in the various input cell collections were successfully combined together
 * to form the combined cells. NOTE: This should not be run during normal usage. It will reduce performance and
 * provide little benefit to the user.
 *
 * @param[in] inputCaloCells vector of calo cells which were combined to create the combined cells.
 */
void AliEmcalCorrectionCellCombineCollections::VerifyCombinedCells(std::vector <AliVCaloCells *> inputCaloCells)
{
  // Input cell arrays properties
  Short_t cellNumber;
  Double_t amplitude, time, eFrac;
  Int_t mcLabel;
  Bool_t cellHighGain;
  Bool_t getCellResult = kFALSE;
  // Combined cell array properties
  Short_t cellNumberCombined;
  Double_t amplitudeCombined, timeCombined, eFracCombined;
  Int_t mcLabelCombined;
  Bool_t cellHighGainCombined;
  Bool_t getCellResultCombined = kFALSE;

  unsigned int combinedCellsIndex = 0;
  for (auto caloCells : inputCaloCells)
  {
    AliDebugStream(2) << "Verifying cells from cell collection " << caloCells->GetName() << "\n";
    for (unsigned int i = 0; i < static_cast<unsigned int>(caloCells->GetNumberOfCells()); i++)
    {
      // Get input calo cell
      getCellResult = caloCells->GetCell(i, cellNumber, amplitude, time, mcLabel, eFrac);
      if (!getCellResult) {
        AliWarning(TString::Format("Could not get cell %i from cell collection %s", i, caloCells->GetName()));
      }
      // Get high gain attribute in addition to cell
      cellHighGain = caloCells->GetHighGain(i);
      // Get combined calo cell
      getCellResultCombined = fCombinedCells->GetCell(combinedCellsIndex, cellNumberCombined, amplitudeCombined, timeCombined, mcLabelCombined, eFracCombined);
      if (!getCellResultCombined) {
        AliWarning(TString::Format("Could not get cell %i from combined cell collection %s", combinedCellsIndex, fCombinedCells->GetName()));
      }
      // Get high gain attribute in addition to cell
      cellHighGainCombined = fCombinedCells->GetHighGain(combinedCellsIndex);

      if (cellNumberCombined != cellNumber ||
        amplitudeCombined != amplitude ||
        timeCombined != time ||
        mcLabelCombined != mcLabel ||
        eFracCombined != eFrac ||
        cellHighGainCombined != cellHighGain)
      {
        std::stringstream errorMessage;
        errorMessage << std::boolalpha;
        errorMessage << "Cell mistmatch at index " << i << " of cell collection " << caloCells->GetName() << ", combined cells " << combinedCellsIndex << " of combined cell collection " << fCombinedCells->GetName() << "\n";
        errorMessage << "Input cell: cell number: " << cellNumber << ", amplitude: " << amplitude << ", time: " << time << ", mcLabel: " << mcLabel << ", eFrac: " << eFrac << ", cellHighGain: " << cellHighGain << "\n";
        errorMessage << "Combined cell: cell number: " << cellNumberCombined << ", amplitude: " << amplitudeCombined << ", time: " << timeCombined << ", mcLabel: " << mcLabelCombined << ", eFrac: " << eFracCombined << ", cellHighGain: " << cellHighGainCombined << "\n";
        AliFatal(errorMessage.str().c_str());
      }

      combinedCellsIndex++;
    }
  }

  AliDebugStream(1) << "Input cells were successfully combined into the combined cells.\n";
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

