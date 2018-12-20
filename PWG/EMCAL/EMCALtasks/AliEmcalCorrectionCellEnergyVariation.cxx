// AliEmcalCorrectionCellEnergyVariation
//

#include "AliEmcalCorrectionCellEnergyVariation.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellEnergyVariation);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellEnergyVariation> AliEmcalCorrectionCellEnergyVariation::reg("AliEmcalCorrectionCellEnergyVariation");

/**
 * Default constructor
 */
AliEmcalCorrectionCellEnergyVariation::AliEmcalCorrectionCellEnergyVariation() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellEnergyVariation"),
  fEnergyScaleShift(0.)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellEnergyVariation::~AliEmcalCorrectionCellEnergyVariation()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellEnergyVariation::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  GetProperty("energyScaleShift", fEnergyScaleShift);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellEnergyVariation::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();
}

/**
 * Called before the first event to initialize the correction.
 */
void AliEmcalCorrectionCellEnergyVariation::ExecOnce()
{
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellEnergyVariation::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  Short_t  absId  =-1;
  Double_t ecell = 0;
  Double_t tcell = 0;
  Double_t efrac = 0;
  Int_t  mclabel = -1;
  
  // Loop over all EMCal cells
  for (Int_t iCell = 0; iCell < fCaloCells->GetNumberOfCells(); iCell++){
    
    // Get cell
    Bool_t getCellResult = fCaloCells->GetCell(iCell, absId, ecell, tcell, mclabel, efrac);
    if (!getCellResult) {
      AliWarning(TString::Format("Could not get cell %i from cell collection %s", iCell, fCaloCells->GetName()));
    }
    
    // Get high gain attribute in addition to cell
    // NOTE: GetCellHighGain() uses the cell position, not cell index, and thus should _NOT_ be used!
    Bool_t cellHighGain = fCaloCells->GetHighGain(iCell);
    
    // Scale cell energy
    if (TMath::Abs(fEnergyScaleShift) > 1e-6) {
      
      ecell *= (1 + fEnergyScaleShift);
      
      if (ecell > 0.) {
        fCaloCells->SetCell(iCell, absId, ecell, tcell, mclabel, efrac, cellHighGain);
      }
      
    }
    
  }

  return kTRUE;
}
