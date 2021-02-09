// AliEmcalCorrectionCellEnergyVariation
//

#include <TGrid.h>
#include <TFile.h>

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
  fMinCellE(0.),
  fMaxCellE(150.),
  fEnergyScaleFactorConstant(1.),
  fEnergyScaleFunction(nullptr)
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
  
  GetProperty("minCellE", fMinCellE);
  GetProperty("maxCellE", fMaxCellE);
  GetProperty("energyScaleFactorConstant", fEnergyScaleFactorConstant);
  
  std::string path;
  std::string name;
  GetProperty("energyScaleFunctionPath", path);
  GetProperty("energyScaleFunctionName", name);
  
  // Load TF1 from file, if provided
  if (!path.empty() && !name.empty()) {
    LoadEnergyScaleFunction(path, name);
  }

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
  // If a constant scale factor is supplied, create a TF1 that is constant in E_cell
  if (TMath::Abs(1 - fEnergyScaleFactorConstant) > 1e-6) {
    
    // If a TF1 was already loaded, throw an error
    if (fEnergyScaleFunction) {
      AliError(Form("%s: fEnergyScaleFunction was already loaded! Do not apply multiple scale factors.", GetName()));
    }
    
    fEnergyScaleFunction = new TF1("energyScaleFunction", "[0]", fMinCellE, fMaxCellE);
    fEnergyScaleFunction->SetParameter(0, fEnergyScaleFactorConstant);
  }
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
    
    // Scale cell energy by TF1, if supplied
    if (fEnergyScaleFunction && ecell > fMinCellE && ecell < fMaxCellE) {
      
      ecell *= fEnergyScaleFunction->Eval(ecell);
      
      if (ecell > 0.) {
        fCaloCells->SetCell(iCell, absId, ecell, tcell, mclabel, efrac, cellHighGain);
      }
    }
    
  }

  return kTRUE;
}

/**
 * Load the energy scale function TF1 from a file into the member fEnergyScaleFunction
 * @param path Path to the file containing the TF1
 * @param name Name of the TF1
 */
void AliEmcalCorrectionCellEnergyVariation::LoadEnergyScaleFunction(const std::string & path, const std::string & name)
{
  TString fname(path);
  if (fname.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }
  
  TFile* file = TFile::Open(path.data());
  
  if (!file || file->IsZombie()) {
    AliErrorStream() << "Could not open energy scale function file\n";
    return;
  }
  
  TF1* energyScaleFunction = dynamic_cast<TF1*>(file->Get(name.data()));
  
  if (energyScaleFunction) {
    AliInfoStream() << Form("Cell energy scale function %s loaded from file %s.", name.data(), path.data()) << "\n";
  }
  else {
    AliErrorStream() << Form("Cell energy scale function %s not found in file %s.", name.data(), path.data()) << "\n";
    return;
  }
  
  fEnergyScaleFunction = static_cast<TF1*>(energyScaleFunction->Clone());
  
  file->Close();
  delete file;
}
