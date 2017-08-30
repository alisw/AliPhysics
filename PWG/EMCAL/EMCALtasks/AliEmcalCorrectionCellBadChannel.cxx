// AliEmcalCorrectionCellBadChannel
//

#include <TObjArray.h>
#include <TFile.h>
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"

#include "AliEmcalCorrectionCellBadChannel.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellBadChannel);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellBadChannel> AliEmcalCorrectionCellBadChannel::reg("AliEmcalCorrectionCellBadChannel");

/**
 * Default constructor
 */
AliEmcalCorrectionCellBadChannel::AliEmcalCorrectionCellBadChannel() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellBadChannel")
  ,fCellEnergyDistBefore(0)
  ,fCellEnergyDistAfter(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellBadChannel::~AliEmcalCorrectionCellBadChannel()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellBadChannel::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  AliWarning("Init EMCAL cell bad channel removal");
  
  GetProperty("createHistos", fCreateHisto);

  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;

  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellBadChannel::UserCreateOutputObjects()
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fCellEnergyDistBefore = new TH1F("hCellEnergyDistBefore","hCellEnergyDistBefore;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistBefore);
    fCellEnergyDistAfter = new TH1F("hCellEnergyDistAfter","hCellEnergyDistAfter;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistAfter);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellBadChannel::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  if (!fEventManager.InputEvent()) {
    AliError("Event ptr = 0, returning");
    return kFALSE;
  }
  
  CheckIfRunChanged();
  
  // CONFIGURE THE RECO UTILS -------------------------------------------------

  fRecoUtils->SwitchOnBadChannelsRemoval();
  
  // START PROCESSING ---------------------------------------------------------
  // Test if cells present
  if (fCaloCells->GetNumberOfCells()<=0)
  {
    AliWarning(Form("Number of EMCAL cells = %d, returning", fCaloCells->GetNumberOfCells()));
    return kFALSE;
  }
  
  // mark the cells not recalibrated
  fRecoUtils->ResetCellsCalibrated();

  if(fCreateHisto)
    FillCellQA(fCellEnergyDistBefore); // "before" QA
  
  // CELL RECALIBRATION -------------------------------------------------------
  // update cell objects
  UpdateCells();
  
  if(fCreateHisto)
    FillCellQA(fCellEnergyDistAfter); // "after" QA

  return kTRUE;
}

/**
 * This function is called if the run changes (it inherits from the base component),
 * to load a new bad channel and fill relevant variables.
 */
Bool_t AliEmcalCorrectionCellBadChannel::CheckIfRunChanged()
{
  Bool_t runChanged = AliEmcalCorrectionComponent::CheckIfRunChanged();
  
  if (runChanged) {
    // init bad channels
    Int_t fInitBC = InitBadChannels();
    if (fInitBC==0) {
      AliError("InitBadChannels returned false, returning");
    }
    if (fInitBC==1) {
      AliWarning("InitBadChannels OK");
    }
    if (fInitBC>1) {
      AliWarning(Form("No external hot channel set: %d - %s", fEventManager.InputEvent()->GetRunNumber(), fFilepass.Data()));
    }
  }
  return runChanged;
}
