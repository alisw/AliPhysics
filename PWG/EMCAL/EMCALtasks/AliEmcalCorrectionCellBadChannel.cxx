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

//________________________________________________________________________
AliEmcalCorrectionCellBadChannel::AliEmcalCorrectionCellBadChannel() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellBadChannel")
  ,fCellEnergyDistBefore(0)
  ,fCellEnergyDistAfter(0)
{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
}

//________________________________________________________________________
AliEmcalCorrectionCellBadChannel::~AliEmcalCorrectionCellBadChannel()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionCellBadChannel::Initialize()
{
  // Initialization
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Initialize();
  // Do base class initializations and if it fails -> bail out
  //AliAnalysisTaskEmcal::ExecOnce();
  //if (!fInitialized) return;
  
  AliWarning("Init EMCAL cell bad channel removal");
  
  GetProperty("createHistos", fCreateHisto);

  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;

  // missalignment function -- TODO: do we need this?
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  return kTRUE;
}

//________________________________________________________________________
void AliEmcalCorrectionCellBadChannel::UserCreateOutputObjects()
{
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fCellEnergyDistBefore = new TH1F("hCellEnergyDistBefore","hCellEnergyDistBefore;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistBefore);
    fCellEnergyDistAfter = new TH1F("hCellEnergyDistAfter","hCellEnergyDistAfter;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistAfter);
  }
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionCellBadChannel::Run()
{
  // Run
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Run();
  
  if (!fEvent) {
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

//________________________________________________________________________
Bool_t AliEmcalCorrectionCellBadChannel::CheckIfRunChanged()
{
  Bool_t runChanged = AliEmcalCorrectionComponent::CheckIfRunChanged();
  
  if (runChanged) {
    // init bad channels
    Int_t fInitBC = InitBadChannels();
    if (fInitBC==0)
    AliError("InitBadChannels returned false, returning");
    if (fInitBC==1)
    AliWarning("InitBadChannels OK");
    if (fInitBC>1)
    AliWarning(Form("No external hot channel set: %d - %s", fEvent->GetRunNumber(), fFilepass.Data()));
  }
  return runChanged;
}
