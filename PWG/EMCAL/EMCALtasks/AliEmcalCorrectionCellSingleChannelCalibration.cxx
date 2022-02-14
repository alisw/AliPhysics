// AliEmcalCorrectionCellSingleChannelCalibration
//

#include <TObjArray.h>
#include <TFile.h>
#include "AliEMCALGeometry.h"
#include "AliOADBContainer.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"
#include "AliDataFile.h"

#include "AliEmcalCorrectionCellSingleChannelCalibration.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellSingleChannelCalibration);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellSingleChannelCalibration> AliEmcalCorrectionCellSingleChannelCalibration::reg("AliEmcalCorrectionCellSingleChannelCalibration");

/**
 * Default constructor
 */
AliEmcalCorrectionCellSingleChannelCalibration::AliEmcalCorrectionCellSingleChannelCalibration() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellSingleChannelCalibration")
  ,fCellSingleChannelEnergyDistBefore(0)
  ,fCellSingleChannelEnergyDistAfter(0)
  ,fUseAutomaticRecalib(1)
  ,fCustomRecalibFilePath("")
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellSingleChannelCalibration::~AliEmcalCorrectionCellSingleChannelCalibration()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellSingleChannelCalibration::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  AliWarning("Init EMCAL cell single channel recalibration");
  
  if(fFilepass.Contains("LHC14a1a")) fUseAutomaticRecalib = kTRUE;

  // check the YAML configuration if a custom energy calibration is requested (default is empty string "")
  GetProperty("customRecalibFilePath",fCustomRecalibFilePath);

  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
    
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellSingleChannelCalibration::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fCellSingleChannelEnergyDistBefore = new TH1F("hCellSingleChannelEnergyDistBefore","hCellSingleChannelEnergyDistBefore;E_{cell} (GeV)",1000,0,10);
    fOutput->Add(fCellSingleChannelEnergyDistBefore);
    fCellSingleChannelEnergyDistAfter = new TH1F("hCellSingleChannelEnergyDistAfter","hCellSingleChannelEnergyDistAfter;E_{cell} (GeV)",1000,0,10);
    fOutput->Add(fCellSingleChannelEnergyDistAfter);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellSingleChannelCalibration::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  if (!fEventManager.InputEvent()) {
    AliError("Event ptr = 0, returning");
    return kFALSE;
  }
  
  CheckIfRunChanged();
  
  // CONFIGURE THE RECO UTILS -------------------------------------------------
  fRecoUtils->SwitchOnRecalibration();
  
  // START PROCESSING ---------------------------------------------------------
  // Test if cells present
  if (fCaloCells->GetNumberOfCells()<=0)
  {
    AliDebug(2, Form("Number of EMCAL cells = %d, returning", fCaloCells->GetNumberOfCells()));
    return kFALSE;
  }
  
  // mark the cells not recalibrated
  fRecoUtils->ResetCellsCalibrated();
  
  
  if(fCreateHisto)
    FillCellQA(fCellSingleChannelEnergyDistBefore); // "before" QA
  
  // CELL RECALIBRATION -------------------------------------------------------
  // update cell objects
  UpdateCells();
  
  if(fCreateHisto)
    FillCellQA(fCellSingleChannelEnergyDistAfter); // "after" QA
  
  // switch off recalibrations so those are not done multiple times
  // this is just for safety, the recalibrated flag of cell object
  // should not allow for farther processing anyways
  fRecoUtils->SwitchOffRecalibration();

  return kTRUE;
}

/**
 * Initialize the energy calibration.
 */
Int_t AliEmcalCorrectionCellSingleChannelCalibration::InitRecalib()
{
  if (!fEventManager.InputEvent())
    return 0;
  
  AliInfo("Initialising recalibration factors");
  
  // init default maps first
  if (!fRecoUtils->GetEMCALRecalibrationFactorsArray())
    fRecoUtils->InitEMCALRecalibrationFactors() ;
  
  Int_t runRC = fEventManager.InputEvent()->GetRunNumber();
  
  std::unique_ptr<AliOADBContainer> contRF;
  std::unique_ptr<TFile> recalibFile;
  if (fBasePath!="")
  { //if fBasePath specified
    AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));

    recalibFile = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALSingleChannelCalibrations.root",fBasePath.Data()),"read"));
    if (!recalibFile || recalibFile->IsZombie())
    {
      AliFatal(Form("EMCALSingleChannelCalibrations.root not found in %s",fBasePath.Data()));
      return 0;
    }
    
    contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(recalibFile->Get("AliEMCALSingleChannelCoefficient")));
  }
  else if (fCustomRecalibFilePath!="")
  { //if custom recalib requested
    AliInfo(Form("Loading custom Recalib OADB from given path %s",fCustomRecalibFilePath.Data()));
        
    recalibFile = std::unique_ptr<TFile>(TFile::Open(Form("%s",fCustomRecalibFilePath.Data()),"read"));
    if (!recalibFile || recalibFile->IsZombie())
    {
      AliFatal(Form("Recalibration file not found. Provided path was: %s",fCustomRecalibFilePath.Data()));
      return 0;
    }
    
    contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(recalibFile->Get("AliEMCALSingleChannelCoefficient")));
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading Recalib OADB from OADB/EMCAL");
    
    recalibFile = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALSingleChannelCalibrations.root").data(),"read"));
    if (!recalibFile || recalibFile->IsZombie())
    {
      AliFatal("OADB/EMCAL/EMCALSingleChannelCalibrations.root was not found");
      return 0;
    }
    
    contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(recalibFile->Get("AliEMCALSingleChannelCoefficient")));
  }
  if(!contRF) {
    AliFatal("No OADB container found");
    return 0;
  }
  contRF->SetOwner(true);
  
  TObjArray *recal=(TObjArray*)contRF->GetObject(runRC);
  if (!recal)
  {
    AliFatal(Form("No Objects for run: %d",runRC));
    return 2;
  }
  
  //AliDebug(1, recalib->Print());
  
  Int_t sms = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
  for (Int_t i=0; i<sms; ++i)
  {
    TH2F *h = fRecoUtils->GetEMCALSingleChannelRecalibrationFactors(i);
    if (h)
      delete h;
    h = (TH2F*)recal->FindObject(Form("EMCALSCCalibMap_Mod%d",i));
    if (!h)
    {
      AliFatal(Form("Could not load EMCALSCCalibMap_Mod%d",i));
      continue;
    }
    h->SetDirectory(0);
    fRecoUtils->SetEMCALSingleChannelRecalibrationFactors(i,h);
  }
  
  return 1;
}

/**
 * This function is called if the run changes (it inherits from the base component),
 * to load a new energy calibration and fill relevant variables.
 */
Bool_t AliEmcalCorrectionCellSingleChannelCalibration::CheckIfRunChanged()
{
  Bool_t runChanged = AliEmcalCorrectionComponent::CheckIfRunChanged();
  
  if (runChanged) {
    // init recalibration factors
    if(fUseAutomaticRecalib)
    {
      Int_t fInitRecalib = InitRecalib();
      if (fInitRecalib==0) {
        AliError("InitRecalib returned false, returning");
      }
      if (fInitRecalib==1) {
        AliWarning("InitRecalib OK");
      }
      if (fInitRecalib>1) {
        AliWarning(Form("No recalibration available: %d - %s", fEventManager.InputEvent()->GetRunNumber(), fFilepass.Data()));
      }
    }
  }
  return runChanged;
}
