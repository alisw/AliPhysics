// AliEmcalCorrectionCellTimeCalib
//

#include <TObjArray.h>
#include <TFile.h>
#include "AliEMCALGeometry.h"
#include "AliOADBContainer.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"

#include "AliEmcalCorrectionCellTimeCalib.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellTimeCalib);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellTimeCalib> AliEmcalCorrectionCellTimeCalib::reg("AliEmcalCorrectionCellTimeCalib");

/**
 * Default constructor
 */
AliEmcalCorrectionCellTimeCalib::AliEmcalCorrectionCellTimeCalib() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellTimeCalib")
  ,fCalibrateTime(kFALSE)
  ,fCalibrateTimeL1Phase(kFALSE)
  ,fUseAutomaticTimeCalib(1)
  ,fCellTimeDistBefore(0)
  ,fCellTimeDistAfter(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellTimeCalib::~AliEmcalCorrectionCellTimeCalib()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellTimeCalib::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  AliWarning("Init EMCAL time calibration");
  
  GetProperty("createHistos", fCreateHisto);
  
  fCalibrateTime = kTRUE;

  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
    
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellTimeCalib::UserCreateOutputObjects()
{
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fCellTimeDistBefore = new TH1F("hCellTimeDistBefore","hCellTimeDistBefore;t_cell",1000,-10e-6,10e-6);
    fOutput->Add(fCellTimeDistBefore);
    fCellTimeDistAfter = new TH1F("hCellTimeDistAfter","hCellTimeDistAfter;t_cell",1000,-10e-6,10e-6);
    fOutput->Add(fCellTimeDistAfter);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellTimeCalib::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  if (!fEventManager.InputEvent()) {
    AliError("Event ptr = 0, returning");
    return kFALSE;
  }

  CheckIfRunChanged();
  
  // CONFIGURE THE RECO UTILS -------------------------------------------------
  
  // allows time calibration
  if (fCalibrateTime)
    fRecoUtils->SwitchOnTimeRecalibration();
  else
    fRecoUtils->SwitchOffTimeRecalibration();
  
  // allows time calibration with L1 phase
  if (fCalibrateTimeL1Phase)
    fRecoUtils->SwitchOnL1PhaseInTimeRecalibration();
  else
    fRecoUtils->SwitchOffL1PhaseInTimeRecalibration();
  
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
    FillCellQA(fCellTimeDistBefore); // "before" QA
  
  // CELL RECALIBRATION -------------------------------------------------------
  // cell objects will be updated
  UpdateCells();
  
  if(fCreateHisto)
    FillCellQA(fCellTimeDistAfter); // "after" QA
  
  return kTRUE;
}

/**
 * Initialize the time calibration.
 */
Int_t AliEmcalCorrectionCellTimeCalib::InitTimeCalibration()
{
  // Initialising bad channel maps
  
  if (!fEventManager.InputEvent())
    return 0;

  AliInfo("Initialising time calibration map");
  
  // init default maps first
  if (!fRecoUtils->GetEMCALTimeRecalibrationFactorsArray())
    fRecoUtils->InitEMCALTimeRecalibrationFactors() ;
  
  Int_t runBC = fEventManager.InputEvent()->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
    
    TFile *fbad=new TFile(Form("%s/EMCALTimeCalib.root",fBasePath.Data()),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("EMCALTimeCalib.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile(Form("%s/EMCALTimeCalib.root",fBasePath.Data()),"AliEMCALTimeCalib");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading time calibration OADB from $ALICE_PHYSICS/OADB/EMCAL");
    
    TFile *fbad=new TFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root","read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root was not found");
      return 0;
    }
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root","AliEMCALTimeCalib");
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external time calibration set for run number: %d", runBC));
    delete contBC;
    return 2;
  }
  
  // The calibration object is accessed by specifying a pass
  // For run 1, the actual pass is used (fFilepass, as determined in AliEmcalCorrectionComponent::GetPass())
  // For run 2, the pass is always set to pass1 (as a convention)
  // Other exceptions are hard-coded below
  
  TString pass = fFilepass;
  if (fFilepass=="spc_calo") pass = "pass3";
  if (fRun > 209121) pass = "pass1";
  
  TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(pass);
  if (!arrayBCpass)
  {
    AliError(Form("No external time calibration set for: %d -%s", runBC,pass.Data()));
    delete contBC;
    return 2;
  }
  
  arrayBCpass->Print();
  
  for(Int_t i = 0; i < 4; i++)
  {
    TH1F *h = fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(i);
    if (h)
      delete h;
    
    h = (TH1F*)arrayBCpass->FindObject(Form("hAllTimeAvBC%d",i));
    
    if (!h)
    {
      AliError(Form("Can not get hAllTimeAvBC%d",i));
      continue;
    }
    h->SetDirectory(0);
    fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(i,h);
  }
  
  delete contBC;
  
  return 1;
}

/**
 * Initialize the L1 phase time calibration.
 */
Int_t AliEmcalCorrectionCellTimeCalib::InitTimeCalibrationL1Phase()
{
  // Initialising run-by-run L1 phase in time calibration maps
  
  if (!fEventManager.InputEvent())
    return 0;

  AliInfo("Initialising run-by-run L1 phase in time calibration map");
  
  // init default maps first
  if (!fRecoUtils->GetEMCALL1PhaseInTimeRecalibrationArray())
    fRecoUtils->InitEMCALL1PhaseInTimeRecalibration() ;
  
  Int_t runBC = fEventManager.InputEvent()->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
    
    TFile *timeFile=new TFile(Form("%s/EMCALTimeL1PhaseCalib.root",fBasePath.Data()),"read");
    if (!timeFile || timeFile->IsZombie())
    {
      AliFatal(Form("EMCALTimeL1PhaseCalib.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }
    
    if (timeFile) delete timeFile;
    
    contBC->InitFromFile(Form("%s/EMCALTimeL1PhaseCalib.root",fBasePath.Data()),"AliEMCALTimeL1PhaseCalib");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading L1 phase in time calibration OADB from $ALICE_PHYSICS/OADB/EMCAL");
    
    TFile *timeFile=new TFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeL1PhaseCalib.root","read");
    if (!timeFile || timeFile->IsZombie())
    {
      AliFatal("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeL1PhaseCalib.root was not found");
      return 0;
    }
    
    if (timeFile) delete timeFile;
    
    contBC->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeL1PhaseCalib.root","AliEMCALTimeL1PhaseCalib");
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external L1 phase in time calibration set for run number: %d", runBC));
    delete contBC;
    return 2;
  }
  
  // The calibration object is accessed by specifying a pass
  // For run 2 (which is the only time L1-phase is implemented), the pass is always set to pass1 (as a convention)
  // Other exceptions are hard-coded below
  
  TString pass = "pass1";

  if ( fFilepass=="muon_calo_pass1" && fRun > 209121 && fRun < 244284 )
    pass = "pass0";//period LHC15a-m
  
  TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(pass);
  if (!arrayBCpass)
  {
    AliError(Form("No external L1 phase in time calibration set for: %d -%s", runBC,pass.Data()));
    delete contBC;
    return 2;
  }
  
  arrayBCpass->Print();
  
  
  TH1C *h = fRecoUtils->GetEMCALL1PhaseInTimeRecalibrationForAllSM();
  if (h) delete h;
  
  h = (TH1C*)arrayBCpass->FindObject(Form("h%d",runBC));
  
  if (!h) {
    AliFatal(Form("There is no calibration histogram h%d for this run",runBC));
  }
  h->SetDirectory(0);
  fRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(h);
  
  delete contBC;
  
  return 1;
}

/**
 * This function is called if the run changes (it inherits from the base component),
 * to load a new time calibration and fill relevant variables.
 */
Bool_t AliEmcalCorrectionCellTimeCalib::CheckIfRunChanged()
{
  Bool_t runChanged = AliEmcalCorrectionComponent::CheckIfRunChanged();
  
  if (runChanged) {
 
    // define what recalib parameters are needed for various switches
    // this is based on implementation in AliEMCALRecoUtils
    Bool_t needTimecalib   = fCalibrateTime;
    if(fRun>209121) fCalibrateTimeL1Phase = kTRUE;
    Bool_t needTimecalibL1Phase = fCalibrateTime & fCalibrateTimeL1Phase;
    
    // init time calibration
    if (needTimecalib && fUseAutomaticTimeCalib) {
      Int_t initTC = InitTimeCalibration();
      if (!initTC) {
        AliError("InitTimeCalibration returned false, returning");
      }
      if (initTC==1) {
        AliWarning("InitTimeCalib OK");
      }
      if (initTC > 1) {
        AliWarning(Form("No external time calibration available: %d - %s", fEventManager.InputEvent()->GetRunNumber(), fFilepass.Data()));
      }
    }
    
    // init time calibration with L1 phase
    if (needTimecalibL1Phase && fUseAutomaticTimeCalib) {
      Int_t initTCL1Phase = InitTimeCalibrationL1Phase();
      if (!initTCL1Phase) {
        AliError("InitTimeCalibrationL1Phase returned false, returning");
      }
      if (initTCL1Phase==1) {
        AliWarning("InitTimeCalibL1Phase OK");
      }
      if (initTCL1Phase > 1) {
        AliWarning(Form("No external time calibration L1 phase available: %d - %s", fEventManager.InputEvent()->GetRunNumber(), fFilepass.Data()));
      }
    }
  }
  return runChanged;
}
