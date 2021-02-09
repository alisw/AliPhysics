// AliEmcalCorrectionCellEnergy
//

#include <TObjArray.h>
#include <TFile.h>
#include "AliEMCALGeometry.h"
#include "AliOADBContainer.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"
#include "AliDataFile.h"

#include "AliEmcalCorrectionCellEnergy.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellEnergy);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellEnergy> AliEmcalCorrectionCellEnergy::reg("AliEmcalCorrectionCellEnergy");

/**
 * Default constructor
 */
AliEmcalCorrectionCellEnergy::AliEmcalCorrectionCellEnergy() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellEnergy")
  ,fCellEnergyDistBefore(0)
  ,fCellEnergyDistAfter(0)
  ,fUseAutomaticRecalib(1)
  ,fUseAutomaticRunDepRecalib(1)
  ,fUseNewRunDepTempCalib(0)
  ,fDisableTempCalib(0)
  ,fUseShaperCorrection(0)
  ,fCustomRecalibFilePath("")
  ,fLoad1DRecalibFactors(0)
{
}

/**
 * Destructor
 */
AliEmcalCorrectionCellEnergy::~AliEmcalCorrectionCellEnergy()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionCellEnergy::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  AliWarning("Init EMCAL cell recalibration");
  
  if(fFilepass.Contains("LHC14a1a")) fUseAutomaticRecalib = kTRUE;

  // check the YAML configuration if the Run2 calibration is requested (default is false)
  GetProperty("enableNewTempCalib",fUseNewRunDepTempCalib);

  // check the YAML configuration if temperature calibration is disabled explicitly (default is false)
  GetProperty("disableTempCalib",fDisableTempCalib);

  // check the YAML configuration if the shaper nonlinearity correction is requested (default is false)
  GetProperty("enableShaperCorrection",fUseShaperCorrection);

  // check the YAML configuration if a custom energy calibration is requested (default is empty string "")
  GetProperty("customRecalibFilePath",fCustomRecalibFilePath);

  //
  GetProperty("load1DRecalibFactors",fLoad1DRecalibFactors);

  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;

  fRecoUtils->SetUse1DRecalibration(fLoad1DRecalibFactors);
    
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionCellEnergy::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  if (fCreateHisto){
    fCellEnergyDistBefore = new TH1F("hCellEnergyDistBefore","hCellEnergyDistBefore;E_{cell} (GeV)",7000,0,70);
    fOutput->Add(fCellEnergyDistBefore);
    fCellEnergyDistAfter = new TH1F("hCellEnergyDistAfter","hCellEnergyDistAfter;E_{cell} (GeV)",7000,0,70);
    fOutput->Add(fCellEnergyDistAfter);
  }
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionCellEnergy::Run()
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
    FillCellQA(fCellEnergyDistBefore); // "before" QA
  
  // CELL RECALIBRATION -------------------------------------------------------
  // update cell objects
  UpdateCells();
  
  if(fCreateHisto)
    FillCellQA(fCellEnergyDistAfter); // "after" QA
  
  // switch off recalibrations so those are not done multiple times
  // this is just for safety, the recalibrated flag of cell object
  // should not allow for farther processing anyways
  fRecoUtils->SwitchOffRecalibration();

  return kTRUE;
}

/**
 * Initialize the energy calibration.
 */
Int_t AliEmcalCorrectionCellEnergy::InitRecalib()
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

    recalibFile = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALRecalib%s.root",fBasePath.Data(), fLoad1DRecalibFactors ? "_1D" : "" ),"read"));
    if (!recalibFile || recalibFile->IsZombie())
    {
      AliFatal(Form("EMCALRecalib%s.root not found in %s", fLoad1DRecalibFactors ? "_1D" : "", fBasePath.Data()));
      return 0;
    }
    
    contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(recalibFile->Get("AliEMCALRecalib")));
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
    
    contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(recalibFile->Get("AliEMCALRecalib")));
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading Recalib OADB from OADB/EMCAL");
    
    recalibFile = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB(Form("EMCAL/EMCALRecalib%s.root", fLoad1DRecalibFactors ? "_1D" : "")).data(),"read"));
    if (!recalibFile || recalibFile->IsZombie())
    {
      AliFatal(Form("OADB/EMCAL/EMCALRecalib%s.root was not found", fLoad1DRecalibFactors ? "_1D" : ""));
      return 0;
    }
    
    contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(recalibFile->Get("AliEMCALRecalib")));
  }
  if(!contRF) {
    AliError("No OADB container found");
    return 0;
  }
  contRF->SetOwner(true);
  
  TObjArray *recal=(TObjArray*)contRF->GetObject(runRC);
  if (!recal)
  {
    AliError(Form("No Objects for run: %d",runRC));
    return 2;
  }
  
  TObjArray *recalpass=(TObjArray*)recal->FindObject(fFilepass);
  if (!recalpass)
  {
    AliError(Form("No Objects for run: %d - %s",runRC,fFilepass.Data()));
    return 2;
  }
  
  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  if (!recalib)
  {
    AliError(Form("No Recalib histos found for  %d - %s",runRC,fFilepass.Data()));
    return 2;
  }
  
  //AliDebug(1, recalib->Print());
 

  if(fLoad1DRecalibFactors){
    TH1S *h = fRecoUtils->GetEMCALChannelRecalibrationFactors1D();
    if (h)
      delete h;
    h=(TH1S*)recalib->FindObject("EMCALRecalFactors");
      
    if (!h)
    {
      AliError("Can not get EMCALRecalFactors");
    }
    h->SetDirectory(0);
    fRecoUtils->SetEMCALChannelRecalibrationFactors1D(h);
  }else{
    Int_t sms = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
    for (Int_t i=0; i<sms; ++i)
    {
      TH2F *h = fRecoUtils->GetEMCALChannelRecalibrationFactors(i);
      if (h)
        delete h;
      h = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
      if (!h)
      {
        AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
        continue;
      }
      h->SetDirectory(0);
      fRecoUtils->SetEMCALChannelRecalibrationFactors(i,h);
    }
  }
  
  return 1;
}

/**
 * Initialize the temperature calibration.
 */
Int_t AliEmcalCorrectionCellEnergy::InitRunDepRecalib()
{
  if (!fEventManager.InputEvent())
    return 0;

  Int_t runRC = fEventManager.InputEvent()->GetRunNumber();

  // init default maps first
  if (!fRecoUtils->GetEMCALRecalibrationFactorsArray())
    fRecoUtils->InitEMCALRecalibrationFactors();

  // Treat new temp. calibration differently. Loading of two OADB objects required for calibration
  // Calibration can be switched between std. run 1 and run2 calibration using: enableNewTempCalib: true in the YAML configuration
  // Calibration can be completely switched off using: disableTempCalib: true in the YAML configuration (this should only be done by experts!)
  if (fDisableTempCalib){
      AliInfo("Temperature calibration switched off. This mode is only meant for experts, please be sure you really wanted to switch it off.");
      return 2;
  }
  
  if(fUseNewRunDepTempCalib){
    AliInfo("Initialising New recalibration factors");

    // two files and two OADB containers are needed for the correction factor
    std::unique_ptr<AliOADBContainer> contTemperature;
    std::unique_ptr<AliOADBContainer> contParams;
    std::unique_ptr<TFile> runDepTemperatureFile;
    std::unique_ptr<TFile> temperatureCalibParamFile;

    if (fBasePath!="")
    { //if fBasePath specified in the ->SetBasePath()
      runDepTemperatureFile = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALTemperatureCalibSM.root",fBasePath.Data()),"read"));
      if (!runDepTemperatureFile || runDepTemperatureFile->IsZombie()) {
        AliFatal(Form("EMCALTemperatureCalibSM.root not found in %s",fBasePath.Data()));
        return 0;
      }

      temperatureCalibParamFile = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALTemperatureCalibParam.root",fBasePath.Data()),"read"));
      if (!temperatureCalibParamFile || temperatureCalibParamFile->IsZombie()) {
        AliFatal(Form("EMCALTemperatureCalibParam.root not found in %s",fBasePath.Data()));
        return 0;
      }

      contTemperature = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(runDepTemperatureFile->Get("AliEMCALTemperatureCalibSM")));
      contParams = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(temperatureCalibParamFile->Get("AliEMCALTemperatureCalibParam")));
    }
    else
    { // Else choose the one in the $ALICE_PHYSICS directory or on EOS via the wrapper function
      runDepTemperatureFile = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCalibSM.root").data(),"read"));
      if (!runDepTemperatureFile || runDepTemperatureFile->IsZombie()) {
        AliFatal("OADB/EMCAL/EMCALTemperatureCalibSM.root was not found");
        return 0;
      }

      temperatureCalibParamFile = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCalibParam.root").data(),"read"));
      if (!temperatureCalibParamFile || temperatureCalibParamFile->IsZombie()) {
        AliFatal("OADB/EMCAL/EMCALTemperatureCalibParam.root was not found");
        return 0;
      }

      contTemperature = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(runDepTemperatureFile->Get("AliEMCALTemperatureCalibSM")));
      contParams = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(temperatureCalibParamFile->Get("AliEMCALTemperatureCalibParam")));
    }

    if(!contTemperature || !contParams) {
      AliError("Temperature or parametrization OADB container not found");
      return 0;
    }
    contTemperature->SetOwner(true);
    contParams->SetOwner(true);

    TObjArray *arrayParams=(TObjArray*)contParams->GetObject(runRC);
    if (!arrayParams)
    {
      AliError(Form("No temperature calibration parameters can be found for run number: %d", runRC));
      return 0;
    }
    TH1D *hRundepTemp = (TH1D*)contTemperature->GetObject(runRC);
    TH1F *hSlopeParam = (TH1F*)arrayParams->FindObject("hParamSlope");
    TH1F *hA0Param = (TH1F*)arrayParams->FindObject("hParamA0");

    if (!hRundepTemp || !hSlopeParam || !hA0Param)
    {
      AliError(Form("Histogram missing for Run2 temperature calibration for run number: %d", runRC));
      return 0;
    }

    Int_t nSM = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
    for (Int_t ism=0; ism<nSM; ++ism)
    {
      Double_t temperature = hRundepTemp->GetBinContent(ism+1);
      for (Int_t icol=0; icol<48; ++icol)
      {
        for (Int_t irow=0; irow<24; ++irow)
        {
          Int_t absID     = fGeom->GetAbsCellIdFromCellIndexes(ism, irow, icol); // original calibration factor
          Float_t factor  = fRecoUtils->GetEMCALChannelRecalibrationFactor(ism,icol,irow);
          Float_t slope   = 0;
          Float_t offset  = 0;
          slope           = hSlopeParam->GetBinContent(absID+1);
          offset          = hA0Param->GetBinContent(absID+1);
          // Correction is the inverse of the calculated factor
          if(slope || offset)
            factor *= 1 / (offset + (slope * temperature) ); // correction dependent on T
          fRecoUtils->SetEMCALChannelRecalibrationFactor(ism,icol,irow,factor);
        } // columns
      } // rows
    } // SM loop

    return 1;

  // standard treatment of Run1 data
  } else {
    if(runRC > 197692){
      AliInfo("Temperature calibration could not be loaded. Please use enableNewTempCalib: true in your configuration file for Run2 data!");
      return 0;
    }
    AliInfo("Initialising Run1 recalibration factors");

    std::unique_ptr<AliOADBContainer> contRF;
    std::unique_ptr<TFile> runDepRecalibFile;
    if (fBasePath!="")
    { //if fBasePath specified in the ->SetBasePath()
      AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));

      runDepRecalibFile = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALTemperatureCorrCalib.root",fBasePath.Data()),"read"));
      if (!runDepRecalibFile || runDepRecalibFile->IsZombie())
      {
        AliFatal(Form("EMCALTemperatureCorrCalib.root not found in %s",fBasePath.Data()));
        return 0;
      }

      contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(runDepRecalibFile->Get("AliEMCALRunDepTempCalibCorrections")));
    }
    else
    { // Else choose the one in the $ALICE_PHYSICS directory or on EOS via the wrapper function
      AliInfo("Loading Recalib OADB from OADB/EMCAL");

      runDepRecalibFile = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root").data(),"read"));
      if (!runDepRecalibFile || runDepRecalibFile->IsZombie())
      {
        AliFatal("OADB/EMCAL/EMCALTemperatureCorrCalib.root was not found");
        return 0;
      }

      contRF = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(runDepRecalibFile->Get("AliEMCALRunDepTempCalibCorrections")));
    }
    if(!contRF) {
      AliError("No OADB container found");
      return 0;
    }
    contRF->SetOwner(true);

    TH1S *rundeprecal=(TH1S*)contRF->GetObject(runRC);

    if (!rundeprecal)
    {
      AliWarning(Form("No TemperatureCorrCalib Objects for run: %d",runRC));
      // let's get the closest runnumber instead then..
      Int_t lower = 0;
      Int_t ic = 0;
      Int_t maxEntry = contRF->GetNumberOfEntries();

      while ((ic < maxEntry) && (contRF->UpperLimit(ic) < runRC)) {
        lower = ic;
        ic++;
      }

      Int_t closest = lower;
      if ((ic<maxEntry) &&
          (contRF->LowerLimit(ic)-runRC) < (runRC - contRF->UpperLimit(lower))) {
        closest = ic;
      }

      AliWarning(Form("TemperatureCorrCalib Objects found closest id %d from run: %d", closest, contRF->LowerLimit(closest)));
      rundeprecal = (TH1S*) contRF->GetObjectByIndex(closest);
    }

    Int_t nSM = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();

    //AliDebug(1, rundeprecal->Print());

    for (Int_t ism=0; ism<nSM; ++ism)
    {
      for (Int_t icol=0; icol<48; ++icol)
      {
        for (Int_t irow=0; irow<24; ++irow)
        {
          Float_t factor = fRecoUtils->GetEMCALChannelRecalibrationFactor(ism,icol,irow);

          Int_t absID = fGeom->GetAbsCellIdFromCellIndexes(ism, irow, icol); // original calibration factor
          factor *= rundeprecal->GetBinContent(absID) / 10000. ; // correction dependent on T

          fRecoUtils->SetEMCALChannelRecalibrationFactor(ism,icol,irow,factor);
        } // columns
      } // rows
    } // SM loop

    return 1;
  }
}

/**
 * This function is called if the run changes (it inherits from the base component),
 * to load a new energy calibration and fill relevant variables.
 */
Bool_t AliEmcalCorrectionCellEnergy::CheckIfRunChanged()
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
    
    if(fUseAutomaticRunDepRecalib)
    {
      Int_t fInitRunDepRecalib = InitRunDepRecalib();
      if (fInitRunDepRecalib==0) {
        AliError("InitrunDepRecalib returned false, returning");
      }
      if (fInitRunDepRecalib==1) {
        AliWarning("InitRecalib OK");
      }
      if (fInitRunDepRecalib>1) {
        AliWarning(Form("No Temperature recalibration available: %d - %s", fEventManager.InputEvent()->GetRunNumber(), fFilepass.Data()));
      }
    }
  }
  if(fUseShaperCorrection)
  {
    fRecoUtils->SetUseTowerShaperNonlinarityCorrection(kTRUE);
  }
  return runChanged;
}
