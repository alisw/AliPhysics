// AliEmcalCorrectionCellTimeCalib
//

#include <memory>

#include <TObjArray.h>
#include <TFile.h>
#include "AliEMCALGeometry.h"
#include "AliOADBContainer.h"
#include "AliEMCALRecoUtils.h"
#include "AliAODEvent.h"
#include "AliDataFile.h"

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
  ,fCellTimeDistBefore(0)
  ,fCellTimeDistAfter(0)
  ,fCalibrateTime(kFALSE)
  ,fCalibrateTimeL1Phase(kFALSE)
  ,fDoMergedBCs(kFALSE)
  ,fDoCalibrateLowGain(kFALSE)
  ,fDoCalibMergedLG(kFALSE)
  ,fUseAutomaticTimeCalib(1)
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
  
  fCalibrateTime = kTRUE;

  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;

  GetProperty("doMergedBCs", fDoMergedBCs);    

  if (fDoMergedBCs)
    fRecoUtils->SetUseOneHistForAllBCs(fDoMergedBCs);

  GetProperty("doCalibrateLowGain", fDoCalibrateLowGain);

  GetProperty("doCalibMergedLG", fDoCalibMergedLG);    

  if (fDoCalibrateLowGain || fDoCalibMergedLG)
    fRecoUtils->SwitchOnLG();
  else
    fRecoUtils->SwitchOffLG();

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
    fCellTimeDistBefore = new TH1F("hCellTimeDistBefore","hCellTimeDistBefore;t_{cell} (s)",1000,-10e-6,10e-6);
    fOutput->Add(fCellTimeDistBefore);
    fCellTimeDistAfter = new TH1F("hCellTimeDistAfter","hCellTimeDistAfter;t_{cell} (s)",1000,-10e-6,10e-6);
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
  
  std::unique_ptr<AliOADBContainer> contTimeCalib;
  std::unique_ptr<TFile> timeCalibFile;
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
    
    timeCalibFile = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALTimeCalib%s.root",fBasePath.Data(), fDoMergedBCs ? "MergedBCs" : "" ),"read"));
    if (!timeCalibFile || timeCalibFile->IsZombie())
    {
      AliFatal(Form("EMCALTimeCalib%s.root was not found in the path provided: %s", fDoMergedBCs ? "MergedBCs" : "" ,fBasePath.Data()));
      return 0;
    }
    
    contTimeCalib = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(timeCalibFile->Get("AliEMCALTimeCalib")));
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading time calibration OADB from $ALICE_PHYSICS/OADB/EMCAL");
    
    timeCalibFile = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB(Form("EMCAL/EMCALTimeCalib%s.root", fDoMergedBCs ? "MergedBCs" : "")).data(),"read"));
    if (!timeCalibFile || timeCalibFile->IsZombie())
    {
      AliFatal(Form("OADB/EMCAL/EMCALTimeCalib%s.root was not found", fDoMergedBCs ? "MergedBCs" : ""));
      return 0;
    }
    
    contTimeCalib = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(timeCalibFile->Get("AliEMCALTimeCalib")));
  }
  if(!contTimeCalib){
    AliError("No OADB container found");
    return 0;
  }
  contTimeCalib->SetOwner(true);
  
  TObjArray *arrayBC=(TObjArray*)contTimeCalib->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external time calibration set for run number: %d", runBC));
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
    return 2;
  }
  
  arrayBCpass->Print();
  
  if(!fDoMergedBCs){
    for(Int_t i = 0; i < 4; i++)
    {
      TH1F *h = (TH1F*)fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(i);
      if (h)
        delete h;
    
      h = (TH1F*)arrayBCpass->FindObject(Form("hAllTimeAvBC%d",i));
    
      if (!h)
      {
        AliError(Form("Can not get hAllTimeAvBC%d",i));
        continue;
      }
    
      // Shift parameters for bc0 and bc1 in this pass
      if ( pass=="spc_calo" && (i==0 || i==1) ) 
      {
        for(Int_t icell = 0; icell < h->GetNbinsX(); icell++) 
          h->SetBinContent(icell,h->GetBinContent(icell)-100);
      }
    
      h->SetDirectory(0);
      fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(i,h);

      if(fDoCalibrateLowGain){
        TH1F *hLG = (TH1F*)fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(i+4);
        if (hLG)
          delete hLG;
      
        hLG = (TH1F*)arrayBCpass->FindObject(Form("hAllTimeAvLGBC%d",i));
      
        if (!hLG)
        {
          AliError(Form("Can not get hAllTimeAvLGBC%d",i));
          continue;
        }
      
        // Shift parameters for bc0 and bc1 in this pass
        if ( pass=="spc_calo" && (i==0 || i==1) ) 
        {
          for(Int_t icell = 0; icell < hLG->GetNbinsX(); icell++) 
            hLG->SetBinContent(icell,hLG->GetBinContent(icell)-100);
        }
      
        hLG->SetDirectory(0);
        fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(i+4,hLG);
      }
    }
  }else{
  
    TH1S *h = (TH1S*)fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(0);//HG cells
    if (h)
      delete h;
  
    h = (TH1S*)arrayBCpass->FindObject("hAllTimeAv");
  
    if (!h)
      AliError("Can not get hAllTimeAv");
    
    h->SetDirectory(0);
    fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(0,h);//HG cells

    if(fDoCalibrateLowGain && !fDoCalibMergedLG){
      TH1S *hLG = (TH1S*)fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(1);//LG cells
      if (hLG)
        delete hLG;
    
      hLG = (TH1S*)arrayBCpass->FindObject("hAllTimeAvLG");
    
      if (!hLG)
        AliError("Can not get hAllTimeAvLG");
      
      hLG->SetDirectory(0);
      fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(1,hLG);//LG cells
    }
  }

  if(fDoCalibMergedLG){

    std::unique_ptr<AliOADBContainer> contTimeCalibLG;
    std::unique_ptr<TFile> timeCalibFileLG;
    if (fBasePath!="")
    { //if fBasePath specified in the ->SetBasePath()
      AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
      
      timeCalibFileLG = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALTimeCalibMergedBCsLG.root",fBasePath.Data()),"read"));
      if (!timeCalibFileLG || timeCalibFileLG->IsZombie())
      {
        AliFatal(Form("EMCALTimeCalibMergedBCsLG.root was not found in the path provided: %s",fBasePath.Data()));
        return 0;
      }
      
      contTimeCalibLG = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(timeCalibFileLG->Get("AliEMCALTimeCalib")));
    }
    else
    { // Else choose the one in the $ALICE_PHYSICS directory
      AliInfo("Loading time calibration OADB from $ALICE_PHYSICS/OADB/EMCAL");
      
      timeCalibFileLG = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeCalibMergedBCsLG.root").data(),"read"));
      if (!timeCalibFileLG || timeCalibFileLG->IsZombie())
      {
        AliFatal("OADB/EMCAL/EMCALTimeCalibMergedBCsLG.root was not found");
        return 0;
      }
      
      contTimeCalibLG = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(timeCalibFileLG->Get("AliEMCALTimeCalib")));
    }
    if(!contTimeCalibLG){
      AliError("No OADB container found");
      return 0;
    }
    contTimeCalibLG->SetOwner(true);
    
    TObjArray *arrayBCLG=(TObjArray*)contTimeCalibLG->GetObject(runBC);
    if (!arrayBCLG)
    {
      AliError(Form("No external time calibration set for run number: %d", runBC));
      return 2;
    }
     
    TObjArray *arrayBCpassLG=(TObjArray*)arrayBCLG->FindObject(pass);
    if (!arrayBCpassLG)
    {
      AliError(Form("No external time calibration set for: %d -%s", runBC,pass.Data()));
      return 2;
    }
    
    arrayBCpassLG->Print();

    TH1S *hLG = (TH1S*)fRecoUtils->GetEMCALChannelTimeRecalibrationFactors(1);//LG cells
    if (hLG)
      delete hLG;
     
    hLG = (TH1S*)arrayBCpassLG->FindObject("hAllTimeAvLG");
      
    if (!hLG)
      AliError("Can not get hAllTimeAvLG");
        
    hLG->SetDirectory(0);
    fRecoUtils->SetEMCALChannelTimeRecalibrationFactors(1,hLG);//LG cells

  }
  
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
  
  std::unique_ptr<AliOADBContainer> contTimeCalib;
  std::unique_ptr<TFile> timeFile;
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    AliInfo(Form("Loading time calibration OADB from given path %s",fBasePath.Data()));
    
    timeFile = std::unique_ptr<TFile>(TFile::Open(Form("%s/EMCALTimeL1PhaseCalib.root",fBasePath.Data()),"read"));
    if (!timeFile || timeFile->IsZombie())
    {
      AliFatal(Form("EMCALTimeL1PhaseCalib.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }
    
    contTimeCalib = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(timeFile->Get("AliEMCALTimeL1PhaseCalib")));
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading L1 phase in time calibration OADB from OADB/EMCAL");
    
    timeFile = std::unique_ptr<TFile>(TFile::Open(AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeL1PhaseCalib.root").data(),"read"));
    if (!timeFile || timeFile->IsZombie())
    {
      AliFatal("OADB/EMCAL/EMCALTimeL1PhaseCalib.root was not found");
      return 0;
    }
    
    contTimeCalib = std::unique_ptr<AliOADBContainer>(static_cast<AliOADBContainer *>(timeFile->Get("AliEMCALTimeL1PhaseCalib")));
  }
  if(!contTimeCalib){
    AliError("No OADB container found");
    return 0;
  }
  contTimeCalib->SetOwner(true);
  
  TObjArray *arrayBC=(TObjArray*)contTimeCalib->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external L1 phase in time calibration set for run number: %d", runBC));
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
    return 2;
  }
  
  arrayBCpass->Print();
  
  
  TH1C *h = fRecoUtils->GetEMCALL1PhaseInTimeRecalibrationForAllSM(0);
  if (h) delete h;
  
  h = (TH1C*)arrayBCpass->FindObject(Form("h%d",runBC));
  
  if (!h) {
    AliFatal(Form("There is no calibration histogram h%d for this run",runBC));
  }
  h->SetDirectory(0);
  fRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(h,0);

  //Now special case for PAR runs
  fRecoUtils->SwitchOffParRun();
  //access tree from OADB file
  TTree *tGID = (TTree*)arrayBCpass->FindObject(Form("h%d_GID",runBC));
  if(tGID){//check whether present = PAR run
    fRecoUtils->SwitchOnParRun();
    //access tree branch with PARs
    ULong64_t parGlobalBCs;
    tGID->SetBranchAddress("GID",&parGlobalBCs);
    //set number of PARs in run
    Short_t nPars = (Short_t) tGID->GetEntries();
    fRecoUtils->SetNPars((Short_t)nPars);
    //set global ID for each PAR
    for (Short_t iParNumber = 0; iParNumber < nPars; ++iParNumber) {
      tGID->GetEntry(iParNumber);
      fRecoUtils->SetGlobalIDPar(parGlobalBCs,iParNumber);
    }//loop over entries  

    //access GlobalID hiostograms for each PAR
    for(Short_t iParNumber=1; iParNumber< fRecoUtils->GetNPars()+1;iParNumber++){
      TH1C *hPar = (TH1C*)arrayBCpass->FindObject( Form("h%d_%llu",runBC,fRecoUtils->GetGlobalIDPar(iParNumber-1) ) );
      if (!hPar) AliError( Form("Could not load h%d_%llu",runBC,fRecoUtils->GetGlobalIDPar(iParNumber-1) ) );
      hPar->SetDirectory(0);
      fRecoUtils->SetEMCALL1PhaseInTimeRecalibrationForAllSM(hPar,iParNumber);
    }//loop over PARs
  }//end if tGID present  
  
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

    if(fRun<225000 && fCalibrateTime && fDoMergedBCs){
      AliWarning("The merged BC histograms don't exist for Run1, falling back to the 4 BC histograms. For question contact constantin.loizides@cern.ch");
      fDoMergedBCs=kFALSE;
      fRecoUtils->SetUseOneHistForAllBCs(kFALSE);
      if(fDoCalibrateLowGain && !fDoCalibMergedLG){
        AliWarning("The merged BC LG histograms don't exist for Run1, using the all period merged LG histograms");
        fDoCalibMergedLG=kTRUE;
      }
    }

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
