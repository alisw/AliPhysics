// AliEmcalCorrectionCellEnergy
//
//  Based on elements of AliEMCALTenderSupply:
/// \author Deepa Thomas (Utrecht University)
/// \author Jiri Kral (University of Jyvaskyla), mods/rewrite
/// \author Salvatore Aiola, make it work for AODs
/// \author C. Loizides, make it work for AODs
/// \author Gustavo Conesa, LPSC-Grenoble, several mods.

#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TInterpreter.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TTree.h>
#include "AliAODMCParticle.h"
#include "AliEMCALAfterBurnerUF.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecParam.h"
#include "AliMagF.h"
#include "AliOADBContainer.h"

#include "AliEmcalCorrectionCellEnergy.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionCellEnergy);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionCellEnergy> AliEmcalCorrectionCellEnergy::reg("AliEmcalCorrectionCellEnergy");

//________________________________________________________________________
AliEmcalCorrectionCellEnergy::AliEmcalCorrectionCellEnergy() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionCellEnergy")
  ,fUseAutomaticRecalib(1)
  ,fUseAutomaticRunDepRecalib(1)
  ,fCellEnergyDistBefore(0)
  ,fCellEnergyDistAfter(0)
{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
}

//________________________________________________________________________
AliEmcalCorrectionCellEnergy::~AliEmcalCorrectionCellEnergy()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionCellEnergy::Initialize()
{
  // Initialization
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Initialize();
  // Do base class initializations and if it fails -> bail out
  //AliAnalysisTaskEmcal::ExecOnce();
  //if (!fInitialized) return;
  
  AliWarning("Init EMCAL cell recalibration");
  
  GetProperty("createHistos", fCreateHisto);
  
  if (fCreateHisto){
    fCellEnergyDistBefore = new TH1F("hCellEnergyDistBefore","hCellEnergyDistBefore;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistBefore);
    fCellEnergyDistAfter = new TH1F("hCellEnergyDistAfter","hCellEnergyDistAfter;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistAfter);
  }

  if(fFilepass.Contains("LHC14a1a")) fUseAutomaticRecalib = kTRUE;
  
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;
    
  // missalignment function -- TO DO: not sure if we need this or not
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionCellEnergy::Run()
{
  // Run
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Run();
  
  if (!fEvent) {
    AliError("Event ptr = 0, returning");
    return kFALSE;
  }
  
  // Initialising parameters once per run number
  
  if (RunChanged())
  {
    fRun = fEvent->GetRunNumber();
    AliWarning(Form("Run changed, initializing parameters for %d", fRun));
    if (dynamic_cast<AliAODEvent*>(fEvent)) {
      AliWarning("=============================================================");
      AliWarning("===  Running on AOD is not equivalent to running on ESD!  ===");
      AliWarning("=============================================================");
    }
    
    // init geometry if not already done
    fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(fRun);
    
    if (!fGeom)
    {
      AliFatal("Can not create geometry");
      return kFALSE;
    }
    
    // init recalibration factors
    if(fUseAutomaticRecalib)
    {
      Int_t fInitRecalib = InitRecalib();
      if (fInitRecalib==0)
        AliError("InitRecalib returned false, returning");
      if (fInitRecalib==1)
        AliWarning("InitRecalib OK");
      if (fInitRecalib>1)
        AliWarning(Form("No recalibration available: %d - %s", fEvent->GetRunNumber(), fFilepass.Data()));
    }
    
    if(fUseAutomaticRunDepRecalib)
    {
      Int_t fInitRunDepRecalib = InitRunDepRecalib();
      if (fInitRunDepRecalib==0)
        AliError("InitrunDepRecalib returned false, returning");
      if (fInitRunDepRecalib==1)
        AliWarning("InitRecalib OK");
      if (fInitRunDepRecalib>1)
        AliWarning(Form("No Temperature recalibration available: %d - %s", fEvent->GetRunNumber(), fFilepass.Data()));
    }
    
    //AliDebug(2, Form("%s", fRecoUtils->Print("")));
  }
  
  // CONFIGURE THE RECO UTILS -------------------------------------------------
  fRecoUtils->SwitchOnRecalibration();
  
  // START PROCESSING ---------------------------------------------------------
  // Test if cells present
  if (fCaloCells->GetNumberOfCells()<=0)
  {
    AliDebug(2, Form("Number of EMCAL cells = %d, returning", fCaloCells->GetNumberOfCells()));
    return kTRUE;
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

//_____________________________________________________
Int_t AliEmcalCorrectionCellEnergy::InitRecalib()
{
  // Initialising recalibration factors.
  
  if (!fEvent)
    return 0;
  
  AliInfo("Initialising recalibration factors");
  
  // init default maps first
  if (!fRecoUtils->GetEMCALRecalibrationFactorsArray())
    fRecoUtils->InitEMCALRecalibrationFactors() ;
  
  Int_t runRC = fEvent->GetRunNumber();
  
  AliOADBContainer *contRF=new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified
    AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));
    
    TFile *fRecalib= new TFile(Form("%s/EMCALRecalib.root",fBasePath.Data()),"read");
    if (!fRecalib || fRecalib->IsZombie())
    {
      AliFatal(Form("EMCALRecalib.root not found in %s",fBasePath.Data()));
      return 0;
    }
    
    if (fRecalib) delete fRecalib;
    
    contRF->InitFromFile(Form("%s/EMCALRecalib.root",fBasePath.Data()),"AliEMCALRecalib");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading Recalib OADB from $ALICE_PHYSICS/OADB/EMCAL");
    
    TFile *fRecalib= new TFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALRecalib.root","read");
    if (!fRecalib || fRecalib->IsZombie())
    {
      AliFatal("$ALICE_PHYSICS/OADB/EMCAL/EMCALRecalib.root was not found");
      return 0;
    }
    
    if (fRecalib) delete fRecalib;
    
    contRF->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALRecalib.root","AliEMCALRecalib");
  }
  
  TObjArray *recal=(TObjArray*)contRF->GetObject(runRC);
  if (!recal)
  {
    AliError(Form("No Objects for run: %d",runRC));
    delete contRF;
    return 2;
  }
  
  TObjArray *recalpass=(TObjArray*)recal->FindObject(fFilepass);
  if (!recalpass)
  {
    AliError(Form("No Objects for run: %d - %s",runRC,fFilepass.Data()));
    delete contRF;
    return 2;
  }
  
  TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
  if (!recalib)
  {
    AliError(Form("No Recalib histos found for  %d - %s",runRC,fFilepass.Data()));
    delete contRF;
    return 2;
  }
  
  //AliDebug(1, recalib->Print());
  
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
  
  delete contRF;
  
  return 1;
}

//_____________________________________________________
Int_t AliEmcalCorrectionCellEnergy::InitRunDepRecalib()
{
  // Initialising recalibration factors.
  
  if (!fEvent)
    return 0;
  
  AliInfo("Initialising recalibration factors");
  
  // init default maps first
  if (!fRecoUtils->GetEMCALRecalibrationFactorsArray())
    fRecoUtils->InitEMCALRecalibrationFactors() ;
  
  Int_t runRC = fEvent->GetRunNumber();
  
  AliOADBContainer *contRF=new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    AliInfo(Form("Loading Recalib OADB from given path %s",fBasePath.Data()));
    
    TFile *fRunDepRecalib= new TFile(Form("%s/EMCALTemperatureCorrCalib.root",fBasePath.Data()),"read");
    if (!fRunDepRecalib || fRunDepRecalib->IsZombie())
    {
      AliFatal(Form("EMCALTemperatureCorrCalib.root not found in %s",fBasePath.Data()));
      return 0;
    }
    
    if (fRunDepRecalib) delete fRunDepRecalib;
    
    contRF->InitFromFile(Form("%s/EMCALTemperatureCorrCalib.root",fBasePath.Data()),"AliEMCALRunDepTempCalibCorrections");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading Recalib OADB from $ALICE_PHYSICS/OADB/EMCAL");
    
    TFile *fRunDepRecalib= new TFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTemperatureCorrCalib.root","read");
    if (!fRunDepRecalib || fRunDepRecalib->IsZombie())
    {
      AliFatal("$ALICE_PHYSICS/OADB/EMCAL/EMCALTemperatureCorrCalib.root was not found");
      return 0;
    }
    
    if (fRunDepRecalib) delete fRunDepRecalib;
    
    contRF->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALTemperatureCorrCalib.root","AliEMCALRunDepTempCalibCorrections");
  }
  
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
  Int_t nbins = rundeprecal->GetNbinsX();
  
  // Avoid use of Run1 param for Run2
  if(nSM > 12 && nbins < 12288)
  {
    AliError(Form("Total SM is %d but T corrections available for %d channels, skip Init of T recalibration factors",nSM,nbins));
    
    delete contRF;
    
    return 2;
  }
  
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
  
  delete contRF;
  
  return 1;
}
