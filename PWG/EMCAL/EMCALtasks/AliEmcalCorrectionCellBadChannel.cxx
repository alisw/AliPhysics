// AliEmcalCorrectionCellBadChannel
//

#include <TObjArray.h>
#include <TFile.h>
#include "AliEMCALGeometry.h"
#include "AliOADBContainer.h"
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
 
  if (fCreateHisto){
    fCellEnergyDistBefore = new TH1F("hCellEnergyDistBefore","hCellEnergyDistBefore;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistBefore);
    fCellEnergyDistAfter = new TH1F("hCellEnergyDistAfter","hCellEnergyDistAfter;E_cell",1000,0,10);
    fOutput->Add(fCellEnergyDistAfter);
  }
  
  // init reco utils
  if (!fRecoUtils)
    fRecoUtils  = new AliEMCALRecoUtils;

  // missalignment function -- TO DO: do we need this?
  fRecoUtils->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);

  return kTRUE;
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
    
    // init bad channels
    Int_t fInitBC = InitBadChannels();
    if (fInitBC==0)
      AliError("InitBadChannels returned false, returning");
    if (fInitBC==1)
      AliWarning("InitBadChannels OK");
    if (fInitBC>1)
      AliWarning(Form("No external hot channel set: %d - %s", fEvent->GetRunNumber(), fFilepass.Data()));
    
    //AliDebug(1, Form("%s",fRecoUtils->Print("")));
  }
  
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

//_____________________________________________________
Int_t AliEmcalCorrectionCellBadChannel::InitBadChannels()
{
  // Initialising bad channel maps
  
  if (!fEvent)
    return 0;
  
  AliInfo("Initialising Bad channel map");
  
  // init default maps first
  if (!fRecoUtils->GetEMCALBadChannelStatusMapArray())
    fRecoUtils->InitEMCALBadChannelStatusMap() ;
  
  Int_t runBC = fEvent->GetRunNumber();
  
  AliOADBContainer *contBC = new AliOADBContainer("");
  if (fBasePath!="")
  { //if fBasePath specified in the ->SetBasePath()
    AliInfo(Form("Loading Bad Channels OADB from given path %s",fBasePath.Data()));
    
    TFile *fbad=new TFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal(Form("EMCALBadChannels.root was not found in the path provided: %s",fBasePath.Data()));
      return 0;
    }
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile(Form("%s/EMCALBadChannels.root",fBasePath.Data()),"AliEMCALBadChannels");
  }
  else
  { // Else choose the one in the $ALICE_PHYSICS directory
    AliInfo("Loading Bad Channels OADB from $ALICE_PHYSICS/OADB/EMCAL");
    
    TFile *fbad=new TFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALBadChannels.root","read");
    if (!fbad || fbad->IsZombie())
    {
      AliFatal("$ALICE_PHYSICS/OADB/EMCAL/EMCALBadChannels.root was not found");
      return 0;
    }
    
    if (fbad) delete fbad;
    
    contBC->InitFromFile("$ALICE_PHYSICS/OADB/EMCAL/EMCALBadChannels.root","AliEMCALBadChannels");
  }
  
  TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runBC);
  if (!arrayBC)
  {
    AliError(Form("No external hot channel set for run number: %d", runBC));
    delete contBC;
    return 2;
  }
  
  Int_t sms = fGeom->GetEMCGeometry()->GetNumberOfSuperModules();
  for (Int_t i=0; i<sms; ++i)
  {
    TH2I *h = fRecoUtils->GetEMCALChannelStatusMap(i);
    if (h)
      delete h;
    h=(TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
    
    if (!h)
    {
      AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
      continue;
    }
    h->SetDirectory(0);
    fRecoUtils->SetEMCALChannelStatusMap(i,h);
  }
  
  delete contBC;
  
  return 1;
}
