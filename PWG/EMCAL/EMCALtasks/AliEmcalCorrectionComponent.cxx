// AliEmcalCorrectionComponent
//

#include "AliEmcalCorrectionComponent.h"

#include <TFile.h>
#include <TH1.h>

#include "AliEmcalList.h"
#include "AliEMCALRecoUtils.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliClusterContainer.h"
#include "AliTrackContainer.h"
#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliOADBContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionComponent);
/// \endcond

// Must create an instance of the map, since it is static
AliEmcalCorrectionComponentFactory::map_type * AliEmcalCorrectionComponentFactory::componentMap = new AliEmcalCorrectionComponentFactory::map_type;

/**
 * Default constructor
 */
AliEmcalCorrectionComponent::AliEmcalCorrectionComponent() :
  TNamed("AliEmcalCorrectionComponent", "AliEmcalCorrectionComponent"),
  fUserConfiguration(),
  fDefaultConfiguration(),
  fCreateHisto(kTRUE),
  fRun(-1),
  fFilepass(""),
  fGetPassFromFileName(kTRUE),
  fEventManager(),
  fEsdMode(0),
  fMCEvent(0),
  fCent(0),
  fNcentBins(4),
  fCentBin(0),
  fNbins(250),
  fMinBinPt(0),
  fMaxBinPt(250),
  fGeom(0),
  fIsEmbedded(kFALSE),
  fMinMCLabel(0),
  fClusterCollArray(),
  fParticleCollArray(),
  fCaloCells(0),
  fRecoUtils(0),
  fOutput(0),
  fBasePath("")

{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

/**
 * Standard constructor
 */
AliEmcalCorrectionComponent::AliEmcalCorrectionComponent(const char * name) :
  TNamed(name, name),
  fUserConfiguration(),
  fDefaultConfiguration(),
  fCreateHisto(kTRUE),
  fRun(-1),
  fFilepass(""),
  fGetPassFromFileName(kTRUE),
  fEventManager(),
  fEsdMode(0),
  fMCEvent(0),
  fCent(0),
  fNcentBins(4),
  fCentBin(0),
  fNbins(250),
  fMinBinPt(0),
  fMaxBinPt(250),
  fGeom(0),
  fIsEmbedded(kFALSE),
  fMinMCLabel(0),
  fClusterCollArray(),
  fParticleCollArray(),
  fCaloCells(0),
  fRecoUtils(0),
  fOutput(0),
  fBasePath("")
{
  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

/**
 * Destructor
 */
AliEmcalCorrectionComponent::~AliEmcalCorrectionComponent()
{
}

/**
 * Initialize basic variables in the correction component from the configuration file.
 */
Bool_t AliEmcalCorrectionComponent::Initialize()
{
  // Read in pass. If it is empty, set flag to automatically find the pass from the filename.
  std::string tempString = "";
  GetProperty("pass", tempString);
  fFilepass = tempString.c_str();
  if (fFilepass != "") {
    fGetPassFromFileName = kFALSE;
    // Handle the "default" value used in MC
    if (fFilepass == "default") {
      AliError("Received \"default\" as pass value. Defaulting to \"pass1\"! In the case of MC, the user should set the proper pass value in their configuration file! For data, empty quotes should be set so that the pass is automatically set.");
      fFilepass = "pass1";
    }
  }

  return kTRUE;
}

/**
 * Create output objects for the analysis. Similar to UserCreateOutputObjects() in
 * AliAnalysisTaskSE
 */
void AliEmcalCorrectionComponent::UserCreateOutputObjects()
{
  // Setup Output
  fOutput = new AliEmcalList();
  fOutput->SetOwner();
}

/**
 * Execute once for the first event to initialize the analysis. Similar to ExecOnce() in
 * AliAnalysisTaskEmcal
 */
void AliEmcalCorrectionComponent::ExecOnce()
{
}

/**
 * Run every event, where the user implements their main analysis. Similar to Run() in
 * AliAnalysisTaskEmcal
 */
Bool_t AliEmcalCorrectionComponent::Run()
{
  AliDebugStream(3) << ": fEventManager.UseEmbeddingEvent(): " << fEventManager.UseEmbeddingEvent() << ", "
           << "fEventManager.InputEvent(): " << fEventManager.InputEvent() << ", "
           << "fEventManager address: " << &fEventManager << "\n";
  if(fGetPassFromFileName)
    GetPass();
  
  return kTRUE;
}

/**
 * Notifying the user that the input data file has
 * changed and performing steps needed to be done.
 */
Bool_t AliEmcalCorrectionComponent::UserNotify()
{
  return kTRUE;
}

/**
 * Calculate \f$\phi\f$ and \f$\eta\f$ difference between a track (t) and a cluster (c). The
 * position of the track is obtained on the EMCAL surface
 * @param[in] t Track to check
 * @param[in] v Cluster to check
 * @param[out] phidiff Distance in \f$\phi\f$ between cluster and track
 * @param[out] etadiff Distance in \f$\eta\f$ between cluster and track
 */
void AliEmcalCorrectionComponent::GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
  phidiff = 999;
  etadiff = 999;
  
  if (!t||!v) return;
  
  Double_t veta = t->GetTrackEtaOnEMCal();
  Double_t vphi = t->GetTrackPhiOnEMCal();
  
  Float_t pos[3] = {0};
  v->GetPosition(pos);
  TVector3 cpos(pos);
  Double_t ceta     = cpos.Eta();
  Double_t cphi     = cpos.Phi();
  etadiff=veta-ceta;
  phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}

/**
 * Remove bad cells from the cell list
 * Recalibrate energy and time cells
 */
void AliEmcalCorrectionComponent::UpdateCells()
{
  
  if (!fEventManager.InputEvent()) return ;
  
  Int_t bunchCrossNo = fEventManager.InputEvent()->GetBunchCrossNumber();
  
  if (fRecoUtils)
    fRecoUtils->RecalibrateCells(fCaloCells, bunchCrossNo);
  
  fCaloCells->Sort();
}

/**
 * Check whether the run changed.
 */
Bool_t AliEmcalCorrectionComponent::CheckIfRunChanged()
{
  // Get run number.
  Bool_t runChanged = fRun != fEventManager.InputEvent()->GetRunNumber();
  
  if (runChanged) {
    fRun = fEventManager.InputEvent()->GetRunNumber();
    AliWarning(Form("Run changed, initializing parameters for %d", fRun));
    
    // init geometry if not already done
    fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(fRun);
    if (!fGeom)
    {
      AliFatal("Can not create geometry");
      return kFALSE;
    }
  }
  return runChanged;
}

/**
 * Check if value is a shared parameter, meaning we should look
 * at another node. Also edits the input string to remove "sharedParameters:"
 * if it exists, making it ready for use.
 *
 * @param[in] value String containing the string value return by the parameter.
 *
 * @return True if the value is shared.
 */
bool AliEmcalCorrectionComponent::IsSharedValue(std::string & value)
{
  std::size_t sharedParameterLocation = value.find("sharedParameters:");
  if (sharedParameterLocation != std::string::npos)
  {
    // "sharedParameters:" is 17 characters long
    value.erase(sharedParameterLocation, sharedParameterLocation + 17);
    return true;
  }
  // Return false otherwise
  return false;
}

/**
 * Get pass from filename. Sets pass in fFilepass.
 */
void AliEmcalCorrectionComponent::GetPass()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TTree *inputTree = mgr->GetTree();
  
  if (!inputTree)
  {
    AliError("Pointer to tree = 0, returning");
    return;
  }
  
  TFile *inputFile = inputTree->GetCurrentFile();
  if (!inputFile) {
    AliError("Null pointer input file, returning");
    return;
  }
  
  TString fname(inputFile->GetName());
  if      (fname.Contains("pass1_pidfix"))                fFilepass = TString("pass1_pidfix");
  else if (fname.Contains("pass3_lowIR_pidfix"))          fFilepass = TString("pass3_lowIR_pidfix");
  else if (fname.Contains("pass4_lowIR_pidfix_cookdedx")) fFilepass = TString("pass4_lowIR_pidfix_cookdedx");
  else if (fname.Contains("pass1")) fFilepass = TString("pass1");
  else if (fname.Contains("pass2")) fFilepass = TString("pass2");
  else if (fname.Contains("pass3")) fFilepass = TString("pass3");
  else if (fname.Contains("pass4")) fFilepass = TString("pass4");
  else if (fname.Contains("pass5")) fFilepass = TString("pass5");
  else if (fname.Contains("LHC11c") && fname.Contains("spc_calo")) fFilepass = TString("spc_calo");
  else if (fname.Contains("calo") || fname.Contains("high_lumi"))
  {
    Printf("%s: Path contains <calo> or <high-lumi>, set as <pass1>", GetName());
    fFilepass = TString("pass1");
  }
  else if (fname.Contains("LHC14a1a"))
  {
    AliInfo("Energy calibration activated for this MC production!");
    fFilepass = TString("LHC14a1a");
  }
  else
  {
    AliFatal(Form("Pass number string not found: %s. Please set the pass number in the configuration!", fname.Data()));
    return;
  }
}

/**
 * Fills the Cell QA histograms
 */
void AliEmcalCorrectionComponent::FillCellQA(TH1F* h){
  TString name = h->GetName();
  
  Short_t  absId  =-1;
  Double_t ecell = 0;
  Double_t tcell = 0;
  Double_t efrac = 0;
  Int_t  mclabel = -1;
  
  for (Int_t iCell = 0; iCell < fCaloCells->GetNumberOfCells(); iCell++){
    
    fCaloCells->GetCell(iCell, absId, ecell, tcell, mclabel, efrac);
    if(name.Contains("Energy")){
      h->Fill(ecell);
    }
    else if(name.Contains("Time")){
      h->Fill(tcell);
    }
    
  }
  
}

/**
 * Initialize the bad channel map.
 */
Int_t AliEmcalCorrectionComponent::InitBadChannels()
{
  if (!fEventManager.InputEvent())
    return 0;
  
  AliInfo("Initialising Bad channel map");
  
  // init default maps first
  if (!fRecoUtils->GetEMCALBadChannelStatusMapArray())
    fRecoUtils->InitEMCALBadChannelStatusMap() ;
  
  Int_t runBC = fEventManager.InputEvent()->GetRunNumber();
  
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

