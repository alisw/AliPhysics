#include "AliMultDepSpecAnalysisTaskUE.h"
/// \cond CLASSIMP
ClassImp(AliMultDepSpecAnalysisTaskUE);
/// \endcond

using std::string;
using std::vector;
using std::array;

 //****************************************************************************************
 /**
  * ROOT I/O Constructor.
  */
 //****************************************************************************************
AliMultDepSpecAnalysisTaskUE::AliMultDepSpecAnalysisTaskUE() : AliMultDepSpecAnalysisTask(),
fIsUE(true),
fPtLeadCut(0.),
// Histograms
///fAxes(),
fHistLeadPt(),
fHistLeadPhi(),
fHistPtLeadCutLoss(),
fHistMCResoPtLead(),
fHistMCResoPhiLead(),
fHistDiffToMCPtLead(),
fHistDiffToMCPhiLead(),
fHistPlateau(),
// Event properties
fPtLead(0),
fPhiLead(0),
fMCPtLead(0),
fMCPhiLead(0),
fMCPtOfLead(0),
fMCPhiOfLead(0)
{
  // ROOT IO constructor, don't allocate memory here!
}

//****************************************************************************************
/**
 * Constructor.
 */
//****************************************************************************************
AliMultDepSpecAnalysisTaskUE::AliMultDepSpecAnalysisTaskUE(const char* name) : AliMultDepSpecAnalysisTask(name),
fIsUE(true),
fPtLeadCut(0.),
// Histograms
///fAxes(),
fHistLeadPt(),
fHistLeadPhi(),
fHistPtLeadCutLoss(),
fHistMCResoPtLead(),
fHistMCResoPhiLead(),
fHistDiffToMCPtLead(),
fHistDiffToMCPhiLead(),
fHistPlateau(),
// Event properties
fPtLead(0),
fPhiLead(0),
fMCPtLead(0),
fMCPhiLead(0),
fMCPtOfLead(0),
fMCPhiOfLead(0)
{
}

//****************************************************************************************
/**
 * Destructor
 */
//****************************************************************************************
AliMultDepSpecAnalysisTaskUE::~AliMultDepSpecAnalysisTaskUE(){
  //AliMultDepSpecAnalysisTask::~AliMultDepSpecAnalysisTask();
}

//****************************************************************************************
/**
 * Define default axis properties .
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::DefineDefaultAxes(int maxMult)
{
  // add dimensions used in base class
  AliMultDepSpecAnalysisTask::DefineDefaultAxes(maxMult);

  // additional dimensions for UE studies
  std::vector<double> ptBins = {0.1, 1.0, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0};

  SetAxis(pt_lead_meas, "pt_lead_meas", "#it{p}_{T,Lead} (GeV/#it{c})", ptBins);
  SetAxis(phi_lead_meas, "phi_lead_meas", "#it{#it{#phi}_{Lead}}", {0, 2*TMath::Pi()}, 180);
  SetAxis(delta_pt_lead, "delta_pt_lead", "#Delta(#it{p}_{T,Lead})", {0., 2.}, 20);
  SetAxis(delta_phi_lead, "delta_phi_lead", "#Delta(#it{#phi}_{Lead})", {0, 2*TMath::Pi()}, 180);

}

//****************************************************************************************
/**
 * Function executed once before the event loop. Create histograms here.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::BookHistograms()
{
  // book all default histograms
  AliMultDepSpecAnalysisTask::BookHistograms();

  // book UE specific histograms
  BookHistogram(fHistLeadPt, "fHistLeadPt", {pt_lead_meas});
  BookHistogram(fHistLeadPhi, "fHistLeadPhi", {phi_lead_meas});
  BookHistogram(fHistPtLeadCutLoss, "fHistPtLeadCutLoss", {pt_meas});
  BookHistogram(fHistPlateau, "fHistPlateau", {pt_lead_meas, mult_meas});

  if(fIsMC)
  {
    BookHistogram(fHistMCResoPtLead, "fHistMCResoPtLead", {delta_pt_lead});
    BookHistogram(fHistMCResoPhiLead, "fHistMCResoPhiLead", {delta_phi_lead});
    BookHistogram(fHistDiffToMCPtLead, "fHistDiffToMCPtLead", {delta_pt_lead});
    BookHistogram(fHistDiffToMCPhiLead, "fHistDiffToMCPhiLead", {delta_phi_lead});
  }

  // check additional required memory
  double requiredMemory =
    fHistLeadPt.GetSize() +
    fHistLeadPhi.GetSize() +
    fHistPtLeadCutLoss.GetSize() +
    fHistPlateau.GetSize() +
    fHistMCResoPtLead.GetSize() +
    fHistMCResoPhiLead.GetSize() +
    fHistDiffToMCPtLead.GetSize() +
    fHistDiffToMCPhiLead.GetSize();

  AliError(Form("Estimated additional memory usage of histograms from UE analysis: %.2f MiB.", requiredMemory/1048576));

}

//****************************************************************************************
/**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitEvent()
{
  fEvent = InputEvent();
  if (!fEvent) {AliError("fEvent not available\n"); return false;}
  if(fIsMC){
    fMCEvent = MCEvent();
    if (!fMCEvent) {AliError("fMCEvent not available\n"); return false;}
  }
  FindLeadingTrack(); // sets fPhiLead, fMCPhilead, fPtLead, fMCPtLead
  if (fPtLead < fPtLeadCut){
    if (fMultMeas != 0) fHistPtLeadCutLoss.Fill(fPtLead);
    return false;
  }
  if (fIsMC && (fPtLead < fPtLeadCut)){
    if (fMultMeas != 0) fHistPtLeadCutLoss.Fill(fPtLead);
    return false;
  }
  return AliMultDepSpecAnalysisTask::InitEvent();
}

//****************************************************************************************
/**
 * Finds the leading track.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FindLeadingTrack()
{
  fPtLead = 0; fPhiLead = 0;
  if (fIsMC) {fMCPtLead = 0; fMCPhiLead = 0; fMCPtOfLead = 0; fMCPhiOfLead = 0;}

  AliVTrack* track = nullptr;
  for (int i = 0; i < fEvent->GetNumberOfTracks(); i++){
    track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(i));
    // Set fPt, fEta, fSigmapt; Check if track in kin range and has good quality

    if(!InitTrack(track)) continue;
    // initialize particle properties
    if(fIsMC)
    {
      // mc lable corresponding to measured track
      // negative lable indicates bad quality track - use it anyway
      int mcLable = TMath::Abs(track->GetLabel());

      // Set fMCPt, fMCEta, fMCIsChargedPrimary; Stop computation if it is not charged primary or secondary
      if(fIsESD){
        if(!InitParticle((AliMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }else{
        if(!InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(mcLable))) continue;
      }
    }
    // Find track with highest Pt
    if (fPt > fPtLead) {
      fPtLead = fPt;
      fPhiLead = fPhi;
      if (fIsMC) {
        fMCPtOfLead = fMCPt;
        fMCPhiOfLead = fMCPhi;
      }
    }

    if ((fIsMC) && (fMCPt > fMCPtLead)){
      fMCPtLead = fMCPt;
      fMCPhiLead = fMCPhi;
    }
  }
}

//****************************************************************************************
/**
 * Fill event histograms. Event related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillEventHistos()
{
  AliMultDepSpecAnalysisTask::FillEventHistos();

  if (fMultMeas == 0) return; // ensure that number of tracks is non-zero
  fHistLeadPt.Fill(fPtLead);
  fHistLeadPhi.Fill(fPhiLead);
  fHistPlateau.Fill(fPtLead, fMultMeas);

  if(fIsMC){
    // resolution of measurement:
    fHistMCResoPtLead.Fill(fabs(fPtLead - fMCPtOfLead));
    fHistMCResoPhiLead.Fill(fabs(fPhiLead - fMCPhiOfLead));

    // maybe one could also determine MC truth particle with highest pt and compare phi
    fHistDiffToMCPtLead.Fill(fabs(fMCPtLead - fPtLead));
    fHistDiffToMCPhiLead.Fill(fabs(fMCPhiLead - fPhiLead));
  }
}

//****************************************************************************************
/**
 * Fill track histograms. Track related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillMeasTrackHistos()
{
  AliMultDepSpecAnalysisTask::FillMeasTrackHistos();
}

//****************************************************************************************
/**
 * Fill measured particle histograms. Track and mc particle related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillMeasParticleHistos()
{
  AliMultDepSpecAnalysisTask::FillMeasParticleHistos();
}

//****************************************************************************************
/**
 * Fill generated particle histograms. MC particle related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillTrueParticleHistos()
{
  AliMultDepSpecAnalysisTask::FillTrueParticleHistos();
}

//****************************************************************************************
/**
 * Initializes track properties and returns false if track bad or out of range.
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitTrack(AliVTrack* track)
{
  bool isValidTrack = AliMultDepSpecAnalysisTask::InitTrack(track);
  // extract additional UE track infos
  // fDeltaLeadingPt = fLeadingPt - fPt;
  // fDeltaLeadingPhi = fLeadingPhi - fPhi;
  // fIsInTransverseRegion = ....

  return isValidTrack;
}

//****************************************************************************************
/**
 * Initializes particle properties and returns false if out of range.
 * Different functions are called for AliMCParticles (ESD) and AliAODMCParticles (AOD).
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitParticle(AliMCParticle* particle)
{
  bool isValidParticle = AliMultDepSpecAnalysisTask::InitParticle(particle);

  return isValidParticle;
}

bool AliMultDepSpecAnalysisTaskUE::InitParticle(AliAODMCParticle* particle)
{
  bool isValidParticle = AliMultDepSpecAnalysisTask::InitParticle(particle);

  return isValidParticle;
}

//****************************************************************************************
/**
 * User defined track selection criteria beyond good quality and fulfilling acceptance.that are applied before counting multiplicities and filling histograms (after InitTrack() was called).
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::SelectTrack()
{
  return fIsUE ? ((-0.5 <= TMath::Cos(fPhiLead-fPhi)) && (TMath::Cos(fPhiLead-fPhi) <= 0.5 )) : true;
}

//****************************************************************************************
/**
 * User defined particle selection criteria that are applied before counting multiplicities and filling histograms (after InitParticle() was called).
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::SelectParticle()
{
  return fIsUE ? ((-0.5 <= TMath::Cos(fMCPhiLead-fMCPhi)) && (TMath::Cos(fMCPhiLead-fMCPhi) <= 0.5 )) : true;
}

//****************************************************************************************
/**
 * Function to hang an instance of this task in a LEGO train.
 */
//****************************************************************************************
AliMultDepSpecAnalysisTaskUE* AliMultDepSpecAnalysisTaskUE::AddTaskMultDepSpecUE(const string& dataSet, int cutModeLow, int cutModeHigh, TString options, bool isMC, bool isUE, double ptLeadMIN)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMultDepSpecUE", "No analysis manager found.");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMultDepSpecUE", "No input event handler found.");
    return nullptr;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  bool isAOD = false;
  if(type.Contains("AOD")){
    isAOD = true;
  }
  else{
    // for ESDs isMC can be determined automatically
    isMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != nullptr);
  }
  string mode = (isMC) ? "MC" : "Data";
  string phiRange = (isUE) ? "UE" : "MB";

  AliMultDepSpecAnalysisTaskUE* returnTask = nullptr;
  char taskName[100] = "";
  for(int cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++){
    sprintf(taskName, "%s_%s_cutMode_%d_%s", dataSet.data(), mode.data(), cutMode, phiRange.data());
    AliMultDepSpecAnalysisTaskUE* task = new AliMultDepSpecAnalysisTaskUE(taskName);
    if(!task->InitTask(isUE, isMC, isAOD, dataSet, options, cutMode, ptLeadMIN))
    {
      delete task;
      task = nullptr;
      continue;
    }
    if(!returnTask) returnTask = task; // return one of the tasks
    task->SaveTrainMetadata();

    // hang task in train
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(taskName, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
  }
  return returnTask;
}

//****************************************************************************************
/**
 * Initialize task for specific cut mode. Creates ESD track cuts and particle compositon objects wit correct settings.
 * Call this funktion in the AddTask function. The resulting settings will be streamed.
 * In AODs the track selections are based on the best possible conversion of the specified esd cuts (see AliESDTrackCuts::AcceptVTrack()).
 * Due to the nature of the AOD format it is not possible to have an exact 1-to-1conversion for all of the cuts (e.g. geom length cut).
 * In addition, some cut variables are not available in early AOD productions: e.g. the golden chi2 cut can only be applied since August 2016.
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitTask(bool isUE, bool isMC, bool isAOD, string dataSet, TString options, int cutMode, float ptLeadMIN)
{
  SetIsUE(isUE);
  SetPtLeadCut(ptLeadMIN);
  return AliMultDepSpecAnalysisTask::InitTask(isMC, isAOD, dataSet, options, cutMode);
}
