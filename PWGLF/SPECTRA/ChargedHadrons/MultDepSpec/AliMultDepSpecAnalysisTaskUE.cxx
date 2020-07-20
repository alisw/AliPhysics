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
AliMultDepSpecAnalysisTaskUE::AliMultDepSpecAnalysisTaskUE() : AliMultDepSpecAnalysisTask()
{
  // ROOT IO constructor, don't allocate memory here!
}

//****************************************************************************************
/**
 * Constructor.
 */
//****************************************************************************************
AliMultDepSpecAnalysisTaskUE::AliMultDepSpecAnalysisTaskUE(const char* name) : AliMultDepSpecAnalysisTask(name)
{
}

//****************************************************************************************
/**
 * Destructor
 */
//****************************************************************************************
AliMultDepSpecAnalysisTaskUE::~AliMultDepSpecAnalysisTaskUE(){
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
  
  // fHistLeadingPt...

  
}

//****************************************************************************************
/**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::InitEvent()
{
  // InitLeadingTrack(); // sets fLeadingPhi, fMCLeadingPhi, fLeadingPt, fMCLeadingPt...
  return AliMultDepSpecAnalysisTask::InitEvent();
}

//****************************************************************************************
/**
 * Fill event histograms. Event related members are set.
 */
//****************************************************************************************
void AliMultDepSpecAnalysisTaskUE::FillEventHistos()
{
  AliMultDepSpecAnalysisTask::FillEventHistos();
  
  // fHistLeadingPt.Fill(fLeadingPt);
  // resolution of measurement:
  // fHistMCLeadingDeltaPt.Fill(abs(fLeadingPt - fMCLeadingPt));
  // fHistMCLeadingDeltaPhi.Fill(abs(fLeadingPhi - fMCLeadingPhi));
  
  // maybe one could also determine MC truth particle with highest pt and compare phi
  // fHistMCLeadingTrackBias();
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
  //return fIsTrackInUE;
  return true;
}

//****************************************************************************************
/**
 * User defined particle selection criteria that are applied before counting multiplicities and filling histograms (after InitParticle() was called).
 */
//****************************************************************************************
bool AliMultDepSpecAnalysisTaskUE::SelectParticle()
{
  // return fIsParticleInUE;
  return true;
}

//****************************************************************************************
/**
 * Function to hang an instance of this task in a LEGO train.
 */
//****************************************************************************************
AliMultDepSpecAnalysisTaskUE* AddTaskMultDepSpecUE(const string& dataSet, int cutModeLow, int cutModeHigh, TString options, bool isMC)
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
  
  AliMultDepSpecAnalysisTaskUE* returnTask = nullptr;
  char taskName[100] = "";
  for(int cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++){
    sprintf(taskName, "%s_%s_cutMode_%d", dataSet.data(), mode.data(), cutMode);
    AliMultDepSpecAnalysisTaskUE* task = new AliMultDepSpecAnalysisTaskUE(taskName);
    if(!task->InitTask(isMC, isAOD, dataSet, options, cutMode))
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
