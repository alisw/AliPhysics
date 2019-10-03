#include "AliFMDMCTrackInspectorTask.h"
#include "AliForwardCorrectionManager.h"
#include "AliAODForwardMult.h"
#include <TROOT.h>
#include <TTree.h>

//____________________________________________________________________
AliFMDMCTrackInspectorTask::AliFMDMCTrackInspectorTask(const char* name, 
						       Bool_t      useTree)
  : AliBaseESDTask(name, "", &(AliForwardCorrectionManager::Instance())),
    fEventInspector("event"), 
    fTrackInspector("fitter")
{
  fCloneList = true;
  fTrackInspector.GetTracker().SetUseTree(useTree);
  if (useTree) DefineOutput(3, TTree::Class());
}

//____________________________________________________________________
TAxis*
AliFMDMCTrackInspectorTask::DefaultEtaAxis() const
{
  static TAxis* a = new TAxis(0, 0, 0);
  return a;
}
//____________________________________________________________________
TAxis*
AliFMDMCTrackInspectorTask::DefaultVertexAxis() const
{
  static TAxis* a = new TAxis(10, -10, 10);
  return a;
}

//____________________________________________________________________
Bool_t 
AliFMDMCTrackInspectorTask::Setup()
{
  DGUARD(fDebug,1,"Setting up the MC hit energy loss task");
  fTrackInspector.Init(); 
  return true; 
}

//____________________________________________________________________
Bool_t 
AliFMDMCTrackInspectorTask::Book()
{
  DGUARD(fDebug,1,"Booking histograms for the MC hit energy loss task");
  fNeededCorrections = 0;
  fExtraCorrections  = 0;

  fTrackInspector.CreateOutputObjects(fList);  
  return true;
}

//____________________________________________________________________
Bool_t 
AliFMDMCTrackInspectorTask::PreData(const TAxis& /*ipz*/, const TAxis& eta)
{
  DGUARD(fDebug,2,"Final setup of the MC hit energy loss task");
  fTrackInspector.SetupForData(eta);
  if (fTrackInspector.GetTracker().GetTree()) 
    PostData(3, fTrackInspector.GetTracker().GetTree());
  return true;
}
//____________________________________________________________________
Bool_t 
AliFMDMCTrackInspectorTask::Event(AliESDEvent& esd)
{
  // --- Read in the data --------------------------------------------
  LoadBranches();

  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = 0;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(&esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  if (found & AliFMDEventInspector::kNoEvent)    return false;
  if (found & AliFMDEventInspector::kNoTriggers) return false;
  if (found & AliFMDEventInspector::kNoSPD)      return false;
  if (found & AliFMDEventInspector::kNoFMD)      return false;
  if (found & AliFMDEventInspector::kNoVertex)   return false;
  if (found & AliFMDEventInspector::kBadVertex)  return false;

  // do not process pile-up, A, C, and E events 
  if (triggers & AliAODForwardMult::kPileUp)     return false;
  if (triggers & AliAODForwardMult::kA)          return false;
  if (triggers & AliAODForwardMult::kC)          return false;
  if (triggers & AliAODForwardMult::kE)          return false;
  
  // We want only the events found by off-line 
  if (!(triggers & AliAODForwardMult::kOffline)) return false;

  // Perhaps we should also insist on MB only 
  // if (fOnlyMB && (!(triggers & AliAODForwardMult::kInel))) return false;

  AliMCEvent* mc = MCEvent();
  if (!mc) return false;
  
  Bool_t ret = fTrackInspector.Event(esd, *mc, cent);

  if (fTrackInspector.GetTracker().GetTree()) 
    PostData(3, fTrackInspector.GetTracker().GetTree());

  return ret;
}
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)

//____________________________________________________________________
void   
AliFMDMCTrackInspectorTask::Print(Option_t* option) const
{
  AliBaseESDTask::Print(option);
  gROOT->IncreaseDirLevel();
  fTrackInspector.Print(option);
  gROOT->DecreaseDirLevel();
}

//
// EOF
// 


