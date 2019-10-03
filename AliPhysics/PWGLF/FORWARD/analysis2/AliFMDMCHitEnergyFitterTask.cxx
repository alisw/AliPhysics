#include "AliFMDMCHitEnergyFitterTask.h"
#include "AliForwardCorrectionManager.h"
#include "AliAODForwardMult.h"
#include "AliMCAuxHandler.h"
#include <TROOT.h>
#include <TNtuple.h>

//____________________________________________________________________
AliFMDMCHitEnergyFitterTask::AliFMDMCHitEnergyFitterTask(const char* name, 
							 Bool_t      useTuple)
  : AliBaseESDTask(name, "", &(AliForwardCorrectionManager::Instance())),
    fEventInspector("event"), 
    fEnergyFitter("fitter", useTuple), 
    fHitHandler(0)
{
  fCloneList = true;
  if (useTuple) DefineOutput(3, TNtuple::Class());
}

//____________________________________________________________________
TAxis*
AliFMDMCHitEnergyFitterTask::DefaultEtaAxis() const
{
  static TAxis* a = new TAxis(0, 0, 0);
  return a;
}
//____________________________________________________________________
TAxis*
AliFMDMCHitEnergyFitterTask::DefaultVertexAxis() const
{
  static TAxis* a = new TAxis(10, -10, 10);
  return a;
}

//____________________________________________________________________
Bool_t 
AliFMDMCHitEnergyFitterTask::Setup()
{
  DGUARD(fDebug,1,"Setting up the MC hit energy loss task");
  fEnergyFitter.Init(); 
  return true; 
}

//____________________________________________________________________
Bool_t 
AliFMDMCHitEnergyFitterTask::Book()
{
  DGUARD(fDebug,1,"Booking histograms for the MC hit energy loss task");
  fNeededCorrections = 0;
  fExtraCorrections  = 0;

  fEnergyFitter.CreateOutputObjects(fList);  
  // We have to add the handler here, since the subsiduary handlers of
  // AliMCEventHandler are not streamed with the parent object - sigh!
  // Hence, we add it at init time of the slaves. 
  fHitHandler = AliMCAuxHandler::Create("FMD", "AliFMDHit");
  return true;
}

//____________________________________________________________________
Bool_t 
AliFMDMCHitEnergyFitterTask::PreData(const TAxis& /*ipz*/, const TAxis& eta)
{
  DGUARD(fDebug,2,"Final setup of the MC hit energy loss task");
  fEnergyFitter.SetupForData(eta);
  if (fEnergyFitter.GetTuple()) 
    PostData(3, fEnergyFitter.GetTuple());
  return true;
}
//____________________________________________________________________
Bool_t 
AliFMDMCHitEnergyFitterTask::Event(AliESDEvent& esd)
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
  
  Bool_t ret = fEnergyFitter.Event(esd, *mc, *fHitHandler);

  if (fEnergyFitter.GetTuple()) 
    PostData(3, fEnergyFitter.GetTuple());

  return ret;
}
#define PFB(N,FLAG)				\
  do {									\
    AliForwardUtil::PrintName(N);					\
    std::cout << std::boolalpha << (FLAG) << std::noboolalpha << std::endl; \
  } while(false)

//____________________________________________________________________
void   
AliFMDMCHitEnergyFitterTask::Print(Option_t* option) const
{
  AliBaseESDTask::Print(option);
  gROOT->IncreaseDirLevel();
  fEnergyFitter.Print(option);
  gROOT->DecreaseDirLevel();
}

//
// EOF
// 


