#include "AliFMDEnergyFitterTask.h"
#include "AliLog.h"
// #include "AliFMDAnaParameters.h"
#include "AliESDEvent.h"
#include "AliAODForwardMult.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>

//====================================================================
AliFMDEnergyFitterTask::AliFMDEnergyFitterTask()
  : AliAnalysisTaskSE(),
    fFirstEvent(true),
    fEventInspector(),
    fEnergyFitter(),
    fList(0)
{
}

//____________________________________________________________________
AliFMDEnergyFitterTask::AliFMDEnergyFitterTask(const char* name)
  : AliAnalysisTaskSE(name), 
    fFirstEvent(true),
    fEventInspector("event"),
    fEnergyFitter("energy"),
    fList(0)
{
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________
AliFMDEnergyFitterTask::AliFMDEnergyFitterTask(const AliFMDEnergyFitterTask& o)
  : AliAnalysisTaskSE(o),
    fFirstEvent(o.fFirstEvent),
    fEventInspector(o.fEventInspector),
    fEnergyFitter(o.fEnergyFitter),
    fList(o.fList) 
{
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________
AliFMDEnergyFitterTask&
AliFMDEnergyFitterTask::operator=(const AliFMDEnergyFitterTask& o)
{
  AliAnalysisTaskSE::operator=(o);

  fFirstEvent        = o.fFirstEvent;
  fEventInspector    = o.fEventInspector;
  fEnergyFitter      = o.fEnergyFitter;
  fList              = o.fList;

  return *this;
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::SetDebug(Int_t dbg)
{
  fEventInspector.SetDebug(dbg);
  fEnergyFitter.SetDebug(dbg);
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::Init()
{
  fFirstEvent = true;
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::InitializeSubs()
{

  // AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  // pars->Init(kTRUE);

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  fcm.Init(fEventInspector.GetCollisionSystem(), 
	   fEventInspector.GetEnergy(),
	   fEventInspector.GetField(), 0);
  TAxis eAxis(0,0,0);
  TAxis vAxis(10,-10,10);
  fEnergyFitter.Init(eAxis);
  fEventInspector.Init(vAxis);
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::UserCreateOutputObjects()
{
  fList = new TList;

  fEventInspector.DefineOutput(fList);
  fEnergyFitter.DefineOutput(fList);

  PostData(1, fList);
}
//____________________________________________________________________
void
AliFMDEnergyFitterTask::UserExec(Option_t*)
{
  // static Int_t cnt = 0;
  // cnt++;
  // Get the input data 
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  // AliInfo(Form("Event # %6d (esd=%p)", cnt, esd));
  if (!esd) { 
    AliWarning("No ESD event found for input event");
    return;
  }

  // On the first event, initialize the parameters 
  if (fFirstEvent && esd->GetESDRun()) { 
    fEventInspector.ReadRunDetails(esd);
    
    AliInfo(Form("Initializing with parameters from the ESD:\n"
		 "         AliESDEvent::GetBeamEnergy()   ->%f\n"
		 "         AliESDEvent::GetBeamType()     ->%s\n"
		 "         AliESDEvent::GetCurrentL3()    ->%f\n"
		 "         AliESDEvent::GetMagneticField()->%f\n"
		 "         AliESDEvent::GetRunNumber()    ->%d\n",
		 esd->GetBeamEnergy(), 
		 esd->GetBeamType(),
		 esd->GetCurrentL3(), 
		 esd->GetMagneticField(),
		 esd->GetRunNumber()));

	      

    // AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
    // pars->SetParametersFromESD(esd);
    // pars->PrintStatus();
    fFirstEvent = false;

    InitializeSubs();
  }
  Bool_t   lowFlux  = kFALSE;
  UInt_t   triggers = 0;
  UShort_t ivz      = 0;
  Double_t vz       = 0;
  UInt_t   found    = fEventInspector.Process(esd, triggers, lowFlux, ivz, vz);
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;
  if (found & AliFMDEventInspector::kNoSPD)     return;
  if (found & AliFMDEventInspector::kNoFMD)     return;
  if (found & AliFMDEventInspector::kNoVertex)  return;
  if (found & AliFMDEventInspector::kBadVertex) return;

  // Get FMD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  // Do the energy stuff 
  if (!fEnergyFitter.Accumulate(*esdFMD, triggers & AliAODForwardMult::kEmpty)){
    AliWarning("Energy fitter failed");
    return;
  }
  PostData(1, fList);
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::Terminate(Option_t*)
{
  TList* list = dynamic_cast<TList*>(GetOutputData(1));
  if (!list) {
    AliError(Form("No output list defined (%p)", GetOutputData(1)));
    if (GetOutputData(1)) GetOutputData(1)->Print();
    return;
  }
  
  // Get our histograms from the container 
  TH1I* hEventsTr    = 0;//static_cast<TH1I*>(list->FindObject("nEventsTr"));
  TH1I* hEventsTrVtx = 0;//static_cast<TH1I*>(list->FindObject("nEventsTrVtx"));
  TH1I* hTriggers    = 0;
  if (!fEventInspector.FetchHistograms(list, hEventsTr, 
				       hEventsTrVtx, hTriggers)) { 
    AliError(Form("Didn't get histograms from event selector "
		  "(hEventsTr=%p,hEventsTrVtx=%p)", 
		  hEventsTr, hEventsTrVtx));
    list->ls();
    return;
  }
  fEnergyFitter.Fit(list);
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::Print(Option_t*) const
{
}

//
// EOF
//
