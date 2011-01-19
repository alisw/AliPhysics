// 
// Histogram and fit the energy loss distributions for the FMD
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - None
// 
// Histograms:
//   
// Corrections used:
//   - None
// 
// 
//
#include "AliFMDEnergyFitterTask.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODForwardMult.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TFile.h>

//====================================================================
AliFMDEnergyFitterTask::AliFMDEnergyFitterTask()
  : AliAnalysisTaskSE(),
    fFirstEvent(true),
    fEventInspector(),
    fEnergyFitter(),
    fList(0)
{
  // 
  // Constructor
  //
}

//____________________________________________________________________
AliFMDEnergyFitterTask::AliFMDEnergyFitterTask(const char* name)
  : AliAnalysisTaskSE(name), 
    fFirstEvent(true),
    fEventInspector("event"),
    fEnergyFitter("energy"),
    fList(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
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
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________
AliFMDEnergyFitterTask&
AliFMDEnergyFitterTask::operator=(const AliFMDEnergyFitterTask& o)
{
  // 
  // Assignment operator 
  // 
  // Parameters:
  //    o Object to assign from 
  // 
  // Return:
  //    Reference to this object 
  //
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
  // 
  // Set the debug level 
  // 
  // Parameters:
  //    dbg Debug level
  //
  fEventInspector.SetDebug(dbg);
  fEnergyFitter.SetDebug(dbg);
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::Init()
{
  // 
  // Initialize the task 
  // 
  //
  fFirstEvent = true;
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::InitializeSubs()
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
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
  // 
  // Create output objects 
  // 
  //
  fList = new TList;

  fEventInspector.DefineOutput(fList);
  fEnergyFitter.DefineOutput(fList);

  PostData(1, fList);
}
//____________________________________________________________________
void
AliFMDEnergyFitterTask::UserExec(Option_t*)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  

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
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  AliInfo(Form("Running terminate of %s", GetName()));
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
  AliInfo("Fitting energy loss spectra");
  fEnergyFitter.Fit(list);

  // Investigate output slot 
  AliAnalysisDataSlot* oslot = GetOutputSlot(1);
  if (oslot) { 
    AliAnalysisDataContainer* ocont = oslot->GetContainer();
    if (ocont) { 
      TFile* ofile = ocont->GetFile();
      if (ofile) {
	AliInfo(Form("Output file %s opened with option %s (%s)", 
		     ofile->GetName(), ofile->GetOption(), 
		     ofile->IsWritable() ? "read-write" : "read-only"));
	ofile->Flush();
      }
      else 
	AliWarning("No output file associated with data container");
    }
    else 
      AliWarning("No container associated with slot 1");
  }
  else 
    AliWarning("Slot number 1 not defined");
  
  // Temporary code to save output to a file - should be disabled once 
  // the plugin stuff works 
  list->ls();
  TDirectory* savdir = gDirectory;
  TFile* tmp = TFile::Open("elossfits.root", "RECREATE");
  list->Write(list->GetName(), TObject::kSingleKey);
  tmp->Write();
  tmp->Close();
  savdir->cd();
}

//____________________________________________________________________
void
AliFMDEnergyFitterTask::Print(Option_t*) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
}

//
// EOF
//
