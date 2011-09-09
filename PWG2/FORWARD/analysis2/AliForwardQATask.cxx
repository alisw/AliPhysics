// 
// Calculate the multiplicity in the forward regions event-by-event 
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODForwardMult 
// 
// Histograms 
//   
// Corrections used 
//
#include "AliForwardQATask.h"
#include "AliForwardUtil.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include "AliAODForwardMult.h"
#include <TH1.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStopwatch.h>

//====================================================================
AliForwardQATask::AliForwardQATask()
  : AliAnalysisTaskSE(),
    fEnableLowFlux(false), 
    fFirstEvent(true),
    fCorrManager(0),
    fESDFMD(),
    fHistos(),
    fEventInspector(),
    fEnergyFitter(),
    fSharingFilter(),
    fDensityCalculator(),
    fList(0)
{
  // 
  // Constructor
  //
}

//____________________________________________________________________
AliForwardQATask::AliForwardQATask(const char* name)
  : AliAnalysisTaskSE(name),
    fEnableLowFlux(false), 
    fFirstEvent(true),
    fCorrManager(0),
    fESDFMD(),
    fHistos(),
    fEventInspector("event"),
    fEnergyFitter("energy"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fList(0)
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  fCorrManager = &AliForwardCorrectionManager::Instance(); 
  fEnergyFitter.SetNParticles(1); // Just find the 1st peak 
  fEnergyFitter.SetDoMakeObject(false); 
  fEnergyFitter.SetUseIncreasingBins(true);
  fEnergyFitter.SetDoFits(kTRUE);
  fEnergyFitter.SetLowCut(0.4);
  fEnergyFitter.SetFitRangeBinWidth(4);
  fEnergyFitter.SetMinEntries(1000);
}

//____________________________________________________________________
AliForwardQATask::AliForwardQATask(const AliForwardQATask& o)
  : AliAnalysisTaskSE(o),
    fEnableLowFlux(o.fEnableLowFlux), 
    fFirstEvent(o.fFirstEvent),
    fCorrManager(o.fCorrManager),
    fESDFMD(o.fESDFMD),
    fHistos(o.fHistos),
    fEventInspector(o.fEventInspector),
    fEnergyFitter(o.fEnergyFitter),
    fSharingFilter(o.fSharingFilter),
    fDensityCalculator(o.fDensityCalculator),
    fList(o.fList) 
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________
AliForwardQATask&
AliForwardQATask::operator=(const AliForwardQATask& o)
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

  fEnableLowFlux     = o.fEnableLowFlux;
  fFirstEvent        = o.fFirstEvent;
  fCorrManager       = o.fCorrManager;
  fEventInspector    = o.fEventInspector;
  fEnergyFitter      = o.fEnergyFitter;
  fSharingFilter     = o.fSharingFilter;
  fDensityCalculator = o.fDensityCalculator;
  fHistos            = o.fHistos;
  fList              = o.fList;

  return *this;
}

//____________________________________________________________________
void
AliForwardQATask::SetDebug(Int_t dbg)
{
  // 
  // Set debug level 
  // 
  // Parameters:
  //    dbg Debug level
  //
  fEventInspector.SetDebug(dbg);
  fEnergyFitter.SetDebug(dbg);
  fSharingFilter.SetDebug(dbg);
  fDensityCalculator.SetDebug(dbg);
}

//____________________________________________________________________
Bool_t 
AliForwardQATask::CheckCorrections(UInt_t what) const
{
  // 
  // Check if all needed corrections are there and accounted for.  If not,
  // do a Fatal exit 
  // 
  // Parameters:
  //    what Which corrections is needed
  // 
  // Return:
  //    true if all present, false otherwise
  //  

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  // Check that we have the energy loss fits, needed by 
  //   AliFMDSharingFilter 
  //   AliFMDDensityCalculator 
  if (what & AliForwardCorrectionManager::kELossFits && !fcm.GetELossFit()) { 
    AliFatal(Form("No energy loss fits"));
    return false;
  }
  return true;
}

//____________________________________________________________________
Bool_t
AliForwardQATask::ReadCorrections(const TAxis*& pe, 
				  const TAxis*& pv, 
				  Bool_t        mc)
{
  //
  // Read corrections
  //
  //
  UInt_t what = AliForwardCorrectionManager::kAll;
  what ^= AliForwardCorrectionManager::kDoubleHit;
  what ^= AliForwardCorrectionManager::kVertexBias;
  what ^= AliForwardCorrectionManager::kAcceptance;
  what ^= AliForwardCorrectionManager::kMergingEfficiency;

  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  if (!fcm.Init(GetEventInspector().GetCollisionSystem(),
		GetEventInspector().GetEnergy(),
		GetEventInspector().GetField(),
		mc,
		what)) return false;
  if (!CheckCorrections(what)) return false;

  // Sett our persistency pointer 
  // fCorrManager = &fcm;

  // Get the eta axis from the secondary maps - if read in
  if (!pe) {
    pe = fcm.GetEtaAxis();
    if (!pe) AliFatal("No eta axis defined");
  }
  // Get the vertex axis from the secondary maps - if read in
  if (!pv) {
    pv = fcm.GetVertexAxis();
    if (!pv) AliFatal("No vertex axis defined");
  }

  return true;
}

//____________________________________________________________________
AliESDEvent*
AliForwardQATask::GetESDEvent()
{
  //
  // Get the ESD event. IF this is the first event, initialise
  //
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return 0;
  }

  // On the first event, initialize the parameters
  if (fFirstEvent && esd->GetESDRun()) {
    GetEventInspector().ReadRunDetails(esd);

    AliInfoF("Initializing with parameters from the ESD:\n"
	     "         AliESDEvent::GetBeamEnergy()   ->%f\n"
	     "         AliESDEvent::GetBeamType()     ->%s\n"
	     "         AliESDEvent::GetCurrentL3()    ->%f\n"
	     "         AliESDEvent::GetMagneticField()->%f\n"
	     "         AliESDEvent::GetRunNumber()    ->%d\n",
	     esd->GetBeamEnergy(),
	     esd->GetBeamType(),
	     esd->GetCurrentL3(),
	     esd->GetMagneticField(),
	     esd->GetRunNumber());

    fFirstEvent = false;

    if (!InitializeSubs()) {
      AliWarning("Initialisation of sub algorithms failed!");
      return 0;
    }
  }
  return esd;
}
//____________________________________________________________________
Bool_t
AliForwardQATask::InitializeSubs()
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  const TAxis* pe = 0;
  const TAxis* pv = 0;

  if (!ReadCorrections(pe,pv)) return false;

  fHistos.Init(*pe);

  fEventInspector.Init(*pv);
  fEnergyFitter.Init(*pe);
  fSharingFilter.Init();
  fDensityCalculator.Init(*pe);

  this->Print();

  return true;
}

//____________________________________________________________________
void
AliForwardQATask::UserCreateOutputObjects()
{
  // 
  // Create output objects 
  // 
  //
  fList = new TList;
  fList->SetOwner();
  
  fEventInspector.DefineOutput(fList);
  fEnergyFitter.DefineOutput(fList);
  fSharingFilter.DefineOutput(fList);
  fDensityCalculator.DefineOutput(fList);

  PostData(1, fList);
}
//____________________________________________________________________
void
AliForwardQATask::UserExec(Option_t*)
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
  AliESDEvent* esd = GetESDEvent();
  if (!esd) { 
    AliWarning("Got no ESD event");
    return;
  }

  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  Double_t vz        = 0;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(esd, triggers, lowFlux, 
					       ivz, vz, cent, nClusters);
  
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;
  if (found & AliFMDEventInspector::kNoSPD)      return;
  if (found & AliFMDEventInspector::kNoFMD)      return;
  if (found & AliFMDEventInspector::kNoVertex)   return;
  if (triggers & AliAODForwardMult::kPileUp)     return;
  if (found & AliFMDEventInspector::kBadVertex)  return;

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;

  // Get FMD data 
  AliESDFMD* esdFMD = esd->GetFMDData();
  
  // Run the energy loss fitter 
  if (!fEnergyFitter.Accumulate(*esdFMD, cent, 
				triggers & AliAODForwardMult::kEmpty)) {
    AliWarning("Energy fitter failed");
    return;
  }
    
  //  // Apply the sharing filter (or hit merging or clustering if you like)
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD)) { 
    AliWarning("Sharing filter failed!");
    return;
  }

  // Calculate the inclusive charged particle density 
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, ivz, lowFlux)) { 
    // if (!fDensityCalculator.Calculate(*esdFMD, fHistos, ivz, lowFlux)) { 
    AliWarning("Density calculator failed!");
    return;
  }
  PostData(1, fList);
}

//____________________________________________________________________
void
AliForwardQATask::Terminate(Option_t*)
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  TStopwatch swt;
  swt.Start();

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
    return;
  }

  TStopwatch swf;
  swf.Start();
  fEnergyFitter.Fit(list);
  swf.Stop();
  AliInfoF("Fitting took %d real-time seconds, and %f CPU seconds", 
	   Int_t(swf.RealTime()), swf.CpuTime());

  fSharingFilter.ScaleHistograms(list,Int_t(hEventsTr->Integral()));
  fDensityCalculator.ScaleHistograms(list,Int_t(hEventsTrVtx->Integral()));

  // Make a deep copy and post that as output 2 
  TList* list2 = static_cast<TList*>(list->Clone(Form("%sResults", 
						      list->GetName())));
  list2->SetOwner();
  PostData(2, list2);

  swt.Stop();
  AliInfoF("Terminate took %d real-time seconds, and %f CPU seconds", 
	   Int_t(swt.RealTime()), swt.CpuTime());

}

//____________________________________________________________________
void
AliForwardQATask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  
  std::cout << ClassName() << ": " << GetName() << "\n" 
	    << "  Enable low flux code:   " << (fEnableLowFlux ? "yes" : "no") 
	    << "\n"
	    << "  Off-line trigger mask:  0x" 
	    << std::hex     << std::setfill('0') 
	    << std::setw (8) << fOfflineTriggerMask 
	    << std::dec     << std::setfill (' ') << std::endl;
  gROOT->IncreaseDirLevel();
  if (fCorrManager) fCorrManager->Print();
  else  
    std::cout << "  Correction manager not set yet" << std::endl;
  GetEventInspector()   .Print(option);
  GetEnergyFitter()     .Print(option);
  GetSharingFilter()    .Print(option);
  gROOT->DecreaseDirLevel();
}

//
// EOF
//
