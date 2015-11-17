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
  : AliBaseESDTask(),
    fEnableLowFlux(false), 
    fESDFMD(),
    fHistos(),
    fEventInspector(),
    fESDFixer(),
    fEnergyFitter(),
    fSharingFilter(),
    fDensityCalculator()
{
  // 
  // Constructor
  //
  fCloneList = true;
}

//____________________________________________________________________
AliForwardQATask::AliForwardQATask(const char* name)
  : AliBaseESDTask(name, "", &(AliForwardCorrectionManager::Instance())),
    fEnableLowFlux(false), 
    fESDFMD(),
    fHistos(),
    fEventInspector("event"),
    fESDFixer("fixer"),
    fEnergyFitter("energy"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density")
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  fEnergyFitter.SetNParticles(1); // Just find the 1st peak 
  fEnergyFitter.SetDoMakeObject(false); 
  fEnergyFitter.SetUseIncreasingBins(true);
  fEnergyFitter.SetDoFits(kTRUE);
  fEnergyFitter.SetLowCut(0.4);
  fEnergyFitter.SetFitRangeBinWidth(4);
  fEnergyFitter.SetMinEntries(1000);
  fCloneList = true;

  // For the QA we always enable fall-back 
  AliForwardCorrectionManager::Instance().SetEnableFallBack(true);
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
  AliBaseESDTask::        SetDebug(dbg);
  GetEnergyFitter()     .SetDebug(dbg);
  GetSharingFilter()    .SetDebug(dbg);
  GetDensityCalculator().SetDebug(dbg);
}

//____________________________________________________________________
TAxis*
AliForwardQATask::DefaultEtaAxis() const
{
  static TAxis* a = new TAxis(240, -6, 6);
  return a;
}
//____________________________________________________________________
TAxis*
AliForwardQATask::DefaultVertexAxis() const
{
  static TAxis* a = new TAxis(10, -10, 10);
  return a;
}

//____________________________________________________________________
Bool_t
AliForwardQATask::Setup()
{
  fEnergyFitter.Init();
  return true;
}

//____________________________________________________________________
Bool_t
AliForwardQATask::Book()
{
  // 
  // Create output objects 
  // 
  //
  UInt_t what = AliForwardCorrectionManager::kAll;
  what ^= AliForwardCorrectionManager::kDoubleHit;
  what ^= AliForwardCorrectionManager::kVertexBias;
  what ^= AliForwardCorrectionManager::kAcceptance;
  what ^= AliForwardCorrectionManager::kMergingEfficiency;
  what ^= AliForwardCorrectionManager::kELossFits;
  fNeededCorrections = what;
  fExtraCorrections  = AliForwardCorrectionManager::kELossFits;
  
  fESDFixer         .CreateOutputObjects(fList);
  fEnergyFitter     .CreateOutputObjects(fList);
  fSharingFilter    .CreateOutputObjects(fList);
  fDensityCalculator.CreateOutputObjects(fList);
  
  return true;
}
//____________________________________________________________________
Bool_t
AliForwardQATask::PreData(const TAxis& /*vertex*/, const TAxis& eta)
{
  // 
  // Initialise the sub objects and stuff.  Called on first event 
  // 
  //
  AliForwardCorrectionManager& fcm = AliForwardCorrectionManager::Instance();
  // We allow fall-back queries so that we may proceed in case we have no 
  // valid corrections 
  if (!fcm.GetELossFit()) { 
    AliWarning("No energy loss fits");
    
    // Fall-back values if we do not have the energy loss fits 
    AliFMDMultCuts& sfLCuts = GetSharingFilter().GetLCuts();
    if (sfLCuts.GetMethod() != AliFMDMultCuts::kFixed) { 
      Double_t cut = 0.15;
      AliWarningF("Using fixed cut @ %f for the lower bound "
		  "of the sharing filter", cut);
      sfLCuts.SetMultCuts(cut);
    }
    AliFMDMultCuts& sfHCuts = GetSharingFilter().GetHCuts();
    if (sfHCuts.GetMethod() != AliFMDMultCuts::kFixed) { 
      Double_t cut = 0.45;
      AliWarningF("Using fixed cut @ %f for the upper bound "
		  "of the sharing filter", cut);
      sfHCuts.SetMultCuts(cut);
    }
    AliFMDMultCuts& dcCuts  = GetDensityCalculator().GetCuts();
    if (dcCuts.GetMethod() != AliFMDMultCuts::kFixed) { 
      Double_t cut = 0.45;
      AliWarningF("Using fixed cut @ %f for the lower bound "
		  "of the density calculator", cut);
      dcCuts.SetMultCuts(cut);
    }
  }
  else 
    fcm.GetELossFit()->CacheBins(GetDensityCalculator().GetMinQuality());

  fHistos.Init(eta);

  // GetEventInspector().SetupForData(vertex);
  GetEnergyFitter()     .SetupForData(eta,
				      GetEventInspector().GetCollisionSystem());
  GetSharingFilter()    .SetupForData(eta);
  GetDensityCalculator().SetupForData(eta);

  return true;
}

//____________________________________________________________________
Bool_t
AliForwardQATask::PreEvent()
{
  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  return true;
}
//____________________________________________________________________
Bool_t
AliForwardQATask::Event(AliESDEvent& esd)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  DGUARD(fDebug,1,"Process the input event");

  if (fFirstEvent) { 
    // If the first event flag wasn't cleared in the above call to
    // GetESDEvent, we should not do anything, since nothing has been
    // initialised yet, so we opt out here (with a warning) 
    AliWarning("Nothing has been initialized yet, opt'ing out");
    return false;
  }

  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(&esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  
  Bool_t ok = true;
  if (found & AliFMDEventInspector::kNoEvent)    ok = false;
  if (found & AliFMDEventInspector::kNoTriggers) ok = false;
  if (found & AliFMDEventInspector::kNoSPD)      ok = false;
  if (found & AliFMDEventInspector::kNoFMD)      ok = false;
  if (found & AliFMDEventInspector::kNoVertex)   ok = false;
  if (triggers & AliAODForwardMult::kPileUp)     ok = false;
  if (triggers & AliAODForwardMult::kA)          ok = false;
  if (triggers & AliAODForwardMult::kC)          ok = false;
  if (triggers & AliAODForwardMult::kE)          ok = false;
  if (!(triggers & AliAODForwardMult::kOffline)) ok = false;
  if (found & AliFMDEventInspector::kBadVertex)  ok = false;
  if (!ok) { 
    DMSG(fDebug,2,"Event failed selection: %s", 
	 AliFMDEventInspector::CodeString(found));
    return false;
  }
  DMSG(fDebug,2,"Event triggers: %s", 
       AliAODForwardMult::GetTriggerString(triggers));

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;

  // Get FMD data 
  AliESDFMD* esdFMD = esd.GetFMDData();

  // Fix up the the ESD 
  GetESDFixer().Fix(*esdFMD, ip);
  
  // Run the energy loss fitter 
  if (!fEnergyFitter.Accumulate(*esdFMD, cent, 
				triggers & AliAODForwardMult::kEmpty)) {
    AliWarning("Energy fitter failed");
    return false;
  }
  
  //  // Apply the sharing filter (or hit merging or clustering if you like)
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD, ip.Z())) { 
    AliWarning("Sharing filter failed!");
    return false;
  }
 
  // Calculate the inclusive charged particle density 
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, lowFlux, cent, ip)) { 
    // if (!fDensityCalculator.Calculate(*esdFMD, fHistos, ivz, lowFlux)) { 
    AliWarning("Density calculator failed!");
    return false;
  }
  
  return true;
}

//____________________________________________________________________
Bool_t
AliForwardQATask::Finalize()
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  if (fDebug) AliInfo("In Forwards terminate");
  TStopwatch swt;
  swt.Start();

  // Get our histograms from the container 
  TH1I* hEventsTr    = 0;
  TH1I* hEventsTrVtx = 0;
  TH1I* hEventsAcc   = 0;
  TH1I* hTriggers    = 0;
  if (!fEventInspector.FetchHistograms(fList, 
				       hEventsTr, 
				       hEventsTrVtx, 
				       hEventsAcc,
				       hTriggers)) { 
    AliErrorF("Didn't get histograms from event selector "
	      "(hEventsTr=%p,hEventsTrVtx=%p,hEventsAcc=%p)", 
	      hEventsTr, hEventsTrVtx,hEventsAcc);
    return false;
  }

  TStopwatch swf;
  swf.Start();
  fEnergyFitter.Fit(fResults);
  swf.Stop();
  AliInfoF("Fitting took %d real-time seconds, and %f CPU seconds", 
	   Int_t(swf.RealTime()), swf.CpuTime());

  fSharingFilter.Terminate(fList,fResults,Int_t(hEventsTr->Integral()));
  fDensityCalculator.Terminate(fList,fResults,Int_t(hEventsTrVtx->Integral()));

  if (fDebug) AliInfoF("Posting post processing results to %s", 
		       fResults->GetName());
  swt.Stop();
  AliInfoF("Finalize took %d real-time seconds, and %f CPU seconds", 
	   Int_t(swt.RealTime()), swt.CpuTime());

  return true;
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
  AliBaseESDTask::Print(option);
  gROOT->IncreaseDirLevel();
  GetESDFixer()         .Print(option);        
  GetEnergyFitter()     .Print(option);
  GetSharingFilter()    .Print(option);
  GetDensityCalculator().Print(option);
  gROOT->DecreaseDirLevel();
}

//
// EOF
//
