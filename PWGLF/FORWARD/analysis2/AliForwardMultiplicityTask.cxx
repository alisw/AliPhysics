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
#include "AliForwardMultiplicityTask.h"
#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODHandler.h"
#include "AliMultiplicity.h"
#include "AliInputEventHandler.h"
#include "AliForwardCorrectionManager.h"
#include "AliAnalysisManager.h"
#include <TH1.h>
#include <TH3D.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TProfile.h>

//====================================================================
AliForwardMultiplicityTask::AliForwardMultiplicityTask()
  : AliForwardMultiplicityBase(),
    fESDFMD(),
    fEventInspector(),
    fSharingFilter(),
    fDensityCalculator(),
    fCorrections(),
    fHistCollector(),
    fEventPlaneFinder()
{
  // 
  // Constructor
  //
  DGUARD(fDebug, 3,"Default CTOR of AliForwardMultiplicityTask");
}

//____________________________________________________________________
AliForwardMultiplicityTask::AliForwardMultiplicityTask(const char* name)
  : AliForwardMultiplicityBase(name),
    fESDFMD(),
    fEventInspector("event"),
    fSharingFilter("sharing"), 
    fDensityCalculator("density"),
    fCorrections("corrections"),
    fHistCollector("collector"),
    fEventPlaneFinder("eventplane")
{
  // 
  // Constructor 
  // 
  // Parameters:
  //    name Name of task 
  //
  DGUARD(fDebug, 3,"named CTOR of AliForwardMultiplicityTask: %s", name);
}

//____________________________________________________________________
AliForwardMultiplicityTask::AliForwardMultiplicityTask(const AliForwardMultiplicityTask& o)
  : AliForwardMultiplicityBase(o),
    fESDFMD(o.fESDFMD),
    fEventInspector(o.fEventInspector),
    fSharingFilter(o.fSharingFilter),
    fDensityCalculator(o.fDensityCalculator),
    fCorrections(o.fCorrections),
    fHistCollector(o.fHistCollector),
    fEventPlaneFinder(o.fEventPlaneFinder)

{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DGUARD(fDebug, 3,"Copy CTOR of AliForwardMultiplicityTask");
}

//____________________________________________________________________
AliForwardMultiplicityTask&
AliForwardMultiplicityTask::operator=(const AliForwardMultiplicityTask& o)
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
  DGUARD(fDebug,3,"Assignment to AliForwardMultiplicityTask");
  if (&o == this) return *this;
  AliForwardMultiplicityBase::operator=(o);

  fEventInspector    = o.fEventInspector;
  fSharingFilter     = o.fSharingFilter;
  fDensityCalculator = o.fDensityCalculator;
  fCorrections       = o.fCorrections;
  fHistCollector     = o.fHistCollector;
  fEventPlaneFinder  = o.fEventPlaneFinder;

  return *this;
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::UserExec(Option_t*)
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  TStopwatch total;
  TStopwatch individual;
  if (fDoTiming) total.Start(true);
  
  DGUARD(fDebug,1,"Process the input event");
  // static Int_t cnt = 0;
  // cnt++;
  // Get the input data 
  AliESDEvent* esd = GetESDEvent();
  if (!esd) return;

  // Clear stuff 
  fHistos.Clear();
  fESDFMD.Clear();
  fAODFMD.Clear();
  fAODEP.Clear();

  // Inspect the event
  if (fDoTiming) individual.Start(true);
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fEventInspector.Process(esd, triggers, lowFlux, 
					       ivz, ip, cent, nClusters);
  if (fDoTiming) fHTiming->Fill(kTimingEventInspector, individual.CpuTime());
  
  if (found & AliFMDEventInspector::kNoEvent)    return;
  if (found & AliFMDEventInspector::kNoTriggers) return;

  // Set trigger bits, and mark this event for storage 
  fAODFMD.SetTriggerBits(triggers);
  fAODFMD.SetSNN(fEventInspector.GetEnergy());
  fAODFMD.SetSystem(fEventInspector.GetCollisionSystem());
  fAODFMD.SetCentrality(cent);
  fAODFMD.SetNClusters(nClusters);
  MarkEventForStore();
 
  if (found & AliFMDEventInspector::kNoSPD)      return;
  if (found & AliFMDEventInspector::kNoFMD)      return;
  if (found & AliFMDEventInspector::kNoVertex)   return;
  
  if (triggers & AliAODForwardMult::kPileUp) return;
  
  fAODFMD.SetIpZ(ip.Z());

  if (found & AliFMDEventInspector::kBadVertex) return;

  // We we do not want to use low flux specific code, we disable it here. 
  if (!fEnableLowFlux) lowFlux = false;

  // Get FMD data 
  AliESDFMD* esdFMD = esd->GetFMDData();  

  // Apply the sharing filter (or hit merging or clustering if you like)
  if (fDoTiming) individual.Start(true);
  if (!fSharingFilter.Filter(*esdFMD, lowFlux, fESDFMD, ip.Z())) { 
    AliWarning("Sharing filter failed!");
    return;
  }
  if (fDoTiming) fHTiming->Fill(kTimingSharingFilter, individual.CpuTime());
  
  // Calculate the inclusive charged particle density 
  if (fDoTiming) individual.Start(true);
  if (!fDensityCalculator.Calculate(fESDFMD, fHistos, lowFlux, cent, ip)) { 
    // if (!fDensityCalculator.Calculate(*esdFMD, fHistos, ivz, lowFlux)) { 
    AliWarning("Density calculator failed!");
    return;
  }
  if (fDoTiming) fHTiming->Fill(kTimingDensityCalculator,individual.CpuTime());

  // Check if we should do the event plane finder
  if (fEventInspector.GetCollisionSystem() == AliFMDEventInspector::kPbPb) {
    if (fDoTiming) individual.Start(true);
    if (!fEventPlaneFinder.FindEventplane(esd, fAODEP, 
					  &(fAODFMD.GetHistogram()), &fHistos))
      AliWarning("Eventplane finder failed!");
    if (fDoTiming) fHTiming->Fill(kTimingEventPlaneFinder,individual.CpuTime());
  }
  
  // Do the secondary and other corrections. 
  if (fDoTiming) individual.Start(true);
  if (!fCorrections.Correct(fHistos, ivz)) { 
    AliWarning("Corrections failed");
    return;
  }
  if (fDoTiming) fHTiming->Fill(kTimingCorrections, individual.CpuTime());

  // Collect our `super' histogram 
  if (fDoTiming) individual.Start(true);
  if (!fHistCollector.Collect(fHistos, fRingSums, 
			      ivz, fAODFMD.GetHistogram(),
			      fAODFMD.GetCentrality())) {
    AliWarning("Histogram collector failed");
    return;
  }
  if (fDoTiming) fHTiming->Fill(kTimingHistCollector, individual.CpuTime());

  if (fAODFMD.IsTriggerBits(AliAODForwardMult::kInel))
    fHData->Add(&(fAODFMD.GetHistogram()));

  PostData(1, fList);

  if (fDoTiming) fHTiming->Fill(kTimingTotal, total.CpuTime());
}

//____________________________________________________________________
void
AliForwardMultiplicityTask::FinishTaskOutput()
{
  if (!fList) 
    Warning("FinishTaskOutput", "No list defined");
  else {
    if (fDebug) 
      fList->ls();
  }
  AliAnalysisTaskSE::FinishTaskOutput();
}


//
// EOF
//
