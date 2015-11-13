//====================================================================
// 
// Base class for classes that calculate the multiplicity in the
// central region event-by-event
// 
// Inputs: 
//   - AliESDEvent 
//
// Outputs: 
//   - AliAODCentralMult 
// 
// Histograms 
//   
// Corrections used 
#include "AliCentralMCMultiplicityTask.h"
#include "AliCentralCorrectionManager.h"
#include "AliCentralCorrAcceptance.h"
#include "AliCentralCorrSecondaryMap.h"
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliFMDEventInspector.h"
#include <AliMCEvent.h>
#include <AliTrackReference.h>
#include <AliStack.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TError.h>
#include <iostream>
#include <iomanip>

//====================================================================
AliCentralMCMultiplicityTask::AliCentralMCMultiplicityTask(const char* name) 
  : AliCentralMultiplicityTask(name),
    fTrackDensity(name),
    fAODMCCentral(kTRUE)
{
  // 
  // Constructor 
  //   
  DGUARD(fDebug,3,"Named CTOR of AliCentralMCMultiplicityTask: %s", 
	 name);
  fBranchNames = 
    "ESD:AliESDRun.,AliESDHeader.,AliMultiplicity.,"
    "SPDVertex.,PrimaryVertex.";
}
//____________________________________________________________________
AliCentralMCMultiplicityTask::AliCentralMCMultiplicityTask() 
  : AliCentralMultiplicityTask(),
    fTrackDensity(),
    fAODMCCentral(kTRUE)
{
  // 
  // Constructor 
  // 
  DGUARD(fDebug, 3,"Default CTOR of AliCentralMCMultiplicityTask");
}
//____________________________________________________________________
void
AliCentralMCMultiplicityTask::CreateBranches(AliAODHandler* ah) 
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user output in AliCentralMCMultiplicityTask");
  AliCentralMultiplicityTask::CreateBranches(ah);

  if (!ah) 
    // AliFatal("No AOD output handler set in analysis manager");
    return;

  
  TObject* obj = &fAODMCCentral;
  ah->AddBranch("AliAODCentralMult", &obj);
}
//____________________________________________________________________
Bool_t
AliCentralMCMultiplicityTask::Book()
{
  AliCentralMultiplicityTask::Book();
  fTrackDensity.CreateOutputObjects(fList);
  return true;
}
//____________________________________________________________________
Bool_t AliCentralMCMultiplicityTask::PreData(const TAxis& v, const TAxis& e)
{
  AliCentralMultiplicityTask::PreData(v, e);
  fAODMCCentral.Init(e);  
  return true;
}

//____________________________________________________________________
Bool_t AliCentralMCMultiplicityTask::PreEvent()
{
  AliCentralMultiplicityTask::PreEvent();  
  fAODMCCentral.Clear("");
  return true;
}
//____________________________________________________________________
Bool_t AliCentralMCMultiplicityTask::Event(AliESDEvent& esd) 
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  DGUARD(fDebug,1,"Process event in AliCentralMCMultiplicityTask");

  fIvz               = 0;
  Bool_t   lowFlux   = kFALSE;
  UInt_t   triggers  = 0;
  UShort_t ivz       = 0;
  TVector3 ip;
  Double_t cent      = -1;
  UShort_t nClusters = 0;
  UInt_t   found     = fInspector.Process(&esd, triggers, lowFlux, 
					  ivz, ip, cent, nClusters);

  // Make sure AOD is filled
  MarkEventForStore();

  // Is this accepted for analysis? 
  Bool_t isAccepted = true;

  // No event or no trigger 
  if (found &  AliFMDEventInspector::kNoEvent)    isAccepted = false;
  if (found &  AliFMDEventInspector::kNoTriggers) isAccepted = false;
  if (found &  AliFMDEventInspector::kNoSPD)      isAccepted = false;
  if (found &  AliFMDEventInspector::kNoVertex)   isAccepted = false;
  if (triggers & AliAODForwardMult::kPileUp)      isAccepted = false;
  if (found &  AliFMDEventInspector::kBadVertex)  isAccepted = false; 

  VtxBin* bin = static_cast<VtxBin*>(fVtxList->At(ivz));
  if (!bin) return false;

  //Doing analysis
  if (isAccepted) {
    const AliMultiplicity* spdmult = esd.GetMultiplicity();
    TH2D&                  aodHist = fAODCentral.GetHistogram();

    ProcessESD(aodHist, spdmult);
    bin->Correct(aodHist, fUseSecondary, fUseAcceptance);
  
    if (triggers & AliAODForwardMult::kInel) 
      fHData->Add(&(fAODCentral.GetHistogram()));
  }

  const AliMCEvent*  mcEvent = MCEvent();
  if (!mcEvent) return false;
  TH2D&              hist    = fAODMCCentral.GetHistogram();


  fTrackDensity.Calculate(*mcEvent, ip, hist, NULL);
  bin->Correct(hist, fUseSecondary, fUseAcceptance, false);

  return true;
}

//____________________________________________________________________
void
AliCentralMCMultiplicityTask::Print(Option_t* option) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  AliCentralMultiplicityTask::Print(option);
  gROOT->IncreaseDirLevel();
  fTrackDensity.Print(option);
  gROOT->DecreaseDirLevel();
}
//
// EOF
//
