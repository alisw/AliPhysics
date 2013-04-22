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
#include "AliForwardCorrectionManager.h"
#include "AliForwardUtil.h"
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
AliCentralMCMultiplicityTask::AliCentralMCMultiplicityTask(const AliCentralMCMultiplicityTask& o)
  : AliCentralMultiplicityTask(o),
    fTrackDensity(o.fTrackDensity),
    fAODMCCentral(o.fAODMCCentral)
{
  //
  // Copy constructor 
  // 
  DGUARD(fDebug, 3,"COPY CTOR of AliCentralMCMultiplicityTask");
}
//____________________________________________________________________
AliCentralMCMultiplicityTask&
AliCentralMCMultiplicityTask::operator=(const AliCentralMCMultiplicityTask& o)
{
  // 
  // Assignment operator 
  //
  DGUARD(fDebug,3,"Assignment of AliCentralMCMultiplicityTask");
  if (&o == this) return *this; 
  AliCentralMultiplicityTask::operator=(o);
  fAODMCCentral     = o.fAODMCCentral;
  fTrackDensity     = o.fTrackDensity;
  return *this;
}
//____________________________________________________________________
void AliCentralMCMultiplicityTask::UserCreateOutputObjects() 
{
  // 
  // Create output objects 
  // 
  //
  DGUARD(fDebug,1,"Create user output in AliCentralMCMultiplicityTask");
  AliCentralMultiplicityTask::UserCreateOutputObjects();

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah)
 {	 
	AliFatal("No AOD output handler set in analysis manager");
  
  
  	TObject* obj = &fAODMCCentral;
  	ah->AddBranch("AliAODCentralMult", &obj);

  }
  fTrackDensity.CreateOutputObjects(fList);

}
//____________________________________________________________________
void AliCentralMCMultiplicityTask::UserExec(Option_t* option) 
{
  // 
  // Process each event 
  // 
  // Parameters:
  //    option Not used
  //  
  DGUARD(fDebug,1,"Process event in AliCentralMCMultiplicityTask");
  fAODMCCentral.Clear("");
  // Call base class 
  AliCentralMultiplicityTask::UserExec(option);
fAODMCCentral.Init(*(fAODCentral.GetHistogram().GetXaxis()));
  // check if we need this event 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah)
{  
  //  AliFatal("No AOD output handler set in analysis manager");

  // if base class did not want this event, then neither to we 
  if (!ah->GetFillAOD() || fIvz <= 0) return;
 } 
  const AliMCEvent*  mcEvent = MCEvent();
  if (!mcEvent) return;
  TH2D&              hist    = fAODMCCentral.GetHistogram();

  Double_t vz = GetManager().GetSecMap()->GetVertexAxis().GetBinCenter(fIvz);

  fTrackDensity.Calculate(*mcEvent, vz, hist, NULL);

  CorrectData(hist, fIvz);

}

//____________________________________________________________________
void AliCentralMCMultiplicityTask::Terminate(Option_t* option) 
{
  // 
  // End of job
  // 
  // Parameters:
  //    option Not used 
  //
  DGUARD(fDebug,1,"Final analysis of merge in AliCentralMCMultiplicityTask");
  AliCentralMultiplicityTask::Terminate(option);
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
