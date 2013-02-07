//
// Calculate the FMD eventplane 
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root
//
/**
 * @file   AliFMDEventPlaneTask.cxx
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:09:40 2013
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_flow
 */
#include <TList.h>
#include <TMath.h>
#include "TH2D.h"
#include "AliLog.h"
#include "TAxis.h"
#include "AliAnalysisManager.h"
#include "AliFMDEventPlaneTask.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODForwardMult.h"
#include "AliAODEvent.h"
#include "AliAODForwardEP.h"

ClassImp(AliFMDEventPlaneTask)
#if 0
; // For emacs 
#endif

AliFMDEventPlaneTask::AliFMDEventPlaneTask()
  : AliAnalysisTaskSE(),
    fSumList(0),	 // Sum list
    fOutputList(0),	 // Output list
    fAOD(0),	         // AOD input event
    fMC(0),              // MC flag
    fEventPlaneFinder(), // EP finder
    fZvertex(1111),	 // Z vertex 
    fCent(-1),		 // Centrality
    fHistCent(),         // Diagnostics histogram
    fHistVertexSel(),    // Diagnostics histogram
    fHistVertexAll()     // Diagnostics histogram
{
  // 
  // Default constructor
  //
  DGUARD(fDebug, 3,"Default CTOR of AliFMDEventPlaneTask");
}
//_____________________________________________________________________
AliFMDEventPlaneTask::AliFMDEventPlaneTask(const char* name) 
  : AliAnalysisTaskSE(name),
    fSumList(0),	             // Sum list
    fOutputList(0),                  // Output list
    fAOD(0),		             // AOD input event
    fMC(0),                          // MC flag
    fEventPlaneFinder("eventPlane"), // EP finder
    fZvertex(1111),	             // Z vertex 
    fCent(-1),                       // Centrality
    fHistCent(0),                    // Diagnostics histogram
    fHistVertexSel(0),               // Diagnostics histogram
    fHistVertexAll(0)                // Diagnostics histogram

{
  // 
  // Constructor
  //
  // Parameters:
  //  name: Name of task
  //
  DGUARD(fDebug, 3,"Named CTOR of AliFMDEventPlaneTask: %s", name);

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

}
//_____________________________________________________________________
AliFMDEventPlaneTask::AliFMDEventPlaneTask(const AliFMDEventPlaneTask& o)
  : AliAnalysisTaskSE(o),
    fSumList(o.fSumList),                    // Sumlist
    fOutputList(o.fOutputList),	             // Output list
    fAOD(o.fAOD),		             // AOD input event
    fMC(o.fMC),                              // MC flag
    fEventPlaneFinder(o.fEventPlaneFinder),  // EP finder
    fZvertex(o.fZvertex),	             // Z vertex 
    fCent(o.fCent),		             // Centrality
    fHistCent(o.fHistCent),                  // Diagnostics histogram
    fHistVertexSel(o.fHistVertexSel),        // Diagnostics histogram
    fHistVertexAll(o.fHistVertexAll)         // Diagnostics histogram
{
  // 
  // Copy constructor 
  // 
  // Parameters:
  //    o Object to copy from 
  //
  DGUARD(fDebug, 3,"Copy CTOR of AliFMDEventPlaneTask");
}
//_____________________________________________________________________
AliFMDEventPlaneTask&
AliFMDEventPlaneTask::operator=(const AliFMDEventPlaneTask& o)
{
  // 
  // Assignment operator 
  //
  DGUARD(fDebug,3,"Assignment of AliFMDEventPlaneTask");
  if (&o == this) return *this;
  fSumList           = o.fSumList;
  fOutputList        = o.fOutputList;
  fAOD               = o.fAOD;
  fMC                = o.fMC;
  fEventPlaneFinder  = o.fEventPlaneFinder;
  fZvertex           = o.fZvertex;
  fCent              = o.fCent;
  fHistCent          = o.fHistCent;
  fHistVertexSel     = o.fHistVertexSel;
  fHistVertexAll     = o.fHistVertexAll;

  return *this;
}
//_____________________________________________________________________
void AliFMDEventPlaneTask::UserCreateOutputObjects()
{
  //
  // Create output objects
  //
  DGUARD(fDebug,1,"Create user objects of AliFMDEventPlaneTask");
  if (!fSumList)
    fSumList = new TList();
  fSumList->SetName("Sums");
  fSumList->SetOwner();

  // Diagnostics histograms
  fHistCent          = new TH1D("hCent", "Centralities", 100, 0, 100);
  fHistVertexSel     = new TH1D("hVertexSel", "Selectec vertices", 40, -20, 20);
  fHistVertexAll     = new TH1D("hVertexAll", "All vertices", 40, -20, 20);

  fSumList->Add(fHistCent);
  fSumList->Add(fHistVertexSel);
  fSumList->Add(fHistVertexAll);

  // Init of EventPlaneFinder
  TAxis* pe = new TAxis(200, -4., 6.);
  fEventPlaneFinder.CreateOutputObjects(fSumList);
  fEventPlaneFinder.SetupForData(*pe);

  PostData(1, fSumList);

}
//_____________________________________________________________________
void AliFMDEventPlaneTask::UserExec(Option_t */*option*/)
{
  // 
  // Called each event
  //
  // Parameters:
  //  option: Not used
  //
  DGUARD(fDebug,1,"Process an event in AliFMDEventPlaneTask");

  // Reset data members
  fCent = -1;
  fZvertex = 1111;

  // Get input event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) return;

  AliAODForwardMult* aodfmult = 
    static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));

  if (!aodfmult) return;
  if (!AODCheck(aodfmult)) return;

  if (fAOD->GetRunNumber() != fEventPlaneFinder.GetRunNumber())
    fEventPlaneFinder.SetRunNumber(fAOD->GetRunNumber());

  AliAODForwardEP aodep;
  TH2D fmdHist = aodfmult->GetHistogram();

  fEventPlaneFinder.FindEventplane(fAOD, aodep, &fmdHist, 0);

  PostData(1, fSumList);

}
//_____________________________________________________________________
void AliFMDEventPlaneTask::Terminate(Option_t */*option*/)
{
  //
  // Terminate - Called after all events
  //
  // Parameters:
  //  option: Not used
  //
  DGUARD(fDebug,1,"Process merged output of AliFMDEventPlaneTask");

  // Reinitiate lists if Terminate is called separately!
  fSumList = dynamic_cast<TList*> (GetOutputData(1));
  if(!fSumList) {
    AliError("Could not retrieve TList fSumList"); 
    return; 
  }
 
  if (!fOutputList)
    fOutputList = new TList();
  fOutputList->SetName("Results");
  fOutputList->SetOwner();

 // Calculations can be done here: Currently there are none

  PostData(2, fOutputList);

}
// _____________________________________________________________________
Bool_t AliFMDEventPlaneTask::AODCheck(const AliAODForwardMult* aodfm) 
{
  // 
  // Function to check that and AOD event meets the cuts
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  DGUARD(fDebug,2,"Check AOD in AliFMDEventPlaneTask");

  if (!aodfm->IsTriggerBits(AliAODForwardMult::kOffline)) return kFALSE;

  fCent = (Double_t)aodfm->GetCentrality();
  if (0. >= fCent || fCent >= 80.) return kFALSE;
  fHistCent->Fill(fCent);

  fZvertex = aodfm->GetIpZ();
  fHistVertexAll->Fill(fZvertex);
  if (TMath::Abs(fZvertex) >= 10.) return kFALSE;
  fHistVertexSel->Fill(fZvertex);

  return kTRUE;

}
//_____________________________________________________________________
//
//
// EOF
