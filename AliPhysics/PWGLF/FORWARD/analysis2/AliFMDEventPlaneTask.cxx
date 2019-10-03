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
  : AliBaseAODTask(),
    fEventPlaneFinder(), // EP finder
    fHistVertexSel()     // Diagnostics histogram
{
  // 
  // Default constructor
  //
  DGUARD(fDebug, 3,"Default CTOR of AliFMDEventPlaneTask");
}
//_____________________________________________________________________
AliFMDEventPlaneTask::AliFMDEventPlaneTask(const char* name) 
  : AliBaseAODTask(name, "AliFMDEventPlaneTask"),
    fEventPlaneFinder("eventPlane"), // EP finder
    fHistVertexSel(0)               // Diagnostics histogram
{
  // 
  // Constructor
  //
  // Parameters:
  //  name: Name of task
  //
  DGUARD(fDebug, 3,"Named CTOR of AliFMDEventPlaneTask: %s", name);
}
//_____________________________________________________________________
Bool_t AliFMDEventPlaneTask::Book()
{
  //
  // Create output objects
  //
  DGUARD(fDebug,1,"Create user objects of AliFMDEventPlaneTask");
  // Diagnostics histograms
  fHistVertexSel     = new TH1D("hVertexSel", "Selected vertices", 40, -20, 20);

  fSums->Add(fHistVertexSel);

  // Init of EventPlaneFinder
  TAxis* pe = new TAxis(200, -4., 6.);
  fEventPlaneFinder.CreateOutputObjects(fSums);
  fEventPlaneFinder.SetupForData(*pe);

  return true;
}
//_____________________________________________________________________
Bool_t AliFMDEventPlaneTask::Event(AliAODEvent& aod)
{
  // 
  // Called each event
  //
  // Parameters:
  //  option: Not used
  //
  DGUARD(fDebug,1,"Process an event in AliFMDEventPlaneTask");

  // Reset data members
  AliAODForwardMult* aodfmult = GetForward(aod);
  fHistVertexSel->Fill(aodfmult->GetIpZ());

  if (aod.GetRunNumber() != fEventPlaneFinder.GetRunNumber())
    fEventPlaneFinder.SetRunNumber(aod.GetRunNumber());

  AliAODForwardEP aodep;
  TH2D& fmdHist = aodfmult->GetHistogram();

  fEventPlaneFinder.FindEventplane(&aod, aodep, &fmdHist, 0);

  return true;
}
//_____________________________________________________________________
Bool_t AliFMDEventPlaneTask::Finalize()
{
  //
  // Terminate - Called after all events
  //
  // Parameters:
  //  option: Not used
  //
  DGUARD(fDebug,1,"Process merged output of AliFMDEventPlaneTask");
  
  // Calculations can be done here: Currently there are none
  // Summed histograms can be found in the list fSums
  // Output should be stored in the output list fResults

  return true;
}
//_____________________________________________________________________
//
//
// EOF
