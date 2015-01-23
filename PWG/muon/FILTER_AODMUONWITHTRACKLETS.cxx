#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskAOD2MuonAOD.h"
#include "AliAODInputHandler.h"
#include "TChain.h"
#include "Riostream.h"

namespace AAF {
  
//
// Convert a full (std) AOD to a muon only AOD
//
//

void FILTER_AODMUONWITHTRACKLETS(const char* from, const char* to)
{
  AliAnalysisManager *mgr = new AliAnalysisManager("AOD2MUONAOD");
  
  AliInputEventHandler* input = new AliAODInputHandler;
  
  mgr->SetInputEventHandler(input);

  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetCreateNonStandardAOD();
  aodHandler->SetOutputFileName(to);
  mgr->SetOutputEventHandler(aodHandler);

  AliAnalysisTask* task = new AliAnalysisTaskAOD2MuonAOD(0,kTRUE);

  mgr->AddTask(task);

  // Connect input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());

  if (!mgr->InitAnalysis())
  {
    std::cout << "Could not InitAnalysis" << std::endl;
    return;
  }

  TChain* chain = new TChain("aodTree");
  
  chain->Add(from);
  
  mgr->StartAnalysis("local",chain);
}

}