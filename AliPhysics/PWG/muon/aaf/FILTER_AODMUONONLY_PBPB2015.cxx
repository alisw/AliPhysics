#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskAOD2MuonAOD.h"
#include "AliAODInputHandler.h"
#include "TChain.h"
#include "Riostream.h"
#include "AliMuonEventCuts.h"

namespace AAF {
  
///
/// Convert a full (std) AOD to a muon only AOD (for real data).
///
/// "muon only" means here that the events kept satisfy at least one of the following
/// conditions :
/// - >= 1 muon track
/// - at least one of the muon L0 inputs (0MSL,0MUL,0MLL,0MSH) present in the event
///
/// This last condition makes this filter period dependent (as we must know the bit pattern
/// of L0s, see below the SetTrigClassPatterns call)
///
/// This version for PbPb does *not* retain the SPD tracklets (in contrast to the pp version
/// AODMUONONLY_PP2015

void FILTER_AODMUONONLY_PBPB2015(const char* from, const char* to)
{
  AliAnalysisManager *mgr = new AliAnalysisManager("AOD2MUONAOD");
  
  AliInputEventHandler* input = new AliAODInputHandler;
  
  mgr->SetInputEventHandler(input);

  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetCreateNonStandardAOD();
  aodHandler->SetOutputFileName(to);
  mgr->SetOutputEventHandler(aodHandler);

  AliMuonEventCuts* eventCuts = new AliMuonEventCuts("L0cutter","");
  
  eventCuts->SetTrigClassPatterns("0MSL|0MUL|0MSH|0MLL","0MSL:17,0MSH:18,0MLL:19,0MUL:20");
  
  AliAnalysisTask* task = new AliAnalysisTaskAOD2MuonAOD(0,kFALSE,eventCuts);

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