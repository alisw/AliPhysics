// WARNING This does not work yet. (Only TPC supported). Please avoid using this. Run the PID response during the filtering, and store the fields you need.

void runOnNano_PIDResponse()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("train");

  AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);

  // PID response
  gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  
  AliAnalysisTaskSE* task = (AliAnalysisTaskSE*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/tutorial/AddTaskSimple.C");
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  
  mgr->InitAnalysis();
  mgr->PrintStatus();
  
  // Create chain of input files
  TChain * chain = new TChain("aodTree");
  chain->Add("AliAOD.NanoAOD.root");

  TStopwatch watch;
  mgr->StartAnalysis("local", chain, 1000);
  watch.Print();
}
