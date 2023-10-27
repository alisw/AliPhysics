///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskGenMcFlattenicity                           //
//            Last update: 28/04/2023                            //
//                                                               //
///////////////////////////////////////////////////////////////////

AliAnalysisTask *AddTaskGenMcFlattenicity(double minpT = 0.5,
                                          double maxEta = 0.8,
                                          bool isPP = kTRUE,
                                          int inGenerator = 0,
                                          TString suffixName = "") {
  // inGenerator: 0 (Monash), 1 (Monash NoCR), 2 (Monash Ropes), 3 (Epos LHC), 4
  // (Herwig), 5 (AMPT), 6 (AMPT no string melting)
  AliAnalysisTaskGenMcFlattenicity *taskUE =
      new AliAnalysisTaskGenMcFlattenicity("taskKno");
  taskUE->SetPtMin(minpT);
  taskUE->SetEtaMax(maxEta);
  taskUE->SetIsPP(isPP);
  taskUE->SetGenerator(inGenerator);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Printf("AliAnalysisTaskSimSpectraLF: No analysis manager to connect to.");
    return 0x0;
  }
  // get the input event handler this handler is part of the managing system and
  // feeds events to your task
  if (!mgr->GetMCtruthEventHandler()) {
    Printf("AliAnalysisTaskSimSpectraLF: This task requires an input MC event "
           "handler.");
    return 0x0;
  }

  mgr->AddTask(taskUE);

  // Create containers for input/output

  TString finDirname = "";
  TString inname = "cinput";
  TString outBasic = "cList";

  finDirname += suffixName.Data();
  inname += finDirname.Data();
  outBasic += finDirname.Data();

  // Input and Output Slots
  //===========================================================================

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWGMM_GenKnoUe_%1.2f", minpT);

  AliAnalysisDataContainer *coutSim = mgr->CreateContainer(
      outBasic, TList::Class(), AliAnalysisManager::kOutputContainer,
      outputfile.Data());

  mgr->ConnectInput(taskUE, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskUE, 1, coutSim);

  return taskUE;
}
