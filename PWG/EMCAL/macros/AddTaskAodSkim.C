// $Id: AddTaskSkim.C 4586 2013-01-16 15:32:16Z loizides $

AliAodSkimTask *AddTaskAodSkim(const Double_t mine=5,
                               const UInt_t trigsel = AliVEvent::kINT7,
                               const char *name=0) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskAodSkim", "No analysis manager found.");
    return 0;
  }

  TString tname(name);
  if (name==0) 
    tname = Form("aodskim_mine%.1f_trigsel%u",mine,trigsel);

  TString fname(tname);
  fname += ".root";

  cout << "Calling AddTaskAodSkim with " << endl;
  cout << "-> mine=" << mine << endl;
  cout << "-> trigsel=" << trigsel << endl;
  if (name) 
    cout << "-> name=" << name << endl;
  else 
    cout << "-> name=0" << endl;
  cout << "Writing skimmed aod tree to " << fname.Data() << endl;
 
  AliAODHandler *output = new AliAODHandler;
  output->SetOutputFileName(fname);
  output->SetFillAOD(1);
  output->SetFillAODforRun(1);
  output->SetFillExtension(0);
  mgr->SetOutputEventHandler(output);

  const char *name = "skimtask";
  AliAodSkimTask *task = new AliAodSkimTask(tname);
  task->SetClusMinE(mine);
  task->SelectCollisionCandidates(trigsel);
  
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  TString contName(tname);
  contName += "_histos";
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(contName.Data(), 
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput );
  mgr->ConnectOutput (task, 1, coutput );

  return task;
}
