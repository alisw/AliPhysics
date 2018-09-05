AliAnalysisTaskNucleiKineCor* AddTaskNucleiKineCor(Double_t trigpt=5, Double_t p0=0, TString contname="cor")
{
  AliAnalysisTaskNucleiKineCor *task = new AliAnalysisTaskNucleiKineCor();
  TString name(Form("task_%s_%.1f_%.1f",contname.Data(),trigpt,p0));
  if (contname.Contains("_")) {
    name = contname;
    TObjArray *obj=contname.Tokenize("_");
    if (obj->GetEntries()<3) {
      cerr << "Contname must be given as string_val_val" << endl;
      return 0;
    }
    TString par1(obj->At(1)->GetName());
    trigpt=0.01*par1.Atoi();
    TString par2(obj->At(2)->GetName());
    p0=0.01*par2.Atoi();
    delete obj;
    cout << "Found parameters: " << trigpt << " " << p0 << endl;
  }
  task->SetPt(trigpt);
  task->SetP0(p0);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cout<<"AliAnalysisTaskNucleiKine: No analysis manager to connect to."<<endl;
    return NULL;
  }
  task->SetName(name);
  mgr->AddTask(task);

  // Create containers for input/output
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("output_%.2f_%.2f",trigpt,p0), 
							   TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
