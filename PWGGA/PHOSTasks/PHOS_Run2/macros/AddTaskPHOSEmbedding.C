AliAnalysisTaskPHOSEmbedding* AddTaskPHOSEmbedding(
    const char* name     = "PHOSEmbedding",
    const UInt_t trigger = AliVEvent::kINT7,
    const TString parname = "Pi0",
    const TString filepath = "Embedding_Pi0_AODMC_246980.txt"
    )
{
  //Add a task AliAnalysisTaskPHOSEmbedding to the analysis train
  //Author: Daiki Sekihata
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSEmbedding", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSEmbedding", "This task requires an input event handler");
    return NULL;
  }

  //AliCDBManager *cdb = AliCDBManager::Instance();
  //cdb->SetDefaultStorage("raw://");
  //cdb->SetRun(246982);//dummy run number
  TString taskname = Form("%s_%s",name,parname.Data());

  AliAnalysisTaskPHOSEmbedding* task = new AliAnalysisTaskPHOSEmbedding(taskname);
  task->SelectCollisionCandidates(trigger);
  task->SetParticle(parname);

  TObjArray *array = new TObjArray();
  TObjString* ostr = 0x0;//TString can not be stored in TObjArray, because TString class is not inherited from TObject.
  TString line = "";
  ifstream fin;

  fin.open(filepath);
  while(1){
    fin >> line;
    if(fin.eof()) break;

    //printf("found AOD MC path = %s.\n", line.Data());

    ostr = new TObjString(line.Data());
    array->Add(ostr);
  }
  fin.close();

  printf("%d AOD MC paths are found in %s.\n", array->GetEntries(),filepath.Data());

  task->SetInputFileArray(array);//array of path to MC AOD.

  //Need MagFeild for tender
  //((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
 
  //TString outputFile = AliAnalysisManager::GetCommonFileName();
  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s","list_test"), THashList::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());//event characerization
  //mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

