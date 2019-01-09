AliAnalysisTaskPHOSEmbedding* AddTaskPHOSEmbedding(
    const char* name     = "PHOSEmbedding",
    const UInt_t trigger = AliVEvent::kINT7,
    const TString parname = "Pi0",
    const TString filepath = "Embedding_Pi0_AODMC_246980.txt",
    const Double_t Ecorrection = 1.02
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
  task->SetSignalCalibration(Ecorrection);

  TF1 *f1nonlin = new TF1("f1nonlin","[2]*(1.+[0]/(1. + TMath::Power(x/[1],2)))",0,100);
  f1nonlin->SetNpx(1000);
  f1nonlin->SetParNames("a","b (GeV)","c");
  f1nonlin->FixParameter(0,-0.06);
  f1nonlin->FixParameter(1,  0.7);
  f1nonlin->FixParameter(2, 1.00);
  task->SetUserNonlinearity(f1nonlin);

  //Need MagFeild for tender
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
 
  TString outputFile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("hist_%s",taskname.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFile.Data(),"PWGGA_PHOSTasks_PHOSRun2"));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

