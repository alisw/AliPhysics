

void    readCDB (TObject *task1);
//_____________________________________________________________________________
AliAnalysisTask  *AddTaskT0Calib(Int_t runNumber)
{
  //
  // add calibration task
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskT0Calib", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskT0Calib", "This task requires an input event handler");
    return NULL;
  }  

  // set TPC OCDB parameters
  //ConfigOCDB(runNumber);

  // setup task
  AliT0CalibOffsetChannelsTask  *task1 = new AliT0CalibOffsetChannelsTask("CalibObjectsTrain1");
  readCDB(task1);
  mgr->AddTask(task1);
  
  //  AliT0AnalysisTaskQA * task2 = new AliT0AnalysisTaskQA("QA task");
  //    mgr->AddTask(task2);

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("T0Calib",TObjArray::Class(), AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root");  

  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,1,coutput1);
  return task1;
}
//_____________________________________________________________________________
void    readCDB (TObject *task1) {

  Float_t zero_timecdb[24]={0};
  Float_t *timecdb = zero_timecdb;
  Float_t cfdvalue[24][5];
  Float_t zero_shiftcdb[4]={0};
  Float_t *shiftcdb = zero_shiftcdb;
  AliT0CalibOffsetChannelsTask *mytask = (AliT0CalibOffsetChannelsTask*)task1;

  //  AliCDBManager* man = AliCDBManager::Instance();
  //  man->SetDefaultStorage("raw://");
  // man->SetRun(158124);
  AliCDBEntry *entryCalib1 = man->Get("T0/Calib/TimeDelay");
  if(!entryCalib1) {
    AliError::(Form("Cannot find any AliCDBEntry for [Calib, TimeDelay]!"));
    return;
  }
  else
    {
      AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entryCalib1->GetObject();
      timecdb = clb->GetTimeEq();
      for(Int_t i=0; i<24; i++) 
	for (Int_t i0=0; i0<5; i0++){
	  cfdvalue[i][i0] = clb->GetCFDvalue(i, i0);
	}
    }
  for (Int_t i=0; i<24; i++) {
    mytask->SetCFDvalue(i,cfdvalue[i][0]);
    mytask->SetTimeEq(i,timecdb[i]);
  } 

  AliCDBEntry *entryCalib2 = man->Get("T0/Calib/TimeAdjust");
  if(!entryCalib2) {
     AliError(Form("Cannot find any AliCDBEntry for [Calib, TimeAdjust]!"));
  }
 else
    {
      AliT0CalibSeasonTimeShift *clb1 = (AliT0CalibSeasonTimeShift*)entryCalib2->GetObject();
      shiftcdb = clb1->GetT0Means();
    }
  
  for (Int_t i=0; i<4; i++)  mytask->SetT0Means(i,shiftcdb[i]);
}
