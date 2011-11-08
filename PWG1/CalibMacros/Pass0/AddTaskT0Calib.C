

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
  readCDB(task1, runNumber);
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
void    readCDB (TObject *task1,  Int_t runNumber) {

  Float_t zero_timecdb[24]={0};
  Float_t *timecdb = zero_timecdb;
  Float_t cfdvalue[24][5];
  for(Int_t i=0; i<24; i++) 
    for (Int_t i0=0; i0<5; i0++)
      cfdvalue[i][i0] = 0;
      
  Float_t zero_shiftcdb[4]={0};
  Float_t *shiftcdb = zero_shiftcdb;
  AliT0CalibOffsetChannelsTask *mytask = (AliT0CalibOffsetChannelsTask*)task1;

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(runNumber);
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);
 
  AliCDBEntry *entry4 = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");
  if (!entry4) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phase = (AliLHCClockPhase*)entry4->GetObject();
  Float_t fGRPdelays = l1Delay - phase->GetMeanPhase();

  AliCDBEntry *entryCalib0 = man->Get("T0/Calib/Latency");
  if(!entryCalib0) {
    AliError::(Form("Cannot find any AliCDBEntry for [Calib, Latency]!"));
    return;
  }
  AliT0CalibLatency *calibda=(AliT0CalibLatency*)entryCalib0->GetObject();
  Float_t fLatencyL1 = calibda->GetLatencyL1();
  Float_t fLatencyHPTDC = calibda->GetLatencyHPTDC();
 
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
    Float_t cfdmean = cfdvalue[i][0];
    if( cfdvalue[i][0] < 500 || cfdvalue[i][0] > 50000) cfdmean = ( 1000.*fLatencyHPTDC - 1000.*fLatencyL1 + 1000.*fGRPdelays)/24.4;
     mytask->SetCFDvalue(i, cfdmean);
    mytask->SetTimeEq(i, timecdb[i]);
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
