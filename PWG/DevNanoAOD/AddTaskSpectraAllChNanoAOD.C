AliAnalysisTaskSpectraAllChNanoAOD* AddTaskSpectraAllChAOD(
						       Bool_t mc=kFALSE,
						       Double_t CentCutMin=0,
						       Double_t CentCutMax=100,
						       Double_t QvecCutMin=0,
						       Double_t QvecCutMax=100,
						       Double_t EtaMin=-0.8,
						       Double_t EtaMax=0.8,
						       Double_t pt=50.,
						       Double_t ptTofMatch=.6,
                                                       UInt_t trkbit=1,
						       Double_t DCA=100000,
						       UInt_t minNclsTPC=70,
						       Double_t nsigmacut=5.,
						       Int_t PIDtype=3,
						       TString opt=""){
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddAliAnalysisTaskSpectraAllChAOD", "No analysis manager to connect to.");
      return NULL;
    }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AliAnalysisTaskSpectraAllChAOD", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD"))
    {
      ::Error("AliAnalysisTaskSpectraAllChAOD", "This task requires to run on AOD");
      return NULL;
    }
  
  AliSpectraAODTrackCuts  * trcuts = new AliSpectraAODTrackCuts(Form("TrackCuts%s",opt.Data()));  
  trcuts->SetDCA(DCA);
  trcuts->SetTrackBits(trkbit);
  trcuts->SetPt(pt);
  trcuts->SetPtTOFMatching(ptTofMatch);   
  trcuts->SetEta(EtaMin,EtaMax);
  trcuts->SetMinTPCcls(minNclsTPC);
  trcuts->PrintCuts();
  
  AliSpectraAODEventCuts * evcuts = new AliSpectraAODEventCuts(Form("EventCuts%s",opt.Data()));
  evcuts->SetQVectorCut(QvecCutMin,QvecCutMax);
  evcuts->SetCentralityCutMax(CentCutMax);  
  evcuts->SetCentralityCutMin(CentCutMin);
  if(mc==1)evcuts->SetIsMC(kTRUE);
  evcuts->PrintCuts();
  
  AliHelperPID *pid=new AliHelperPID();
  pid->SetName(Form("HelperPID%s",opt.Data()));
  pid->SetNSigmaCut(nsigmacut);
  pid->SetPIDType(PIDtype);
  if(PIDtype==3){
    AliPIDCombined *pidc=new AliPIDCombined();
    pidc->SetDefaultTPCPriors();
    pid->SetPIDCombined(pidc);
  }
  
  AliAnalysisTaskSpectraAllChNanoAOD *task = new AliAnalysisTaskSpectraAllChNanoAOD(Form("TaskAODSpectraCent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_TrBit%d%s",	
                                                                                         CentCutMin,
                                                                                         CentCutMax,
                                                                                         QvecCutMin,
                                                                                         QvecCutMax,
                                                                                         EtaMin,
                                                                                         EtaMax,
                                                                                         trkbit,
                                                                                         opt.Data()));
  task->SetEventCuts(evcuts);
  task->SetTrackCuts(trcuts);
  task->SetHelperPID(pid);
  if(mc==1)task->SetIsMC(kTRUE);
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  
  TString typeofdata=mc?"MC":"Data";
  
  outputFileName += Form(":SpectraESE_%s%s",typeofdata.Data(),opt.Data());
  
  cout<<"outputFileName:  "<<outputFileName<<endl;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
  AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("chist%s",opt.Data()),                      TList::Class(),     AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("cvcut%s",opt.Data()), AliSpectraAODEventCuts::Class(),     AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("ctcut%s",opt.Data()), AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputpt4 = mgr->CreateContainer(Form("cpid%s",opt.Data()),               AliHelperPID::Class(),     AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->AddTask(task);
  
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputpt1);
  mgr->ConnectOutput(task, 2, coutputpt2);
  mgr->ConnectOutput(task, 3, coutputpt3);
  mgr->ConnectOutput(task, 4, coutputpt4);
  
  return task;
}
