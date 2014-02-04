AliAnalysisTaskSpectraAOD* AddTaskSpectraAOD(Bool_t mc=kFALSE,
					     Double_t CentCutMin=0,
					     Double_t CentCutMax=100,
					     Double_t QvecCutMin=0,
					     Double_t QvecCutMax=100,
					     Double_t EtaMin=-0.8,
					     Double_t EtaMax=0.8,
					     Double_t Nsigmapid=3.,
					     Double_t pt=5.,
					     Double_t p=5.,
					     Double_t y=.5,
					     Double_t ptTofMatch=.6,
					     UInt_t trkbit=1,
					     UInt_t trkbitQVector=1,
					     Double_t DCA=100000,
					     UInt_t minNclsTPC=70,
					     Int_t nrebin=0,
					     TString opt=""){
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    {
      ::Error("AddAliAnalysisTaskSpectraAOD", "No analysis manager to connect to.");
      return NULL;
    }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) 
    {
      ::Error("AddTaskITSsaTracks", "This task requires an input event handler");
      return NULL;
    }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD"))
    {
      ::Error("AddTaskITSsaTracks", "This task requires to run on AOD");
      return NULL;
    }
  
  using namespace AliSpectraNameSpace;
  
  AliSpectraAODPID *pid = new AliSpectraAODPID(); 
  pid->SetNSigmaCut(Nsigmapid);
  
  AliSpectraAODTrackCuts  * trcuts = new AliSpectraAODTrackCuts("Track Cuts");  
  trcuts->SetDCA(DCA);
  trcuts->SetTrackBits(trkbit);
  trcuts->SetPt(pt);
  trcuts->SetP(p);
  trcuts->SetY(y);
  trcuts->SetPtTOFMatching(ptTofMatch);   
  trcuts->SetEta(EtaMin,EtaMax);
  trcuts->SetMinTPCcls(minNclsTPC);
  trcuts->PrintCuts();
  
  AliSpectraAODEventCuts * evcuts = new AliSpectraAODEventCuts("Event Cuts");
  evcuts->SetQVectorCut(QvecCutMin,QvecCutMax);
  evcuts->SetCentralityCutMax(CentCutMax);  
  evcuts->SetCentralityCutMin(CentCutMin);
  evcuts->SetTrackBits(trkbitQVector);
  if(mc==1)evcuts->SetIsMC(kTRUE);
  evcuts->PrintCuts();
  
  AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD(Form("TaskAODSpectraCent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d%s",	
  								       CentCutMin,
  								       CentCutMax,
  								       QvecCutMin,
  								       QvecCutMax,
  								       EtaMin,
  								       EtaMax,
  								       Nsigmapid,
  								       trkbit,
								       opt.Data()));
  task->SetPID(pid);  
  task->SetEventCuts(evcuts);
  task->SetTrackCuts(trcuts);
  task->SetNRebin(nrebin);
  if(mc==1)task->SetIsMC(kTRUE);
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  
  TString typeofdata=mc?"MC":"Data";
  
  outputFileName += Form(":OutputAODSpectraTask_%s_Cent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d%s",typeofdata.Data(),evcuts->GetCentralityMin(),evcuts->GetCentralityMax(),evcuts->GetQVectorCutMin(), evcuts->GetQVectorCutMax(),trcuts->GetEtaMin(),trcuts->GetEtaMax(),pid->GetNSigmaCut(),trcuts->GetTrackType(),opt.Data());
  
  cout<<"outputFileName:  "<<outputFileName<<endl;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
  AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("chistpt", AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("cvcutpt", AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("ctcutpt", AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputpt4 = mgr->CreateContainer("cpidpt",  AliSpectraAODPID::Class(),     AliAnalysisManager::kOutputContainer,outputFileName);

  mgr->AddTask(task);
  
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputpt1);
  mgr->ConnectOutput(task, 2, coutputpt2);
  mgr->ConnectOutput(task, 3, coutputpt3);
  mgr->ConnectOutput(task, 4, coutputpt4);

  return task;
}
