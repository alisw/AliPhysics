AliAnalysisTaskSpectraBoth* AddTaskSpectraBoth(Bool_t mc=kFALSE,
					     Double_t CentCutMin=0,
					     Double_t CentCutMax=100,
					     Double_t QvecCutMin=0,
					     Double_t QvecCutMax=100,
					     Double_t EtaMin=-0.8,
					     Double_t EtaMax=0.8,
					     Double_t Nsigmapid=3.,
					     Double_t pt=5.,
					     Double_t p=5.,
					     Double_t ymin=-0.5,
					     Double_t ymax=.5, 	
					     Double_t ptTofMatch=.6,
					     UInt_t trkbit=1,
					     UInt_t trkbitQVector=1,
					     Bool_t UseCentPatchAOD049=kFALSE,
					     Double_t DCA=100000,
					     UInt_t minNclsTPC=70,
					     Int_t nrebin=0,
					     TString centestimator="V0M",
					     Int_t pidmethod=3,
					     TString taskname="")
 {
  
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
  
  //TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  //if(type.Contains("ESD"))
   // {
    //  ::Error("AddTaskITSsaTracks", "This task requires to run on AOD");
    //  return NULL;
 //   }
  
  using namespace AliSpectraNameSpaceBoth;
  
  TString opt=Form("Cent%.3fto%.3f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%dEst_%s_Pid_%d_Y%.1fto%.1f",CentCutMin,CentCutMax,QvecCutMin,QvecCutMax,EtaMin,EtaMax,Nsigmapid,trkbit,centestimator.Data(),pidmethod,ymin,ymax);
  if(taskname.Length()>0)
	opt=taskname;	

	
  AliSpectraBothPID *pid = new AliSpectraBothPID(); 
  pid->SetNSigmaCut(Nsigmapid);
  if(pidmethod==0)
	pid->SetPIDtype(AliSpectraBothPID::kNSigmaTPC);
  else if (pidmethod==1)
	pid->SetPIDtype(AliSpectraBothPID::kNSigmaTOF);
  else if  (pidmethod==2)		
  	pid->SetPIDtype(AliSpectraBothPID::kNSigmacircleTPCTOF);
  else if (pidmethod==3)
	pid->SetPIDtype(AliSpectraBothPID::kNSigmasquareTPCTOF);	  
  else 
	pid->SetPIDtype(AliSpectraBothPID::kNSigmaTPCorTOF);
  AliSpectraBothTrackCuts  * trcuts = new AliSpectraBothTrackCuts("Track Cuts");  
  trcuts->SetDCA(DCA);
  trcuts->SetTrackBits(trkbit);
  trcuts->SetPt(pt);
  trcuts->SetP(p);
  trcuts->SetY(ymax,ymin);
  trcuts->SetPtTOFMatching(ptTofMatch);   
  trcuts->SetEta(EtaMin,EtaMax);
  trcuts->SetMinTPCcls(minNclsTPC);
  trcuts->PrintCuts();
  
  AliSpectraBothEventCuts * evcuts = new AliSpectraBothEventCuts("Event Cuts");
  evcuts->SetQVectorCut(QvecCutMin,QvecCutMax);
  evcuts->SetCentralityCutMax(CentCutMax);  
  evcuts->SetCentralityCutMin(CentCutMin);
  evcuts->SetTrackBits(trkbitQVector);
  evcuts->SetCentEstimator(centestimator);	
  if(mc==1)evcuts->SetIsMC(kTRUE);
  evcuts->PrintCuts();
  
  AliAnalysisTaskSpectraBoth *task = new AliAnalysisTaskSpectraBoth(Form("TaskBothSpectra%s",opt.Data()));
  task->SetPID(pid);  
  task->SetEventCuts(evcuts);
  task->SetTrackCuts(trcuts);
  task->SetNRebin(nrebin);
  if(mc==1)task->SetIsMC(kTRUE);
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  
  TString typeofdata=mc?"MC":"Data";
  outputFileName += Form(":OutputBothSpectraTask_%s_%s",typeofdata.Data(),opt.Data());
  
  TString tmpstring= Form("OutputBothSpectraTask_%s_%s",typeofdata.Data(),opt.Data());
  
  cout<<"outputFileName:  "<<outputFileName<<endl;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();      
  AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer(Form("%schistpt",tmpstring.Data()), AliSpectraBothHistoManager::Class(),  AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer(Form("%scvcutpt",tmpstring.Data()), AliSpectraBothEventCuts::Class(),    AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer(Form("%sctcutpt",tmpstring.Data()), AliSpectraBothTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *coutputpt4 = mgr->CreateContainer(Form("%scpidpt",tmpstring.Data()),  AliSpectraBothPID::Class(),     AliAnalysisManager::kOutputContainer,outputFileName);
  
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputpt1);
  mgr->ConnectOutput(task, 2, coutputpt2);
  mgr->ConnectOutput(task, 3, coutputpt3);
  mgr->ConnectOutput(task, 4, coutputpt4);
  
  mgr->AddTask(task);
  return task;
}
