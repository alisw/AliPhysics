void AddTaskPIDflowQA( Float_t centrMin=0.,
                       Float_t centrMax=100.,
                       TString fileNameBase="outputPIDQA"  )
{   
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("event cuts");
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  //cutsEvent->SetRefMultMethod(AliFlowEventCuts::kTPConly);
  //cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kSPD1tracklets);
  cutsEvent->SetNContributorsRange(2);
  cutsEvent->SetPrimaryVertexZrange(-10.,10.);
  cutsEvent->SetCutSPDvertexerAnomaly(); //"Francesco's cut"
  cutsEvent->SetCutZDCtiming();
  
	AliESDtrackCuts*  cutsTrack = new AliESDtrackCuts("cutsTrack","cutsTrack");
	cutsTrack->SetPtRange(0.2,5.);
 	cutsTrack->SetEtaRange(-0.8,0.8);
 	cutsTrack->SetMinNClustersTPC(70);
 	cutsTrack->SetMaxChi2PerClusterTPC(4.0);
 	cutsTrack->SetMaxDCAToVertexXY(0.3);
 	cutsTrack->SetMaxDCAToVertexZ(0.3);
 	cutsTrack->SetDCAToVertex2D(kTRUE);
 	cutsTrack->SetAcceptKinkDaughters(kFALSE);
	cutsTrack->SetRequireTPCRefit(kTRUE);
	cutsTrack->SetRequireITSRefit(kTRUE);
	cutsTrack->SetMinNClustersITS(2);
	
	//task1->SetNsigmaDCAcut(5.0,5.0);
	//task1->SetMCOn();	

  TString centralityName("");
  centralityName+=Form("%.0f",centrMin);
  centralityName+="-";
  centralityName+=Form("%.0f",centrMax);

  TString fileName(fileNameBase);
  fileName.Append(Form("%.0f%.0f",centrMin,centrMax));
  fileName.Append(".root");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowEvent", "No analysis manager to connect to.");
    return;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFlowEvent", "This task requires an input event handler");
    return;
  }  

  AliAnalysisTaskPIDflowQA* task = new AliAnalysisTaskPIDflowQA(Form("taskPIDQA%.0f%.0f",centrMin,centrMax));
	task->SetAliESDtrackCuts(cutsTrack);
  task->SetEventCuts(cutsEvent);
  task->SelectCollisionCandidates(AliVEvent::kMB);
  //old
  task->GetESDpid()->GetTPCResponse().SetBetheBlochParameters(0.0283086,
                                                    2.63394e+01,
                                                    5.04114e-11,
                                                    2.12543e+00,
                                                    4.88663e+00 );
  //new
  //task->GetESDpid()->GetTPCResponse().SetBetheBlochParameters(1.28949/50.,
  //                                                  2.74095e+01,
  //                                                  TMath::Exp(-3.21763e+01),
  //                                                  2.44026,
  //                                                  6.58800);


  mgr->AddTask(task);

  AliAnalysisDataContainer* coutputQAtask = mgr->CreateContainer(Form("output%s",centralityName.Data()),
                                                 TList::Class(),
                                                 AliAnalysisManager::kOutputContainer,
                                                 fileName.Data());
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutputQAtask);
}
