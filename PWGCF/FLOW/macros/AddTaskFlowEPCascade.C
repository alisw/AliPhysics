AliAnalysisTask *AddTaskFlowEPCascade(Float_t centrMin=0, Float_t centrMax=100,
				      TString fileName, TString suffixName="",
				      Double_t vtxCut = 10., 
				      Double_t etaCut = 0.9, 
				      Int_t nTPCcls = 70,
				      Double_t mSigma = 0.0024)
{
  /*
  //-E-V-E-N-T- -c-u-t-s-----------------------------------------------------
  AliFlowEventCuts* cutsEvent 
    = new AliFlowEventCuts(Form("event_cuts_%s", suffixName.Data()));
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kV0);
  cutsEvent->SetNContributorsRange(2);
  cutsEvent->SetPrimaryVertexZrange( -10., 10. );
  cutsEvent->SetCutSPDvertexerAnomaly();
  cutsEvent->SetCutZDCtiming();
  cutsEvent->SetCutTPCmultiplicityOutliers();

  AliFlowTrackCuts* cutsRP 
    = (AliFlowTrackCuts*) AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
  cutsRP->SetPtRange(0.15, 10.);  //added to extend pt range
  */

  AliFlowTrackCuts * cutsDaughter 
    = new AliFlowTrackCuts(Form("daughter_cuts_%s",suffixName.Data()));
  cutsDaughter->SetPtRange(0.15,10.0);
  cutsDaughter->SetEtaRange(-etaCut,etaCut);
  cutsDaughter->SetMinNClustersTPC(nTPCcls);
  cutsDaughter->SetMaxChi2PerClusterTPC(4.0);
  cutsDaughter->SetAODfilterBit(1); //TPC track only
  cutsDaughter->SetRequireITSRefit(kFALSE);
  cutsDaughter->SetRequireTPCRefit(kTRUE);
  cutsDaughter->SetMinNClustersITS(0);
  cutsDaughter->SetAcceptKinkDaughters(kFALSE);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //____________________________________________//                           
  double pdg = 1.32171; double hwid = mSigma;
  double sep = 6*hwid; double XiStart = pdg-sep; 
  double XiBands[3][2];
  XiBands[0][0]=(XiStart+6*hwid*0)-2*hwid; 
  XiBands[0][1]=(XiStart+6*hwid*0)+2*hwid;
  XiBands[1][0]=(XiStart+6*hwid*1)-2*hwid;   
  XiBands[1][1]=(XiStart+6*hwid*1)+2*hwid;
  XiBands[2][0]=(XiStart+6*hwid*2)-2*hwid; 
  XiBands[2][1]=(XiStart+6*hwid*2)+2*hwid;

  double OmegaPDG = 1.67245;   
  double Omegasg = mSigma; // 1 sigma?
  double OmegaStart = OmegaPDG-6*Omegasg; 
  double OmegaBands[3][2];
  OmegaBands[0][0]=(OmegaStart+6*Omegasg*0)-2*Omegasg;   
  OmegaBands[0][1]=(OmegaStart+6*Omegasg*0)+2*Omegasg;
  OmegaBands[1][0]=(OmegaStart+6*Omegasg*1)-2*Omegasg;     
  OmegaBands[1][1]=(OmegaStart+6*Omegasg*1)+2*Omegasg;
  OmegaBands[2][0]=(OmegaStart+6*Omegasg*2)-2*Omegasg;   
  OmegaBands[2][1]=(OmegaStart+6*Omegasg*2)+2*Omegasg;

  AliAnalysisTaskFlowEPCascade * task
    = new AliAnalysisTaskFlowEPCascade(Form("v2EP_%s", 
					    suffixName.Data()), 
				       centrMin, centrMax,
				       XiBands, OmegaBands);
  // task->SetFlowEventCuts(cutsEvent);
  //task->SetFlowTrackCuts(cutsRP);
  task->SetVertexCut(vtxCut);
  task->SetFlowDauTrackCuts(cutsDaughter);
  task->SelectCollisionCandidates(AliVEvent::kMB);
  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("hist_%s", suffixName.Data()), 
			 TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s.root",fileName.Data()));
  // Connect input/output                                                      
  mgr->ConnectInput(task, 0, cinput);
  //mgr->ConnectInput(task, 1, cinput1);
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}

