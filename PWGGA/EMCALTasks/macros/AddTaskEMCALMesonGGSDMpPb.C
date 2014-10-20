void AddTaskEMCALMesonGGSDMpPb(Int_t calibration = 0) 
{
	
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAliEMCALMesonGGSDMpPb", "No analysis manager to connect to.");
    return NULL;
  }  

  //#####################################################
  // Private Recalibrator:
     Int_t GoodTasks  [4] = {1, 0,0,0};
     Int_t RecalScheme[4] = {calibration, 3,5,7};
  
  // Int_t calibration: 
  // 0:  no recalibration! 
  // 1:  constant scale factor.
  // 2:  Symmetric Decay Method
  // 3:  J's fit to LHC12f1a MC single photons, 4 Aug 2013
  // 4:  J's fit to the test beam data, 4 Aug 2013
  // 5:  Based on kSDM/kTBCv3 (for MC)
  // 6:  kBeamTestCorrectedv2 - in AliROOT! 
  // 7:  kPi0MCv3 - in AliROOT! 
  // 8:  kSDMv5 - based on J's fit to the noNL MC/data and kPi0MCv5 - 28 Oct 2013 (call it kSDMv5)
  // 9:  kPi0MCv5 - J's fit to LHC12f1a/b MC single photons, 28 Oct 2013 (call it kPi0MCv5)
  // 10: kBTCv6 - J played with test beam data - 19 Nov 2013
  // 11: kPi0MCv6 - J played with test beam mc - 19 Nov 2013
  //
  //#####################################################  

  AliAnalysisTaskEMCALMesonGGSDMpPb* task[4];
  AliAnalysisDataContainer*  coutput[4];  
  char saythis[500];

  for(int i=0; i<4; i++){
    if(GoodTasks[i]==0)
      continue;
    
    sprintf(saythis,"EMCALMesonGGSDMpPbTask_%d",i);
    task[i] = new AliAnalysisTaskEMCALMesonGGSDMpPb(saythis);
    //task[i]->SelectCollisionCandidates(AliVEvent::kAny);
    //task[i]->SelectCollisionCandidates(AliVEvent::kMB);//LHC11a
    task[i]->SelectCollisionCandidates(AliVEvent::kINT7);//LHC13b/c
    task[i]->SetRecalScheme(RecalScheme[i]);

    task[i]->SetdRmin_ClustTrack(0.025);//not used at the moment! (Gustavo uses 0.04)
    task[i]->SetFidPhiMinMax(1.39626, 3.15); // full emcal. 
    //task[i]->SetFidPhiMinMax(1.39626, 2.10); //pPb no TRD
    //task[i]->SetFidPhiMinMax(2.10, 3.15); //pPb with TRD
    task[i]->SetFidEtaMinMax(-0.65, 0.65);
    
    mgr->AddTask(task[i]);
    sprintf(saythis,"cont_AliAnalysisTaskEMCALMesonGGSDMpPb_%d",RecalScheme[i]);
    coutput[i] = 
      mgr->CreateContainer(saythis,
			   TList::Class(),
			   AliAnalysisManager::kOutputContainer,
			   "EMCALMesonGGSDMpPbTask.root");
    mgr->ConnectInput (task[i],0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[i],1,coutput[i]); 
    //AliMCEventHandler* handler = new AliMCEventHandler;
    //handler->SetReadTR(kFALSE);
    //mgr->SetMCtruthEventHandler(handler);
  }

}
