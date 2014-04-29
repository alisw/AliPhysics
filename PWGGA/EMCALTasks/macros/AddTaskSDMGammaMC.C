void AddTaskAliAnalysisTaskGammaMC() 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAliGammaMC", "No analysis manager to connect to.");
    return NULL;
  }  

  //#####################################################
  // Private Recalibrator:
     Int_t GoodTasks  [4] = {1,0,0,0};
     Int_t RecalScheme[4] = {0,3,5,7};
  //
  // 0: no recalibration! 
  // 1: constant scale factor.
  // 2: Symmetric Decay Method
  // 3: J's fit to LHC12f1a MC single photons, 4 Aug 2013
  // 4: J's fit to the test beam data, 4 Aug 2013
  // 5: Based on kSDM/kTBCv3 (for MC)
  // 6: kBeamTestCorrectedv2 - in AliROOT! 
  // 7: kPi0MCv3 - in AliROOT! 
  //
  //#####################################################  

  AliAnalysisTaskGammaMC* task[4];
  AliAnalysisDataContainer*  coutput[4];  
  char saythis[500];

  for(int i=0; i<4; i++){
    if(GoodTasks[i]==0)
      continue;
    
    sprintf(saythis,"GammaMCTask_%d",i);
    task[i] = new AliAnalysisTaskGammaMC(saythis);
    task[i]->SelectCollisionCandidates(AliVEvent::kINT7);//LHC13b/c
    //task[i]->SelectCollisionCandidates(AliVEvent::kMB);
    task[i]->SetRecalScheme(RecalScheme[i]);
    
    task[i]->SetFidPhiMinMax(1.39626, 3.15);// my defaults: 1.39626, 3.15
    //task[i]->SetFidPhiMinMax(1.39626, 2.10); //pPb no TRD
    //task[i]->SetFidPhiMinMax(2.10, 3.15); //pPb with TRD    
    task[i]->SetFidEtaMinMax(-0.65, 0.65);

    mgr->AddTask(task[i]);
    sprintf(saythis,"cont_AliAnalysisTaskGammaMC_%d",RecalScheme[i]);
    coutput[i] = 
      mgr->CreateContainer(saythis,
			   TList::Class(),
			   AliAnalysisManager::kOutputContainer,
			   "GammaMCTask.root");
    mgr->ConnectInput (task[i],0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[i],1,coutput[i]);
        
    RequestMemory(task[i],1000*1024); // request 1.0GB memory for task
    
    //AliMCEventHandler* handler = new AliMCEventHandler;
    //handler->SetReadTR(kFALSE);
    //mgr->SetMCtruthEventHandler(handler);
  }
  

}
