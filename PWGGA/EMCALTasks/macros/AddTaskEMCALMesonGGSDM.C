void AddTaskEMCALMesonGGSDM(Int_t calibration = 0) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAliEMCALMesonGGSDM", "No analysis manager to connect to.");
    return NULL;
  }  

  //#####################################################
  // Private Recalibrator:
     Int_t GoodTasks  [4] = {1, 0,0,0};
     Int_t RecalScheme[4] = {calibration, 3,5,7};
  //
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

  //#####################################################  
  // MC Generator Types: 
  //
     char MyMCType[50] = "";
  // doesn't work yet... set below!!!!!! 
  //
  // for all Primary particles, assign "ANY"
  // 
  // lhc13e7: 
  //         hijing_0
  //         pi0_1    - flat pT, number varies with centrality
  //         eta_2    - flat pT, number varies with centrality
  //         pi0EMC_3 - flat pT, 1 pi0 into EMCal
  //         pi0PHS_4 - flat pT, 1 pi0 into PHOS
  //         etaEMC_5 - flat pT, 1 eta into EMCal
  //         etaPHS_6 - flat pT, 1 eta into PHOS 
  //
  // lhc12i3: 
  //         Pythia
  //         BOX      - flat pT, 10 pi0 into 2pi
  //         BOX      - flat pT, 10 eta into 2pi
  //         PARAM    - flat pT, 1 pi0 into EMCal
  //         PARAM    - flat pT, 1 pi0 into PHOS 
  //         PARAM    - flat pT, 1 eta into EMCal
  //         PARAM    - flat pT, 1 eta into PHOS 
  //
  // dangerous... right now, it will only take the first BOX that it finds. 
  //
  //#####################################################  

  AliAnalysisTaskEMCALMesonGGSDM* task[4];
  AliAnalysisDataContainer*  coutput[4];  
  char saythis[500];

  for(int i=0; i<4; i++){
    if(GoodTasks[i]==0)
      continue;
    
    sprintf(saythis,"EMCALMesonGGSDM_%d",i);
    task[i] = new AliAnalysisTaskEMCALMesonGGSDM(saythis);
    task[i]->SelectCollisionCandidates(AliVEvent::kMB);//LHC11a
    //task[i]->SelectCollisionCandidates(AliVEvent::kINT7);//LHC13b/c, LHC11c?
    //task[i]->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task[i]->SetRecalScheme(RecalScheme[i]);

    task[i]->SetdRmin_ClustTrack(0.025);//not used at the moment! (Gustavo uses 0.04)
    task[i]->SetFidPhiMinMax(1.39626, 3.15);// from: emc->GetArm1PhiMin()*TMath::DegToRad()=1.39626
    //task[i]->SetFidPhiMinMax(1.39626, 2.10); //pPb no TRD
    //task[i]->SetFidPhiMinMax(2.10, 3.15); //pPb with TRD
    task[i]->SetFidEtaMinMax(-0.65, 0.65);

    //task[i]->SetMyMCType(**MyMCType);
    //task[i]->SetMyMCType("hijing_0");
    task[i]->SetMyMCType("Pythia");
    
    mgr->AddTask(task[i]);
    sprintf(saythis,"cont_AliAnalysisTaskEMCALMesonGGSDM_%d",RecalScheme[i]);
    coutput[i] = 
      mgr->CreateContainer(saythis,
			   TList::Class(),
			   AliAnalysisManager::kOutputContainer,
			   "EMCALMesonGGSDM.root");
    mgr->ConnectInput (task[i],0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[i],1,coutput[i]);
        
    //AliMCEventHandler* handler = new AliMCEventHandler;
    //handler->SetReadTR(kFALSE);
    //mgr->SetMCtruthEventHandler(handler);
  }
  
}
