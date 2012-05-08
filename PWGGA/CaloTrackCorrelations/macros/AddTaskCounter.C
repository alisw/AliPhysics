//_______________________________________________
void AddTaskCounter(const TString trigger = "MB")
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  AliAnalysisTaskCounter * counter =  new AliAnalysisTaskCounter(Form("Counter%s",trigger.Data()));
  //if(kRun > 140000 && kRun < 146900) counter ->RejectFastCluster();
  //if     (kCollision=="pp"  )   counter->SetZVertexCut(10.);  //Open cut
  //else if(kCollision=="PbPb")   counter->SetZVertexCut(10.);  //Centrality defined in this range.
  
  if(trigger=="EMC7")
  {
    printf("counter trigger EMC7\n");
    counter->SelectCollisionCandidates(AliVEvent::kEMC7);
  }
  else if (trigger=="INT7")
  {
    printf("counter trigger INT7\n");
    counter->SelectCollisionCandidates(AliVEvent::kINT7);
  }
  if(trigger=="EMC1")
  {
    printf("counter trigger EMC1\n");
    counter->SelectCollisionCandidates(AliVEvent::kEMC1);
  }
  else if(trigger=="MB")
  {
    printf("counter trigger MB\n");
    counter->SelectCollisionCandidates(AliVEvent::kMB);
  }
  else if(trigger=="PHOS")
  {
    printf("counter trigger PHOS\n");
    counter->SelectCollisionCandidates(AliVEvent::kPHI7);
  }
  else if(trigger=="PHOSPb")
  {
    printf("counter trigger PHOSPb\n");
    counter->SelectCollisionCandidates(AliVEvent::kPHOSPb);
  }
  else if(trigger=="AnyINT")
  {
    printf("counter trigger AnyINT\n");
    counter->SelectCollisionCandidates(AliVEvent::kAnyINT);
  }  
  else if(trigger=="INT")
  {
    printf("counter trigger AnyINT\n");
    counter->SelectCollisionCandidates(AliVEvent::kAny);
  }
  else if(trigger=="EMCEGA")
  {
    printf("counter trigger EMC Gamma\n");
    counter->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  } 
  else if(trigger=="EMCEJE")
  {
    printf("counter trigger EMC Jet\n");
    counter->SelectCollisionCandidates(AliVEvent::kEMCEJE);
  }
  else if(trigger=="Central")
  {
    printf("counter trigger Central\n");
    counter->SelectCollisionCandidates(AliVEvent::kCentral);
  } 
  else if(trigger=="SemiCentral")
  {
    printf("counter trigger SemiCentral\n");
    counter->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  }
  
  TString outputFile = AliAnalysisManager::GetCommonFileName(); 
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutput = 
  mgr->CreateContainer(Form("Counter%s",trigger.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,  outputFile.Data());
  mgr->AddTask(counter);
  mgr->ConnectInput  (counter, 0, cinput1);
  mgr->ConnectOutput (counter, 1, coutput);
  
}

