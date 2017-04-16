/// \file AddTaskCounter.C
/// \ingroup CaloTrackCorrMacros
/// \brief Counting task configurations.
///
/// Example of configuration of AliAnalysisTaskCounter,
/// simple task counting events
/// Gets the cross section from file pyxsec if requested.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)


///
/// Main method for the configuration.
/// It creates a Counter task, configures it and adds it to the analysis manager.
///
AliAnalysisTaskCounter * AddTaskCounter(const TString trigger = "",
                                        Bool_t xsOn = kFALSE)
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  AliAnalysisTaskCounter * counter =  new AliAnalysisTaskCounter(Form("Counter%s",trigger.Data()));
  //if(kRun > 140000 && kRun < 146900) counter ->RejectFastCluster();
  //if     (kCollision=="pp"  )   counter->SetZVertexCut(10.);  //Open cut
  //else if(kCollision=="PbPb")   counter->SetZVertexCut(10.);  //Centrality defined in this range.
  
  if(xsOn) counter->SwitchOnMCCrossSectionCalculation();
  else     counter->SwitchOffMCCrossSectionCalculation();
  
  if(trigger=="EMC8")
  {
    printf("counter trigger EMC8\n");
    counter->SelectCollisionCandidates(AliVEvent::kEMC8);
  }
  else if(trigger=="EMC7")
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
  else if(trigger=="SemiOrCentral")
  {
    printf("counter trigger SemiCentral Or Central\n");
    counter->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kCentral);
  }
  
  TString outputFile = AliAnalysisManager::GetCommonFileName(); 
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutput = 
  mgr->CreateContainer(Form("Counter%s",trigger.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,  outputFile.Data());
  mgr->AddTask(counter);
  mgr->ConnectInput  (counter, 0, cinput1);
  mgr->ConnectOutput (counter, 1, coutput);
  
  return counter;
  
}

