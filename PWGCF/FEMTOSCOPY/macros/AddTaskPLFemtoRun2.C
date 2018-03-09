#include <stdio.h>
#include <string.h>

AliAnalysisTaskSE *AddTaskPLFemtoRun2(Int_t system=0/*0=pp,1=PbPb*/,
					       Bool_t theMCon=kFALSE,Bool_t OnlineSearchV0=kFALSE,
                 TString whichV0="Lambda",TString whichV0region="signal", TString trigger="kINT7",
                 bool isRun1=true, const char *cutVariation = "0")
{
  TString suffix;
  suffix.Form("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) 
    {
      ::Error("AddTaskPLFemtoSpectra", "No analysis manager to connect to.");
      return NULL;
    }

  if (theMCon) {
      // IMPORTANT - SET WHEN USING DIFFERENT PASS
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
            gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/"
                                       "AddTaskPIDResponse.C (kTRUE, kTRUE, "
                                       "kTRUE, \"1\")"));
    } else {
      AliAnalysisTaskPIDResponse *pidResponse =
          reinterpret_cast<AliAnalysisTaskPIDResponse *>(
            gInterpreter->ExecuteMacro(
              "$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C)"));
    }

  printf("CREATE TASK\n");
  // create the task

  AliAnalysisTaskPLFemto *task = new AliAnalysisTaskPLFemto(Form("AliAnalysisPLFemto_%s",cutVariation) ,OnlineSearchV0,whichV0,whichV0region);
  task->SetRun1(isRun1);
  if(trigger == "kINT7") task->SetTrigger(AliVEvent::kINT7);
  else if(trigger == "kHighMultV0") task->SetTrigger(AliVEvent::kHighMultV0);
  else if(trigger == "kMB") task->SetTrigger(AliVEvent::kMB);
  else {
      std::cout << "Requested trigger type unknown \n";
      return nullptr;
    }
  task->SelectCollisionCandidates(task->GetTrigger());
  if (task->GetTrigger() == AliVEvent::kHighMultV0) task->SetV0Percentile(0.1);

  if(suffix == "0") task->SetSystematics(AliFemtoCutValues::kDefault);
  else if(suffix == "1") {
    task->SetSystematics(AliFemtoCutValues::kProtonVariationLowerPtThresholdDown);
    task->SetIsLightweight(true);
  }
  else if(suffix == "2") {
    task->SetSystematics(AliFemtoCutValues::kProtonVariationLowerPtThresholdUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "3") {
    task->SetSystematics(AliFemtoCutValues::kProtonVariationEtaRangeUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "4") {
    task->SetSystematics(AliFemtoCutValues::kProtonVariationEtaRangeDown);
    task->SetIsLightweight(true);
  }
  else if(suffix == "5") {
    task->SetSystematics(AliFemtoCutValues::kProtonVariationNsigmaUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "6") {
    task->SetSystematics(AliFemtoCutValues::kProtonVariationNsigmaDown);
    task->SetIsLightweight(true);
  }
  else if(suffix == "7") {
    task->SetSystematics(AliFemtoCutValues::kProtonTPCClusterUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "8") {
    task->SetSystematics(AliFemtoCutValues::kProtonVariationFilterBitGlobal);
    task->SetIsLightweight(true);
  }
  else if(suffix == "9") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationLowerPtThresholdDown);
    task->SetIsLightweight(true);
  }
  else if(suffix == "10") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationLowerPtThresholdUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "11") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationCosinePointingUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "12") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationNsigmaDown);
    task->SetIsLightweight(true);
  }
  else if(suffix == "13") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationTPCClusterUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "14") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationEtaRangeUp);
    task->SetIsLightweight(true);
  }
  else if(suffix == "15") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationEtaRangeDown);
    task->SetIsLightweight(true);
  }
  else if(suffix == "16") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationDCAatV0DecayDown);
    task->SetIsLightweight(true);
  }
  else if(suffix == "17") {
    task->SetSystematics(AliFemtoCutValues::kV0VariationDCADaughtersToPVUp);
    task->SetIsLightweight(true);
  }
  else {
      std::cout << "Unknown cut variation " << suffix << "\n";
      return nullptr;
    }

  task->SetDebugLevel(0);
  task->SetMC(theMCon);

  mgr->AddTask(task);

  // Create and connect containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if(trigger == "kHighMultV0") suffix += "_HM";
  outputfile += ":PWGCF_PLFemto_";
  outputfile += suffix;

  // ------ input data ------
  TString input = "infemto";
  input += suffix;
  TString output1 = "Evtinfo_";
  output1 += suffix;
  TString outputEvt = "AliEventCuts_";
  outputEvt += suffix;
  TString output2 = "SPdir_";//directory for single particle quantities
  output2 += suffix;
  TString output3 = "PIDdir_";//directory for PID quantities
  output3 += suffix;
  TString output4 = "TPdir_";//directory for two particle quantities
  output4 += suffix;

  AliAnalysisDataContainer *cinput0  =  mgr->CreateContainer(input,TChain::Class(),AliAnalysisManager::kInputContainer);

 // ----- output data -----

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output1,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output2,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output3,TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(output4,TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(outputEvt,TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);
  mgr->ConnectOutput(task,5,coutput5);

  return task;
}



