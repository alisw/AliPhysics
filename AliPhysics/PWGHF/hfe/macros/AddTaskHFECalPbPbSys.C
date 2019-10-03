AliAnalysisTask *AddTaskHFECalPbPbSys(Bool_t MassConst, Bool_t MassWidthCut, Bool_t MassCal, Bool_t MassNonlinear ,Double_t asspTCut, Double_t angleCut, Double_t MassCut, Double_t NsigCut,TString ID="phoSys0")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEECalPbPbSys", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFECalPbPbSys", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFECalPbPbSys", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  
  //analysis task 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/PbPb/ConfigHFECal.C");
  AliAnalysisTaskHFECal *hfetask = ConfigHFECal(MCthere,MassConst,MassWidthCut,MassCal,MassNonlinear,asspTCut,angleCut,MassCut,NsigCut,0);
  mgr->AddTask(hfetask);
  // semi-central
  hfetask->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  TString containerName3 = mgr->GetCommonFileName();
  containerName3 += ":PWGHF_hfeCalSemiCentral";
  containerName3 += ID;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFE_Results_SemiCentral_%s",ID.Data()), TList::Class(),AliAnalysisManager::kOutputContainer, containerName3.Data());
  mgr->ConnectInput(hfetask, 0, cinput);
  mgr->ConnectOutput(hfetask, 1, coutput1);
 
  // central
  AliAnalysisTaskHFECal *hfetask1 = ConfigHFECal(MCthere,MassConst,MassWidthCut,MassCal,MassNonlinear,asspTCut,angleCut,MassCut,NsigCut,0);
  mgr->AddTask(hfetask1);
  hfetask1->SelectCollisionCandidates(AliVEvent::kCentral);
  TString containerName1 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_hfeCalCentral";
  containerName1 += ID;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFE_Results_Central_%s",ID.Data()), TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  mgr->ConnectInput(hfetask1, 0, cinput);
  mgr->ConnectOutput(hfetask1, 1, coutput1);

  //trigger
  AliAnalysisTaskHFECal *hfetaskTrig = ConfigHFECal(MCthere,MassConst,MassWidthCut,MassCal,MassNonlinear,asspTCut,angleCut,MassCut,NsigCut,0);
  mgr->AddTask(hfetaskTrig);
  hfetaskTrig->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeCalkTrig";
  containerName2 += ID;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFE_Results_EMCalTrig_%s",ID.Data()), TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(hfetaskTrig, 0, cinput);
  mgr->ConnectOutput(hfetaskTrig, 1, coutput1);
 

  //MB trigger
  AliAnalysisTaskHFECal *hfetaskMB = ConfigHFECal(MCthere,MassConst,MassWidthCut,MassCal,MassNonlinear,asspTCut,angleCut,MassCut,NsigCut,0);
  mgr->AddTask(hfetaskMB);
  hfetaskMB->SelectCollisionCandidates(AliVEvent::kMB);
  TString containerName4 = mgr->GetCommonFileName();
  containerName4 += ":PWGHF_hfeCalkMB";
  containerName4 += ID;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("HFE_Results_EMCalMB_%s",ID.Data()), TList::Class(),AliAnalysisManager::kOutputContainer, containerName4.Data());
  mgr->ConnectInput(hfetaskMB, 0, cinput);
  mgr->ConnectOutput(hfetaskMB, 1, coutput1);
  

  return NULL;
}
