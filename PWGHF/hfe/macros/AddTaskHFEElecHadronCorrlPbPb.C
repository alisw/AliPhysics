AliAnalysisTaskElecHadronCorrel *AddTaskHFEElecHadronCorrlPbPb(Bool_t TrigSelCen = kTRUE,
                                                               Double_t CentrMin = 0,
                                                               Double_t CentMax = 7,
                                                               Double_t TPCNsigMinE = -2,
                                                               Double_t TPCNsigMaxE = 2,
                                                               Double_t TPCNsigMinH = -10,
                                                               Double_t TPCNsigMaxH = -3.5,
                                                               Double_t SSM02Min = 0.03,
                                                               Double_t SSM02Max = 0.5,
                                                               Double_t SSM20Min = 0.03,
                                                               Double_t SSM20Max = 0.3,
                                                               Double_t Disp = 1,
                                                               Double_t EovPMin = 0.8,
                                                               Double_t EovPMax = 1.2,
                                                               Double_t InvM = 0.1,
                                                               TString ContNameExt = "Central", 
                                                               TString TaskName = "hfeCorrl"
                                                               )

{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEElecHadronCorrlPbPb", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFEElecHadronCorrlPbPb", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
/*  if (type=="AOD"){
    ::Error("AddTaskHFEElecHadronCorrlPbPb", "The tasks exits because AODs are in input");
    return NULL;
  }
  */
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }

  //analysis task 
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEElecHadronCorrelPbPb.C");
  AliAnalysisTaskElecHadronCorrel *taskHFEeh = ConfigHFEElecHadronCorrelPbPb(MCthere,TrigSelCen,CentrMin,CentMax,TPCNsigMinE,TPCNsigMaxE,TPCNsigMinH,TPCNsigMaxH,SSM02Min,SSM02Max,SSM20Min,SSM20Max,Disp,EovPMin,EovPMax,InvM,TaskName);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeCalPbPbeh";
  containerName += ContNameExt;
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ccontainer0",TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data());

  mgr->ConnectInput(taskHFEeh,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHFEeh,1,coutput3);

  mgr->AddTask(taskHFEeh);

  return taskHFEeh;
}
