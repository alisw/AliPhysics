AliAnalysisTask *AddTaskHFECalPbPbSys(int trigtype, int TPCclust, int Nits, double nSigMim, double Mimeop, double Maxeop, int PIDorder){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFE", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFE", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFE", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  cout<<"AddTaskHFE - MC config is: "<<MCthere<<endl;

  //============= Set Task Name ===================
  //TString taskName=("AliAnalysisTaskHFE.cxx+");
  //===============================================
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFECalSys_PbPb.C");

  AliAnalysisTaskHFE *hfetask = ConfigHFECalSys_PbPb(MCthere,TPCclust,Nits,nSigMim,Mimeop,Maxeop,PIDorder);
  //RequestMemory(hfetask, 250*1024);
  if(trigtype==0)hfetask->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  if(trigtype==1)hfetask->SelectCollisionCandidates(AliVEvent::kCentral);
  hfetask->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  mgr->AddTask(hfetask);

  int inSig = (int)(nSigMim*10.0);
  int iMimeop = (int)(Mimeop*10.0);
  int iMaxeop = (int)(Maxeop*10.0);

  //find input container
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString containerName = mgr->GetCommonFileName();
  TString appendix(TString::Format(":PWGHF_hfeCalPbPbEGATPC%dITS%dnSig%dMimeop%dMaxeop%dPID%d",TPCclust,Nits,inSig,iMimeop,iMaxeop,PIDorder));
  if(trigtype==1)appendix(TString::Format(":PWGHF_hfeCentTPC%dITS%dnSig%dMimeop%dMaxeop%dPID%d",TPCclust,Nits,inSig,iMimeop,iMaxeop,PIDorder));
  //containerName += ":PWGHF_hfeCalPbPbEGA";
  containerName += appendix;
  
  hfetask->ConnectOutput(1, mgr->CreateContainer("HFE_Results_EMCALSys", TList::Class(),
						 AliAnalysisManager::kOutputContainer, containerName.Data()));
  hfetask->ConnectOutput(2, mgr->CreateContainer("HFE_QA_EMCALSys", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));

  mgr->ConnectInput  (hfetask,  0, cinput );
  
  /*
  //<--- task2 for central trigger

  AliAnalysisTaskHFE *hfetask2 = ConfigHFECalstandard_PbPb(MCthere,TPCclust);
  hfetask2->SelectCollisionCandidates(AliVEvent::kCentral);
  mgr->AddTask(hfetask2);

  //find input container
  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeCalPbPbCent";
  
  hfetask2->ConnectOutput(1, mgr->CreateContainer("HFE_Results_EMCALCent", TList::Class(),
						 AliAnalysisManager::kOutputContainer, containerName2.Data()));
  hfetask2->ConnectOutput(2, mgr->CreateContainer("HFE_QA_EMCALCent", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName2.Data()));

  mgr->ConnectInput  (hfetask2,  0, cinput );
  */ 

/*
  AliAnalysisTaskHFE *trdtask = ConfigHFEtrd(MCthere);

  //----------------------
  //create data containers
  //----------------------
 
  trdtask->ConnectOutput(1, mgr->CreateContainer("HFE_Results", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  trdtask->ConnectOutput(2, mgr->CreateContainer("HFE_QA", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput  (trdtask,  0, cinput );
*/

  //return hfetask;
  return NULL;
}
