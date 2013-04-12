//=========================================================================//
//                                                                         //
//           Analysis AddTask for Particle Ratio Fluctuation Study         //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch                   //
//                       Fri Apr 11 13:31:07 CEST 2013
//                                                                         //
//=========================================================================//

void AddAliEbyEParticleRatioFluctuationTask(Double_t vz=10,Double_t ptl=0.5, Double_t pth=5, Int_t AODfilterBit = 768, Int_t pidtype = 2, Int_t requestTofPid = 1, Double_t nSigmaCut = 3., Int_t ikey = 0, TString analdata = "AOD", TString analtype = "PbPb", TString centralityEstimator = "V0M", TString ctaskname = "2011") {

  Double_t vx = 3.; Double_t vy = 3.;

  TString taskname = "EbyECF_";
  taskname += ctaskname;
  taskname += "_";
  taskname += analdata;
  taskname += "_";
  taskname += analtype;
  taskname += "_";
  taskname += centralityEstimator;
  taskname += "_";
  taskname += Form("%d",AODfilterBit);
  taskname += "_";
  taskname += Form("pt_%.1f_%.1f", ptl, pth);
  taskname += "_";

  Bool_t isMC = 0;  
  if(analdata == "AODMC" || analdata == "MC") 
    isMC = 1;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); 

 
  TString basefilename = AliAnalysisManager::GetCommonFileName();
  
  
  AliHelperPID* help = new AliHelperPID();
  help->SetNSigmaCut(nSigmaCut);
  help->SetPIDType(pidtype);        
  help->SetfRequestTOFPID(requestTofPid);
  help->SetfPtTOFPID(ptl);
  help->SetisMC(isMC);
  
  if(ikey == 0 ) {
  AliEbyEParticleRatioFluctuationTask *task[8];
  AliAnalysisDataContainer *cout[8];
  for(Int_t i = 0; i < 8 ; i ++) {
    Double_t eta = 0.1 + 0.1*i;
    TString taskname1 = taskname;
    taskname1 += Form("eta_%.2f",eta);
    task[i] = new AliEbyEParticleRatioFluctuationTask(taskname1.Data());
    task[i]->SetVertexDiamond(vx,vy,vz);
    task[i]->SetAODtrackCutBit(AODfilterBit);
    task[i]->SetKinematicsCuts(ptl,pth,eta);
    task[i]->SetHelperPID(help);
    mgr->AddTask(task[i]);
    
    cout[i] = mgr->CreateContainer(Form("%s",taskname1.Data()),TList::Class(), 
				   AliAnalysisManager::kOutputContainer,
				   Form("%s:CFEbyE_PR",basefilename.Data()));
    mgr->ConnectInput(task[i], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[i], 1, cout[i]);
    
  }
   
  AliEbyEParticleRatioFluctuationTask *taskqa;
  taskqa = new AliEbyEParticleRatioFluctuationTask("QACFEbyEPR");
  //  taskqa->SetPIDMethod(1);
  taskqa->RunQA();
  taskqa->SetVertexDiamond(vx,vy,vz);
  taskqa->SetAODtrackCutBit(AODfilterBit);
  taskqa->SetKinematicsCuts(ptl,pth,0.8);
  taskqa->SetHelperPID(help);
  mgr->AddTask(taskqa);
  
  AliAnalysisDataContainer *coutqa 
    = mgr->CreateContainer(Form("QA_%s",taskname1.Data()),TList::Class(),
			   AliAnalysisManager::kOutputContainer,
			   Form("%s:CFEbyE_PR",basefilename.Data()));
  mgr->ConnectInput(taskqa, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskqa, 1, coutqa);
  }
  else if(ikey > 0 && ikey < 9) {
    Double_t eta = 0.1*ikey;
    TString taskname1 = taskname;
    taskname1 += Form("eta_%.2f",eta);
    AliEbyEParticleRatioFluctuationTask *task1 
      = new AliEbyEParticleRatioFluctuationTask(taskname1.Data());
    task1->SetVertexDiamond(vx,vy,vz);
    task1->SetAODtrackCutBit(AODfilterBit);
    task1->SetKinematicsCuts(ptl,pth,eta);
    task1->SetHelperPID(help);
    mgr->AddTask(task1);
    
    AliAnalysisDataContainer *cout1 
      = mgr->CreateContainer(Form("%s",taskname1.Data()),TList::Class(),
			     AliAnalysisManager::kOutputContainer,
			     Form("%s:CFEbyE_PR",basefilename.Data()));
    mgr->ConnectInput(task1, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task1, 1, cout1);
    
   
  AliEbyEParticleRatioFluctuationTask *taskqa1;
  taskqa1 = new AliEbyEParticleRatioFluctuationTask("QACFEbyEPR");
  taskqa1->RunQA();
  taskqa1->SetVertexDiamond(vx,vy,vz);
  taskqa1->SetAODtrackCutBit(AODfilterBit);
  taskqa1->SetKinematicsCuts(ptl,pth,eta);
  taskqa1->SetHelperPID(help);
  mgr->AddTask(taskqa1);
  
  AliAnalysisDataContainer *coutqa1 
    = mgr->CreateContainer(Form("QA_%s",taskname1.Data()),TList::Class(), 
			   AliAnalysisManager::kOutputContainer,
			   Form("%s:CFEbyE_PR",basefilename.Data()));
  mgr->ConnectInput(taskqa1, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskqa1, 1, coutqa1);
  
  } else return;

  return;
}
