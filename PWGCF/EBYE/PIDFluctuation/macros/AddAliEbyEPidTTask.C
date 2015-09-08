//=========================================================================//
//                                                                         //
//           Analysis AddTask for Particle Ratio Fluctuation Study         //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch 
//                        Thu Jun 19 11:44:51 CEST 2014
//                                                                         //
//=========================================================================//


AliEbyEPidTTask *AddAliEbyEPidTTask(Bool_t isModeAOD    = 0,
				    Int_t aodFilterBit  = 768, 
				    Int_t cuttype       = 9,
				    Bool_t IsKMbOnly    = 0,
				    Bool_t IisMC        = 0,

				    Int_t pidtype       = 2, 
				    Int_t requestTofPid = 1,
				    Double_t nSigmaCut  = 3.,
				    Double_t lptfortof  = 0.9,
				    
				    Double_t ptl        = 0.2, 
				    Double_t pth        = 5.0, 
				    Double_t gEta       = 0.8,
				    
				    Double_t dcaxy      = 2.4,
				    Double_t dcaz       = 3.2,
				    
				    Double_t vz         = 15.,
				    
				    TString ctaskname = "D2011") {
  
 
  Double_t vx = 5.; Double_t vy = 5.;
 
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }
  
  Bool_t isMC = IisMC;

 // Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);
  if (isMC)
    Info("AddTaskNetParticle", "This task has MC.");
  
   TString taskname = ctaskname;


  AliEbyEPidTTask *task = new AliEbyEPidTTask(taskname.Data());
  if (!task) {
    Error("EbyEPidRatio", "Task could not be created.");
    return NULL;
  }
  
  if (isMC) task->SetIsMC();
  
  if (isModeAOD) {
    task->SetIsAOD();                       
    task->SetAODtrackCutBit(aodFilterBit);
  }
  task->SetVertexDiamond(vx,vy,vz);
  task->SetKinematicsCuts(ptl,pth,gEta);
  task->SetDca(dcaxy,dcaz);
  if (IsKMbOnly) task->SetIsKMb();
  if (!isModeAOD) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/PIDFluctuation/macros/configureNetChargeTrackCut.C"); 
    AliESDtrackCuts *cuts = configureNetChargeTrackCut(taskname,cuttype,10001006, gEta, dcaxy,dcaz); 
    task->SetAnalysisCutObject(cuts);
  }
  
  AliHelperPID* help = new AliHelperPID();
  help->SetNSigmaCut(nSigmaCut);
  help->SetPIDType(pidtype); // kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF
  if (requestTofPid) {
    help->SetfRequestTOFPID(requestTofPid);
    if (ptl != 0 ) help->SetfPtTOFPID(lptfortof);
  }
  
//  if (isMC) help->SetisMC(1); 
//  else help->SetisMC(0);
  
  if (pidtype == 3){
    AliPIDCombined *pidc=new AliPIDCombined();
    pidc->SetDefaultTPCPriors();
    help->SetPIDCombined(pidc);
  }
  task->SetHelperPID(help);
  mgr->AddTask(task);
  
  TString basefilename = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *coutqa = mgr->CreateContainer(Form("%s_QA",taskname.Data()),TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s",basefilename.Data()));
  AliAnalysisDataContainer *coutt = mgr->CreateContainer("Event",TTree::Class(),AliAnalysisManager::kOutputContainer, Form("%s",basefilename.Data()));

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutqa);
  mgr->ConnectOutput(task, 2, coutt);
  return task;
}

