//=========================================================================//
//                                                                         //
//           Analysis AddTask for Particle Ratio Fluctuation Study         //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch                   //
//                       Thu Dec 19 09:09:38 CET 2013
//                            Aaded on 2015/02/07
//                                                                         //
//=========================================================================//

AliEbyENetChargeFluctuationTask *AddAliEbyENetChargeFluctuationTaskSW(const Char_t *centralityEstimator = "V0M",
					  Bool_t isModeAOD    = 1,
					  Int_t  aodFilterBit = 768, 
					  Int_t  sysType      = 0,   // 0-pp,1-pA,2-AA, 3,4,5 mc
					  Int_t  cuttype      = 9,   // esd cut type
					  
					  Int_t pidtype       = 2, 
					  Int_t requestTofPid = 1,
					  Double_t nSigmaCut  = 3.,
					  Double_t lptfortof  = 0.5,
					
					  Double_t ptl        = 0.5, 
					  Double_t pth        = 5.0, 
					  
					  Double_t gEta       = 0.8,
					  Double_t gRap       = 0.5,

					  Double_t dcaxy     = 2.4,
					  Double_t dcaz      = 3.2,
					  
					  Double_t vz         = 10.,
					  Int_t nSample       = 25,
					  Int_t analType      = 1,
					  Int_t isBasic       = 0,
					  const Char_t *taskname="ND") {
  
  Double_t vx = 3.; Double_t vy = 3.; 
 

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }

  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);
  if (isMC)
    Info("AddTaskNetParticle", "This task has MC.");

  AliEbyENetChargeFluctuationTask *task = new AliEbyENetChargeFluctuationTask(taskname);
  if (!task) {
    Error("EbyEPidRatio", "Task could not be created.");
    return NULL;
  }


  Printf("============================== I am here very great ================================");

  if (isMC) task->SetIsMC();

  Int_t sysii = sysType;
  if (sysType > 2) {task->SetIsMC(); sysii = sysType - 3;}

  
  if (isModeAOD) {
    task->SetIsAOD();                       
    task->SetTrackFilterBit(aodFilterBit);
  }
  task->SetSystemType(sysii);
  task->SetEventSelectionBit(AliVEvent::kMB);
  task->SetCentralityEstimator(centralityEstimator);
  task->SetVertexDiamond(vx,vy,vz);
  task->SetKinematicsCuts(ptl,pth,gEta,gRap);
  task->SetNSubSamples(nSample);
  task->SetDca(dcaxy,dcaz);
  if (!isModeAOD) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/PIDFluctuation/macros/configureNetChargeTrackCut.C"); 
    // gROOT->LoadMacro("./configureNetChargeTrackCut.C"); 
    AliESDtrackCuts *cuts = configureNetChargeTrackCut(taskname,cuttype,10001006, gEta, 5.,5.); 
     task->SetAnalysisCutObject(cuts);
  }

  AliHelperPID* help = new AliHelperPID();
  help->SetNSigmaCut(nSigmaCut);
  help->SetPIDType(pidtype); // kNSigmaTPC,kNSigmaTOF, kNSigmaTPCTOF
  if (requestTofPid) {
    help->SetfRequestTOFPID(requestTofPid);
    if (ptl != 0 ) help->SetfPtTOFPID(lptfortof);
  }

  if (sysType > 2) help->SetisMC(1); else help->SetisMC(0);

  if (pidtype == 3){
    AliPIDCombined *pidc=new AliPIDCombined();
    pidc->SetDefaultTPCPriors();
    help->SetPIDCombined(pidc);
  }
  task->SetHelperPID(help);
  task->SetAnal(analType,isBasic);
  mgr->AddTask(task);

    
  TString commonname   = Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname);
  
  AliAnalysisDataContainer *cout1 
    = mgr->CreateContainer(Form("%s_phy",taskname), TList::Class(),  
			   AliAnalysisManager::kOutputContainer, commonname);
  AliAnalysisDataContainer *cout2 
    = mgr->CreateContainer(Form("%s_qa",taskname), TList::Class(),  
			   AliAnalysisManager::kOutputContainer, commonname);
  AliAnalysisDataContainer *cout3 
    = mgr->CreateContainer(Form("%s_dca",taskname), TList::Class(),  
			   AliAnalysisManager::kOutputContainer, commonname);
  AliAnalysisDataContainer *cout4 
    = mgr->CreateContainer(Form("%s_eff",taskname), TList::Class(),  
			   AliAnalysisManager::kOutputContainer, commonname);
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout1);
  mgr->ConnectOutput(task, 2, cout2);
  mgr->ConnectOutput(task, 3, cout3);
  mgr->ConnectOutput(task, 4, cout4);
  
  return task;
}
