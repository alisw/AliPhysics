//=========================================================================//
//                                                                         //
//           Analysis AddTask for Particle Ratio Fluctuation Study         //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch                   //
//                       Thu Dec 19 09:09:38 CET 2013
//                                                                         //
//=========================================================================//

void AddAliEbyENetChargeFluctuationBinTask(const Char_t *taskname="TOFTPC",
					   const Char_t *centralityEstimator = "V0M",
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
					   Int_t isBasic       = 0) {
  
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


  AliHelperPID* help = new AliHelperPID();
help->SetNameTitle(taskname,taskname);
help->Print();
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

  

  AliEbyENetChargeFluctuationTask *task[8];
  AliAnalysisDataContainer *cout1[8];
  AliAnalysisDataContainer *cout2[8];
  AliAnalysisDataContainer *cout3[8];
  AliAnalysisDataContainer *cout4[8];
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/PIDFluctuation/macros/configureNetChargeTrackCut.C"); 
  AliESDtrackCuts *cuts[8];

  for (Int_t iEta = 0; iEta < 8; iEta++ ) {

    TString taskname_ctr = Form("%s_E%d",taskname,iEta);

    Double_t eEta = 0.1 + iEta * 0.1;

    task[iEta] = new AliEbyENetChargeFluctuationTask(Form("%s_%s",taskname,taskname_ctr.Data()));
    if (!task[iEta]) {
      Error("EbyEPidRatio", "Task could not be created.");
      return NULL;
    }
    
    Printf("=========== Run For the Eta Window %s ==========",taskname_ctr.Data());
    
    if (isMC) task[iEta]->SetIsMC();
    
    Int_t sysii = sysType;
    if (sysType > 2) {task[iEta]->SetIsMC(); sysii = sysType - 3;}
    
    if (isModeAOD) {
      task[iEta]->SetIsAOD();                       
      task[iEta]->SetTrackFilterBit(aodFilterBit);
    }
    task[iEta]->SetSystemType(sysii);
    task[iEta]->SetEventSelectionBit(AliVEvent::kMB);
    task[iEta]->SetCentralityEstimator(centralityEstimator);
    task[iEta]->SetVertexDiamond(vx,vy,vz);
    task[iEta]->SetKinematicsCuts(ptl,pth,eEta,gRap);
    task[iEta]->SetNSubSamples(nSample);
    task[iEta]->SetDca(dcaxy,dcaz);
    if (iEta == 7)  task[iEta]->DoBasicQA();
    
    if (!isModeAOD) {
        cuts[iEta] = configureNetChargeTrackCut(taskname_ctr.Data(),cuttype,10001006, eEta, 5.,5.); 
        task[iEta]->SetAnalysisCutObject(cuts[iEta]);
	Printf("=========== Cut Applied %s ==========",taskname_ctr.Data());
    }
    
    
    task[iEta]->SetHelperPID(help);
    task[iEta]->SetAnal(analType,isBasic);
    mgr->AddTask(task[iEta]);
    
     
    TString commonname   = Form("%s:%s_E%d", AliAnalysisManager::GetCommonFileName(),taskname, iEta);
    
    cout1[iEta] = mgr->CreateContainer(Form("%s_phy",taskname_ctr.Data()), TList::Class(),  
				       AliAnalysisManager::kOutputContainer, commonname);
    cout2[iEta] = mgr->CreateContainer(Form("%s_qa",taskname_ctr.Data()), TList::Class(),  
				       AliAnalysisManager::kOutputContainer, commonname);
    cout3[iEta] = mgr->CreateContainer(Form("%s_dca",taskname_ctr.Data()), TList::Class(),  
				       AliAnalysisManager::kOutputContainer, commonname);
    cout4[iEta] = mgr->CreateContainer(Form("%s_eff",taskname_ctr.Data()), TList::Class(),  
				       AliAnalysisManager::kOutputContainer, commonname);
    
    mgr->ConnectInput(task[iEta],  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[iEta], 1, cout1[iEta]);
    mgr->ConnectOutput(task[iEta], 2, cout2[iEta]);
    mgr->ConnectOutput(task[iEta], 3, cout3[iEta]);
    mgr->ConnectOutput(task[iEta], 4, cout4[iEta]);
  }
  
  return;
}
