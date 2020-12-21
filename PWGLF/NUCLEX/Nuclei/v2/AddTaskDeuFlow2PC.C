AliAnalysisTaskDeuFlow2PC *AddTaskDeuFlow2PC( Bool_t krunMCtruth  = kFALSE,
					      Int_t year = 2017,
					      //Bool_t krunMCtruth  = kTRUE,
					      Bool_t kusecontainer = kFALSE, 
					      //TString CollidingSystem = "PbPb",
					      TString CollidingSystem = "pPb",
					      //TString CollidingSystem = "pPb",
					      Float_t  centMin = 0.,
					      Float_t  centMax = 100.,
					      //Int_t  centMax = 80,
					      const char* outfoldersuffix = "HH",
					      const char* outlistsuffix = "HH",
					      //Bool_t bCentralTrigger = kTRUE,
					      Int_t  fistpart = 1, //2 pi 3 ka //4 p 
					      Int_t  secondpart = 5,
					      
					      Float_t nsigmatpcpidfirst     = 3.,//3.,
					      Float_t nsigmatpctofpidfirst  = 3.,//3.,
					      Float_t nsigmatpcpidsec       = 3., 
					      Float_t nsigmatpctofpidsecond = 3.,
					      Int_t  trackbuffersize = 20200,
					      Int_t  multfirstpart = 3000, // 3000 for pions, 1000 ok for protons
					      Int_t  multsecondpart= 1000,
					      
					      //Bool_t kapplyttc = kTRUE,
					      Bool_t kapplyttc = kFALSE,
					      Float_t dphismin = 0.2,    // 0.04 TPC 0.06 global //paper k = 0.017 p = 0.045// 0.045 for pi paper, 0.045 for protons, 0.015 for pi my cut
					      Float_t detasmin = 0.02,     // 0.1 TPC  0.15 global //paper k = 0.02 p = 0.01
					      Float_t radius = 1.2,       // radius whre phistar is evaluated
					      
					      short nevmixing = 7,//10, //7, //nev mix
					      
					      Float_t momemtumlimitforTOFPIDfirst = 0.4, //0.4 for k //paper 0.5
					      Float_t momemtumlimitforTOFPIDsecond= 0.8, // 0.8 for p
					      Bool_t useHMtrigger = kFALSE,
					      Int_t filterBit = 4,
					      Bool_t kusecrrfindratiocut = kFALSE,
					      //Bool_t kusetpcip = kFALSE,
					      Bool_t kusetpcip = kTRUE,
					      Float_t cutipxyprim = 2.4,      // TPC 2.4 protons 1. pions GLOBAL 0.1 protons pions 0.1 
					      Float_t cutipzprim = 3.2,       // TPC 3.2 protons 1. pions GLOBAL 0.15 protons pions 0.15
					      Float_t cutipxysec = 2.4,       // TPC 2.4 protons 1. pions GLOBAL 0.1 protons pions and kaon
					      Float_t cutipzsec  = 3.2,       // TPC 3.2 protons 1. pions GLOBAL 0.15 protons pions and kaon
					      Float_t minptforprim = 0.15,    // 0.4 protons  0.14 pions 0.15 ka
					      Float_t maxptforprim = 30,     // 3.  for protons; 1.pions ;2. for Ka //sharp cuts 1.4 for K and 1.8 for p
					      Float_t minptforsec = 0.2,      //
					      Float_t maxptforsec = 7.0,  
					      Bool_t kpropagateglobal = kTRUE,
					      Bool_t kDoSphericity = kFALSE,
					      //Bool_t kDoSphericity = kTRUE,
					      const char* outlistsubwagon = ""  // "pXi",
					      ) {
  


  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskhDFemto", "No analysis manager found.");
    return 0;
  }


  // Check the analysis type using the event handlers connected to the analysis
  // manager. The availability of MC handler cann also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskhDFemto", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  //cout << "Found " <<type << " event handler" << endl;


  // Create and configure the task, add it to manager.
  //===========================================================================
  TString combinedName;
  combinedName.Form("TaskhDFemto%s%s", outlistsuffix, outlistsubwagon);
  AliAnalysisTaskDeuFlow2PC *task = new AliAnalysisTaskDeuFlow2PC(combinedName);

  // // FIXME : migliora e pp option
  // if(bCentralTrigger)
  //task->SelectCollisionCandidates(AliVEvent::kINT7 |AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  // else
  //task->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kMB);

  //  task->SelectCollisionCandidates(AliVEvent::kMB);
  //  task->SelectCollisionCandidates(AliVEvent::kINT7);
  // task->SetTriggerMask(AliVEvent::kINT7);

  task->SetAnalysisType (type);
  task->SetReadMCTruth(krunMCtruth);
  task->SetUseContainer(kusecontainer);
  task->SetCollidingSystem(CollidingSystem);
  task->SetCentrality(centMin,centMax);
  task->SetFirstParticle(fistpart);
  task->SetSecondParticle(secondpart);
  task->SetTrackBufferSize(trackbuffersize);
  task->SetFirstPartMaxMult(multfirstpart);
  task->SetSecondPartMaxMult(multsecondpart);
  task->SetnSigmaTPCPIDfirstParticle(nsigmatpcpidfirst);
  task->SetnSigmaTPCTOFPIDfirstParticle(nsigmatpctofpidfirst);
  task->SetnSigmaTPCPIDsecondParticle(nsigmatpcpidsec);
  task->SetnSigmaTPCTOFPIDsecondParticle(nsigmatpctofpidsecond);
  task->SetApplyTtc(kapplyttc);
  task->SetDphisMin(dphismin);
  task->SetDetasMin(detasmin);
  task->SetNEventsToMix(nevmixing);
  task->SetMomentumLimitForTOFPIDfirst(momemtumlimitforTOFPIDfirst);
  task->SetMomentumLimitForTOFPIDsecond(momemtumlimitforTOFPIDsecond);
  task->SetApplyRatioCrRnFindCut(kusecrrfindratiocut);
  task->SetCutOnTPCIP(kusetpcip); 
  task->SetIPCutxyPrim(cutipxyprim); 
  task->SetIPCutzPrim(cutipzprim); 
  task->SetIPCutxySec(cutipxysec); 
  task->SetIPCutzSec(cutipzsec); 
  task->SetMinPtPrim(minptforprim);
  task->SetMaxPtPrim(maxptforprim);
  task->SetMinPtSec(minptforsec);
  task->SetMaxPtSec(maxptforsec);
  task->SetPropagateGlobal(kpropagateglobal);
  task->DoSphirocity(kDoSphericity);
  task->SetRadius(radius);
  task->SetYear(year);

  task->SetHMTrigger(useHMtrigger);
  task->SetFilterBit(filterBit);
 
  mgr->AddTask(task);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();

  outputfile += Form(":PWGCFFEMTO_outputHHTask_%s",outfoldersuffix); 
  AliAnalysisDataContainer *cout_HH  = mgr->CreateContainer(combinedName,  TList::Class(),
							     AliAnalysisManager::kOutputContainer,outputfile);
  
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("Signal", TTree::Class(),AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("Bkg"   , TTree::Class(),AliAnalysisManager::kOutputContainer, outputfile);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout_HH);
  mgr->ConnectOutput (task,  2, coutput2);
  mgr->ConnectOutput (task,  3, coutput3);

  if(!task) {
    Error("AddTaskhDFemto","AliAnalysisTaskhDFemto not created!");
    return;
  }

  return task;


}
