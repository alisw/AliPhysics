AliAnalysisTaskhCascadeFemto *AddTaskHCascadeFemto ( Bool_t krunMCtruth  = kFALSE,
                                                     Bool_t kusecontainer = kFALSE,
                                                     const char* outfoldersuffix = "pXi",
                                                     const char* outlistsuffix = "wttc",
                                                     Bool_t bCentralTrigger = kTRUE,
                                                     Int_t  fistpart = 1,
                                                     Int_t  secondpart = 2,
                                                     Int_t  xicuts = 0, // kloose, kdefault, ktight, if none of those then reco cuts
                                                     Bool_t kcascadesidebands = kFALSE,
                                                     Float_t masswincasc = 0.003, // 0.003 2s 0.008 4s, this always for side bands!!!
                                                     Float_t nsigmatpcpidfirst = 3.,
                                                     Float_t nsigmatpctofpidfirst = 3.,
                                                     Float_t nsigmatpcpidsecdau = 4.,
                                                     Int_t  trackbuffersize = 20200,
                                                     Int_t  multfirstpart = 3000, // 3000 for pions, 1000 ok for protons
                                                     Int_t  multsecondpart= 20,

                                                     Bool_t kapplyttc = kFALSE,
                                                     Float_t dphismin = 0.06,    // 0.04 TPC 0.06 global // 0.045 for pi paper, 0.045 for protons
                                                     Float_t detasmin = 0.15,    // 0.1 TPC  0.15 global

                                                     short nevmixing = 5,

						     Float_t momemtumlimitforTOFPID = 0.75,
                                                     Bool_t kusecrrfindratiocut = kFALSE,
                                                     Bool_t kusecrrowcut = kFALSE,
                                                     Bool_t kusetpcip = kFALSE,
                                                     Float_t cutipxy = 0.1,      // TPC 2.4 protons 1. pions GLOBAL 0.1 protons pions
                                                     Float_t cutipz = 0.15,      // TPC 3.2 protons 1. pions GLOBAL 0.15 protons pions
                                                     Float_t minptforprim = 0.7, // 0.7 protons  0.14 pions
                                                     Float_t maxptforprim = 4.,  // 4.  proton s 2.   pions

                                                     Float_t minptforcasc = 0.8, 
                                                     Float_t maxptforcasc = 10.,  
                                                     Float_t cutipbac = 0.1,     //0.03,
                                                     Bool_t kapplyycutcasc = kFALSE,
                                                     Bool_t kpropagateglobal = kTRUE,
                                                     Bool_t kpropagatefixedr = kFALSE,
                                                     Bool_t kcutonttcprop = kFALSE,
                                                     Bool_t kapplyshareddaughtercut = kTRUE, 
                                                     const char* outlistsubwagon = ""  // "pXi"
                                                   ) {


  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHCascadeFemto", "No analysis manager found.");
    return 0;
  }


  // Check the analysis type using the event handlers connected to the analysis
  // manager. The availability of MC handler cann also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHCascadeFemto", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  //cout << "Found " <<type << " event handler" << endl;


  // Create and configure the task, add it to manager.
  //===========================================================================
  TString combinedName;
  combinedName.Form("TaskhCascadeFemto%s%s", outlistsuffix, outlistsubwagon);
  AliAnalysisTaskhCascadeFemto *task = new AliAnalysisTaskhCascadeFemto(combinedName);

  if(bCentralTrigger)
      task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else
      task->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kMB);

  task->SetAnalysisType (type);
  task->SetReadMCTruth(krunMCtruth);
  task->SetUseContainer(kusecontainer);

  task->SetMassHist(500, 1.,1.5);
  if (secondpart == 2) task->SetMassHistGrandmother (350,1.1,1.8);
  else if (secondpart == 3) task->SetMassHistGrandmother (120,1.62, 1.74);
 
  task->SetCentrality(0., 90.);

  task->SetFirstParticle(fistpart);
  task->SetSecondParticle(secondpart);
  task->SetTrackBufferSize(trackbuffersize);
  task->SetFirstPartMaxMult(multfirstpart);
  task->SetSecondPartMaxMult(multsecondpart);
  task->SetMassWindowCascades(masswincasc);
  task->SetnSigmaTPCPIDfirstParticle(nsigmatpcpidfirst);
  task->SetnSigmaTPCTOFPIDfirstParticle(nsigmatpctofpidfirst);
  task->SetnSigmaTPCPIDsecondParticleDau(nsigmatpcpidsecdau);
  task->SetApplyTtc(kapplyttc);
  task->SetDphisMin(dphismin);
  task->SetDetasMin(detasmin);
  task->SetCascadeSideBands(kcascadesidebands);
  task->SetXiCuts(xicuts);
  task->SetNEventsToMix(nevmixing);
  task->SetMomentumLimitForTOFPID(momemtumlimitforTOFPID);
  task->SetApplyRatioCrRnFindCut(kusecrrfindratiocut);
  task->SetApplyCrossedRowCut(kusecrrowcut);
  task->SetCutOnTPCIP(kusetpcip); 
  task->SetIPCutxy(cutipxy); 
  task->SetIPCutz(cutipz); 
  task->SetMinPtPrim(minptforprim);
  task->SetMaxPtPrim(maxptforprim);
  task->SetMinPtCasc(minptforcasc);
  task->SetMaxPtCasc(maxptforcasc);
  task->SetIPCutBac(cutipbac);
  task->SetApplyYcutCasc(kapplyycutcasc);
  task->SetPropagateGlobal(kpropagateglobal);
  task->SetPropagateAtFixedR(kpropagatefixedr);
  task->SetCutOnttcProp(kcutonttcprop);
  task->SetApplySharedDaughterCutXi(kapplyshareddaughtercut);

  mgr->AddTask(task);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();

  outputfile += Form(":PWGCFFEMTO_outputHCascadeTask_%s",outfoldersuffix); 
  AliAnalysisDataContainer *cout_pXi  = mgr->CreateContainer(Form("%s_%s",combinedName.Data(),outfoldersuffix),  TList::Class(),
                                                               AliAnalysisManager::kOutputContainer,outputfile);

//  AliAnalysisDataContainer *cout_pXi2 = mgr->CreateContainer(Form("cfcontCutsXi_%s_%s",outfoldersuffix,combinedName.Data()),
//                                                              AliCFContainer::Class(),
//                                                              AliAnalysisManager::kOutputContainer,
//                                                              outputfile );


  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout_pXi);
//  mgr->ConnectOutput(task, 2, cout_pXi2);

  if(!task) {
            Error("AddTaskhCascadeFemto","AliAnalysisTaskhCascadeFemto not created!");
            return;
  }

  return task;


}
