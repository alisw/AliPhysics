AliAnalysisTaskSEHFCJqa* AddTaskSEHFCJQA(TString fileout="standard",Int_t readmc=kFALSE,TString cutfile="HFJetVertexCuts.root",TString containerprefix="c",Float_t minC=0., Float_t maxC=7.5)
{  
  //
  // Configuration macro for the task to analyze the fraction of prompt charm
  // using the D0 impact parameter
  // andrea.rossi@ts.infn.it
  //
  //==========================================================================
  //+ all electrons -jet study: pt, eta, phi, deltaphi


  Int_t last=0;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCharmFraction", "No analysis manager to connect to.");
    return NULL;
  }   
  
  TString str,containername;
  if(fileout=="standard"){
    fileout=AliAnalysisManager::GetCommonFileName();
    fileout+=":PWG3_HFCJ_";
    fileout+="QA";
    if(containerprefix!="c")fileout+=containerprefix;
    str="HFCJqa";
  }
  else {
    str=fileout;
    str.ReplaceAll(".root","");
  }
  str.Prepend("_");

  AliAnalysisTaskSEHFCJqa *hfTask;
  
  AliRDHFCutsD0toKpi *cuts;

  if(!gSystem->AccessPathName(cutfile.Data(),kFileExists)){
    TFile *f=TFile::Open(cutfile.Data());
    cuts= (AliRDHFCutsD0toKpi*)f->Get("EventTrackCuts");
    cuts->PrintAll();
    hfTask = new AliAnalysisTaskSEHFCJqa("AliAnalysisTaskSEHFCJqa");
    hfTask->SetCutObject(cuts);
  }
  else {// possible (!not standard!) selection for pp2012 data triggered with EMCAL
    cuts=new AliRDHFCutsD0toKpi();
    hfTask = new AliAnalysisTaskSEHFCJqa("AliAnalysisTaskSEHFCJqa");
   
    // AliRDHFJetsCutsVertex *cuts2=new AliRDHFJetsCutsVertex("jetCuts");

    //cuts for jets
    //    cuts->SetJetRadius(0.4);
    //    cuts->SetMaxEtaJet(0.5);//0.9-R
    //    cuts->SetMinPtJet(0);
    //    cuts->SetMaxPtJet(200);
    //    cuts->ResetMaskAndEnableMBTrigger();
    //cuts->SetUseAnyTrigger();
    //cuts->SetTriggerMask(0);
    cuts->SetTriggerMask(AliVEvent::kEMC1 | AliVEvent::kEMC7 | AliVEvent::kEMC8);
    cuts->SetTriggerClass("");

    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetMinNClustersITS(2);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetPtRange(1,1.e10);
    

    if(minC>0&&minC<maxC){
      // Pb-Pb
      cuts->SetTriggerClass("");
      cuts->ResetMaskAndEnableMBTrigger();
      cuts->EnableCentralTrigger();
      cuts->EnableSemiCentralTrigger();
      cuts->SetUseCentrality(AliRDHFCuts::kCentV0M);
      cuts->SetMinCentrality(minC);
      cuts->SetMaxCentrality(maxC);
      
    }

    cuts->AddTrackCuts(esdTrackCuts);
    //cuts for vertexing

    //    ::Error("AddTaskSEHFJets","No Cut Object");
    
  }
  hfTask->SetFilterBit(AliAODTrack::kTrkGlobalNoDCA);
  //  hfTask->SetLoadJet(2,"clustersAOD_ANTIKT04_B0_Filter00768_Cut00150_Skip00");
  //clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip02");
 hfTask->SetLoadJet(2,"clustersAOD_ANTIKT04_B0_Filter00768_Cut00150_Skip02");


  hfTask->SetCutObject(cuts);


  //  hfTask->SetJetCuts(cuts);
  //  hfTask->SetElectronCuts(cuts);
  //  hfTask->SetFilterBitElectron(AliAODTrack::kTrkGlobalNoDCA);
  //  hfTask->SetReadMC(readmc);

  
  mgr->AddTask(hfTask);
 
 
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput =   mgr->GetCommonInputContainer();
  //mgr->CreateContainer("cinput",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(hfTask,0,cinput);
  

  //Now container for general properties histograms
  containername="outputNentries";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer(containername.Data(),TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,1,coutputNentries);


  containername="listTrackProperties";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistTrackProp = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,2,clistTrackProp);

  containername="listJetsProperties";
  containername.Prepend(containerprefix.Data());
  containername.Append(str.Data());
  AliAnalysisDataContainer *clistJetsProp = mgr->CreateContainer(containername.Data(),TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   fileout.Data());
  
  mgr->ConnectOutput(hfTask,3,clistJetsProp);

  return hfTask;
}
