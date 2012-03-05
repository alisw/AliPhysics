//---------------------------------------------------------------------//
// Macro to add the task for jet spectrum and jet shape stduies in pp. //
//---------------------------------------------------------------------// 

AliAnalysisTaskJetProperties *AddTaskJetProperties(Char_t* bJet="clustersAOD",
						   Char_t* jetFinder="ANTIKT",
						   Float_t radius=0.4,
						   UInt_t filterMask=256,
						   Float_t ptTrackCut = 0.15){
  
  //*********************************************************************************************//
  //bJet can take the following name (make sure the same branch is ON for JetCluster task)       //
  //      "clustersAOD"                                                                          //
  //      "clustersAODMC2"                                                                       //
  //      "clustersAODMC"                                                                        //
  //      "clustersMCKINE2"                                                                      //
  //      "clustersMCKINE"                                                                       //
  //                                                                                             //
  //      "jetsAOD"                                                                              //
  //      "jetsAODMC2"                                                                           //
  //      "jetsAODMC"                                                                            //
  //      "jetsMCKINE2"                                                                          //
  //      "jetsMCKINE"                                                                           //
  //---------------------------------------------------------------------------------------------//
  // Example to add this task in AnalysisTrainPWG4Jets.C:                                        //
  //---------------------------------------------------------------------------------------------//
  // AliAnalysisTaskJetProperties *taskJetProp = 0;                                              //
  // taskJetProp = AddTaskJetPropertiesPP("clustersAOD","ANTIKT", 0.4, kHighPtFilterMask, 0.15); //
  //*********************************************************************************************//
  

  Int_t debug = -1; // debug level, -1: not set here
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetPropertiesPP", "No analysis manager to connect to.");
    return NULL;
  }
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetPropertiesPP", "This task requires an input event handler");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Printf("########## AddTaskJetProerties: Data Type: %s", type.Data());
  
  TString JetBranch(bJet);   
  TString Method(jetFinder);
  Method.ToUpper();
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskJetProperties *task = new 
    AliAnalysisTaskJetProperties(Form("JetProperties_%s_%s",JetBranch.Data(), Method.Data()));
  
  
  if(debug>=0) task->SetDebugLevel(debug);
  
  // attach the filter mask and options
  Float_t PtTrackMin = ptTrackCut;
  TString cAdd = "";
  cAdd += Form("_%s",Method.Data());
  cAdd += Form("%02d",(int)((radius+0.01)*10.));
  cAdd += Form("_B0");
  cAdd += Form("_Filter%05d",filterMask);
  cAdd += Form("_Cut%05d",(int)((1000.*PtTrackMin)));
  cAdd += Form("_Skip00");
  if(JetBranch.Length())JetBranch += cAdd;
  Printf("########## AddTaskJetProperties: JetBranch '%s'", JetBranch.Data());
  Printf("########## AddTaskJetProperties: JetFinder '%s'", Method.Data());  
  
  if(JetBranch.Length()) task->SetJetBranch(JetBranch);
  //seting the track type by looking jet branch
  //three type of tracks imlmented: kTrackAOD, kTrackKine, and kTrackAODMC 
  //otherwise undifined track type kTrackUndef 
  
  if(JetBranch.Contains("AODMC"))      task->SetTrackType(AliAnalysisTaskJetProperties::kTrackAODMC);
  else if(JetBranch.Contains("MCKINE"))task->SetTrackType(AliAnalysisTaskJetProperties::kTrackKine);
  else if(JetBranch.Contains("AOD"))   task->SetTrackType(AliAnalysisTaskJetProperties::kTrackAOD);
  else task->SetTrackType(AliAnalysisTaskJetProperties::kTrackUndef);//undefined track type
  
  //setting the jet rejection
  //two types implemented:kNoReject and kReject1Trk
  TString contName="";
  Bool_t Is1TrackJetReject = kFALSE;//=kTRUE if you want to reject single track jet
  if(Is1TrackJetReject){//by default no rejection
    task->SetJetRejectType(AliAnalysisTaskJetProperties::kReject1Track);
    contName="_No1TrackJet";
  }
  
  task->SetEventCuts(8.0,2);//VtxZ=+-8cm, nContributors>2
  task->SetFilterMask(filterMask);
  task->SetTrackCuts();// default : pt > 0.150 GeV, |eta|<0.9, full phi acc
  Float_t minJetPt  = 10.0; 
  Float_t minJetEta = -0.5; 
  Float_t maxJetEta =  0.5;
  Int_t tmpR = (Int_t)10*radius;
  if(tmpR==2){minJetEta=-0.7;maxJetEta=0.7;}
  if(tmpR==6){minJetEta=-0.3;maxJetEta=0.3;}
  task->SetJetCuts(minJetPt,minJetEta,maxJetEta);
  
  mgr->AddTask(task);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  
  AliAnalysisDataContainer *coutput_JetProperties = 
    mgr->CreateContainer(Form("jetprop%s_%s", contName.Data(),JetBranch.Data()),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("%s:PWG4_JetProp%s_%s", 
			      AliAnalysisManager::GetCommonFileName(), contName.Data(),JetBranch.Data()));
  
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, coutput_JetProperties);
  
  return task;
}
