//Bool_t AddTrackCutsLHC10h(AliAnalysisTaskESDfilter* esdFilter);
void CreateAOD()
{
  // Main
  
  // LoadLibraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libMatrix.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");


  //Data
  TChain *chain       = new TChain("esdTree") ;
  chain->Add("AliESDs.root");
  
  if(chain){
    AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
    
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
    
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
    mgr->SetMCtruthEventHandler(mcHandler);
    
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("AliAOD.root");
    ////aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);
    
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdHandler);
    
    //mgr->SetDebugLevel(-1); // For debugging, do not uncomment if you want no messages.
    
    //-------------------------------------------------------------------------
    //Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1;
    
    coutput1 = mgr->GetCommonOutputContainer();


   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Barrel tracks filter
   AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
//   esdfilter->SetTimeZeroType(AliESDpid::kTOF_T0);
//   esdfilter->DisableCascades();
//   esdfilter->DisableKinks();
  
   mgr->AddTask(esdfilter);
  
   
   // Filtering of MC particles (decays conversions etc)
   // this task has to go AFTER all other filter tasks
   // since it fills the AODMC array with all
   // selected MC Particles, only this way we have the 
   // AODMCparticle information available for following tasks
   AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
   mgr->AddTask(kinefilter);

   
  Printf("%s%d: Creating Track Cuts for LH10h",(char*)__FILE__,__LINE__);

  // Cuts on primary tracks
  AliESDtrackCuts* esdTrackCutsL = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  
  // ITS stand-alone tracks
  AliESDtrackCuts* esdTrackCutsITSsa = new AliESDtrackCuts("ITS stand-alone Track Cuts", "ESD Track Cuts");
  esdTrackCutsITSsa->SetRequireITSStandAlone(kTRUE);
  
  // Pixel OR necessary for the electrons
  AliESDtrackCuts *itsStrong = new AliESDtrackCuts("ITSorSPD", "pixel requirement for ITS");
  itsStrong->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  
  
  // PID for the electrons
  AliESDpidCuts *electronID = new AliESDpidCuts("Electrons", "Electron PID cuts");
  electronID->SetTPCnSigmaCut(AliPID::kElectron, 3.);
  
  // tighter cuts on primary particles for high pT tracks
  // take the standard cuts, which include already 
  // ITSrefit and use only primaries...
  
  // ITS cuts for new jet analysis 
  //  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
  //  AliESDtrackCuts* esdTrackCutsHG0 = CreateTrackCutsPWGJE(10001006);

  AliESDtrackCuts *jetCuts1006 = new AliESDtrackCuts("AliESDtrackCuts"); 

  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  jetCuts1006->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  jetCuts1006->SetMinNClustersTPC(70);
  jetCuts1006->SetMaxChi2PerClusterTPC(4);
  jetCuts1006->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
  jetCuts1006->SetAcceptKinkDaughters(kFALSE);
  jetCuts1006->SetRequireTPCRefit(kTRUE);
  jetCuts1006->SetMaxFractionSharedTPCClusters(0.4);
  // ITS
  jetCuts1006->SetRequireITSRefit(kTRUE);
  //accept secondaries
  jetCuts1006->SetMaxDCAToVertexXY(2.4);
  jetCuts1006->SetMaxDCAToVertexZ(3.2);
  jetCuts1006->SetDCAToVertex2D(kTRUE);
  //reject fakes
  jetCuts1006->SetMaxChi2PerClusterITS(36);
  jetCuts1006->SetMaxChi2TPCConstrainedGlobal(36);

  jetCuts1006->SetRequireSigmaToVertex(kFALSE);

  jetCuts1006->SetEtaRange(-0.9,0.9);
  jetCuts1006->SetPtRange(0.15, 1E+15.);

  AliESDtrackCuts* esdTrackCutsHG0 = jetCuts1006->Clone("JetCuts10001006");
  esdTrackCutsHG0->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);


  // throw out tracks with too low number of clusters in
  // the first pass (be consistent with TPC only tracks)
  // N.B. the number off crossed rows still acts on the tracks after
  // all iterations if we require tpc standalone, number of clusters
  // and chi2 TPC cuts act on track after the first iteration
  //   esdTrackCutsH0->SetRequireTPCStandAlone(kTRUE);
  //   esdTrackCutsH0->SetMinNClustersTPC(80); // <--- first pass
  
  
  // the complement to the one with SPD requirement
  //  AliESDtrackCuts* esdTrackCutsHG1 = CreateTrackCutsPWGJE(10011006);
  AliESDtrackCuts* esdTrackCutsHG1 = jetCuts1006->Clone("JetCuts10011006");
  esdTrackCutsHG1->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);

  // the tracks that must not be taken pass this cut and
  // non HGC1 and HG
  //  AliESDtrackCuts* esdTrackCutsHG2 = CreateTrackCutsPWGJE(10021006);
  AliESDtrackCuts* esdTrackCutsHG2 = jetCuts1006->Clone("JetCuts10021006");
  esdTrackCutsHG2->SetMaxChi2PerClusterITS(1E10);


  // standard cuts also used in R_AA analysis
  //   "Global track RAA analysis QM2011 + Chi2ITS<36";
  //  AliESDtrackCuts* esdTrackCutsH2 = CreateTrackCutsPWGJE(1000);
  AliESDtrackCuts* esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  esdTrackCutsH2->SetMinNCrossedRowsTPC(120);
  esdTrackCutsH2->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCutsH2->SetMaxChi2PerClusterITS(36);
  esdTrackCutsH2->SetMaxFractionSharedTPCClusters(0.4);
  esdTrackCutsH2->SetMaxChi2TPCConstrainedGlobal(36);

  esdTrackCutsH2->SetEtaRange(-0.9,0.9);
  esdTrackCutsH2->SetPtRange(0.15, 1e10);


  //  AliESDtrackCuts* esdTrackCutsGCOnly = CreateTrackCutsPWGJE(10041006);
  AliESDtrackCuts* esdTrackCutsGCOnly = jetCuts1006->Clone("JetCuts10041006");
  esdTrackCutsGCOnly->SetRequireITSRefit(kFALSE);



  // TPC only tracks
  AliESDtrackCuts* esdTrackCutsTPCCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  esdTrackCutsTPCCOnly->SetMinNClustersTPC(70);
  
  // Compose the filter
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  // 1, 1<<0
  trackFilter->AddCuts(esdTrackCutsL);
  // 2 1<<1
  trackFilter->AddCuts(esdTrackCutsITSsa);
  // 4 1<<2
  trackFilter->AddCuts(itsStrong);
  itsStrong->SetFilterMask(1);        // AND with Standard track cuts 
  // 8 1<<3
  trackFilter->AddCuts(electronID);
  electronID->SetFilterMask(4);       // AND with Pixel Cuts
  // 16 1<<4
  trackFilter->AddCuts(esdTrackCutsHG0);
  // 32 1<<5
  trackFilter->AddCuts(esdTrackCutsHG1);
  // 64 1<<6
  trackFilter->AddCuts(esdTrackCutsHG2);
  // 128 1<<7
  trackFilter->AddCuts(esdTrackCutsTPCCOnly); // add QM TPC only track cuts
  esdfilter->SetTPCOnlyFilterMask(128);
  // 256 1<<8
  trackFilter->AddCuts(esdTrackCutsGCOnly);
  // 512 1<<9                         
  AliESDtrackCuts* esdTrackCutsHG1_tmp = new AliESDtrackCuts(*esdTrackCutsHG1); // avoid double delete
  trackFilter->AddCuts(esdTrackCutsHG1_tmp); // add once more for tpc only tracks
  // 1024 1<<10                        
  trackFilter->AddCuts(esdTrackCutsH2); // add r_aa cuts

  
  
  esdfilter->SetGlobalConstrainedFilterMask(1<<8|1<<9); // these tracks are written out as global constrained tracks
  esdfilter->SetHybridFilterMaskGlobalConstrainedGlobal((1<<4)); // these normal global tracks will be marked as hybrid
  esdfilter->SetWriteHybridGlobalConstrainedOnly(kTRUE); // write only the complement
  //     esdfilter->SetTPCConstrainedFilterMask(1<<11); // these tracks are written out as tpc constrained tracks

  esdfilter->SetTrackFilter(trackFilter);
      
   

   // Filter with cuts on V0s
   AliESDv0Cuts*   esdV0Cuts = new AliESDv0Cuts("Standard V0 Cuts pp", "ESD V0 Cuts");
   esdV0Cuts->SetMinRadius(0.2);
   esdV0Cuts->SetMaxRadius(200);
   esdV0Cuts->SetMinDcaPosToVertex(0.05);
   esdV0Cuts->SetMinDcaNegToVertex(0.05);
   esdV0Cuts->SetMaxDcaV0Daughters(1.5);
   esdV0Cuts->SetMinCosinePointingAngle(0.99);
   AliAnalysisFilter* v0Filter = new AliAnalysisFilter("v0Filter");
   v0Filter->AddCuts(esdV0Cuts);

   esdfilter->SetV0Filter(v0Filter);

 
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());
   mgr->ConnectInput  (kinefilter,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (kinefilter,  0, mgr->GetCommonOutputContainer());
   AliAnalysisDataContainer *coutputEx = mgr->CreateContainer("cFilterList", TList::Class(),
								   AliAnalysisManager::kOutputContainer,"pyxsec_hists.root");
   mgr->ConnectOutput (kinefilter,  1,coutputEx);

    
    //-----------------------
    // Run the analysis
    //-----------------------    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
    
    cout <<" Analysis ended sucessfully "<< endl ;
    
  }
  else cout << "Chain was not produced ! "<<endl;
  
}
