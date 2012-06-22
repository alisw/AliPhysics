void runAODProof(Int_t c=4, const char * proofMode = "full")
{ //1 data AOD049
  //2 MC AOD048
  //3 data AOD086
  //4 MC AOD090
  
  //  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  
  // Load common libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so"); 
  gSystem->Load("libGui.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libProof.so");
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libSTEER.so");
  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gSystem->Load("libPWGLFspectra.so");
  
  AliAnalysisAlien * handler = new AliAnalysisAlien("test");
  handler->SetOverwriteMode();
  handler->SetRunMode(proofMode);
  handler->SetProofReset(0);
  handler->SetAliROOTVersion("v5-03-33-AN");
  
  //handler->SetNproofWorkers(80);
  //handler->SetNproofWorkersPerSlave(4);
  handler->SetProofCluster(Form("%s@alice-caf.cern.ch", gSystem->Getenv("CAFUSER")));
  //handler->SetProofCluster(Form("%s@skaf.saske.sk",gSystem->Getenv("CAFUSER")));
  // Set handler for Real DATA:
  if (c == 1){
    handler->SetProofDataSet("/default/lmilano/LHC10h_000138653_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138662_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138666_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138730_AOD049_p2#aodTree|/default/lmilano/LHC10h_000138732_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139507_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139465_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139437_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139107_AOD049_p2#aodTree|/default/lmilano/LHC10h_000139510_AOD049_p2#aodTree");
  }
  if (c == 2){
    handler->SetProofDataSet("/default/lmilano/LHC11a10a_138653_AOD048#aodTree|/default/lmilano/LHC11a10a_138662_AOD048#aodTree|/default/lmilano/LHC11a10a_138666_AOD048#aodTree|/default/lmilano/LHC11a10a_138730_AOD048#aodTree|/default/lmilano/LHC11a10a_138732_AOD048#aodTree|/default/lmilano/LHC11a10a_139507_AOD048#aodTree|/default/lmilano/LHC11a10a_139465_AOD048#aodTree|/default/lmilano/LHC11a10a_139437_AOD048#aodTree|/default/lmilano/LHC11a10a_139107_AOD048#aodTree|/default/lmilano/LHC11a10a_139510_AOD048#aodTree");      
  }
  if (c == 3){
    handler->SetProofDataSet("/default/lmilano/LHC10h_000138653_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138666_AOD086_p2#aodTree|/default/lmilano/LHC10h_000139107_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138275_AOD086_p2#aodTree|/default/lmilano/LHC10h_000139465_AOD086_p2#aodTree|/default/lmilano/LHC10h_000139437_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138442_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138396_AOD086_p2#aodTree|/default/lmilano/LHC10h_000138364_AOD086_p2#aodTree");      
  }
  if (c == 4){
    handler->SetProofDataSet("/default/lmilano/LHC11a10a_bis_138653_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138666_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_139107_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138275_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_139465_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_139437_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138442_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138396_AOD090#aodTree|/default/lmilano/LHC11a10a_bis_138364_AOD090#aodTree");        }

  handler->SetAliRootMode("default");
  handler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF -I$ALICE_ROOT/PWGLF");
  handler->SetAdditionalLibs("libPWGLFspectra.so");
  // gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
  // gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
  // gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
  // gROOT->LoadMacro("AliSpectraAODPID.cxx+g");
  // gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");
  
  // handler->SetAdditionalLibs("AliSpectraAODHistoManager.cxx AliSpectraAODHistoManager.h AliSpectraAODPID.cxx AliSpectraAODPID.h AliSpectraAODTrackCuts.cxx AliSpectraAODTrackCuts.h AliSpectraAODEventCuts.cxx AliSpectraAODEventCuts.h AliAnalysisTaskSpectraAOD.cxx AliAnalysisTaskSpectraAOD.h");
  // handler->SetAnalysisSource("Histograms.h HistogramNames.h AliSpectraAODHistoManager.cxx+ AliSpectraAODTrackCuts.cxx+ AliSpectraAODEventCuts.cxx+ AliSpectraAODPID.cxx+ AliAnalysisTaskSpectraAOD.cxx+");
  handler->SetClearPackages();
  
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  mgr->SetGridHandler(handler);
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  // Add PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  Bool_t isMC = kFALSE;
  if (c == 2 || c == 4) isMC = kTRUE;   
  Printf("-------------------------------adding in runAOD AddTaskPIDResponse");
  // Add PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse();
  // Printf("OADB PATH:::::%s",taskPID->GetOADBPath());
  // taskPID->SetOADBPath("alien:///alice/cern.ch/user/a/akalweit/ForLeornado/OADB");
  // Printf("OADB PATH:::::%s",taskPID->GetOADBPath());
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
  AliVZEROEPSelectionTask *selTask = AddTaskVZEROEPSelection();
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AddTaskSpectraAOD.C");
  //LOOP OVER SELECTION
  //                            0    1    2    3    4    5
  Bool_t mc=kFALSE;
  Double_t CentCutMin[4]= {     0,  20,  20,  20};
  Double_t CentCutMax[4]= {   100,  50,  50,  50};
  Double_t QvecCutMin[4]=    {  0,   0,   0, 1.5};
  Double_t QvecCutMax[4]=   { 100, 100, 0.4, 100};
  Double_t EtaMin[4]={       -0.8,-0.8,-0.8,-0.8};
  Double_t EtaMax[4]={        0.8, 0.8, 0.8, 0.8};
  Double_t Nsigmapid=3.;
  Double_t pt=5.;
  Double_t p=5.;
  Double_t y=.5;
  Double_t ptTofMatch=.6;
  UInt_t trkbit=1024;
  UInt_t trkbitQVector=1;
  Bool_t UseCentPatchAOD049=kFALSE;
  Double_t DCA=100000;
  UInt_t minNclsTPC=70;
  
  if(c==2||c==4)mc=kTRUE;
  if(c==1)UseCentPatchAOD049=kTRUE;
  for(Int_t icut=0;icut<4;icut++){
    //if(icut!=0)continue;
    AliAnalysisTaskSpectraAOD *taskAOD =AddTaskSpectraAOD(mc,CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC);
    
    taskAOD->GetOutputSlot(1)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(1)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(2)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(2)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(3)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(3)->GetContainer()->GetName(),taskAOD->GetName()));
    taskAOD->GetOutputSlot(4)->GetContainer()->SetName(Form("%s_%s",taskAOD->GetOutputSlot(4)->GetContainer()->GetName(),taskAOD->GetName()));
  
  }
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof");
}
