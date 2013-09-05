AliAnalysisTaskSEDplusCorrelations *AddTaskDplusCorrWT(Bool_t system=kFALSE,				                              					                                     Bool_t readMC=kFALSE,
					     Bool_t ifTrig=kTRUE,	     
					     Bool_t mixing=kFALSE,
					     Int_t select=1,
					     Bool_t reco=kTRUE,		         					     Int_t usedisp=0,
					     TString finDirname="",
					     TString filename="",
					     TString finAnObjname="AnalysisCuts")
{

    //                                                                                                                                    
  // Test macro for the AliAnalysisTaskSE for D+ candidates 

  //Invariant mass histogram and                                                 
  // association with MC truth (using MC info in AOD)                                                                                   
  
  // Get the pointer to the existing analysis manager via the static access method.                                                     
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDplusCorr", "No analysis manager to connect to.");
  }

  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
      filecuts=TFile::Open(filename.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	AliFatal("Input file not found : check your cut object");
      }
  }
  
  TFile* filecuts1=new TFile("AssocPartCuts.root");
	  if(!filecuts1->IsOpen()){
		  cout<<"Input file1 not found: exit"<<endl;
		  return;
  }

  
  //Analysis Task

  
  AliRDHFCutsDplustoKpipi* analysiscuts=new AliRDHFCutsDplustoKpipi();
  if(stdcuts) {
    cout<<analysiscuts->GetNPtBins()<<endl;
    if(system==0) analysiscuts->SetStandardCutsPP2010();
    else if(system==1){
      analysiscuts->SetStandardCutsPbPb2011();
      analysiscuts->SetMinCentrality(minC);
      analysiscuts->SetMaxCentrality(maxC);
      //      analysiscuts->SetUseAOD049(kTRUE);
      analysiscuts->SetUseCentrality(AliRDHFCuts::kCentV0M);
    }
  }
  else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);

 
  
  AliHFAssociatedTrackCuts* corrCuts=new AliHFAssociatedTrackCuts();
  corrCuts = (AliHFAssociatedTrackCuts*)filecuts1->Get("AssociatedCuts");
  corrCuts->SetName("AssociatedCuts");
  corrCuts->PrintAll();
  
    
  AliAnalysisTaskSEDplusCorrelations *dplusTask = new AliAnalysisTaskSEDplusCorrelations("DplusAnalysis",analysiscuts,corrCuts);


  dplusTask->SetReadMC(readMC);
  dplusTask->SetEventMix(mixing);
  dplusTask->SetCorrelator(select);
  dplusTask->SetUseDisplacement(usedisp);
  dplusTask->SetTrigEfficiency(ifTrig);
  dplusTask->SetDebugLevel(0);
  dplusTask->SetMassLimits(0.2);
  dplusTask->SetUseBit(kTRUE);
  dplusTask->SetSystem(kFALSE);
dplusTask->SetUseReconstruction(kTRUE);   

  //  if (system==0) dplusTask->SetDoImpactParameterHistos(kTRUE);

  mgr->AddTask(dplusTask);
  
  // Create containers for input/output 

  TString inname = "cinputDplus";
  TString outname = "coutputDplus";
  TString cutsname = "coutputDplusCuts";
  TString normname = "coutputDplusNorm";
  //TString assocutsname = "coutputAssoCuts";
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  //assocutsname += finDirname.Data();
  TString centr=Form("%.0f%.0f",analysiscuts->GetMinCentrality(),analysiscuts->GetMaxCentrality());
  inname += centr;
  outname += centr;
  cutsname += centr;
  normname += centr;
  // assocutsname += centr;


  AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer(inname,TChain::Class(),
							       AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGHF_CorrHF_Dplus";
  
  AliAnalysisDataContainer *coutputDplusCuts = mgr->CreateContainer(cutsname,TList::Class(),
								    AliAnalysisManager::kOutputContainer,
								    outputfile.Data());
  
  AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer(outname,TList::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  AliAnalysisDataContainer *coutputDplusNorm = mgr->CreateContainer(normname,AliNormalizationCounter::Class(),
								AliAnalysisManager::kOutputContainer,
								outputfile.Data());
  

  //  AliAnalysisDataContainer *coutputAssoCuts = mgr->CreateContainer(assocutsname,TList::Class(),
  //								   AliAnalysisManager::kOutputContainer,
  //								   outputfile.Data());
    

  mgr->ConnectInput(dplusTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(dplusTask,1,coutputDplus);
  
  mgr->ConnectOutput(dplusTask,2,coutputDplusCuts);

  mgr->ConnectOutput(dplusTask,3,coutputDplusNorm);  

  // mgr->ConnectOutput(dplusTask,4,coutputAssoCuts);
  
  return dplusTask;
}
