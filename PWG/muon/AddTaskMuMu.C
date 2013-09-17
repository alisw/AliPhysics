///
/// Configure a task to get invariant mass spectrum of dimuons
///
/// author: L. Aphecetche (Subatech) (laurent.aphecetche - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMu(const char* outputname, 
							 TList* triggerClassesToConsider, 
							 const char* beamYear, 
							 TArrayF* centralities,
							 Bool_t simulations)                                    
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMuMu", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  cout << "inputDataType=" << inputDataType.Data() << endl;
  
  // Configure analysis
  //===========================================================================  
  
  AliAnalysisTaskMuMu* task;
  
  if (simulations && triggerClassesToConsider )
  {
    triggerClassesToConsider->Add(new TObjString("CMULLO-B-NOPF-MUON"));
	triggerClassesToConsider->Add(new TObjString("CMSNGL-B-NOPF-MUON"));
	triggerClassesToConsider->Add(new TObjString("ANY"));
  }
  
  if ( triggerClassesToConsider ) 
  {
    task = new AliAnalysisTaskMuMu((inputDataType=="ESD"),triggerClassesToConsider,beamYear,centralities);    
  }
  else 
  {
    task = new AliAnalysisTaskMuMu((inputDataType=="ESD"),beamYear,centralities);    
  }
  
  task->AddEventCut("ALL",AliAnalysisTaskMuMu::kEventAll);

if (!simulations) 
  {
  	task->AddEventCut("PSALL",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS);  
//  	task->AddEventCut("OFFLINE1",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventOFFLINEMUL1);  
//  	task->AddEventCut("REJECTED",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventREJECTED);  
   }
  
//  task->AddEventCut("PSALLNOTZEROPILEUP",AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventNOTZEROPILEUP);  

//  task->AddEventCut("PSALLMUL1", AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventOFFLINEMUL1);

//  task->AddEventCut("PSALLMUL2", AliAnalysisTaskMuMu::kEventAll | AliAnalysisTaskMuMu::kEventPS | AliAnalysisTaskMuMu::kEventOFFLINEMUL2);

//   task->AddSingleCut("MATCHLOWRABSDCA",
//                      AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kDCA);
// 
   task->AddSingleCut("MATCHLOWRABS",
                      AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs);

   task->AddSingleCut("MATCHLOWRABSETA",
                      AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs);

  task->AddSingleCut("MATCHLOWRABSETADCA",
                     AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kEta|AliAnalysisTaskMuMu::kDCA);


//  task->AddPairCut("ALL",AliAnalysisTaskMuMu::kAll);

   task->AddPairCut("MATCHLOWRABSBOTH",
                    AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs,
                    AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs);

   task->AddPairCut("MATCHLOWRABSETABOTH",
                    AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kEta,
                    AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kEta);

   task->AddPairCut("MATCHLOWRABSETADCABOTH",
                    AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kEta|AliAnalysisTaskMuMu::kDCA,
                    AliAnalysisTaskMuMu::kAll|AliAnalysisTaskMuMu::kMatchedLow|AliAnalysisTaskMuMu::kRabs|AliAnalysisTaskMuMu::kEta|AliAnalysisTaskMuMu::kDCA);

// igor binning

task->AddBin("psi","pt", 0.0, 0.5,"IGOR");
task->AddBin("psi","pt", 0.5, 1.0,"IGOR");
task->AddBin("psi","pt", 1.0, 1.5,"IGOR");
task->AddBin("psi","pt", 1.5, 2.0,"IGOR");
task->AddBin("psi","pt", 2.0, 2.5,"IGOR");
task->AddBin("psi","pt", 2.5, 3.0,"IGOR");
task->AddBin("psi","pt", 3.0, 3.5,"IGOR");
task->AddBin("psi","pt", 3.5, 4.0,"IGOR");
task->AddBin("psi","pt", 4.0, 4.5,"IGOR");
task->AddBin("psi","pt", 4.5, 5.0,"IGOR");
task->AddBin("psi","pt", 5.0, 6.0,"IGOR");
task->AddBin("psi","pt", 6.0, 7.0,"IGOR");
task->AddBin("psi","pt", 7.0, 8.0,"IGOR");
task->AddBin("psi","pt", 8.0, 9.0,"IGOR");
task->AddBin("psi","pt", 9.0,11.0,"IGOR");
task->AddBin("psi","pt",11.0,13.0,"IGOR");
task->AddBin("psi","pt",13.0,15.0,"IGOR");

// roberta binning
task->AddBin("psi","pt",0.0,1.0,"ROBERTA");
task->AddBin("psi","pt",1.0,2.0,"ROBERTA");
task->AddBin("psi","pt",2.0,3.0,"ROBERTA");
task->AddBin("psi","pt",3.0,4.0,"ROBERTA");
task->AddBin("psi","pt",4.0,5.0,"ROBERTA");
task->AddBin("psi","pt",5.0,6.0,"ROBERTA");
task->AddBin("psi","pt",6.0,8.0,"ROBERTA");

/*

  for ( Int_t i = 0; i < 10; ++i ) 
  {
    Double_t xmin = i*1.0;
    
  	task->AddBin("psi","pt",xmin,xmin+1.0);
  }

  for ( Int_t i = 0; i < 5; ++i ) 
  {
    Double_t xmin = 10+i*2.0;
    
  	task->AddBin("psi","pt",xmin,xmin+2.0);
  }

  for ( Int_t i = 0; i < 6; ++i ) 
  {
    Double_t xmin = -4+i*0.25;
    
  	task->AddBin("psi","y",xmin,xmin+0.25);
  }

  Int_t nphiSteps = 18;
  
  for ( Int_t i = 0; i < nphiSteps; ++i ) 
  {
    Double_t xstep = 2*TMath::Pi()/nphiSteps;
    Double_t xmin = -TMath::Pi() + i*xstep;
    
  	task->AddBin("psi","phi",xmin,xmin+xstep);
  }
  
  */
  
  Double_t ymin(-4.0);
  Double_t ymax(-3.43);

  Double_t ymin(-4.0);
  Double_t ymax(-3.43);

  task->AddBin("psi","y",-4,-2.5,"ILAB");
  
  task->AddBin("psi","y",ymin,ymax,"ICOMMON");

  for ( Int_t i = 0; i < 6; ++i ) 
  {
    Double_t y = -4+i*0.25;
    
    task->AddBin("psi","y",y,y+0.25,"6PACK");
  }

/*  
  for ( Int_t i = 0; i < 20; ++i ) 
  {
    Double_t xmin = i*0.25;
    
  	task->AddBin("psi","y vs pt",xmin,xmin+0.25,ymin,ymax);
  }

  for ( Int_t i = 0; i < 4; ++i ) 
  {
    Double_t xmin = 5+i*1.0;
    
  	task->AddBin("psi","y vs pt",xmin,xmin+1.0,ymin,ymax);
  }

*/
/*
  task->CreateMesh("psi","pt","y");
*/

  task->DisableHistograms("^V02D");
  task->DisableHistograms("^dca");
  task->DisableHistograms("^Chi12");
  task->DisableHistograms("^Rabs12");
  
  mgr->AddTask(task);  
  
  static int n(0);
  
  ++n;
  
  if ( n > 1 ) containerName += Form("%d",n);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();  

  AliAnalysisDataContainer *coutputHC = 
  mgr->CreateContainer("HC",AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputname);

  AliAnalysisDataContainer *coutputCC = 
  mgr->CreateContainer("CC",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer,outputname);
  
  AliAnalysisDataContainer* cparam = 
  mgr->CreateContainer("BIN", AliAnalysisMuMuBinning::Class(),AliAnalysisManager::kParamContainer,outputname);
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputHC);
  mgr->ConnectOutput(task, 2, coutputCC);
  mgr->ConnectOutput(task, 3, cparam);
  
  return task;
}

