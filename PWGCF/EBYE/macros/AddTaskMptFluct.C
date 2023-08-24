
/**************************************************************************
 Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *Macro designed for AliAnalysisMeanPtdata class 
  Author: (Tulika Tripathy, Sadhana Dash),   IIT Bombay                                                                      *                                    *
 * Contributors are mentioned in the code where appropriate.              *
 **************************************************************************/

AliAnalysisMeanPtdata * AddTaskMptFluct
(const char * outfilename="meaanpTfluck",const char * listname="PbPb5TeV",
 /*achar *inputPileupCutFileName   = (char*)"alien:///alice/cern.ch/user/t/tutripat/tulika_profile_final.root",*/
 // int    usePileupCut_PbPb5TeV    =  1 ,  // 0: no                   1: yes  
 int    singlesOnly             =  1,
 int    nContributor             =  1,
 int    fb             =  768,
 int    nClusterMin             =  70,
 int    nTPCCrossRows             = 70,
 double dcaZMax                 =  1.0,///tulika 3.2,
 double dcaXYMax                =  1.0,///tulika 2.4,
 double chi2perTPC              = 36,
 double chi2perITS              = 36,
 double VzMin                 =  -10.0,
 double VzMax                 =  10.0,//vz<10 Tulika
 int pileUpEvent              =0) //0:without pileup cut 1:with pileup cut

 {
   double dcaZMin                = -dcaZMax;
   double dcaXYMin               = -dcaXYMax;
   Bool_t trigger                = kFALSE; //will come to this later
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    //    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    ::Error("AddTaskMptFluct", "No analysis manager to connect to.");
    return NULL;
  }  
  

  //  TString listName;
  TString taskName;
  //  taskName=

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if (inputDataType != "AOD") {
    Printf("ERROR! This task can only run on AODs!");
  }

  
  
  //  char taskName[15];
  //  sprintf(taskName,"ptflucttask");
  

  //  int iTask = 0;
  //  AliAnalysisDataContainer *cinput0;
  //  AliAnalysisDataContainer *coutput1;
  AliAnalysisMeanPtdata *task;
  AliAnalysisDataContainer *taskInputContainer;
  AliAnalysisDataContainer *taskOutputContainer;
  //  TFile  * inputFilePileup  = 0;
   
  // if ( usePileupCut_PbPb5TeV )
  //   {
     /*     TGrid::Connect("alien:");
     inputFilePileup = TFile::Open(inputPileupCutFileName,"OLD");
     if (inputFilePileup)
       {
	 cout << "\n\n\n\n ====================================================" << endl;
	 cout << " Requested file:" << inputPileupCutFileName << " was opened." << endl;
	 cout << "\n\n\n\n ====================================================" << endl;
       }
     else if (!inputFilePileup)
       {
	 cout << "\n\n\n\n ====================================================" << endl;
	 cout << " Requested file:" << inputPileupCutFileName << " was not opened. ABORT." << endl;
	 cout << "\n\n\n\n ====================================================" << endl;
	 return 0;*/
       


     //     TProfile *hProfPileupCut = (TProfile*)inputFilePileup->Get("hProf");
     //     TProfile *hProfPileupCut = (TProfile*)inputFilePileup->Get("_profV0MvsTPCout");

     //     task = new AliAnalysisMeanPtdata("meanptdataPbPb5.02");
  //  task = new AliAnalysisMeanPtdata("meanptdataPbPb5.02",/*usePileupCut_PbPb5TeV,*/singlesOnly,nContributor,nClusterMin,dcaZMax,dcaXYMax,chi2perTPC,chi2perITS);
  task = new AliAnalysisMeanPtdata("PbPb5TeV");
  //     task->SetUsePileupCut_PbPb5TeV( usePileupCut_PbPb5TeV);
     //     task->SetPileupCut_PbPb5TeV( hProfPileupCut);
     task->SetNClusterMin(         nClusterMin     );
     task->SetNTPCCrossRows(         nTPCCrossRows     );
     task->SetFilterBit(         fb     );
     task->SetSinglesOnly(         singlesOnly     );
     task->SetNContributors(       nContributor     );
     task->SetDcaZMin(             dcaZMin         );
     task->SetDcaZMax(             dcaZMax         );
     task->SetDcaXYMin(            dcaXYMin        );
     task->SetDcaXYMax(            dcaXYMax        );
     task->SetChi2PerTPCCluster(          chi2perTPC      ); 
     task->SetChi2PerITSCluster(          chi2perITS      ); 
     task->SetVzMin(             VzMin        );
     task->SetVzMax(             VzMax        );
     task->SetPileUpEvent(    pileUpEvent     );


     //     TString outputfile = AliAnalysisManager::GetCommonFileName();
     //     outputfile += ":PWG2FEMTO";    
     
     //     cinput0 = mgr->GetCommonInputContainer();
   
     //     coutput1 = mgr->CreateContainer("containerName", TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:Meanpt", AliAnalysisManager::GetCommonFileName()));

     //     coutput1  = mgr->CreateContainer(containerName,  TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);


  // AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
   
  //  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outfilename, TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:Meanpt", AliAnalysisManager::GetCommonFileName()));


     taskOutputContainer = mgr->CreateContainer(listname, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", AliAnalysisManager::GetCommonFileName(),outfilename));



     mgr->AddTask(task);
     mgr->ConnectInput (task,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(task,1,taskOutputContainer);

     //  mgr->ConnectInput (task, 0, cinput0);
     //  mgr->ConnectOutput(task,1,coutput1);

     //  iTask++;
 
  return task;
}   


