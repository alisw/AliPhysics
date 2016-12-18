AliAnalysisTask *AddTask_taku_lowmassTree(Bool_t getFromAlien=kFALSE,
					  TString cFileName = "Config_taku_lowmassT.C", 
					  Char_t* outputFileName="LMEE.root",
					  ULong64_t triggerMask = AliVEvent::kINT7,
					  TString PhotonCuts = "00000062007300008500200000",
					  TString PeriodName = "LHC13c"
					  )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_taku_lowmass", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train

  //TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  TString configBasePath= "$PWD/";
  TString cFileName = "Config_taku_lowmassT.C";

  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/t/tgunji/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
    std::cout << " Config from Alien:  " << std::endl;
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);
  


  std::cout << "Configpath:  " << configFilePath << std::endl;
  
  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);  
  std::cout << "hasMC :  " << hasMC << std::endl;
  

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();


  //=========  Add V0ReaderV1 task ================================             
  //=========  Add My main task ================================             
  gROOT->LoadMacro(configFilePath.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
    AliV0ReaderV1 *fV0ReaderV1 = AddTask_V0Reader(PhotonCuts, PeriodName);
    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);
  }    

  AliAnalysisTaskTGReducedTree *taskpt = AddTask_LowMassTree(triggerMask, PhotonCuts);
  if (!hasMC) taskpt->UsePhysicsSelection(kTRUE);
  // Add task(s)
  mgr->AddTask(taskpt);


  AliAnalysisDataContainer *coutputpt = mgr->CreateContainer("electron", 
							     TList::Class(), 
							     //TTree::Class(), 
                                                             AliAnalysisManager::kOutputContainer, outputFileName);

  /*
  AliAnalysisDataContainer *cReducedEvent =
    mgr->CreateContainer("ReducedEventDQ",
			 AliDielectronTGReducedInfo::Class(),
			 AliAnalysisManager::kExchangeContainer,
			 "reducedEvent");
  */

  // Connect input/output                                                                                                                                 
  mgr->ConnectInput(taskpt,  0, mgr->GetCommonInputContainer());
  mgr->ConnectInput(taskpt, 0, cinput);
  //mgr->ConnectOutput(taskpt, 1, cReducedEvent);
  mgr->ConnectOutput(taskpt, 1, coutputpt);



  return taskpt;

}
