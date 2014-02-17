AliAnalysisTask *AddTaskNuclexFilter( Bool_t onlyMuon = kTRUE,
				      Bool_t keepAllEvents = kTRUE,
				      Int_t mcMode = 0,
				      Int_t nsigmaTrk1 = 3,
				      Int_t partType1 = 2,
				      Int_t nsigmaTrk2 = 5, 
				      Int_t partType2 = 7
				      )
{

  //get the current analysis manager

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskNuclexFilter", "No analysis manager found.");
    return 0;
  }
  
  //check for output aod handler
  if (!mgr->GetOutputEventHandler()||mgr->GetOutputEventHandler()->IsA()!=AliAODHandler::Class()) {
    Warning("AddTaskNuclExFilter","No AOD output handler available. Not adding the task!");
    return 0;
  }

  //Do we have an MC handler?
  //  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0)||hasMC_aod;
  /*  
  //Do we run on AOD?
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  //Allow merging of the filtered aods on grid trains
  if(mgr->GetGridHandler()) {
    printf(" SET MERGE FILTERED AODs \n");
    mgr->GetGridHandler()->SetMergeAOD(kTRUE);
  }
  
  if(isAOD) {
    //add options to AliAODHandler to duplicate input event
    AliAODHandler *aodHandler = (AliAODHandler*)mgr->GetOutputEventHandler();
    aodHandler->SetCreateNonStandardAOD();
    }
 */

  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  if ( !isAOD) {
    Printf("ERROR! This task can only run on AODs!");
  }

  
 
  //========= Add task to the ANALYSIS manager =====
  
  cout<<"\n\n=============== addtask Nuclex =================== "<<endl<<endl;

  //AliAnalysisTaskSE *esdnuclexfilter = new AliAnalysisTaskESDNuclExFilter("NuclEx Filter",onlyMuon,keepAllEvents,mcMode,nsigmaTrk1,partType1,nsigmaTrk2,partType2);

  AliAnalysisTaskESDNuclExFilter  *esdnuclexfilter = new AliAnalysisTaskESDNuclExFilter("NuclEx Filter",onlyMuon,keepAllEvents,mcMode,nsigmaTrk1,partType1,nsigmaTrk2,partType2);
  
  esdnuclexfilter->SetWriteMuonAOD(kTRUE);

  //================================================
  //              data containers
  //================================================
  //            find input container

  mgr->ConnectInput  (esdnuclexfilter,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (esdnuclexfilter,  0, mgr->GetCommonOutputContainer());


  return esdnuclexfilter;
}
