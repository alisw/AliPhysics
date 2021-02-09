AliAnalysisTask *AddTaskHFEminiEventCreator(
					    Bool_t IsMCevent = kTRUE,
					    Double_t TPCchi2 = 4.,
					    Int_t MinTPCNcluster = 100,
					    Int_t MinTPCclusterPID = 80,
					    Double_t TPCclusterRatio = 0.6,
					    Int_t MinNclusterITS = 3,
					    Bool_t checkITSLayerstatus = kFALSE,
					    Double_t eta = 0.8,
					    Double_t ptMin = 0.5,
					    Double_t ptMax = 100.,
					    Double_t Vz = 10.,
					    Double_t dcaxy = 1.,
					    Double_t dcaz = 2.,
					    Double_t prodVz = 0.5,
					    Double_t spdResolution = 0.25,
					    Double_t nsigmaTPClow = -1.,
					    Double_t nsigmaTPChigh = 3.,
					    Double_t nsigmaTOF = 3.,
					    TString collisionSystem = "pp",
					    TString taskName = "Test" ){
  
  printf("Adding mini event creator\n");
  
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hfe_HFEAOD", "No analysis manager found.");
    return 0;
  }
  
  if(dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler())) isMC = kTRUE;
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


  //Set the cuts ----
  
  Bool_t isInt7 = kTRUE;
  Bool_t isRemoveFirstEvent = kTRUE;

  TString commontaskName = "MiniTree";
  commontaskName        += taskName;
  
  AliHFEminiEventCreator *miniEventCreator = new AliHFEminiEventCreator(commontaskName);

  miniEventCreator->SetIsMCEvent(IsMCevent);
  miniEventCreator->SetChi2TPCCut( TPCchi2 );
  miniEventCreator->SetMinClusterTPC( MinTPCNcluster );
  miniEventCreator->SetMinClusterTPCPID( MinTPCclusterPID );
  miniEventCreator->SetClusterRatioTPC( TPCclusterRatio );
  miniEventCreator->SetClusterITS( MinNclusterITS );
  miniEventCreator->SetITSLayerStatus( checkITSLayerstatus );
  miniEventCreator->SetKinematicsCuts( ptMin, ptMax, eta );
  miniEventCreator->SetVertexZCut( Vz );
  miniEventCreator->SetDCAcut( dcaxy, dcaz );
  miniEventCreator->SetProductionVzToSPDVz( prodVz );
  miniEventCreator->SetSPDResolution( spdResolution );
  miniEventCreator->SetTPCnSigma( nsigmaTPClow, nsigmaTPChigh );
  miniEventCreator->SetTOFnSigma( nsigmaTOF );  
  miniEventCreator->SetCollisionSystem(collisionSystem);
  if(isInt7) miniEventCreator->SelectCollisionCandidates(AliVEvent::kINT7);
  else miniEventCreator->SelectCollisionCandidates();
  
  if(isRemoveFirstEvent) miniEventCreator->SetRemoveFirstEventFromChunk();
 
  AliHFEpidTPC *tpcpid = miniEventCreator->GetTPCResponse();
  
  mgr->AddTask(miniEventCreator);
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(commontaskName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
  mgr->ConnectInput(miniEventCreator, 0, cinput);
  mgr->ConnectOutput(miniEventCreator, 1, coutput);
  
  return miniEventCreator;
}
