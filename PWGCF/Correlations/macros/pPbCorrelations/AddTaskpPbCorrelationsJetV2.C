AliAnalysisTaskSEpPbCorrelationsJetV2* AddTaskpPbCorrelationsJetV2(
								       TString  fListName      ="pPbCorrelations_1",
								       TString  fListName1     ="Corr_1",
								       TString  fListName2     ="QA_1",
								       TString  fCollisiontype ="pPb",
								       Bool_t fDataType        =kTRUE,//TRUE=real data, FALSE=MC
								       Bool_t frun2            =kTRUE,
								       Bool_t fFMDcut          =kTRUE,
								       TString anamode         ="FMDAFMDC",//TPCTPC, TPCFMDA, TPCFMDC, TPCTPCFMDA, TPCTPCFMDC, FMDAFMDC, Rcp
								       TString anacent         ="V0A",
								       TString assomode        ="hadron",
								       Int_t ffilterbit        =32,
								       Int_t fFMDcutpar        =7,
								       Bool_t fprimTPC         =kFALSE,
                                                                       Bool_t fprimFMD         =kFALSE,
                                                                       Bool_t fcentcalib       =kFALSE,
                                                                       Double_t fReduceDphi    =-1., // 1.5707, 0.9, -1
                                                                       Bool_t fSymmetricFMD    =kFALSE,
                                                                       Bool_t IsLikeSign       =kTRUE,
								       Float_t fminpt          =0.5,
								       Float_t fmaxpt          =20.0,
								       Float_t fAsscoptCut     =0.5,
                                                                       Int_t fMinNTracksInPool =5000,
								       Int_t fMinNEventsInPool =5, 
								       Double_t dCenMin = 0.,
								       Double_t dCenMax = 10.
                                                                      )
{
  // Get the current analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { Error("AddTaskpPbCorrelationsJetV2.C", "No Analysis Manager"); return 0x0;}

  //PVz Binning for pool PP or PbPb
  Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
  Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;

  //Mult Binning for pool pbpb
  Double_t cent_mult_binlimitsPbPb[] = {  0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
					100,110,120,130,140,150,160,170,180,190,
					200,210,220,230,240,250,260,270,280,290,
					300,500,1000,2000};
  Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;

  Double_t cent_mult_binlimitsPP[] = {0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200};
  Int_t cent_mult_bin_numbPP = sizeof(cent_mult_binlimitsPP)/sizeof(Double_t) - 1;

  //Cent Binning for pool pPb
  Double_t cent_mult_binlimitspPb[] = {0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100};
  Int_t cent_mult_bin_numbpPb = sizeof(cent_mult_binlimitspPb)/sizeof(Double_t) - 1;

  //Correlation task
  AliAnalysisTaskSEpPbCorrelationsJetV2 *myTask = new AliAnalysisTaskSEpPbCorrelationsJetV2(fListName.Data());

  myTask->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);
  myTask->SetFilterBit(ffilterbit);
  myTask->SetAnalysisMode(anamode);
  myTask->SetAssociatedTrack(assomode);
  myTask->SetDatatype(fDataType);
  myTask->SetCentCalib(fcentcalib);
  myTask->SetRunType(frun2);
  myTask->SetFMDcut(fFMDcut);
  myTask->SetFMDcutpar(fFMDcutpar);
  myTask->SetReduceDphi(fReduceDphi);
  myTask->SetSymmetricFMD(fSymmetricFMD);
  myTask->SetLikeSign(IsLikeSign);
  myTask->SetPtMin(fminpt);
  myTask->SetPtMax(fmaxpt);
  myTask->SetAssoCut(fAsscoptCut);
  myTask->SetCentrality(dCenMin,dCenMax);
//  myTask->SetTPCTPCList(TPCTPC_Fit);

  myTask->SetMinNTracksInPool(fMinNTracksInPool);
  myTask->SetMinEventsToMix(fMinNEventsInPool);
			    
  myTask->SetAnalysisCent(anacent);
  myTask->SetAnalysisCollisionType(fCollisiontype);
  myTask->SetmcprimFMD(fprimFMD);
  myTask->SetmcprimTPC(fprimTPC);


  myTask->SetPoolCentBinLimits(cent_mult_bin_numbpPb,cent_mult_binlimitspPb);    
  mgr->AddTask(myTask);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(fListName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(fListName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(fListName2.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  //____________________________________________//
  mgr->ConnectInput(myTask,0,cinput);
  mgr->ConnectOutput(myTask,1,coutput);
  mgr->ConnectOutput(myTask,2,coutput2);
  mgr->ConnectOutput(myTask,3,coutput3);
  
  

  return myTask;
}
