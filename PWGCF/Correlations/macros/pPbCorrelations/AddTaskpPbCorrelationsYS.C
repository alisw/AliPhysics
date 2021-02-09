AliAnalysisTaskSEpPbCorrelationsYS* AddTaskpPbCorrelationsYS(
								       TString  fListName      ="pPbCorrelations_1",
								       TString  fListName1     ="Corr_1",
								       TString  fListName2     ="QA_1",
								       TString  fCollisiontype ="HMPPV0",//MBPP,HMPP, pPb,PbPb
								       Bool_t  fDataType       =kTRUE,//TRUE=real data, FALSE=MC
								       Bool_t frun2            =kTRUE,
								       Bool_t fFMDcut          =kTRUE,
								       TString anamode         ="TPCFMD",//TPCTPC, TPCV0A, TPCV0C, V0AV0C,TPCFMD, TPCFMDC, FMDFMD, SE__CA
								       TString anacent         ="Manual",
								       TString assomode        ="hadron",
								       Int_t ffilterbit        =32,
								       Int_t fFMDcutpar        =12,
								       Bool_t fmakehole        =kFALSE,
								       Bool_t fptdiff          =kFALSE,
								       Float_t fmaxpt          =3.0,
								       Int_t fMinNTracksInPool =5000,
								       Int_t fMinNEventsInPool =5,
								       Bool_t fefficalib       =kFALSE,
								       Float_t fminpt=0.2,
								       Float_t fcuthighmult=0.1
								       //Bool_t fFillcorrelation=kTR
								       )
{
  // Get the current analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {Error("AddTaskpPbCorrelationsYS.C", "No Analysis Manager");return 0;}

  //PVz Binning for pool PP or PbPb
  //Double_t pvzbinlimits[] = {-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12};
  Double_t pvzbinlimits[] = {-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12};
  Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;

  //Mult Binning for pool pbpb
  Double_t cent_mult_binlimitsPbPb[] = {  0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
					100,110,120,130,140,150,160,170,180,190,
					200,210,220,230,240,250,260,270,280,290,
					  300,500,1000,1500,2000,3000};
  Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;

  //Cent Binning for pool pp

  /*Double_t cent_mult_binlimitsPP[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,
					 14,15,16,17,18,19,20,21,22,23,24,25,26,27,
					 28,29,30,31,32,33,34,35,36,37,38,39,40,42,
					 44,46,48,50,52,54,56,58,60,65,70,75,80,90};
  */
  Double_t cent_mult_binlimitsPP[] = {0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200};
  Int_t cent_mult_bin_numbPP = sizeof(cent_mult_binlimitsPP)/sizeof(Double_t) - 1;

  //  Double_t cent_mult_binlimitspPbManual[] = {0,10,20,30,40,50,60,70,80,90,100,150,200};
  Double_t cent_mult_binlimitspPbManual[] = {0,15,20,30,40,50,60,70,80,90,100,150,200};
  Int_t cent_mult_bin_numbpPbManual = sizeof(cent_mult_binlimitspPbManual)/sizeof(Double_t) - 1;
  
  //Cent Binning for pool	pPb
  Double_t cent_mult_binlimitspPb[] = {0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100};
  Int_t cent_mult_bin_numbpPb = sizeof(cent_mult_binlimitspPb)/sizeof(Double_t) - 1;

  //  Double_t cent_mult_binlimitsHMPP[] = { 0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  //  Double_t cent_mult_binlimitsHMPP[] = {0,0.001,0.0033,0.01,0.02,0.033,0.05,0.1,0.2,0.5,1,2,5,10,15,20,30,40,50,70,80,90,100};
  Double_t cent_mult_binlimitsHMPP[] = {0.,0.1,0.2,0.5,1,2,5,10,15,20,30,40,50,70,80,90,100};
  Int_t cent_mult_bin_numbHMPP = sizeof(cent_mult_binlimitsHMPP)/sizeof(Double_t) - 1;
  
  //Correlation task
  AliAnalysisTaskSEpPbCorrelationsYS *myTask = new AliAnalysisTaskSEpPbCorrelationsYS(fListName.Data());

  myTask->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);
  myTask->SetFilterBit(ffilterbit);
  myTask->SetAnalysisMode(anamode);
  myTask->SetAssociatedTrack(assomode);
  //myTask->SetPID(fpid);
  myTask->SetDatatype(fDataType);
  myTask->SetRunType(frun2);
  myTask->SetFMDcut(fFMDcut);
  myTask->SetFMDcutpar(fFMDcutpar);
  myTask->Setacceptancehole(fmakehole);
  myTask->SetPtdiff(fptdiff);
  myTask->SetPtMax(fmaxpt);
  myTask->SetPtMin(fminpt);
  myTask-> SetHighmultcut(fcuthighmult);
			  
  //myTask->SetMinNTracksInPool(5000);
  myTask->SetMinNTracksInPool(fMinNTracksInPool);
  myTask->SetMinEventsToMix(fMinNEventsInPool);
			    
  myTask->SetAnalysisCent(anacent);//0:V0A 1:ZNA 2:
  myTask->SetAnalysisCollisionType(fCollisiontype);
  
  //  myTask->SetFillCorrelation(kFALSE);
  myTask->SetEfficiencyCorrection(fefficalib);

  //  if(fCollisiontype=="PP")myTask->SetPoolCentBinLimits(cent_mult_bin_numbPP,cent_mult_binlimitsPP);
  //  if(fCollisiontype=="PbPb"){myTask->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);}
    if(anacent=="Manual"){
      //      if(fCollisiontype.Contains("PP"))myTask->SetPoolCentBinLimits(cent_mult_bin_numbPP,cent_mult_binlimitsPP);
      if(fCollisiontype.Contains("HMPP")||fCollisiontype.Contains("MBPP"))myTask->SetPoolCentBinLimits(cent_mult_bin_numbpPbManual,cent_mult_binlimitspPbManual);
      else if(fCollisiontype.Contains("pPb"))myTask->SetPoolCentBinLimits(cent_mult_bin_numbpPbManual,cent_mult_binlimitspPbManual);
      else if(fCollisiontype.Contains("PbPb"))myTask->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);
    }else{
      if(fCollisiontype=="pPb" ||fCollisiontype=="PbPb"){myTask->SetPoolCentBinLimits(cent_mult_bin_numbpPb,cent_mult_binlimitspPb);}    
      else if(fCollisiontype.Contains("HMPP")|| fCollisiontype=="PP"||fCollisiontype=="MBPP") myTask->SetPoolCentBinLimits(cent_mult_bin_numbHMPP,cent_mult_binlimitsHMPP);
    }
    
  mgr->AddTask(myTask);

  //cout<<"hogehoge"<<endl;
  //  gSystem->Exec("alien_cp alien:///alice/cern.ch/user/y/ysekiguc/correction.root ./");
  //cout<<"hogehoge"<<endl;


  // Create containers for input/output
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  ///  TString output1name="Corr";
  //TString output2name="QA";

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(fListName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(fListName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(fListName2.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  //____________________________________________//
  mgr->ConnectInput(myTask,0,cinput);
  
  
  mgr->ConnectOutput(myTask,1,coutput);
  mgr->ConnectOutput(myTask,2,coutput2);
  mgr->ConnectOutput(myTask,3,coutput3);

  
  //  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  //  mgr->ConnectInput(myTask,1,valid);
  
  return myTask;
}
