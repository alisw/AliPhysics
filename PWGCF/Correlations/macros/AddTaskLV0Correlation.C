AliLeadingV0Correlation* AddTaskLV0Correlation(TString  fListName                = "LV0Correlation",
											   TString  fCollisiontype			 = "PP",
											   Bool_t	fAnalysisMC              = 1,
											   Int_t    fCase                    = 2,
											   Bool_t   fRemoveAutoCorr          = 0,
											   Double_t fPVzCut                  = 10,
											   Int_t    fFilterBit               = 256,
											   Double_t fRapidityCut             = 0.75,
											   Int_t	fmaxEventsinPool         = 2000,
											   Int_t	fminTracksinPool         = 100,
											   Int_t    fMinEventsToMix          = 5,
											   Double_t fV0radius			     = 0.3,
											   Double_t fV0PostoPVz			     = 0.060,
											   Double_t fV0NegtoPVz			     = 0.060,
											   Double_t fDCAV0Daughters		     = 1.00,
											   Double_t fCPAK0				     = 0.960,
											   Double_t fCPALam				     = 0.997,
											   Double_t fRejectLamK0		     = 0.005,
											   Double_t fRejectK0Lam		     = 0.010,
											   Double_t fSigmaPID                = 3.0,
											   Double_t fCutCTK0				 = 20.0,
											   Double_t fCutCTLa				 = 30.0,
											   Double_t fMassCutK0               = 0.0105,             
											   Double_t fMassCutLa               = 0.0105,
											   Double_t fTriglow			     = 6.0,
											   Double_t fTrighigh                = 12.0,
											   Double_t fTPCClusters			 = 70,
											   Double_t fTPCfindratio            = 0.8) 
{
	// Get the current analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {Error("AddTaskLV0Correlation.C", "No Analysis Manager ");return 0;}
	
	//PVz Binning for pool PP or PbPb	
	Double_t pvzbinlimits[] = {-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12};
	Int_t pvzbinnumb = sizeof(pvzbinlimits)/sizeof(Double_t) - 1;
	
	//Mult Binning for pool	pp 
	Double_t cent_mult_binlimitsPP[] = {  0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
		100,110,120,130,140,150,160,170,180,190,
		200,210,220,230,240,250,260,270,280,290,
		300,500,1000,2000};
	
	Int_t cent_mult_bin_numbPP = sizeof(cent_mult_binlimitsPP)/sizeof(Double_t) - 1;
	
	//Cent Binning for pool	PbPb 
	Double_t cent_mult_binlimitsPbPb[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,
		14,15,16,17,18,19,20,21,22,23,24,25,26,27,
		28,29,30,31,32,33,34,35,36,37,38,39,40,42,
		44,46,48,50,52,54,56,58,60,65,70,75,80,90};
	
	Int_t cent_mult_bin_numbPbPb = sizeof(cent_mult_binlimitsPbPb)/sizeof(Double_t) - 1;
	//Correlation task
	AliLeadingV0Correlation *myTask = new AliLeadingV0Correlation(fListName.Data());
	
	myTask->SetCollidingSystem(fCollisiontype);
	myTask->SetMCAnalysis(fAnalysisMC);
	myTask->SetCase(fCase);
	myTask->SetRemoveAutoCorr(fRemoveAutoCorr);
	myTask->SetMaxNEventsInPool(fmaxEventsinPool);
	myTask->SetMinNTracksInPool(fminTracksinPool);
	myTask->SetMinEventsToMix(fMinEventsToMix);
	myTask->SetPoolPVzBinLimits(pvzbinnumb,pvzbinlimits);
	if(fCollisiontype=="PP")myTask->SetPoolCentBinLimits(cent_mult_bin_numbPP,cent_mult_binlimitsPP);
	if(fCollisiontype=="PbPb")myTask->SetPoolCentBinLimits(cent_mult_bin_numbPbPb,cent_mult_binlimitsPbPb);
	myTask->SetPrimeryVertexCut(fPVzCut);
	myTask->SetFilterBit(fFilterBit);
	myTask->SetCutRap(fRapidityCut);
	myTask->SetV0Radius(fV0radius);
	myTask->SetV0PostoPVz(fV0PostoPVz);
	myTask->SetV0NegtoPVz(fV0NegtoPVz);
	myTask->SetDCAV0Daughters(fDCAV0Daughters);
	myTask->SetCPAK0(fCPAK0);
	myTask->SetCPALam(fCPALam);
	myTask->SetRejectLamK0(fRejectLamK0);
	myTask->SetRejectK0Lam(fRejectK0Lam);
	myTask->SetSigmaPID(fSigmaPID);
	myTask->SetCTK0(fCutCTK0);
	myTask->SetCTLa(fCutCTLa);
	myTask->SetMassCutK0(fMassCutK0);
	myTask->SetMassCutLa(fMassCutLa);
	myTask->SetTrigLow(fTriglow);
	myTask->SetTrigHigh(fTrighigh);
	myTask->SetTPCClusters(fTPCClusters);
	myTask->SetTPCFinables(fTPCfindratio);
	
	mgr->AddTask(myTask);
	
	// Create containers for input/output
	TString outputFileName = AliAnalysisManager::GetCommonFileName();
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(fListName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
	
	//____________________________________________//
	mgr->ConnectInput(myTask,0,cinput);	
	mgr->ConnectOutput(myTask,1,coutput);
	
return myTask;
}



