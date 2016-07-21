AliAnalysisTaskSE *AddTaskLambdacTMVA(TString finname,Int_t storeNtuple,Bool_t readMC,Bool_t MCPid,Bool_t realPid,Bool_t resPid,Bool_t keepLcNoQuark,Bool_t isHijing,
Int_t syst=0, Int_t bit=0, TString postname="",Int_t useNtrkWeight = 0,Int_t storeNtupleReco = 0)
{
  //==============================================================================                                                      
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLambdac", "No analysis manager to connect to.");
    return NULL;
  }


  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( finname.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
      filecuts=TFile::Open(finname.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	AliFatal("Input file not found : check your cut object");
      }
  }
  AliRDHFCutsLctopKpi *analysiscuts = new AliRDHFCutsLctopKpi();
  // syst = 0 : pp, syst = 1: PbPb, syst = 2 : pPb
  if(stdcuts) {
   if(syst==0) analysiscuts->SetStandardCutsPP2010();
   if(syst==1) analysiscuts->SetStandardCutsPbPb2011();
   if(syst==2) analysiscuts->SetStandardCutsPPb2013();
  }
  else analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiAnalysisCuts");
  analysiscuts->SetName("LctopKpiAnalysisCuts");
  analysiscuts->SetMinPtCandidate(-1.);
  analysiscuts->SetMaxPtCandidate(10000.);

  // Analysis task                                                                                                                     
  AliAnalysisTaskSELambdacTMVA *lambdacTask = new AliAnalysisTaskSELambdacTMVA("LambdacAnalysis",storeNtuple,storeNtupleReco,analysiscuts);
  //if(storeNtuple<0 || storeNtuple>2) {AliFatal("Invalid storeNtuple argument - check value");}
  lambdacTask->SetReadMC(readMC);
	lambdacTask->SetKeepLcNotFromQuark(keepLcNoQuark);
	lambdacTask->SetCollisionSystem(syst);
  if(MCPid) lambdacTask->SetMCPid();
  if(resPid) lambdacTask->SetResonantPid();
  if(realPid) lambdacTask->SetRealPid();
	if(isHijing) lambdacTask->SetIsHijing();
  lambdacTask->SetAnalysis(kTRUE);

  //bit:0 nocut, 1:LcCut, 2:LcPID, 3: Both
  lambdacTask->SetUseFilterBitCut(bit==1||bit==3?1:0);
  lambdacTask->SetUseFilterBitPID(bit>1?1:0);

  lambdacTask->SetDebugLevel(0);
  mgr->AddTask(lambdacTask);

	if(useNtrkWeight>0){
		TH1F *hNchPrimaries;
		if(useNtrkWeight==1) hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithD");
		else if(useNtrkWeight==2) hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithCand");
		else if(useNtrkWeight==3) hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvSel");
		else if(useNtrkWeight==4) hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrPSSel");
		else {
			AliFatal("useNtrkWeight value not a valid option - choice from 1-4");
			return 0x0;
		}
		if(hNchPrimaries){ 
			lambdacTask->SetUseNchWeight(kTRUE);
			lambdacTask->SetMCNchHisto(hNchPrimaries);
		}
		else {
			AliFatal("Histogram for multiplicity weights not found");
			return 0x0;
		}
	}


  //
  // Create containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassLambdac";

  TString finDirname="pp";
  TString inname = "cinputLc";
  TString outname = "coutputLc";
  TString cutsname = "coutputLcCuts";
  TString normname = "coutputLcNorm";
  TString normnament = "coutputLcNormNt";
  TString ntuplename = "fNtupleLambdac";
  TString ntuplenamereco = "fNtupleLambdacReco";
  TString nev2 = "coutputNev";
  TString outname2 = "coutputLambdacMC";
  TString aPrioriname = "coutputAPriori";
  TString multiplicityname = "coutputMultiplicity";
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  normnament += finDirname.Data();
  nev2 += finDirname.Data();
  outname2 += finDirname.Data();
  aPrioriname += finDirname.Data();
  multiplicityname += finDirname.Data();

  inname +=  postname.Data();
  outname +=  postname.Data();
  cutsname +=  postname.Data();
  normname += postname.Data();
  normnament += postname.Data();
  ntuplename +=  postname.Data();
  ntuplenamereco +=  postname.Data();
  nev2 +=  postname.Data();
  outname2 +=  postname.Data();


  //input container
  AliAnalysisDataContainer *cinputLambdac = mgr->CreateContainer(inname,TChain::Class(),
								 AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(lambdacTask,0,mgr->GetCommonInputContainer());

  
  AliAnalysisDataContainer *coutputLambdacCuts = mgr->CreateContainer(cutsname,TList::Class(),
								      AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,2,coutputLambdacCuts);

  AliAnalysisDataContainer *coutputLambdac = mgr->CreateContainer(outname,TList::Class(),
								  AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,1,coutputLambdac);

  AliAnalysisDataContainer *coutputLambdacNev = mgr->CreateContainer(nev2,TH1F::Class(),
								     AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,3,coutputLambdacNev);

 
  AliAnalysisDataContainer *coutputLambdacNorm = mgr->CreateContainer(normname,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  mgr->ConnectOutput(lambdacTask,4,coutputLambdacNorm);

	TString foutname = Form("InvMassLambdac_nt1.root",postname.Data());
  if (storeNtuple) {
    AliAnalysisDataContainer *coutputLambdac2 = mgr->CreateContainer(ntuplename,TNtuple::Class(),
								     AliAnalysisManager::kOutputContainer,foutname);
    coutputLambdac2->SetSpecialOutput();
    mgr->ConnectOutput(lambdacTask,5,coutputLambdac2);    

    coutputLambdacNorm = mgr->CreateContainer(normnament,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,foutname);
    mgr->ConnectOutput(lambdacTask,4,coutputLambdacNorm);
	}
	if(storeNtupleReco>0){
		AliAnalysisDataContainer *coutputLambdac3 = mgr->CreateContainer(ntuplenamereco,TNtuple::Class(),
				AliAnalysisManager::kOutputContainer,foutname);
		coutputLambdac3->SetSpecialOutput();
		mgr->ConnectOutput(lambdacTask,6,coutputLambdac3);    
	}


  return lambdacTask;
}
