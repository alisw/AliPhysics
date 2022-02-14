 AliAnalysisTaskSELambdac* AddTaskLambdac(TString finname,Bool_t storeNtuple,Bool_t readMC,Bool_t MCPid,Bool_t realPid,Bool_t resPid,Bool_t useKF,
				  Bool_t fillVarHists=kFALSE, Bool_t priorsHists=kFALSE, Bool_t multiplicityHists=kFALSE, Int_t syst=0, Int_t bit=0,  TString postname="")
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
	Printf("FATAL: Input file not found : check your cut object");
      }
  }
  AliRDHFCutsLctopKpi* prodcuts=new AliRDHFCutsLctopKpi();
  // syst = 0 : pp, syst = 1: PbPb, syst = 2 : pPb
  if(stdcuts) {
    if(syst==0) prodcuts->SetStandardCutsPP2010();
    if(syst==1) prodcuts->SetStandardCutsPbPb2011();
    if(syst==2) prodcuts->SetStandardCutsPPb2013();
  }
  else   prodcuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiProdCuts");
  prodcuts->SetName("LctopKpiProdCuts");
  prodcuts->SetMinPtCandidate(-1.);
  prodcuts->SetMaxPtCandidate(10000.);

  AliRDHFCutsLctopKpi *analysiscuts = new AliRDHFCutsLctopKpi();
  if(stdcuts) {
   if(syst==0) analysiscuts->SetStandardCutsPP2010();
   if(syst==1) analysiscuts->SetStandardCutsPbPb2011();
   if(syst==2) analysiscuts->SetStandardCutsPPb2013();
  }
  else analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get("LctopKpiAnalysisCuts");
  analysiscuts->SetName("LctopKpiAnalysisCuts");
  analysiscuts->SetMinPtCandidate(-1.);
  analysiscuts->SetMaxPtCandidate(10000.);

  // Aanalysis task                                                                                                                     
  AliAnalysisTaskSELambdac *lambdacTask = new AliAnalysisTaskSELambdac("LambdacAnalysis",storeNtuple,analysiscuts,prodcuts);
  lambdacTask->SetReadMC(readMC);
  if(MCPid) lambdacTask->SetMCPid();
  if(resPid) lambdacTask->SetResonantPid();
  if(realPid) lambdacTask->SetRealPid();
  lambdacTask->SetFillVarHists(fillVarHists);
  lambdacTask->SetPriorsHists(priorsHists);
  lambdacTask->SetMultiplicityHists(multiplicityHists);
  lambdacTask->SetAnalysis(kTRUE);

  //bit:0 nocut, 1:LcCut, 2:LcPID, 3: Both
  lambdacTask->SetUseFilterBitCut(bit==1||bit==3?1:0);
  lambdacTask->SetUseFilterBitPID(bit>1?1:0);

  lambdacTask->SetDebugLevel(0);
  if(useKF) {
    lambdacTask->SetUseKF();
    Float_t cuts[10]={0.1,0.1,1.5,0.5,0.1,1.5,0.5,0.1,1.5,0.5};
    lambdacTask->SetCutsKF(cuts);
  }
  mgr->AddTask(lambdacTask);

  //
  // Create containers for input/output
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_InvMassLambdac";

  TString finDirname="pp";
  TString inname = "cinputLc";
  TString outname = "coutputLc";
  TString cutsname = "coutputLcCuts";
  TString normname = "coutputLcNorm";
  TString ntuplename = "coutputLc2";
  TString nev2 = "coutputNev";
  TString outname2 = "coutputLambdacMC";
  TString aPrioriname = "coutputAPriori";
  TString multiplicityname = "coutputMultiplicity";
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  ntuplename += finDirname.Data();
  nev2 += finDirname.Data();
  outname2 += finDirname.Data();
  aPrioriname += finDirname.Data();
  multiplicityname += finDirname.Data();

  inname +=  postname.Data();
  outname +=  postname.Data();
  cutsname +=  postname.Data();
  normname += postname.Data();
  ntuplename +=  postname.Data();
  nev2 +=  postname.Data();
  outname2 +=  postname.Data();
  aPrioriname +=  postname.Data();
  multiplicityname +=  postname.Data();


  AliAnalysisDataContainer *cinputLambdac = mgr->CreateContainer(inname,TChain::Class(),
								 AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(lambdacTask,0,mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutputLambdacCuts = mgr->CreateContainer(cutsname,TList::Class(),
								      AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,2,coutputLambdacCuts);

  AliAnalysisDataContainer *coutputLambdac = mgr->CreateContainer(outname,TList::Class(),
								  AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,1,coutputLambdac);

  AliAnalysisDataContainer *coutputLambdacMC = mgr->CreateContainer(outname2,TList::Class(),
								    AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,3,coutputLambdacMC);

  AliAnalysisDataContainer *coutputLambdacNev = mgr->CreateContainer(nev2,TH1F::Class(),
								     AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectOutput(lambdacTask,4,coutputLambdacNev);

  if(priorsHists){
    AliAnalysisDataContainer *coutputAPriori = mgr->CreateContainer(aPrioriname,TList::Class(),
								  AliAnalysisManager::kOutputContainer,outputfile.Data());
    mgr->ConnectOutput(lambdacTask,5,coutputAPriori);
  }
  if(multiplicityHists){
    AliAnalysisDataContainer *coutputMultiplicity = mgr->CreateContainer(multiplicityname,TList::Class(),
								       AliAnalysisManager::kOutputContainer,outputfile.Data());
    mgr->ConnectOutput(lambdacTask,6,coutputMultiplicity);
  }
  AliAnalysisDataContainer *coutputLambdacNorm = mgr->CreateContainer(normname,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

 mgr->ConnectOutput(lambdacTask,7,coutputLambdacNorm);

  if (storeNtuple) {
    AliAnalysisDataContainer *coutputLambdac2 = mgr->CreateContainer(ntuplename,TNtuple::Class(),
								     AliAnalysisManager::kOutputContainer,"InvMassLambdac_nt1.root");
    coutputLambdac2->SetSpecialOutput();
    mgr->ConnectOutput(lambdacTask,7,coutputLambdac2);
  }


  return lambdacTask;
}
