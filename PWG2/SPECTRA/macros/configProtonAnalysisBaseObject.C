AliProtonAnalysisBase *GetProtonAnalysisBaseObject(const char* analysisLevel = "ESD",
						   Bool_t kAnalyzeMC = kTRUE,
						   const char* esdAnalysisType = "Hybrid",
						   const char* pidMode = "Bayesian",
						   Bool_t kUseOnlineTrigger = kFALSE,
						   Bool_t kUseOfflineTrigger = kFALSE,
						   Bool_t kRunQA = kFALSE) {
  //Function to setup the AliProtonAnalysisBase object and return it
  AliProtonAnalysisBase *baseAnalysis = new AliProtonAnalysisBase();
  //baseAnalysis->SetDebugMode();
  if(kRunQA) baseAnalysis->SetRunQA();
  baseAnalysis->SetAnalysisLevel(analysisLevel);
  if(analysisLevel == "ESD") {
    if(kAnalyzeMC)
      baseAnalysis->SetTriggerMode(AliProtonAnalysisBase::kMB2);
    //use the offline trigger
    if(kUseOnlineTrigger) baseAnalysis->UseOnlineTrigger();

    //use the offline trigger
    if(kUseOfflineTrigger) baseAnalysis->OfflineTriggerInit();

    baseAnalysis->SetMinTPCClusters(80);
    baseAnalysis->SetMaxChi2PerTPCCluster(3.5);
    /*baseAnalysis->SetMaxCov11(2.0);
    baseAnalysis->SetMaxCov22(2.0);
    baseAnalysis->SetMaxCov33(0.5);
    baseAnalysis->SetMaxCov44(0.5);
    baseAnalysis->SetMaxCov55(2.0);*/
    baseAnalysis->SetMinTPCdEdxPoints(80);
    switch(esdAnalysisType) {
    case "TPC":
      baseAnalysis->SetAnalysisMode(AliProtonAnalysisBase::kTPC);
      baseAnalysis->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
      baseAnalysis->SetTPCpid();
      baseAnalysis->SetMaxSigmaToVertexTPC(2.0);
      //baseAnalysis->SetMaxDCAXYTPC(1.5);
      //baseAnalysis->SetMaxDCAZTPC(1.5);
      break;
    case "Hybrid":
      baseAnalysis->SetAnalysisMode(AliProtonAnalysisBase::kHybrid);
      baseAnalysis->SetPhaseSpace(9, -0.9, 0.9, 6, 0.45, 1.05);
      //baseAnalysis->SetPhaseSpace(18, -0.9, 0.9, 32, 0.5, 1.3);
      baseAnalysis->SetTPCpid();
      //baseAnalysis->SetMaxSigmaToVertex(3.0);
      //baseAnalysis->SetMaxDCAXY(0.5);
      //baseAnalysis->SetMaxDCAZ(0.7);
      baseAnalysis->SetMaxDCA3D(2.0);
      //baseAnalysis->SetPointOnITSLayer6();
      //baseAnalysis->SetPointOnITSLayer5();
      //baseAnalysis->SetPointOnITSLayer4();
      //baseAnalysis->SetPointOnITSLayer3();
      //baseAnalysis->SetPointOnITSLayer2();
      //baseAnalysis->SetPointOnITSLayer1();
      baseAnalysis->SetPointOnSPDLayers();
      baseAnalysis->SetMinITSClusters(2);
      break;
    case "FullHybrid":
      baseAnalysis->SetAnalysisMode(AliProtonAnalysisBase::kFullHybrid);
      baseAnalysis->SetPhaseSpace(9, -0.9, 0.9, 6, 0.45, 1.05);
      //baseAnalysis->SetPhaseSpace(18, -0.9, 0.9, 32, 0.5, 1.3);
      baseAnalysis->SetTPCpid();
      //baseAnalysis->SetMaxSigmaToVertex(3.0);
      //baseAnalysis->SetMaxDCAXY(0.2);
      //baseAnalysis->SetMaxDCAZ(0.7);
      //baseAnalysis->SetMaxDCA3D(0.2);
      baseAnalysis->SetPtDependentDCAxy(5,2.89575e+02,6.62161e+01,1.99085e+00);
      //baseAnalysis->SetPointOnITSLayer6();
      //baseAnalysis->SetPointOnITSLayer5();
      //baseAnalysis->SetPointOnITSLayer4();
      //baseAnalysis->SetPointOnITSLayer3();
      //baseAnalysis->SetPointOnITSLayer2();
      //baseAnalysis->SetPointOnITSLayer1();
      baseAnalysis->SetPointOnSPDLayers();
      baseAnalysis->SetMinITSClusters(2);
      break;
    case "Global":
      baseAnalysis->SetAnalysisMode(AliProtonAnalysisBase::kGlobal);
      baseAnalysis->SetPhaseSpace(20, -1.0, 1.0, 48, 0.3, 1.5);
      baseAnalysis->SetMaxSigmaToVertex(2.0);
      //baseAnalysis->SetMaxDCAXY(2.0);
      //baseAnalysis->SetMaxDCAZ(2.0);
      baseAnalysis->SetTPCRefit();
      baseAnalysis->SetPointOnITSLayer1();
      baseAnalysis->SetPointOnITSLayer2();
      //baseAnalysis->SetPointOnITSLayer3();
      //baseAnalysis->SetPointOnITSLayer4();
      baseAnalysis->SetPointOnITSLayer5();
      baseAnalysis->SetPointOnITSLayer6();
      baseAnalysis->SetMinITSClusters(5);
      baseAnalysis->SetITSRefit();
      baseAnalysis->SetESDpid();
      baseAnalysis->SetTOFpid();
      break;
    default:
      break;
    }
    baseAnalysis->SetAcceptedVertexDiamond(1.,1.,10.);
    baseAnalysis->SetMinNumOfContributors(0);
    //baseAnalysis->SetEtaMode();
    switch(pidMode) {
    case "Bayesian":
      baseAnalysis->SetPIDMode(AliProtonAnalysisBase::kBayesian);
      //Momentum dependent priors
      /*TFile *f = TFile::Open("$ALICE_ROOT/PWG2/data/PriorProbabilities.root ");
	TF1 *fitElectrons = (TF1 *)f->Get("fitElectrons");
	TF1 *fitMuons = (TF1 *)f->Get("fitMuons");
	TF1 *fitPions = (TF1 *)f->Get("fitPions");
	TF1 *fitKaons = (TF1 *)f->Get("fitKaons");
	TF1 *fitProtons = (TF1 *)f->Get("fitProtons");
	baseAnalysis->SetPriorProbabilityFunctions(fitElectrons,
	fitMuons,
	fitPions,
	fitKaons,
	fitProtons);*/
      //Fixed prior probabilities
      Double_t partFrac[5] = {0.01, 0.01, 0.85, 0.10, 0.05};
      if(!baseAnalysis->IsPriorProbabilityFunctionUsed())
	baseAnalysis->SetPriorProbabilities(partFrac);
      break;
    case "Ratio":
      baseAnalysis->SetPIDMode(AliProtonAnalysisBase::kRatio);
      baseAnalysis->SetRatio(-0.2);
      break;
    case "Sigma":
      baseAnalysis->SetPIDMode(AliProtonAnalysisBase::kSigma1);
      baseAnalysis->SetNSigma(4);
      break;
    default:
      break;
    }//PID mode
  }//ESD
  if(analysisLevel == "MC") 
    baseAnalysis->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
  if(analysisLevel == "AOD")
    baseAnalysis->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);

  return baseAnalysis;
}
