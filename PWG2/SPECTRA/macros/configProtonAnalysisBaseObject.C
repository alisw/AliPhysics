AliProtonAnalysisBase *GetProtonAnalysisBaseObject(const char* analysisLevel = "ESD",
						   const char* esdAnalysisType = "Hybrid",
						   const char* pidMode = "Bayesian") {
						   
  //Function to setup the AliProtonAnalysisBase object and return it
  AliProtonAnalysisBase *baseAnalysis = new AliProtonAnalysisBase();
  //baseAnalysis->SetDebugMode();
  baseAnalysis->SetAnalysisLevel(analysisLevel);
  if(analysisLevel == "ESD") {  
    baseAnalysis->SetTriggerMode(AliProtonAnalysisBase::kMB2);
    baseAnalysis->SetMinTPCClusters(110);
    baseAnalysis->SetMaxChi2PerTPCCluster(2.2);
    baseAnalysis->SetMaxCov11(0.5);
    baseAnalysis->SetMaxCov22(0.5);
    baseAnalysis->SetMaxCov33(0.5);
    baseAnalysis->SetMaxCov44(0.5);
    baseAnalysis->SetMaxCov55(0.5);
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
      baseAnalysis->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
      baseAnalysis->SetTPCpid();
      baseAnalysis->SetMaxSigmaToVertex(2.0);
      /*baseAnalysis->SetMaxDCAXY(1.5);
	baseAnalysis->SetMaxDCAZ(1.5);*/
      baseAnalysis->SetPointOnITSLayer6();
      baseAnalysis->SetPointOnITSLayer5();
      //baseAnalysis->SetPointOnITSLayer4();
      //baseAnalysis->SetPointOnITSLayer3();
      baseAnalysis->SetPointOnITSLayer2();
      baseAnalysis->SetPointOnITSLayer1();
      baseAnalysis->SetMinITSClusters(4);
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
    baseAnalysis->SetAcceptedVertexDiamond(5.,5.,15.);
    baseAnalysis->SetEtaMode();
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
      break;
    case "Sigma1":
      baseAnalysis->SetPIDMode(AliProtonAnalysisBase::kSigma1);
      baseAnalysis->SetNSigma(3);
      baseAnalysis->SetdEdxBandInfo("$ALICE_ROOT/PWG2/data/protonsdEdxInfo.dat");
      break;
    case "Sigma2":
      baseAnalysis->SetPIDMode(AliProtonAnalysisBase::kSigma2);
      baseAnalysis->SetNSigma(3);
      baseAnalysis->SetdEdxBandInfo("$ALICE_ROOT/PWG2/data/protonsdEdxInfo.dat");
      break;
    default:
      break;
    }//PID mode
  }//ESD
  if(analysisLevel == "MC") 
    baseAnalysis->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);

  return baseAnalysis;
}