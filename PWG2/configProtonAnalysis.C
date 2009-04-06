//__________________________________________________//
AliProtonAnalysis *GetProtonAnalysisObject(const char* analysisLevel = "ESD", 
					   const char* esdAnalysisType = "Hybrid", 
					   const char* pidMode = "Bayesian") {
  //Function to setup the AliProtonAnalysis object and return it
  AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(analysisLevel,esdAnalysisType,pidMode);

  AliProtonAnalysis *analysis = new AliProtonAnalysis();
  analysis->SetBaseAnalysis(baseAnalysis);
  //if(analysisBase->GetEtaMode()) analysis->SetEtaMode();
  analysis->InitAnalysisHistograms(baseAnalysis->GetNBinsX(),
				   baseAnalysis->GetMinX(),
				   baseAnalysis->GetMaxX(),
				   baseAnalysis->GetNBinsY(),
				   baseAnalysis->GetMinY(),
				   baseAnalysis->GetMaxY());
    
  return analysis;
}

//__________________________________________________//
AliProtonQAAnalysis *GetProtonQAAnalysisObject(const char* analysisLevel = "ESD",
					       const char* esdAnalysisType = "Hybrid",
					       const char* pidMode = "Bayesian") {
  //Function to setup the AliProtonQAAnalysis object and return it
  AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(analysisLevel,esdAnalysisType,pidMode);

  AliProtonQAAnalysis *analysis = new AliProtonQAAnalysis();
  analysis->SetQAOn();
  analysis->SetRunMCAnalysis();
  //analysis->SetMCProcessId(4);//4: weak decay - 13: hadronic interaction
  //analysis->SetMotherParticlePDGCode(3122);//3122: Lambda
  analysis->SetQAYPtBins(baseAnalysis->GetNBinsX(),
			 baseAnalysis->GetMinX(),
			 baseAnalysis->GetMaxX(),
			 baseAnalysis->GetNBinsY(),
			 baseAnalysis->GetMinY(),
			 baseAnalysis->GetMaxY());

  return analysis;
}

//__________________________________________________//
AliProtonAnalysisBase *GetProtonAnalysisBaseObject(const char* analysisLevel = "ESD",
						   const char* esdAnalysisType = "Hybrid",
						   const char* pidMode = "Bayesian") {
  //Function to setup the AliProtonAnalysisBase object and return it
  AliProtonAnalysisBase *baseAnalysis = new AliProtonAnalysisBase();
  //baseAnalysis->SetDebugMode();
  baseAnalysis->SetAnalysisLevel(analysisLevel);
  if(analysisLevel == "ESD") {  
    baseAnalysis->SetTriggerMode(AliProtonAnalysisBase::kMB2);
    switch(esdAnalysisType) {
    case "TPC":
      baseAnalysis->SetAnalysisMode(AliProtonAnalysisBase::kTPC);
      baseAnalysis->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
      baseAnalysis->SetTPCpid();
      baseAnalysis->SetMinTPCClusters(100);
      baseAnalysis->SetMaxChi2PerTPCCluster(2.2);
      baseAnalysis->SetMaxCov11(0.5);
      baseAnalysis->SetMaxCov22(0.5);
      baseAnalysis->SetMaxCov33(0.5);
      baseAnalysis->SetMaxCov44(0.5);
      baseAnalysis->SetMaxCov55(0.5);
      baseAnalysis->SetMaxSigmaToVertexTPC(2.0);
      //baseAnalysis->SetMaxDCAXYTPC(1.5);
      //baseAnalysis->SetMaxDCAZTPC(1.5);
      break;
    case "Hybrid":
      baseAnalysis->SetAnalysisMode(AliProtonAnalysisBase::kHybrid);
      baseAnalysis->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
      baseAnalysis->SetTPCpid();
      baseAnalysis->SetMinTPCClusters(110);
      baseAnalysis->SetMaxChi2PerTPCCluster(2.2);
      baseAnalysis->SetMaxCov11(0.5);
      baseAnalysis->SetMaxCov22(0.5);
      baseAnalysis->SetMaxCov33(0.5);
      baseAnalysis->SetMaxCov44(0.5);
      baseAnalysis->SetMaxCov55(0.5);
      baseAnalysis->SetMaxSigmaToVertex(2.0);
      /*baseAnalysis->SetMaxDCAXY(1.5);
	baseAnalysis->SetMaxDCAZ(1.5);*/
      baseAnalysis->SetPointOnITSLayer6();
      baseAnalysis->SetPointOnITSLayer5();
      //baseAnalysis->SetPointOnITSLayer4();
      //baseAnalysis->SetPointOnITSLayer3();
      baseAnalysis->SetPointOnITSLayer2();
      baseAnalysis->SetPointOnITSLayer1();
      baseAnalysis->SetMinITSClusters(5);
      break;
    case "Global":
      baseAnalysis->SetAnalysisMode(AliProtonAnalysisBase::kGlobal);
      baseAnalysis->SetPhaseSpace(20, -1.0, 1.0, 48, 0.3, 1.5);
      baseAnalysis->SetMinTPCClusters(110);
      baseAnalysis->SetMaxChi2PerTPCCluster(2.2);
      baseAnalysis->SetMaxCov11(0.5);
      baseAnalysis->SetMaxCov22(0.5);
      baseAnalysis->SetMaxCov33(0.5);
      baseAnalysis->SetMaxCov44(0.5);
      baseAnalysis->SetMaxCov55(0.5);
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
      /*TFile *f = TFile::Open("PriorProb/PriorProbabilities.root ");
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
    case "Sigma":
      baseAnalysis->SetPIDMode(AliProtonAnalysisBase::kSigma);
      break;
    default:
      break;
    }//PID mode
  }//ESD
  if(analysisLevel == "MC") 
    baseAnalysis->SetPhaseSpace(56, -1.0, 1.0, 16, 0.1, 1.5);

  return baseAnalysis;
}
