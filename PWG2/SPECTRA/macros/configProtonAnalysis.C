//__________________________________________________//
AliProtonAnalysis *GetProtonAnalysisObject(const char* analysisLevel = "ESD",
					   Bool_t kAnalyzeMC = kTRUE,
					   const char* esdAnalysisType = "Hybrid", 
					   const char* pidMode = "Bayesian",
					   Bool_t kUseOnlineTrigger = kFALSE,
					   Bool_t kUseOfflineTrigger = kFALSE,
					   Bool_t kRunQA = kFALSE) {
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C");  
  //Function to setup the AliProtonAnalysis object and return it
  AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(analysisLevel,kAnalyzeMC,esdAnalysisType,pidMode,kUseOnlineTrigger,kUseOfflineTrigger,kRunQA);

  AliProtonAnalysis *analysis = new AliProtonAnalysis();
  analysis->SetBaseAnalysis(baseAnalysis);
  if(baseAnalysis->GetAnalysisMode() == AliProtonAnalysisBase::kGlobal) {
    Double_t gY[17] = {-0.9,-0.75,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.75,0.9};
    Double_t gPt[18] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.1,2.4,2.7,3.1};
    analysis->InitAnalysisHistograms(16,gY,17,gPt);
    //Double_t gPt[27] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.1,2.4,2.7,3.0,3.4,3.8,4.2,4.7,5.1,6.0,7.0,9.0,12.0};
    //analysis->InitAnalysisHistograms(16,gY,26,gPt);
  }
  else
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
					       const char* pidMode = "Bayesian",
					       Bool_t kUseOnlineTrigger = kFALSE,
					       Bool_t kUseOfflineTrigger = kFALSE) {
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C"); 
  //Function to setup the AliProtonQAAnalysis object and return it
  AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(analysisLevel,kTRUE,esdAnalysisType,pidMode,kUseOnlineTrigger,kUseOfflineTrigger);

  AliProtonQAAnalysis *analysis = new AliProtonQAAnalysis();
  analysis->SetBaseAnalysis(baseAnalysis);
  analysis->SetRunMCAnalysis();
  //analysis->SetMCProcessId(4);//4: weak decay - 13: hadronic interaction
  //analysis->SetMotherParticlePDGCode(3122);//3122: Lambda
  analysis->SetRunEfficiencyAnalysis(kFALSE);//use cuts in the eff. analysis
  analysis->SetQAYPtBins(baseAnalysis->GetNBinsX(),
			 baseAnalysis->GetMinX(),
			 baseAnalysis->GetMaxX(),
			 baseAnalysis->GetNBinsY(),
			 baseAnalysis->GetMinY(),
			 baseAnalysis->GetMaxY());

  return analysis;
}

//__________________________________________________//

