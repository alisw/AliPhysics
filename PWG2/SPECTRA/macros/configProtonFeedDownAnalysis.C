//__________________________________________________________

Double_t test1(Double_t *x,Double_t *p)
{
	/*if(x[0]<1.5) 
		return -1.0*x[0]+2.5; 
	else 
		return 1.0;*/
	if(x[0]<3.0) 
		return -0.2*x[0]+1.6; 
	else 
		return 1.0;	
}
//_____________________________________________________________
Double_t test2(Double_t *x,Double_t *p)
{
	if(x[0]<1.5) 
		return -1.0*x[0]+2.5; 
	else 
		return 1.0;
}

//__________________________________________________//
AliProtonFeedDownAnalysis *GetProtonFeedDownAnalysisObject(const char* analysisLevel = "ESD", 
					   const char* esdAnalysisType = "TPC", 
					   const char* pidMode = "Bayesian") {
gROOT->LoadMacro(" $ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C"); 
  //Function to setup the AliProtonAnalysis object and return it
  AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(analysisLevel,esdAnalysisType,pidMode);

  AliProtonFeedDownAnalysis *analysis = new AliProtonFeedDownAnalysis();
  analysis->SetBaseAnalysis(baseAnalysis);
  //if(analysisBase->GetEtaMode()) analysis->SetEtaMode();
  analysis->InitAnalysisHistograms(baseAnalysis->GetNBinsX(),
				   baseAnalysis->GetMinX(),
				   baseAnalysis->GetMaxX(),
				   baseAnalysis->GetNBinsY(),
				   baseAnalysis->GetMinY(),
				   baseAnalysis->GetMaxY());
TF1* weightfunction=new TF1("weightfunction","1");
//TF1* weightfunction=new TF1("weightfunction",test1,0.5,4.0,0);	
//TF1* weightfunction=new TF1("weightfunction",test2,0.5,4.0,0);			   
 analysis->SetWeightFunction(weightfunction);	   
    
  return analysis;
}


