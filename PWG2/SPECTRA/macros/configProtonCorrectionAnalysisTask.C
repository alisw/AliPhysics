//______________________________________________________
AliProtonCorrectionAnalysisTask* GetAliProtonCorrectionAnalysisTask(const char* mode = "ESD",const char* analysisType  = "Hybrid",const char* pidMode = "Bayesian",Bool_t fIsOn_AliProtonAbsorptionCorrection=kTRUE, Bool_t fIsOn_AliProtonFeedDownAnalysis=kTRUE,Bool_t fIsOn_AliProtonSpectraCorrection=kTRUE) {
  AliProtonCorrectionAnalysisTask *taskProtons = new AliProtonCorrectionAnalysisTask("TaskProtonsProtonCorrection");
  if(fIsOn_AliProtonAbsorptionCorrection||fIsOn_AliProtonFeedDownAnalysis||fIsOn_AliProtonSpectraCorrection) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonAnalysisBaseObject.C"); 
    AliProtonAnalysisBase *baseAnalysis = GetProtonAnalysisBaseObject(mode,analysisType,pidMode);
    taskProtons->SetBaseAnalysis(baseAnalysis);
  }	
  else
    return 0x0;
  if(fIsOn_AliProtonAbsorptionCorrection) {
    AliProtonAbsorptionCorrection* absorptioncorrection=new AliProtonAbsorptionCorrection();
    taskProtons->SetAnalysisObjectAbsorptionCorrection(absorptioncorrection);
  }
  if(fIsOn_AliProtonFeedDownAnalysis) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonFeedDownAnalysis.C");
    AliProtonFeedDownAnalysis* analysisFeedDown = GetProtonFeedDownAnalysisObject();
    taskProtons->SetAnalysisObjectFeedDown(analysisFeedDown);
  }	
  if(fIsOn_AliProtonSpectraCorrection) {
    AliProtonSpectraCorrection* spectracorrection=new AliProtonSpectraCorrection();
    taskProtons->SetAnalysisObjectSpectraCorrection(spectracorrection);
  }
  return taskProtons;
}


