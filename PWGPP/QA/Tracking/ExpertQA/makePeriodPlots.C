makePeriodPlots( const char* mergedQAfile, const char *periodName ){
//what do I do?
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_ROOT/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx++");
  gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++");

  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis();
    a->ConnectGenericHistos( mergedQAfile );
    a->SetPeriodName( periodName );
    a->SetMakeFitPerfomancePlots(kTRUE);
    a->RunPeriod();
}
