makePeriodPlots( const char *ch, const char *periodName ){

//  gSystem->Load("/hera/alice/tbroeker/testground/trackDump/highPt/AliHighPtTreeAnalysis_C.so");
  gSystem->Load("/hera/alice/tbroeker/highPt/AliHighPtTreeAnalysis_C.so");

  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis();
    a->ConnectGenericHistos( ch );
    a->SetPeriodName( periodName );
    a->SetMakeFitPerfomancePlots(kTRUE);
    a->RunPeriod();

}
