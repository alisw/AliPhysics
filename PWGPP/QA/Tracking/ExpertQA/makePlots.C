makePlots( const char *ch, const char *codeDir ){

//  gSystem->Load("/hera/alice/tbroeker/testground/trackDump/highPt/AliHighPtTreeAnalysis_C.so");
  gSystem->Load(Form("%s/AliHighPtTreeAnalysis_C.so",codeDir));

  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis( ch );

     a->Loop();

}
