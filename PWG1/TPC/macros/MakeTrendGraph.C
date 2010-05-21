//
// macro to extract TPC Performance values using AliTPCPerformanceSummary class
// only for one run!
//

void MakeTrendGraph(const char *infilelist, const char* outfile) {

  gSystem->AddIncludePath("-I$ALICE_ROOT/PWG1/TPC");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
  gSystem->Load("libTPCcalib");  
  
  AliTPCPerformanceSummary::ProduceTrends(infilelist, outfile);
}
