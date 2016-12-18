int runtest(const TString &testname) {
  TestTHistManager::THistManagerTestSuite tester;
  if(testname == "build_simple") return tester.TestBuildSimpleHistograms();
  else if(testname == "build_grouped") return tester.TestBuildGroupedHistograms();
  else if(testname == "fill_simple") return tester.TestFillSimpleHistograms();
  else if(testname == "fill_grouped") return tester.TestFillGroupedHistograms();
  else return 1;
}
