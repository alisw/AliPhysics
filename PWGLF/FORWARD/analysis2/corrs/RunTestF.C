void
RunTestF(UShort_t single=0, Bool_t old=false)
{
  gSystem->AddIncludePath(Form("%s-DTEST_SHIFT -DTEST_FITTER -I$ANA_SRC -I.", 
			       (old ? "-DNO_SIGMA_SHIFT " : "")));
  gROOT->SetMacroPath(Form("%s:%s",gROOT->GetMacroPath(),"$ANA_SRC"));
  gROOT->LoadMacro("AliLandauGaus.h++g");
  gROOT->LoadMacro("AliLandauGausFitter.h+g");
  gROOT->LoadMacro("TestF.C++g");
  
  Double_t vv[] = { 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4 };
  TArrayD v(9, vv);

  // TestShift ts;
  // ts.ScanOne(true, v, 8);
  // ts.ScanTwo(v, 8);

  TestFit tf;
  tf.ScanOne(true, 0, v.GetArray(), 5);

  if (gROOT->IsBatch()) return;
  TFile::Open("shiftSigmaXi.root");
  TFile::Open("fitSigma.root");
  new TBrowser();
}

  
