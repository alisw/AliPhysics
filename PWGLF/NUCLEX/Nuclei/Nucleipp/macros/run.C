Int_t run(const TString& config)
{
//
// load libraries for the analysis
// and run the config macro
//
	TStopwatch timer;
	timer.Start();
	
	gROOT->SetMacroPath(Form("%s:$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/", gROOT->GetMacroPath()));
	
	gSystem->AddIncludePath("-I\"$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/\"");
	
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnM2.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnGaussianExpTail.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnSecondaries.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnEfficiency.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnCorr.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnPt.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnRatio.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnSpectra.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnBA.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/AliLnDriver.cxx+g");
	gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/Nuclei/B2/macros/RatioMult.C+g");
	
	gROOT->ProcessLine(Form(".x %s", config.Data()));
	
	timer.Stop();
	timer.Print();
	
	return 0;
}
