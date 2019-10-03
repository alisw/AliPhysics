void MakeTrend()
{
  const char* fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2";
  gROOT->SetMacroPath(Form("%s:%s:%s/scripts:$(ANA_SRC)/scan", 
			   gROOT->GetMacroPath(), 
			   fwd, fwd));
  gSystem->AddIncludePath(Form("-I%s", fwd));
  gSystem->AddIncludePath(Form("-I%s/scripts", fwd));

  gROOT->LoadMacro("Trend.C++g");

  Trend t;
  t.AddDCCut("mpv", "0.6 0.7 0.8 0.82 0.85");
  t.AddSLCut("fix", "0.05 0.15 0.2 0.3 0.4");
  t.AddSLCut("mpv", "0.01 0.05 0.1");
  t.AddSLCut("sig", "1 2 3");
  t.AddSHCut("mpv", "0.6 0.7 0.8 0.82 0.85");
  t.AddSHCut("sig", "1 2");
  t.AddSHCut("xi",  "1 2");
  t.AddRun("137848");
  t.AddRun("138190");
  t.AddCentrality( 0, 5);
  t.AddCentrality( 5,10);
  t.AddCentrality(10,20);
  t.AddCentrality(20,30);
  t.SetOrder("dc sl sh");
  
  t.Run();

  TFile* f = TFile::Open("trending.root", "READ");
  new TBrowser("b", f);
}

