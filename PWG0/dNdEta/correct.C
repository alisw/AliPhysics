void correct()
{
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/dNdEta/run.C");
  FinishAnalysisAll();
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/dNdEta/drawPlots.C");
  dNdEta();
}


