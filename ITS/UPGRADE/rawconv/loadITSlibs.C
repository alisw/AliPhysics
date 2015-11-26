{
  gROOT->ProcessLine(".L PixConv.cxx+g");
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
}
