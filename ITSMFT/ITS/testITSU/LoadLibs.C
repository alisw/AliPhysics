void LoadLibs(Bool_t anlibs=kTRUE)
{
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");
  //
  if (anlibs) {
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
  }
}

