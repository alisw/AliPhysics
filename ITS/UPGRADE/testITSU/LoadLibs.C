void LoadLibs(Bool_t anlibs=kTRUE)
{
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  gSystem->Load("libITSUpgradeRec.so");
  //
  if (anlibs) {
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
  }
}

