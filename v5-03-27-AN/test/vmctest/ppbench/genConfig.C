// $Id$
//
// Configuration macro for primary event generation for ppbench test (in vmctest).
//
// By I. Hrivnacova, IPN Orsay

void Config()
{
  cout << "Running genConfig.C ... " << endl;

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/commonConfig.C");
  commonConfig(kFALSE);

  // Load Geant3 + Geant3 VMC libraries
  //
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif

  // Create TGeant3
  //  
  new  TGeant3TGeo("C++ Interface to Geant3");

  // AliRoot event generator
  // (it has to be created after MC, as it may use decayer via VMC)
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/genPPbenchConfig.C");
  genGunConfig();

  cout << "Running genConfig.C finished ... " << endl;
}  
