// $Id$
//
// Configuration macro for running aliroot with Geant3
// with primary events read from external file.
//
// By I. Hrivnacova, IPN Orsay

void Config(const TString& det)
{
  cout << "Running g3Config.C ... " << endl;

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/lego/commonConfig.C");
  commonConfig(det);

  // Load Geant3 + Geant3 VMC libraries
  //
#if defined(__CINT__)
    gSystem->Load("libgeant321");
#endif

  // Create TGeant3
  //  
  new  TGeant3TGeo("C++ Interface to Geant3");

  cout << "Running g3Config.C finished ... " << endl;
}  
