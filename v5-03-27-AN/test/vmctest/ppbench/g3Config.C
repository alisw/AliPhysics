// $Id$
//
// Configuration macro for running aliroot with Geant3
// with primary events read from external file.
//
// By I. Hrivnacova, IPN Orsay

void Config()
{
  cout << "Running g3Config.C ... " << endl;

  // AliRoot setup
  //
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/commonConfig.C");
  commonConfig();

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
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/genPPbenchConfig.C");
  gROOT->LoadMacro("$ALICE_ROOT/test/vmctest/ppbench/genExtFileConfig.C");
            // The event generator selection (srun) is done in genPPbenchConfig.C
            // that´s why we have to load it too
  genExtFileConfig(srun);

  cout << "Running g3Config.C finished ... " << endl;
}  
