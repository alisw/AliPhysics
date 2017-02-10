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
  commonConfig();

  // Create TGeant3
  //  
  new  TGeant3TGeo("C++ Interface to Geant3");

  // AliRoot event generator
  // (it has to be created after MC, as it may use decayer via VMC)
  //
  genExtFileConfig();

  cout << "Running g3Config.C finished ... " << endl;
}  
