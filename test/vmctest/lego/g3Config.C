// $Id$
//
// Configuration macro for running aliroot with Geant3
// with primary events read from external file.
//
// By I. Hrivnacova, IPN Orsay

void Config(const char * det)
{
  cout << "Running g3Config.C ... " << endl;

  // AliRoot setup
  //
  commonConfig(det);

  // Create TGeant3
  //  
  new  TGeant3TGeo("C++ Interface to Geant3");

  cout << "Running g3Config.C finished ... " << endl;
}  
