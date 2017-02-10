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

            // The event generator selection (srun) is done in genPPbenchConfig.C
            // that´s why we have to load it too
  genExtFileConfig(srun);

  cout << "Running g3Config.C finished ... " << endl;
}  
