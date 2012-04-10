/*

 Macro to initialize: 
 - the OCDB (run number required as input argument)
 - the geometry (expected to be in the current directory)
 to run the Calibration train.
 
 Example:
 .L $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/ConfigCalibTrain.C
 ConfigCalibTrain(129160,"raw://");

*/

void ConfigCalibTrain(Int_t run, const char *ocdb="raw://"){

  // OCDB
  printf("setting run to %d\n",run);
  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(run); 
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {  
    AliCDBManager::Instance()->InitFromSnapshot("OCDB.root");
  }
  // magnetic field
  if ( !TGeoGlobalMagField::Instance()->GetField() ) {
    printf("Loading field map...\n");
    AliGRPManager grpMan;
    if( !grpMan.ReadGRPEntry() ) { 
      printf("Cannot get GRP entry\n"); 
    }
    if( !grpMan.SetMagField() ) { 
      printf("Problem with magnetic field setup\n"); 
    }
  }

  // geometry
  printf("Loading geometry...\n");
  AliGeomManager::LoadGeometry();
  if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD TOF HMPID") ) {
    printf("Problem with align objects\n"); 
  }

}
