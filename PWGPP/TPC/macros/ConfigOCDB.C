/*

 Macro to initialize: 
 - the OCDB (run number required as input argument)
 - the geometry and mag. field initialized from GRP
 
 Example:
 .L $ALICE_PHYSICS/PWGPP/TPC/macros/ConfigOCDB.C
 ConfigOCDB(129160,"raw://");

*/

void ConfigOCDB(Int_t run, const char *ocdb="raw://") {

  // OCDB
  printf("setting run to %d\n",run);
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    Printf("ConfigOCDB: using OCDB snapshot");
    AliCDBManager::Instance()->SetSnapshotMode("OCDB.root");
  }
  else {
    Printf("ConfigOCDB: NOT using OCDB snapshot");
  }
  Printf("Default storage is %s", ocdb);

  AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  AliCDBManager::Instance()->SetRun(run); 

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
  if ( !TGeoGlobalMagField::Instance()->GetField()){
    AliMagF::BMap_t smag = AliMagF::k5kG;
    Double_t bzfac = 1;
    AliMagF* magF= new AliMagF("Maps","Maps", bzfac, 1., smag);
    TGeoGlobalMagField::Instance()->SetField(magF);
  }

  // geometry
  printf("Loading geometry...\n");
  AliGeomManager::LoadGeometry();
  if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD TOF HMPID") ) {
    printf("Problem with align objects\n"); 
  }

   if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {
    Printf("ConfigOCDB: localOCDBaccessConfig detected\n");
     gROOT->LoadMacro("localOCDBaccessConfig.C");
    gInterpreter->ProcessLine("localOCDBaccessConfig();");
   }
 
}
