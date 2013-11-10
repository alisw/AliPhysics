/*

 Macro to initialize: 
 - the OCDB (run number required as input argument)
 - the geometry and mag. field initialized from GRP
 
 Example:
 .L $ALICE_ROOT/PWGPP/TPC/macros/ConfigOCDB.C
 ConfigOCDB(129160,"raw://");

*/

void ConfigOCDB(Int_t run, const char *ocdb="raw://") {

  // OCDB
  printf("setting run to %d\n",run);
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
  if( !AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC") ) {
    printf("Problem with align objects\n"); 
  }
  
  if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {        
    printf("loading localOCDBaccessConfig.C\n");
    gROOT->LoadMacro("localOCDBaccessConfig.C");
    localOCDBaccessConfig();
  }
}
