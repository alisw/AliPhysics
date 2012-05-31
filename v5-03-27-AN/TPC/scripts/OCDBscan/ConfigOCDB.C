//
// Macro to Setup OCDB for calibration scan
// By default - all entries taken from the AliEn OCDB storage 
// This is just example macro
// Responsible: marian.ivanov@cern.ch

 

void ConfigOCDB(Int_t crun=-1){
  // 
  printf("SETUP OCBD for TPC\n");
  //
  AliCDBManager::Instance()->SetDefaultStorage("raw://");  
  Int_t run =crun;
  if (run<0) run =0;
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
  AliTPCcalibDB::Instance()->SetRun(run); 
}



