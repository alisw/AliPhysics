//
// This macro puts the gains in the database
//
void AliPMDPutCalibration()
{
  AliPMDCalibData *calibda=new AliPMDCalibData();
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBId id("PMD/Calib/Data",0,0);
  AliCDBMetaData *md=new AliCDBMetaData();
  man->GetDefaultStorage()->Put(calibda,id,md);
}
