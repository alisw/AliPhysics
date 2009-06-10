//
// This macro puts the gains in the database
//
void AliPMDPutHot()
{
  AliPMDHotData *hotda=new AliPMDHotData();
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBId id("PMD/Calib/Hot",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *md=new AliCDBMetaData();
  man->GetDefaultStorage()->Put(hotda,id,md);
}
