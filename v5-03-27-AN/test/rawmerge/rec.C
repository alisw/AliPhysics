void rec()
{
  AliReconstruction rec;
  AliCDBManager *cdbm = AliCDBManager::Instance();
  //cdbm->SetRun(run);
  cdbm->SetDefaultStorage("local://$ALICE_ROOT/OCDB");     
  //cdbm->SetDefaultStorage("raw://");     
  cdbm->SetSpecificStorage("GRP/GRP/Data",Form("local://%s","$ALICE_ROOT/test/rawmerge"));  
  rec.SetRunQA(":");
  rec.SetRunReconstruction("ITS TPC");
  rec.Run();
}
