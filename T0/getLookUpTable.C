 void getLookUpTable(Int_t run)
{
 AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage("local:///home/alla/alice/aliroot/master/inst/OCDB/");
  man->SetRun(run);
  AliCDBEntry *entry = man->Get("T0/Calib//LookUp_Table");

    cout<<" AliT0CalibData ::GetLookUp :: "<<entry<<endl;
  AliT0CalibData *clb = (AliT0CalibData*)entry->GetObject();
  cout<<" AliT0CalibData *clb "<<clb <<endl;
  clb->PrintLookup();

}
