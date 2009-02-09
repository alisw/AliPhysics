void AliTRDreadRaw(const char *fname = "raw.root")
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  AliTRDdigitsManager manR;
  manR.CreateArrays();

  AliRawReaderRoot reader(fname, 0);
  reader.Select("TRD");

  while (reader.NextEvent())
    {
      AliTRDrawStreamBase *pstream = AliTRDrawStreamBase::GetRawStream(&reader);
      AliTRDrawStreamBase &stream = *pstream;
      
      Int_t ichambers = 0;
      while (stream.NextChamber(&manR) >= 0)
	ichambers++;
      
      cout << "Chambers loaded " << stream.IsA()->GetName() << " " 
	   << ichambers << endl;
      
      delete pstream
    }
}

void readRaw2(const char *fname = "raw.root")
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  AliRawReaderRoot reader(fname, 0);
  reader.Select("TRD");

  AliTRDrawData rawData;
  AliTRDdigitsManager *man = rawData.Raw2Digits(&reader);
}
