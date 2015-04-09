void MakeOCDBConfigPreprocessor(const char* storageUri="local://$ALICE_ROOT/../AliRoot/OCDB", Int_t firstRun=0, Int_t lastRun=999999999)
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(storageUri);
  // create the object and set some values:
  AliCDBEntry *cdbEntry = new AliCDBEntry();
  TEnv *confEnv = new TEnv();

  confEnv->SetValue("Pedestal","DAQ");
  confEnv->SetValue("Signal","DAQ");
  confEnv->SetValue("Temperature","ON");
  confEnv->SetValue("ErrorHandling","ON");

  cdbEntry->SetObject(confEnv);

  // done; now save it..; add some metadata business etc.
  Int_t firstRun   =  0;
  Int_t lastRun    =  999999999;
  Int_t beamPeriod =  1;

  AliCDBMetaData md;
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("David Silvermyr");
  
  AliCDBId id("EMCAL/Config/Preprocessor", firstRun, lastRun);

  cdbEntry->SetId(id);
  cdbEntry->SetMetaData(&md);

  cdb->Put(cdbEntry);

}
