
AliCDBStorage      *storLoc;
AliCDBStorage      *storGrid;
AliCDBEntry        *entry     = 0;
AliTRDCalPad    *calobj    = 0;
AliTRDCalROC *calobjROC = 0;

void AliTRDDBAccess()
{

  // Single instance of AliCDBManager. 
  AliCDBManager *man = AliCDBManager::Instance();

  printf("\nTRD :: Activating Local storage...\n");
  storLoc = man->GetStorage("local://DBLocal");

  // Create the new database object
  AliTRDCalPad    *calVdrift    = new AliTRDCalPad("Vdrift","TRD drift velocities");
  AliTRDCalROC *calROCVdrift = calVdrift->GetCalROC(0,2,0);
  calROCVdrift->SetValue(1,1,1.50);
  calROCVdrift->SetValue(2,2,1.55);
  calROCVdrift->SetValue(3,3,1.60);

  // Id of the object: AliCDBId("name", firstRun, lastRun)
  AliCDBId id1("TRD/Calib/Vdrift",0,10); 

  // MetaData describing the object
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName("AliTRDCalPad");
  md1->SetResponsible("Christoph Blume");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("05-04-00"); //root version
  md1->SetComment("This is a Vdrift test");
  TObjString *strVdrift = new TObjString("Drift velocities for 540 chambers");
  md1->SetProperty("key1",strVdrift);

  // Store the object into local storage
  // Filename: DBLocal/TRD/Calib/Vdrift/Run0_10_v0_s0.root
  printf("\nTRD :: Storing object into local storage...\n");
  storLoc->Put(calVdrift,id1,md1); 

  // Read, update, store again
  printf("\nTRD :: Retrieve object from local storage...\n");
  entry = storLoc->Get("TRD/Calib/Vdrift",5);

  calobj    = (AliTRDCalPad *) entry->GetObject();
  calobjROC = calobj->GetCalROC(0,2,0);
  printf("\nTRD :: Drift velocities (1): %f (1,1), %f (2,2), %f (3.3)\n"
        ,calobjROC->GetValue(1,1)
        ,calobjROC->GetValue(2,2)
        ,calobjROC->GetValue(3,3));

  // Update object
  calobjROC->SetValue(1,1,1.60);
  calobjROC->SetValue(2,2,1.35);
  calobjROC->SetValue(3,3,1.10);

  // Store into local: filename = Run0_10_v0_s1.root
  printf("\nTRD :: Storing object into local storage...\n");
  storLoc->Put(entry); 

  // Read, update, store again
  printf("\nTRD :: Retrieve object from local storage...\n");
  entry = storLoc->Get("TRD/Calib/Vdrift",5);
  entry = storLoc->Get("TRD/Calib/Vdrift",5);
  entry = storLoc->Get("TRD/Calib/Vdrift",5);

  calobj    = (AliTRDCalPad *) entry->GetObject();
  calobjROC = calobj->GetCalROC(0,2,0);
  printf("\nTRD :: Drift velocities (2): %f (1,1), %f (2,2), %f (3.3)\n"
        ,calobjROC->GetValue(1,1)
        ,calobjROC->GetValue(2,2)
        ,calobjROC->GetValue(3,3));

  printf("\nTRD :: Activating GRID storage...\n");
  storGrid = 0x0;
  storGrid = man->GetStorage("alien://aliendb4.cern.ch:9000;cblume;DBGrid;ALICE::CERN::se01");
  if (!storGrid) {
    printf("\nTRD :: Failure in activating GRID storage...\n");   
  }

  // Store into GRID: filename =  DBGrid/TRD/Calib/Vdrift/Run0_10_v1.root
  if (storGrid) {
    storGrid->Put(entry); 
  }

  entry->Dump();
  cout << entry->IsOwner() << endl;
  
  AliCDBManager::Instance()->Destroy();
  delete entry;
  delete md1;

}
