void DBStorageCalib(){

AliCDBManager *man = AliCDBManager::Instance();

AliCDBStorage *storLoc;
man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");


AliACORDECalibData *calibda = new AliACORDECalibData("Calib");

Float_t Efficiencies[60] = { 
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94
};
Float_t Rates[60] = {

  1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  1.14, 0.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94,
  1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  1.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  0.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94, 0.94


                    };
  
calibda->SetEfficiencies(Efficiencies);
calibda->SetRates(Rates);

// Creation of the object ACORDE Calibration as a MetaData
        
TObjString str("ACORDE Calibration");      // object that will be stored

AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object

AliCDBId id("ACORDE/Calib/Data",0,9999999);

md->SetResponsible("ACORDE team");
md->SetBeamPeriod(0);
md->SetAliRootVersion("v4-05-Release");
md->SetComment("Prototype");
md->PrintMetaData();

storLoc = man->GetDefaultStorage();
storLoc->Put(calibda, id, md);

storLoc->Delete();
delete md;

}

