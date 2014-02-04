void DBStorageCalib(){

AliCDBManager *man = AliCDBManager::Instance();

AliCDBStorage *storLoc;
man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");


AliADCalibData *calibda = new AliADCalibData("Calib");

Float_t Efficiencies[16] = { 
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  0.94, 0.94, 0.94, 0.94, 0.94, 0.94
};
Float_t Rates[16] = {

  1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  1.14, 0.94, 0.94, 0.94, 0.94, 0.94


                    };
  
calibda->SetEfficiencies(Efficiencies);
calibda->SetRates(Rates);

// Creation of the object AD Calibration as a MetaData
        
TObjString str("AD Calibration");      // object that will be stored

AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object

AliCDBId id("AD/Calib/Data",0,9999999);

md->SetResponsible("Mario Rodriguez Cahuantzi");
md->SetBeamPeriod(0);
md->SetAliRootVersion("v5-04-Rev-08");
md->SetComment("Prototype");
md->PrintMetaData();

storLoc = man->GetDefaultStorage();
storLoc->Put(calibda, id, md);

storLoc->Delete();
delete md;

}

