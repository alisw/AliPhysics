
void MakeVZEROCalibEntryPbPb(Int_t run,const char *inputCDB = "raw://"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);
  man->SetRun(run);

  AliCDBEntry *entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibdaorg = (AliVZEROCalibData*)entry->GetObject();
  AliVZEROCalibData *calibda = new AliVZEROCalibData(*calibdaorg);

  const Float_t scale = 70.;
  Float_t p1[64] = {
    0.135569, 0.146405, 0.142425, 0.144278, 0.142307, 0.141648, 0.128477, 0.138239, 
    0.144173, 0.143419, 0.143572, 0.144482, 0.138024, 0.136542, 0.135955, 0.138537, 
    0.148521, 0.141999, 0.139627, 0.130014, 0.134970, 0.135635, 0.139094, 0.140634, 
    0.137971, 0.142080, 0.142793, 0.136054, 0.142778, 0.146045, 0.139133, 0.142080, 
    0.144121, 0.142311, 0.136564, 0.142686, 0.138792, 0.166285, 0.136387, 0.155391, 
    0.176082, 0.140408, 0.164738, 0.144270, 0.142766, 0.147486, 0.141951, 0.138012, 
    0.132394, 0.142849, 0.140477, 0.144592, 0.141558, 0.157646, 0.143758, 0.173385, 
    0.146489, 0.143279, 0.145230, 0.147203, 0.147333, 0.144979, 0.148597, 0.138985 };

  for (Int_t i = 0; i < 64; ++i) {
    Float_t hvorg = calibdaorg->GetMeanHV(i);
    Float_t hvnew = hvorg/(TMath::Power(scale,p1[i]));
    calibda->SetMeanHV(hvnew,i);
  }

  // Creation of the object VZERO Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Pb-Pb VZERO Calibration from RAW OCDB");
  AliCDBId id("VZERO/CalibPbPb/Data",0,AliCDBRunRange::Infinity());

  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage *storLoc = man->GetDefaultStorage();
  storLoc->Put(calibda, id, md);

}
