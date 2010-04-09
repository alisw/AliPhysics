
void MakeVZEROCalibEntryPbPb(Int_t run,const char *inputCDB = "raw://"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);
  man->SetRun(run);

  AliCDBEntry *entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibdaorg = (AliVZEROCalibData*)entry->GetObject();
  AliVZEROCalibData *calibda = new AliVZEROCalibData(*calibdaorg);

  for (Int_t i = 0; i < 64; ++i) {
    calibda->SetTimeOffset(5.0,i);
  }
  for (Int_t i = 0; i < 8; ++i) {
    calibda->SetWidthResolution(2,i);
  }

  Float_t p0[64] = {
    7.094891, 7.124938, 7.089708, 7.098169, 7.094482, 7.147250, 7.170978, 7.183392, 
    7.145760, 7.148096, 7.153840, 7.143544, 7.186069, 7.194580, 7.203516, 7.195176, 
    7.188333, 7.198607, 7.209412, 7.226565, 7.221695, 7.205132, 7.191238, 7.227724, 
    7.232810, 7.252655, 7.230309, 7.140891, 7.273518, 7.242969, 7.252859, 7.252655, 
    7.026802, 7.079913, 7.134147, 7.092387, 7.079561, 7.072848, 7.123192, 7.003141, 
    7.024667, 7.124784, 7.123442, 7.129744, 7.110671, 7.143031, 7.139439, 7.178109, 
    7.247803, 7.139396, 7.293809, 7.094454, 6.992198, 7.206448, 7.244765, 7.056197, 
    7.263595, 7.138569, 7.089582, 7.215683, 7.266183, 7.165123, 7.243276, 7.235135 };
  Float_t p1[64] = {
    0.135569, 0.146405, 0.142425, 0.144278, 0.142307, 0.141648, 0.128477, 0.138239, 
    0.144173, 0.143419, 0.143572, 0.144482, 0.138024, 0.136542, 0.135955, 0.138537, 
    0.148521, 0.141999, 0.139627, 0.130014, 0.134970, 0.135635, 0.139094, 0.140634, 
    0.137971, 0.142080, 0.142793, 0.136054, 0.142778, 0.146045, 0.139133, 0.142080, 
    0.144121, 0.142311, 0.136564, 0.142686, 0.138792, 0.166285, 0.136387, 0.155391, 
    0.176082, 0.140408, 0.164738, 0.144270, 0.142766, 0.147486, 0.141951, 0.138012, 
    0.132394, 0.142849, 0.140477, 0.144592, 0.141558, 0.157646, 0.143758, 0.173385, 
    0.146489, 0.143279, 0.145230, 0.147203, 0.147333, 0.144979, 0.148597, 0.138985 };

  Float_t mipperadc = 0.5;
  for (Int_t i = 0; i < 64; ++i) {
    Float_t hv = TMath::Exp(p1[i]*TMath::Log(0.6/mipperadc)+p0[i]);
    calibda->SetMeanHV(hv,i);
  }

  // Creation of the object VZERO Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Pb-Pb VZERO Calibration from RAW OCDB");
  AliCDBId id("VZERO/Calib/Data",0,AliCDBRunRange::Infinity());

  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB/VZERO/PbPb");
  AliCDBStorage *storLoc = man->GetDefaultStorage();
  storLoc->Put(calibda, id, md);

}
