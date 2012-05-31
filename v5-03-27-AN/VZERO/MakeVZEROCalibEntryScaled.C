
void MakeVZEROCalibEntryScaled(Int_t run, Int_t scale = 1, const char *inputCDB = "raw://", const char *outputCDB = "local://$ALICE_ROOT/OCDB/VZERO/test"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);

  man->SetSpecificStorage("VZERO/Calib/PMGains",
			  "local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("VZERO/Calib/LightYields",
			  "local://$ALICE_ROOT/OCDB");

  man->SetRun(run);

  AliCDBEntry *entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibdaorg = (AliVZEROCalibData*)entry->GetObject();
  AliVZEROCalibData *calibda = new AliVZEROCalibData(*calibdaorg);

  Float_t b[64] = {  7.40,  6.83,  7.02,  6.94,  7.03,  7.04,  7.79,  7.27,
		     6.92,  6.96,  7.01,  6.90,  7.28,  7.38,  7.33,  7.23,
		     6.71,  7.05,  7.17,  7.69,  7.41,  7.38,  7.21,  7.11,
		     7.26,  7.12,  6.98,  7.35,  6.99,  6.79,  7.13,  7.58,
		     6.95,  7.01,  7.33,  7.01,  7.21,  6.01,  7.34,  6.44,
		     5.68,  7.12,  6.07,  6.92,  7.04,  6.82,  7.04,  7.24,
		     7.53,  6.99,  7.10,  6.89,  7.07,  6.35,  6.88,  5.77,
		     6.81,  7.01,  6.89,  6.84,  6.68,  6.95,  6.73,  7.14};

  Float_t lightYieldCorr[64] = {0.01051, 0.00955, 0.00861, 0.00948, 0.01082, 0.00870, 0.01023, 0.01012, 0.01270, 0.01184, 0.01110, 0.01266, 0.00956, 0.00826, 0.00966, 0.00891, 0.01358, 0.01543, 0.01516, 0.01337, 0.01908, 0.01641, 0.01767, 0.01512, 0.01664, 0.01326, 0.01536, 0.00586, 0.01439, 0.01445, 0.01504, 0.01079, 0.00105, 0.00110, 0.00143, 0.00093, 0.00072, 0.01919, 0.00073, 0.02580, 0.02911, 0.00148, 0.03176, 0.00126, 0.00158, 0.00111, 0.02804, 0.00109, 0.00157, 0.00158, 0.00104, 0.00120, 0.00123, 0.00188, 0.00193, 0.03442, 0.00200, 0.00185, 0.00143, 0.00257, 0.00201, 0.00151, 0.00197, 0.00282};

  for (Int_t i = 0; i < 64; ++i) {
    Float_t hv = calibdaorg->GetMeanHV(i)*TMath::Exp(-TMath::Log((Float_t)scale)/b[i]);
    calibda->SetMeanHV(hv,i);
    Float_t mip = (i < 32) ? 6950 : 33690; 
    printf("%d   %.1f %.1f %.1f   %.1f %.1f %.1f  %.1f\n",
	   i,
	   calibdaorg->GetMeanHV(i),
	   1./calibdaorg->GetMIPperADC(i),
	   mip*lightYieldCorr[i]*0.18*TMath::Qe()*calibdaorg->GetGain(i)/0.6e-12,
	   calibda->GetMeanHV(i),
	   1./calibda->GetMIPperADC(i),
	   mip*lightYieldCorr[i]*0.18*TMath::Qe()*calibda->GetGain(i)/0.6e-12,
	   calibda->GetMIPperADC(i)/calibdaorg->GetMIPperADC(i));
  }

  // Creation of the object VZERO Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment(Form("VZERO Calibration from RAW OCDB (HV gains are scaled down by factor of %d)",scale));
  AliCDBId id("VZERO/Calib/Data",run,run);

  man->SetDefaultStorage(outputCDB);
  AliCDBStorage *storLoc = man->GetDefaultStorage();
  storLoc->Put(calibda, id, md);

}
