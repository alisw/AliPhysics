void MakeVZEROSaturationEntry(const char *cdbUri = "local://$ALICE_ROOT/OCDB",
			      Bool_t default = kTRUE,
			      Int_t firstRun = 0,
			      Int_t lastRun = AliCDBRunRange::Infinity())
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbUri);

  // Creation of the signal saturation OCDB object
  // x = Total ADC charge (summed up in 8x25 ns - from -1 to +6)
  // The correction for saturation is:
  // x -> x + alpha * ( x - beta )
  // alpha and beta are calculated for each channel, and are
  // listed below
Double_t alpha[64] = { 7.70e-01 , 7.97e-01 , 7.64e-01 , 7.80e-01 , 6.61e-01 , 7.46e-01 , 6.59e-01 , 6.69e-01 , 6.40e-01 , 6.11e-01 , 6.25e-01 , 4.51e-01 , 5.72e-01 , 3.34e-01 , 6.23e-01 , 5.92e-01 , 5.28e-01 , 4.18e-01 , 3.56e-01 , 3.39e-01 , 4.91e-01 , 3.36e-01 , 9.48e-02 , 1.86e-01 , 1.72e-01 , 2.65e-01 , 2.15e-01 , 1.24e-01 , 2.62e-01 , 2.57e-01 , 2.21e-01 , 1.11e-01 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 1.74e-01 , 0.00e+00 , 1.69e-01 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 1.45e-01 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 5.62e-02 , 0.00e+00 , 0.00e+00 , 1.46e-01 };
Double_t beta[64] = { 1.40e+03 , 1.38e+03 , 1.42e+03 , 1.38e+03 , 1.40e+03 , 1.40e+03 , 1.38e+03 , 1.40e+03 , 1.34e+03 , 1.38e+03 , 1.40e+03 , 1.24e+03 , 1.42e+03 , 1.36e+03 , 1.34e+03 , 1.36e+03 , 1.41e+03 , 1.40e+03 , 1.41e+03 , 1.32e+03 , 1.34e+03 , 1.40e+03 , 1.21e+03 , 1.36e+03 , 1.36e+03 , 1.39e+03 , 1.38e+03 , 1.40e+03 , 1.44e+03 , 1.47e+03 , 1.44e+03 , 1.34e+03 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 1.18e+03 , 0.00e+00 , 1.22e+03 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 1.34e+03 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 1.22e+03 , 0.00e+00 , 0.00e+00 , 1.31e+03 };

  TObjArray *arr = new TObjArray(64);
  arr->SetOwner(1);
  for(Int_t i = 0; i < 64; ++i) {
    TF1 *saturation = new TF1(Form("VZEROSaturationCh_%d",i),"x < [1] ? x : x + [0] * (x - [1])",0.,2500.);
    if (default) {
      saturation->SetParameter(0,0.0);
      saturation->SetParameter(1,0.0);
    }
    else {
      saturation->SetParameter(0,alpha[i]);
      saturation->SetParameter(1,beta[i]);
    }
    arr->AddAt(saturation,i);
  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Signal saturation correction used in reconstruction of data");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/Saturation",firstRun,lastRun);

  storLoc->Put(arr, id, md);

  storLoc->Delete();
  delete md;
}
