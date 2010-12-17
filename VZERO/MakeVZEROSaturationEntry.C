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
  Double_t alpha[64] = { 7.94e-01 , 8.19e-01 , 7.98e-01 , 7.76e-01 , 6.96e-01 , 8.27e-01 , 6.81e-01 , 5.34e-01 , 6.52e-01 , 6.57e-01 , 5.76e-01 , 5.33e-01 , 4.47e-01 , 2.98e-01 , 4.59e-01 , 5.51e-01 , 5.01e-01 , 4.34e-01 , 3.98e-01 , 3.44e-01 , 1.42e-01 , 3.54e-01 , 0.00e+00 , 1.94e-01 , 1.68e-01 , 2.76e-01 , 2.38e-01 , 9.10e-02 , 2.45e-01 , 2.67e-01 , 2.26e-01 , 1.14e-01 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 1.68e-01 , 0.00e+00 , 4.63e-01 , 0.00e+00 , 6.47e-02 , 0.00e+00 , 0.00e+00 , 6.32e-02 };
  Double_t beta[64] = { 1.42e+03 , 1.41e+03 , 1.42e+03 , 1.38e+03 , 1.50e+03 , 1.52e+03 , 1.47e+03 , 1.55e+03 , 1.34e+03 , 1.36e+03 , 1.38e+03 , 1.29e+03 , 1.52e+03 , 1.49e+03 , 1.42e+03 , 1.52e+03 , 1.39e+03 , 1.40e+03 , 1.38e+03 , 1.29e+03 , 1.28e+03 , 1.47e+03 , 0.00e+00 , 1.46e+03 , 1.39e+03 , 1.40e+03 , 1.41e+03 , 1.45e+03 , 1.51e+03 , 1.57e+03 , 1.52e+03 , 1.38e+03 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 0.00e+00 , 1.49e+03 , 0.00e+00 , 1.32e+03 , 0.00e+00 , 1.37e+03 , 0.00e+00 , 0.00e+00 , 1.34e+03 };

  TObjArray *arr = new TObjArray(64);
  arr->SetOwner(1);
  for(Int_t i = 0; i < 64; ++i) {
    TF1 *saturation = new TF1(Form("VZEROSaturationCh_%d",i),"x < [1] ? x : x + [0] * (x - [1])",0.,2000.);
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
