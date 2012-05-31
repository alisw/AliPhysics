void MakeVZEROEqualizationFactorsEntry(Bool_t default = kTRUE, const char *infile = "alpha.dat")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");


  TH1F *eqFactors = new TH1F("VZEROEqualizationFactors","VZERO Equalization Factors for Pb-Pb",64,-0.5,63.5);
  if (default) {
    const Double_t alpha[66] = {0.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
				0.0};
    eqFactors->SetContent(alpha);
  }
  else {
    FILE *falphas;
    if((falphas = fopen(infile,"r")) == NULL){
      printf("Cannot open file %s",infile);
      return;
    }
    Double_t alpha[66], alpha2[66], beta[66];
    alpha[0] = alpha2[0] = beta[0] = alpha[65] = alpha2[65] = beta[65] = 0;
    Int_t tempCh;
    for(Int_t j=0; j<64; ++j) fscanf(falphas,"%d %lf %lf %lf", &tempCh, &beta[j+1], &alpha[j+1], &alpha2[j+1]);
    fclose(falphas);

    eqFactors->SetContent(alpha2);
  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Default entry for VZERO Equalization Factors object");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/EqualizationFactors",0,AliCDBRunRange::Infinity());

  storLoc->Put(eqFactors, id, md);

  storLoc->Delete();
  delete md;

}
