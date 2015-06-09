void MakeADTimeSlewingEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the default time slewing OCDB object - splines at TOF - no correction
  // Other objects to be created per run by DA and Preprocessor
  const Double_t fTOF[4] = {65.30,65.12,56.54,56.71};
  const Double_t fRes = 256/25;
  
  TList *fListSplines = new TList();
  TSpline3 *fTimeSlewingSpline[16];
  for(Int_t i=0; i<16; i++){
  	TH1F *slew = new TH1F("NoSlewing"," ",100,-4,0);
  	for(Int_t j =0; j<100; j++){slew->SetBinContent(j+1,fRes*fTOF[i/4]); slew->SetBinError(j+1,0.1);}
  
	TString TimeSlewingSplineName = "hTimeSlewingSpline";
	TimeSlewingSplineName += i;
	fTimeSlewingSpline[i] = new TSpline3(slew);
	fTimeSlewingSpline[i]->SetName(TimeSlewingSplineName.Data());
	fListSplines->Add(fTimeSlewingSpline[i]);	
	}
	
  TObjString str("AD Time-slewing correction");

  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time-slewing correction used in reconstruction and MC simulation");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("AD/Calib/TimeSlewing",0,AliCDBRunRange::Infinity());

  storLoc->Put(fListSplines, id, md);

  storLoc->Delete();
  delete md;

}
