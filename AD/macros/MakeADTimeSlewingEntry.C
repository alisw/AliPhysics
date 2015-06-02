void MakeADTimeSlewingEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the default time slewing OCDB object - splines at 0 - no correction
  // Other objects to be created per run by DA and Preprocessor
  TH1F *slew = new TH1F("NoSlewing"," ",100,-4,0);
  for(Int_t i =0; i<100; i++){slew->SetBinContent(i+1,0.0); slew->SetBinError(i+1,0.1);}
  
  TList *fListSplines = new TList();
  TSpline3 *fTimeSlewingSpline[16];
  for(Int_t i=0; i<16; i++){
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
