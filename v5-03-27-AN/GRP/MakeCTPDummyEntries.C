void MakeCTPDummyEntries(){

	// Example macro to put in OCDB the dummy entries for CTP configuration and scalers
	// The entries are at present taken from $ALICE_ROOT 
	// Should be used to test the GRP preprocessor 

	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	Char_t * filenameConfig = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.cfg");
	Char_t * filenameScalers = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.cnt");
        Char_t * filenameCTPtime = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.tim");

	AliTriggerConfiguration *runcfg = AliTriggerConfiguration::LoadConfiguration(filenameConfig);
	AliTriggerRunScalers *scalers = AliTriggerRunScalers::ReadScalers(filenameScalers);
	AliCTPTimeParams *ctptime = AliCTPTimeParams::LoadCTPTimeParams(filenameCTPtime);

	// CTP configuration
	AliCDBMetaData* metaconfig = new AliCDBMetaData();
	metaconfig->SetResponsible("Roman Lietava");
	metaconfig->SetComment("Dummy CTP configuration for standalone runs");
	AliCDBId idconfig("GRP/CTP/DummyConfig",0,AliCDBRunRange::Infinity());
	man->Put(runcfg,idconfig, metaconfig);

	// CTP scalers
	AliCDBMetaData* metascalers = new AliCDBMetaData();
	metascalers->SetResponsible("Roman Lietava");
	metascalers->SetComment("Dummy CTP scalers for standalone runs");
	AliCDBId idscalers("GRP/CTP/DummyScalers",0,AliCDBRunRange::Infinity());
	man->Put(scalers,idscalers, metascalers);

	// CTP time parameters
	AliCDBMetaData* metactptime = new AliCDBMetaData();
	metactptime->SetResponsible("Roman Lietava");
	metactptime->SetComment("Dummy CTP time params for standalone runs");
	AliCDBId idctptime("GRP/CTP/DummyCTPtime",0,AliCDBRunRange::Infinity());
	man->Put(ctptime,idctptime, metactptime);

	// CTP LTU configuration

	TObjArray* ltuarray = new TObjArray();
	ltuarray->SetOwner(1);
	AliLTUConfig* ltu;
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("ITSSPD"),14000.,16459.,13);
	ltuarray->AddAtAndExpand(ltu,0);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("ITSSDD"),3126.,16459.,8);
	ltuarray->AddAtAndExpand(ltu,1);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("ITSSSD"),3126.,16459.,17);
	ltuarray->AddAtAndExpand(ltu,2);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("TPC"),3126.,16459.,15);
	ltuarray->AddAtAndExpand(ltu,3);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("TRD"),3126.,16459.,17);
	ltuarray->AddAtAndExpand(ltu,4);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("TOF"),3126.,16459.,14);
	ltuarray->AddAtAndExpand(ltu,5);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("HMPID"),3126.,16459.,19);
	ltuarray->AddAtAndExpand(ltu,6);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("PHOS"),3126.,16459.,19);
	ltuarray->AddAtAndExpand(ltu,7);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("CPV"),3126.,16459.,16);
	ltuarray->AddAtAndExpand(ltu,8);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("PMD"),3126.,16459.,22);
	ltuarray->AddAtAndExpand(ltu,9);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("MUONTRK"),3126.,16459.,8);
	ltuarray->AddAtAndExpand(ltu,10);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("MUONTRG"),3126.,16459.,11);
	ltuarray->AddAtAndExpand(ltu,11);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("FMD"),3126.,16459.,17);
	ltuarray->AddAtAndExpand(ltu,12);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("T0"),3126.,16459.,15);
	ltuarray->AddAtAndExpand(ltu,13);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("VZERO"),2000.,16459.,12);
	ltuarray->AddAtAndExpand(ltu,14);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("ZDC"),3126.,16459.,17);
	ltuarray->AddAtAndExpand(ltu,15);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("ACORDE"),16126.,22459.,18);
	ltuarray->AddAtAndExpand(ltu,16);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("EMCAL"),3126.,16459.,19);
	ltuarray->AddAtAndExpand(ltu,17);
	ltu = new AliLTUConfig((Char_t)AliDAQ::DetectorID("DAQ_TEST"),3126.,16459.,10);
	ltuarray->AddAtAndExpand(ltu,18);

	AliCDBMetaData* md = new AliCDBMetaData();
	md->SetResponsible("Roman Lietava");
	md->SetComment("Example of (dummy -> default settings from 23/11/2010) entry for the detectors' LTU config");
	AliCDBId id("GRP/CTP/DummyLTUConfig",0,AliCDBRunRange::Infinity());
	man->Put(ltuarray,id, md);
	// check if ok
	/*
	man->SetRun(0);
	TObjArray* ltuarrayR = (TObjArray*) man->Get("GRP/CTP/LTUConfig")->GetObject();
	cout << "Array size: " << ltuarrayR->GetEntriesFast() << endl;
        for(Int_t i=0;i<ltuarrayR->GetEntriesFast();i++){
	  if(ltu=(AliLTUConfig*) ltuarrayR->At(i)) ltu->Print();
	  else cout << "--------------------->Empty position " << i << endl;
	}
	*/

	return;
}
