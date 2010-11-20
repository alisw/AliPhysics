void MakeCTPLTUConfigEntry(const char *cdbStorage = "local://$ALICE_ROOT/OCDB"){

	// Example macro to put in OCDB an LTU Config entry
	// Since it is just an example macro, the AliLTUConfig members will be set as 
	// follows for all the detectors:
	// fFineDelay1 = 0;
	// fFineDelay2 = 0;
	// fBCDelaysAdd = 0;

	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage(cdbStorage);

	TObjArray* ltuarray = new TObjArray();
	ltuarray->SetOwner(1);
	AliLTUConfig* ltu;
	for(Int_t i = 0; i<AliDAQ::kNDetectors-2; i++){
		const char* name = AliDAQ::DetectorName(i);
		ltu = new AliLTUConfig((UChar_t)AliDAQ::DetectorID(name),0.,0.,0.);
		ltuarray->AddAtAndExpand(ltu,i);
	}

	AliCDBMetaData* md = new AliCDBMetaData();
	md->SetResponsible("Roman Lietava");
	md->SetComment("Example of (dummy -> everything set to 0) entry for the detectors' LTU config");
	AliCDBId id("GRP/CTP/LTUConfig",0,AliCDBRunRange::Infinity());
	man->Put(ltuarray,id, md);

	return;
}
