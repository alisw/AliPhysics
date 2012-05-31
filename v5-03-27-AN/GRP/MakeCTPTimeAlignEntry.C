void MakeCTPTimeAlignEntry(const char *cdbStorage = "local://$ALICE_ROOT/OCDB",Bool_t storeFill1069 = kFALSE){

	// Example macro to put in OCDB the dummy entries for CTP timing params valid for perioid
        // File *.tip interpretation:
        //  1st column = trigger input name
        //  2nd column = dummy (to keep compatibilyu with *.tim file
        //  3rd column =  time between L0 and reference 
        //  4th column =   as in tim file
        //  5th column =   as in tim file
	AliCDBManager *man = AliCDBManager::Instance();
	man->SetDefaultStorage(cdbStorage);
        Char_t * filenameCTPtimeAlignIdeal = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/stdln.tip");
        Char_t * filenameCTPtimeAlign1069 = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/fill1069.tip");

	AliCTPTimeParams *ctptimealignideal = AliCTPTimeParams::LoadCTPTimeParams(filenameCTPtimeAlignIdeal);
	AliCTPTimeParams *ctptimealign1069 = AliCTPTimeParams::LoadCTPTimeParams(filenameCTPtimeAlign1069);

	AliCDBMetaData* metactptimeideal = new AliCDBMetaData();
	AliCDBMetaData* metactptime1069 = new AliCDBMetaData();
	metactptimeideal->SetResponsible("Roman Lietava");
	metactptime1069->SetResponsible("Roman Lietava");
	metactptimeideal->SetComment("CTP time-alignment params (Ideal)");
	metactptime1069->SetComment("CTP time-alignment params for period of runs corresponding to fill 1069");
	AliCDBId idctptimeideal("GRP/CTP/TimeAlign",0,AliCDBRunRange::Infinity());
	AliCDBId idctptime1069("GRP/CTP/TimeAlign",118556,118780);
	man->Put(ctptimealignideal,idctptimeideal, metactptimeideal);
	if (storeFill1069) man->Put(ctptimealign1069,idctptime1069, metactptime1069);

	return;
}
