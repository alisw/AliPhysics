
/**************************************************************************
	
	Macro created for storing the OCDB Config Rates data of
	ACORDE in $ALICE_ROOT/OCDB/ACORDE/Config/Rates

	From: 

		Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>
		FCFM, BUAP, Puebla, Mexico

	Created:

		March 3rd. 2009 @ CERN


	Further commnents:

		Arturo Fernandez <afernan@mail.cern.ch>

**************************************************************************/
void MakeACORDEOCDBConfigRate()
{

	AliCDBManager *man = AliCDBManager::Instance();

	AliCDBStorage *storLoc;
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");


	AliACORDECalibData *calibda = new AliACORDECalibData("OCDBConfgRates");

        Float_t Rates[60] = 
        {
                1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
                1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
                1.14, 0.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94,
                1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
                1.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94, 0.94,
                0.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94, 0.94
        };

	calibda->SetRates(Rates);
	
	// Creation of the object ACORDE Calibration as a MetaData
        
	TObjString str("ACORDE OCDB Reference Config Rates Data");      // object that will be stored

	AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object

	AliCDBId id("ACORDE/Config/Rates",0,9999999);

	md->SetResponsible("Mario Rodriguez");
	md->SetBeamPeriod(0);
	md->SetAliRootVersion("Current trunk version");
	md->SetComment("Version 1.0 of OCDB Reference Config Data for ACORDE");
	md->PrintMetaData();

	storLoc = man->GetDefaultStorage();
	storLoc->Put(calibda, id, md);

	storLoc->Delete();
	delete md;

}

