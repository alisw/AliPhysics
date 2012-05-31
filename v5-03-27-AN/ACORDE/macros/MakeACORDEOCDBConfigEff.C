
/**************************************************************************
	
	Macro created for storing the OCDB Config Efficiency data of 
	ACORDE in $ALICE_ROOT/OCDB/ACORDE/Config/Efficiency


	From: 

		Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>
		FCFM, BUAP, Puebla, Mexico

	Created:

		March 3rd. 2009 @ CERN


	Further commnents:

		Arturo Fernandez <afernan@mail.cern.ch>

**************************************************************************/
void MakeACORDEOCDBConfigEff()
{

	AliCDBManager *man = AliCDBManager::Instance();

	AliCDBStorage *storLoc;
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");


	AliACORDECalibData *calibda = new AliACORDECalibData("OCDBConfigEff");

 	Float_t Efficiencies[60] = 
	{ 
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94
	};
	calibda->SetEfficiencies(Efficiencies);
	
	// Creation of the object ACORDE Calibration as a MetaData
        
	TObjString str("ACORDE OCDB Config Efficiency Data");      // object that will be stored

	AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object

	AliCDBId id("ACORDE/Config/Efficiency",0,9999999);

	md->SetResponsible("Mario Rodriguez");
	md->SetBeamPeriod(0);
	md->SetAliRootVersion("Current trunk version");
	md->SetComment("Version 1.0 of OCDB Config Efficiency Reference Calib Data for ACORDE");
	md->PrintMetaData();

	storLoc = man->GetDefaultStorage();
	storLoc->Put(calibda, id, md);

	storLoc->Delete();
	delete md;

}

