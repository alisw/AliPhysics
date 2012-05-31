
/**************************************************************************

	Macro created for storing the OCDB Calib data of ACORDE
	in $ALICE_ROOT/OCDB/ACORDE/Calib


	From: 

		Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>
		FCFM, BUAP, Puebla, Mexico

	Created:

		March 3rd. 2009 @ CERN


	Further commnents:

		Arturo Fernandez <afernan@mail.cern.ch>

**************************************************************************/
void MakeACORDEOCDBCalib()
{

	AliCDBManager *man = AliCDBManager::Instance();

	AliCDBStorage *storLoc;
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");


	AliACORDECalibData *calibda = new AliACORDECalibData("OCDBCalib");

	Float_t Efficiencies[60] = 
	{ 
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94
	};
	Float_t Rates[60] = 
	{
  		1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		1.14, 0.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94,
  		1.14, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		1.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94, 0.94,
  		0.94, 0.94, 0.94, 0.94, 1.94, 0.94, 0.94, 0.94, 0.94, 0.94
	};
	Float_t ModulesActivity[60] = 
	{
 		0.92,0.51,0.68,0.76,0.78,0.83,0.00,0.69,0.72,0.86,
		0.86,0.85,0.79,0.75,0.79,0.62,0.82,0.92,0.79,0.78,
		0.00,0.90,0.84,0.95,0.79,0.87,0.91,0.88,0.92,0.82,
		0.80,0.98,1.00,0.89,0.82,0.89,0.85,0.92,0.88,0.91,
		0.86,0.00,0.86,0.92,0.88,0.81,0.45,0.84,0.86,0.60,
		0.84,0.86,0.74,0.24,0.71,0.82,0.56,0.00,0.00,0.79		
	};
  	
	calibda->SetEfficiencies(Efficiencies);
        calibda->SetRates(Rates);
	calibda->SetModulesActivity(ModulesActivity);
	
	// Creation of the object ACORDE Calibration as a MetaData
        
	TObjString str("ACORDE OCDB Reference Calib Data");      // object that will be stored

	AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object

	AliCDBId id("ACORDE/Calib/Data",0,9999999);

	md->SetResponsible("Mario Rodriguez");
	md->SetBeamPeriod(0);
	md->SetAliRootVersion("Current trunk version");
	md->SetComment("Version 1.0 of OCDB Reference Calib Data for ACORDE");
	md->PrintMetaData();

	storLoc = man->GetDefaultStorage();
	storLoc->Put(calibda, id, md);

	storLoc->Delete();
	delete md;

}

