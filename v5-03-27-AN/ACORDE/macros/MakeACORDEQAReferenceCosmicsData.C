
/**************************************************************************

	Macro created to storage in $ALICE_ROOT/QAref/ACORDE/QA/Calib
	the reference data for ACORDE's Quality Assurance

	By now, it is the same information for

	../ACORDE/QA/Calib
	../ACORDE/QA/Cosmics
	../ACORDE/QA/HighMultiplicity
	../ACORDE/QA/LowMultiplicity

	From: 

		Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch>
		FCFM, BUAP, Puebla, Mexico

	Created:

		March 3rd. 2009 @ CERN


	Further commnents:

		Arturo Fernandez <afernan@mail.cern.ch>

**************************************************************************/
void MakeACORDEQAReferenceCosmicsData()
{

	AliCDBManager *man = AliCDBManager::Instance();

	AliCDBStorage *storLoc;
	man->SetDefaultStorage("local://$ALICE_ROOT/QAref");


	AliACORDECalibData *calibda = new AliACORDECalibData("QACosmics");

	Float_t ModulesActivity[60] = 
	{
 		0.92,0.51,0.68,0.76,0.78,0.83,0.00,0.69,0.72,0.86,
		0.86,0.85,0.79,0.75,0.79,0.62,0.82,0.92,0.79,0.78,
		0.00,0.90,0.84,0.95,0.79,0.87,0.91,0.88,0.92,0.82,
		0.80,0.98,1.00,0.89,0.82,0.89,0.85,0.92,0.88,0.91,
		0.86,0.00,0.86,0.92,0.88,0.81,0.45,0.84,0.86,0.60,
		0.84,0.86,0.74,0.24,0.71,0.82,0.56,0.00,0.00,0.79		
	};
  
	calibda->SetModulesActivity(ModulesActivity);
	
	// Creation of the object ACORDE Calibration as a MetaData
        
	TObjString str("ACORDE QA Reference Cosmics Data");      // object that will be stored

	AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object

	AliCDBId id("ACORDE/QA/Cosmic",0,9999999);

	md->SetResponsible("Mario Rodriguez");
	md->SetBeamPeriod(0);
	md->SetAliRootVersion("Current trunk version");
	md->SetComment("Version 1.0 of QA Reference Cosmics Data for ACORDE");
	md->PrintMetaData();

	storLoc = man->GetDefaultStorage();
	storLoc->Put(calibda, id, md);

	storLoc->Delete();
	delete md;

}

