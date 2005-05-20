void DBAccessTutorial(){


/***************************************************************************************************************
 Prelude(0): 	if it doesn't exist, create DB directory: $ALICE_ROOT/DB
 ***************************************************************************************************************/
TString DBFolder="$ALICE_ROOT/DBTutorial";
gSystem->ExpandPathName(DBFolder);

if(!gSystem->OpenDirectory(DBFolder)){
	printf("Warning: folder DBTutorial does not exist, I will create it!");
        TString command = "mkdir "+ DBFolder;
	gSystem->Exec(command.Data());
}

/***************************************************************************************************************
 Prelude(1): 	create a DB object. We'll take AliZDCCalibData 
		(container class for ZDC Calibration data) as an example
 ***************************************************************************************************************/
 
AliZDCCalibData *calibda=new AliZDCCalibData("ZDC");

// Fill array of ZDC energy calibration factors (4 float's) with some numbers
Float_t EnFactorsArray[4]={20.,25.,200.,10.};
calibda->SetEnCalib(EnFactorsArray);

printf("\nI am in Prelude(1). ZDC Energy calibration factors: \n");  
printf("PM0=%f; PM1=%f; PM2=%f; PM3=%f \n",
   calibda->GetEnCalib(0),  calibda->GetEnCalib(1),  calibda->GetEnCalib(2),  calibda->GetEnCalib(3) );  

/***************************************************************************************************************
 Part 1: 	Store the calibration object in the database, using 
		the "Organized" storage system
		(i.e. a well defined directory structure, defined by the 
		object's metadata and set in the DBFodler above specified)
 ***************************************************************************************************************/

// a. Define the object's metadata (set of information which describes the object)
// AliObjectMetaData(char* name, int firstRun, int lastRun,int Period,char* ObjectFormat, char* ResponsibleName, char* extraInfo)
//
// name: Object's name ("Detector/DBType/DetSpecType")
// firstRun: first run for which the object is valid
// lastRun:  last run for which the object is valid
// period: beam period (to be better specified in future)
// objectFormat: a string which describes the object's format
// ResponsibleName: name of the person who created and stored the object
// extraInfo: anything else you may want to know

AliObjectMetaData omd("ZDC/Calib/GainFactors",1,10,1,"AliZDCCalibData: Pedestals (47 floats), En. calib factors (4 floats)","A. Colla", "Tutorial");

// b. Create the specific AliRunDataStorage instance (AliRunDataOrganizedFile in this case)

AliRunDataOrganizedFile *orgf = new AliRunDataOrganizedFile(DBFolder.Data());

// c. Put the object in the database. Use Put method: AliRunDataStorage::Put(TObject object, AliObjectMetaData& metaData)
// 
// the object will be stored in the file DBFolder/ZDC/Calib/GainFactors/Run1-10_v0.root
// Note that the version is automatically determined.

AliRunDataStorage::Instance()->Put(calibda, omd); // we could also type: orgf->Put(calibda, omd)

// delete the AliRunDataStorage instance and the object (optional)
AliRunDataStorage::Instance()->Delete();


/***************************************************************************************************************
 Part 2: 	Now we want to retrieve the object valid, say, for run 5 from the database folder. 
		We will dump it in a local file called "dummyDBTutorail.root", for later usage. 
 ***************************************************************************************************************/

// a. Create the specific AliRunDataStorage instance (AliRunDataOrganizedFile in this case)
AliRunDataOrganizedFile *orgf = new AliRunDataOrganizedFile(DBFolder.Data());

// b. Prepare AliRunDataStorage to save the object in the local file
AliRunDataStorage::Instance()->RecordToFile("dummyDBTutorial.root");

// c. Get the object! 
// TObject *AliRunDataStorage::Get(char* name, int runNumber)
// if no selection criteria are specified, the highest version found is retrieved.

AliZDCCalibData *calibda2 = (AliZDCCalibData*) AliRunDataStorage::Instance()->Get("ZDC/Calib/GainFactors",5);

printf("\nI am in Part2. ZDC Energy calibration factors: \n");  
printf("PM0=%f; PM1=%f; PM2=%f; PM3=%f \n",
   calibda2->GetEnCalib(0),  calibda2->GetEnCalib(1),  calibda2->GetEnCalib(2),  calibda2->GetEnCalib(3) );  

AliRunDataStorage::Instance()->Delete();


/***************************************************************************************************************
 Part 3: 	Now we will retrieve the object from the local file. 
		We will tune the energy calibration factors and we suppose that the new object is valid from run 1 to run 15.
		Finally we will store the object in DBFolder again with the new run range specified in the object's metadata.  
 ***************************************************************************************************************/

// a. Create the specific AliRunDataStorage instance (AliRunDataFile in this case)
AliRunDataFile *df = new AliRunDataFile("dummyDBTutorial.root");
//AliRunDataOrganizedFile *orgf = new AliRunDataOrganizedFile(DBFolder.Data());

// b. Get the object.
AliZDCCalibData *calibda3 = (AliZDCCalibData*) AliRunDataStorage::Instance()->Get("ZDC/Calib/GainFactors",5);

//AliZDCCalibData calibda3copy=*calibda3;

// c. Tune the energy calibration factors.
calibda3->SetEnCalib(21.5, 0); calibda3->SetEnCalib(26.3, 1); 
calibda3->SetEnCalib(201., 2); calibda3->SetEnCalib(9.45, 3);

printf("\nI am in Part3. New ZDC Energy calibration factors: \n");  
printf("PM0=%f; PM1=%f; PM2=%f; PM3=%f \n",
   calibda3->GetEnCalib(0),  calibda3->GetEnCalib(1),  calibda3->GetEnCalib(2),  calibda3->GetEnCalib(3) );  

// d. Get the object's metadata, set the new run range
AliObjectMetaData omd2 = AliRunDataStorage::Instance()->GetObjectMetaData("ZDC/Calib/GainFactors");
omd2.SetRunRange(1,15);

// e. Store the object in the "organized" database.
// the object will be stored in the file DBFolder/ZDC/Calib/GainFactors/Run1-15_v1.root
// Note that, since there is an "overlap" of the new run range with the one specified
// in the previous version, the version is automatically increased.

AliRunDataOrganizedFile *orgf = new AliRunDataOrganizedFile(DBFolder.Data());
AliRunDataStorage::Instance()->Put(calibda3, omd2);

/***************************************************************************************************************
 Part 3: 	Last act. 
		We want to retrieve the object from the "organized" database folder again.
		This time, anyway, we don't want the highest version but version 0 instead.
		Therefore we have to specify the version wanted using the "Select" method of AliRunDataStorage
		and the aliSelectionMetaData object.  
 ***************************************************************************************************************/

// a. Create the specific AliRunDataStorage instance (AliRunDataOrganizedFile in this case)
orgf = new AliRunDataOrganizedFile(DBFolder.Data());

// b. Specify the selection criterion. We want version 0 for the object valid for run range (1,100)
// NOTE(1): TRegexp expression are valid! E.g. we could type "ZDC/Calib/*" if we want to specify version
// for all the "ZDC/Calib" subsamples...
// Note(2): we could also not specify the run range, if we want version 0 for any retrieved object
AliSelectionMetaData smd("ZDC/Calib/GainFactors",1,100,0);
AliRunDataStorage::Instance()->Select(smd);

// c. Get the object
AliZDCCalibData *calibda4 = (AliZDCCalibData*) AliRunDataStorage::Instance()->Get("ZDC/Calib/GainFactors",5);
printf("\nI am in Part4. I've just got \"ZDC/Calib/GainFactors\" object valid for run 5, version 0.\n");
printf("ZDC Energy calibration factors: \n");  
printf("PM0=%f; PM1=%f; PM2=%f; PM3=%f \n",
   calibda4->GetEnCalib(0),  calibda4->GetEnCalib(1),  calibda4->GetEnCalib(2),  calibda4->GetEnCalib(3) );  
   
printf("\nEnd of tutorial. Bye bye!! \n");


}
