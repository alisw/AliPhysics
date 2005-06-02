void AliPHOSCalibDataTest(){


  /***************************************************************************************************************
 Prelude(0): 	if it doesn't exist, create DB directory: $ALICE_ROOT/DB
  ***************************************************************************************************************/
  TString DBFolder="PhosCalibDB";
  gSystem->ExpandPathName(DBFolder);

  if(!gSystem->OpenDirectory(DBFolder)){
    printf("Warning: folder PhosCalibDB does not exist, I will create it!");
    TString command = "mkdir "+ DBFolder;
    gSystem->Exec(command.Data());
  }

  /***************************************************************************************************************
 Prelude(1): 	create a DB object. We'll take AliPHOSCalibData 
		(container class for PHOS Calibration data) as an example
  ***************************************************************************************************************/
 
  AliPHOSCalibData *calibda=new AliPHOSCalibData("PHOS");

  // Fill array of PHOS energy calibration factors (4 float's) with some numbers
  Float_t EnFactorsArray[5][64][56];
  for(Int_t module=0; module<5; module++) {
    for(Int_t column=0; column<64; column++) {
      for(Int_t row=0; row<56; row++) {
	calibda->SetADCchannelEmc(module,column,row,1.1);
	calibda->SetADCpedestalEmc(module,column,row,0.5);
      }
    }
  }

  printf("\nI am in Prelude(1). PHOS Energy calibration factors: \n");  
  calibda->Print("ped gain");

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

  AliObjectMetaData omd("PHOS/Calib/GainFactors",1,10,1,"AliPHOSCalibData: Pedestals and ADC channels (5x64x56 floats)","Yuri Kharlov", "PHOS calibration test");

  // b. Create the specific AliRunDataStorage instance (AliRunDataOrganizedFile in this case)

  AliRunDataOrganizedFile *orgf = new AliRunDataOrganizedFile(DBFolder.Data());

  // c. Put the object in the database. Use Put method: AliRunDataStorage::Put(TObject object, AliObjectMetaData& metaData)
  // 
  // the object will be stored in the file DBFolder/PHOS/Calib/GainFactors/Run1-10_v0.root
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
  AliRunDataStorage::Instance()->RecordToFile("dummyPhosCalibDB.root");

  // c. Get the object! 
  // TObject *AliRunDataStorage::Get(char* name, int runNumber)
  // if no selection criteria are specified, the highest version found is retrieved.

  AliPHOSCalibData *calibda2 = (AliPHOSCalibData*) AliRunDataStorage::Instance()->Get("PHOS/Calib/GainFactors",5);

  printf("\nI am in Part2. PHOS Energy calibration factors: \n");  
  calibda2->Print("gain");

  AliRunDataStorage::Instance()->Delete();


  /***************************************************************************************************************
 Part 3: 	Now we will retrieve the object from the local file. 
		We will tune the energy calibration factors and we suppose that the new object is valid from run 1 to run 15.
		Finally we will store the object in DBFolder again with the new run range specified in the object's metadata.  
  ***************************************************************************************************************/

  // a. Create the specific AliRunDataStorage instance (AliRunDataFile in this case)
  AliRunDataFile *df = new AliRunDataFile("dummyPhosCalibDB.root");
  //AliRunDataOrganizedFile *orgf = new AliRunDataOrganizedFile(DBFolder.Data());

  // b. Get the object.
  AliPHOSCalibData *calibda3 = (AliPHOSCalibData*) AliRunDataStorage::Instance()->Get("PHOS/Calib/GainFactors",5);

  //AliPHOSCalibData calibda3copy=*calibda3;

  // c. Tune the energy calibration factors in the 5th module
  for(Int_t module=4; module<5; module++) {
    for(Int_t column=0; column<64; column++) {
      for(Int_t row=0; row<56; row++) {
	calibda3->SetADCchannelEmc(module,column,row,1.2);
      }
    }
  }

  printf("\nI am in Part3. New PHOS Energy calibration factors: \n");  
  calibda3->Print("gain");

  // d. Get the object's metadata, set the new run range
  AliObjectMetaData omd2 = AliRunDataStorage::Instance()->GetObjectMetaData("PHOS/Calib/GainFactors");
  omd2.SetRunRange(1,15);

  // e. Store the object in the "organized" database.
  // the object will be stored in the file DBFolder/PHOS/Calib/GainFactors/Run1-15_v1.root
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
  // NOTE(1): TRegexp expression are valid! E.g. we could type "PHOS/Calib/*" if we want to specify version
  // for all the "PHOS/Calib" subsamples...
  // Note(2): we could also not specify the run range, if we want version 0 for any retrieved object
  AliSelectionMetaData smd("PHOS/Calib/GainFactors",1,100,0);
  AliRunDataStorage::Instance()->Select(smd);

  // c. Get the object
  AliPHOSCalibData *calibda4 = (AliPHOSCalibData*) AliRunDataStorage::Instance()->Get("PHOS/Calib/GainFactors",5);
  printf("\nI am in Part4. I've just got \"PHOS/Calib/GainFactors\" object valid for run 5, version 0.\n");
  printf("PHOS Energy calibration factors: \n");  
  calibda4->Print("gain");
   
  printf("\nEnd of tutorial. Bye bye!! \n");


}
