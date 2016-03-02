/* $Id: TestPreprocessor.C 30923 2009-02-09 15:47:08Z hristov $ */

// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle

void ADTestPreprocessor()
{
  // load library
  gSystem->Load("$ALICE_ROOT/../src/SHUTTLE/TestShuttle/libTestShuttle");

   // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(666, 0, 123450);

  // TODO if needed, change location of OCDB and Reference test folders
  // by default they are set to $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB and TestReference
  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/OCDB");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB");

  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());


  // TODO(1)
  //
  // The shuttle can read DCS data, if the preprocessor should be tested to process DCS data,
  // some fake data has to be created.
  //
  // The "fake" input data can be taken using either (a) or (b):
  // (a) data from a file: Use ReadDCSAliasMap()
  //     the format of the file is explained in ReadDCSAliasMap()
  //     To use it uncomment the following line:
  //
  //TMap* dcsAliasMap = ReadDCSAliasMap();
  //
  // (b) generated in this macro: Use CreateDCSAliasMap() and its documentation
  //     To use it uncomment the following line:
  //
  TMap* dcsAliasMap = CreateDCSAliasMap();

  // now give the alias map to the shuttle
  shuttle->SetDCSInput(dcsAliasMap);

  // TODO(2)
  //
  // The shuttle can also process files that originate from DCS, DAQ and HLT.
  // To test it, we provide some local files and locations where these would be found when
  // the online machinery would be there.
  // In real life this functions would be produces by the sub-detectors
  // calibration programs in DCS, DAQ or HLT. These files can then be retrieved using the Shuttle.
  //
  // Files are added with the function AliTestShuttle::AddInputFile. The syntax is:
  // AddInputFile(<system>, <detector>, <id>, <source>, <local-file>)
  // In this example we add a file originating from the GDC with the id PEDESTALS
  // Three files originating from different LDCs but with the same id are also added
  // Note that the test preprocessor name is TPC. The name of the detector's preprocessor must follow
  // the "online" naming convention ALICE-INT-2003-039.
//  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "TPC", "PEDESTALS", "GDC0", "file1.root");
//  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "TPC", "DRIFTVELOCITY", "LDC0", "file2a.root");
//  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "TPC", "DRIFTVELOCITY", "LDC1", "file2b.root");
//  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "TPC", "DRIFTVELOCITY", "LDC2", "file2c.root");
//  shuttle->AddInputFile(AliShuttleInterface::kHLT, "TPC", "HLTData", "source1", "hlt_file1.root");
  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "AD0", "AD0da_results", "source1", "/home/mbroz/AD/Data/db/run000000000_DAQ_test_AD0da_results");
  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "AD0", "AD0da_slewing", "source2", "/home/mbroz/AD/Data/db/run000000000_DAQ_test_AD0da_slewing");
//
  // TODO(3)
  //
  // The shuttle can read run type stored in the DAQ logbook.
  // To test it, we must provide the run type manually. They will be retrieved in the preprocessor
  // using GetRunType function.
  shuttle->SetInputRunType("PHYSICS");

  // TODO(4)
  //
  // The shuttle can read run parameters stored in the DAQ run logbook.
  // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
  // using GetRunParameter function.
  shuttle->AddInputRunParameter("totalEvents", "30000");
  shuttle->AddInputRunParameter("NumberOfGDCs", "15");

  // TODO(5)
  //
  // This is for preprocessor that require data from HLT.
  // Since HLT may be switched off, the preprocessor should first query the Run logbook where
  // the HLT status is stored. SHUTTLE implements a shortcut function (GetHLTStatus) that returns
  // a bool directly. 1 = HLT ON, 0 = HLT OFF
  //

  Bool_t hltStatus=kFALSE;
  shuttle->SetInputHLTStatus(hltStatus);


  // TODO(6)
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  AliPreprocessor* adPreprocessor = new AliADPreprocessor(shuttle);

  // Test the preprocessor
  shuttle->Process();

  // TODO(7)
  // In the preprocessor AliShuttleInterface::Store should be called to put the final
  // data to the CDB. To check if all went fine have a look at the files produced in
  // $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/SHUTTLE/Data
  //
  // Check the file which should have been created
  AliCDBEntry* chkEntry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get("AD/Calib/Data", 0);
  if (!chkEntry)
  {
    printf("The file is not there. Something went wrong.\n");
    return;
  }

	
  AliTestDataDCS* output = dynamic_cast<AliTestDataDCS*> (chkEntry->GetObject());
  // If everything went fine, draw the result
  if (output)
    output->Draw();
  //  
}

TMap* CreateDCSAliasMap()
{
  // Creates a DCS structure
  // The structure is the following:
  //   TMap (key --> value)
  //     <DCSAlias> --> <valueList>
  //     <DCSAlias> is a string
  //     <valueList> is a TObjArray of AliDCSValue
  //     An AliDCSValue consists of timestamp and a value in form of a AliSimpleValue

  // In this example 6 aliases exists: DCSAlias1 ... DCSAlias6
  // Each contains 1000 values randomly generated by TRandom::Gaus + 5*nAlias

  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(1);
  TRandom random;

	FILE *fp = fopen("./DCSValues.txt","r");

	char name[50];
	Float_t val;
	while(!(EOF == fscanf(fp,"%s %f",name,&val))){
		TObjArray* valueSet = new TObjArray;
		valueSet->SetOwner(1);

		TString aliasName=name;
		
		//printf("alias: %s\t\t",aliasName.Data());

		int timeStamp=10;
		
		
		if(aliasName.Contains("HV")) {
			for(int i=0;i<20;i++){
				dcsVal = new AliDCSValue((Float_t) (val+random.Gaus(0,val*0.1)), timeStamp+10*i);
				valueSet->Add(dcsVal);
			}
		} else {
			for(int i=0;i<2;i++){
				AliDCSValue* dcsVal = new AliDCSValue((UInt_t) (val), timeStamp+10*i);
				valueSet->Add(dcsVal);
			}
		}
		
		aliasMap->Add(new TObjString(aliasName), valueSet);

	}
	fclose(fp);
  return aliasMap;
}

TMap* ReadDCSAliasMap()
{
  // Open a file that contains DCS input data
  // The CDB framework is used to open the file, this means the file is located
  // in $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/DCS/Data
  // The file contains an AliCDBEntry that contains a TMap with the DCS structure.
  // An explanation of the structure can be found in CreateDCSAliasMap()

  AliCDBEntry *entry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("DET/DCS/Data", 0);
  return dynamic_cast<TMap*> (entry->GetObject());
}

void WriteDCSAliasMap()
{
  // This writes the output from CreateDCSAliasMap to a CDB file

  TMap* dcsAliasMap = CreateDCSAliasMap();

  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Responsible person");
	metaData.SetComment("Test object for TestPreprocessor.C");

  AliCDBId id("DET/DCS/Data", 0, 0);

  // look into AliTestShuttle's CDB main folder

  AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Put(dcsAliasMap, id, &metaData);
}
