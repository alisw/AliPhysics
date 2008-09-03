/* $Id$ */

// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
// CreateDCSAliasMap() creates input that would in the same way come from DCS
// ReadDCSAliasMap() reads from a file
// CreateInputFilesMap() creates a list of local files,
// that can be accessed by the shuttle

void TestPreprocessor()
{
    // load library
    gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");

    // create AliTestShuttle instance
    // The parameters are run, startTime, endTime
    AliTestShuttle* shuttle = new AliTestShuttle(7, 0, 1);


    printf("Test Shuttle temp dir: %s\n", AliShuttleInterface::GetShuttleTempDir());

    // TODO if needed, change location of OCDB and Reference test folders
    // by default they are set to $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB
    // and TestReference

    AliTestShuttle::SetMainCDB("local://TestCDB");
    AliTestShuttle::SetMainRefStorage("local://TestReference");
    
    printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
    printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());
    
    
    // TODO(1)
    //
    // The shuttle can read DCS data, if the preprocessor should be
    // tested to process DCS data,
    // some fake data has to be created.
    //
    // The "fake" input data can be taken using either (a) or (b):
    // (a) data from a file: Use ReadDCSAliasMap()
    //     the format of the file is explained in ReadDCSAliasMap()
    //     To use it uncomment the following line:
    //
    // TMap* dcsAliasMap = ReadDCSAliasMap();
    //
    // (b) generated in this macro: Use CreateDCSAliasMap() and
    //     its documentation
    //     To use it uncomment the following line:
    //

    TMap* dcsAliasMap = CreateDCSAliasMap();

    // now give the alias map to the shuttle

    shuttle->SetDCSInput(dcsAliasMap);

    // TODO(2)
    //
    // The shuttle can also process files that originate from DCS, DAQ and HLT.
    // To test it, we provide some local files and locations where these would
    // be found when
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
    
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "PMD", "PMD_PED.root", "GDC0", "PMD_PED.root");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "PMD", "PMDGAINS.root", "GDC0", "xy.root");
    
    // TODO(3)
    //
    // The shuttle can read run type stored in the DAQ logbook.
    // To test it, we must provide the run type manually. They will be retrieved in the preprocessor
    // using GetRunType function.

    shuttle->SetInputRunType("PEDESTAL_RUN");

    // TODO(4)
    //
    // The shuttle can read run parameters stored in the DAQ run logbook.
    // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
    // using GetRunParameter function.

    shuttle->AddInputRunParameter("totalEvents", "30000");
    shuttle->AddInputRunParameter("NumberOfGDCs", "15");

    // TODO(5)
    //
    // The shuttle can query condition parameters valid from the current run from the OCDB
    // To test it, we must first store the object into the OCDB. It will be retrieved in the preprocessor
    // using GetFromOCDB function.
    
    TObjString obj("This is a condition parameter stored in OCDB");
    AliCDBId id("TPC/Calib/Data", 0, AliCDBRunRange::Infinity());
    AliCDBMetaData md;
    AliCDBEntry entry(&obj, id, &md);

    shuttle->AddInputCDBEntry(&entry);

    
    // TODO(6)
    // Create the preprocessor that should be tested, it registers
    // itself automatically to the shuttle
    // AliPreprocessor* test = new AliTestPreprocessor(shuttle);
    AliPMDPreprocessor* test = new AliPMDPreprocessor(shuttle);

    // Test the preprocessor
    shuttle->Process();
    
    // TODO(7)
    // In the preprocessor AliShuttleInterface::Store should be called to put the final
    // data to the CDB. To check if all went fine have a look at the files produced in
    // $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/SHUTTLE/Data
    //
    // Check the file which should have been created
    AliCDBEntry* chkEntry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
	->Get("PMD/Calib/Ped", 7);
    if (!chkEntry)
    {
	printf("The file is not there. Something went wrong.\n");
	return;
    }
    
    AliPMDPedestal* output = dynamic_cast<AliPMDPedestal*> (chkEntry->GetObject());
    // If everything went fine, draw the result
    if (output)
	output->Print();
}

TMap* CreateDCSAliasMap()
{
    // Creates a DCS structure
    // The structure is the following:
    // TMap (key --> value)
    //     <DCSAlias> --> <valueList>
    //     <DCSAlias> is a string
    //     <valueList> is a TObjArray of AliDCSValue
    //     An AliDCSValue consists of timestamp and a value in form
    //     of a AliSimpleValue
    
    // In this example 6 aliases exists: DCSAlias1 ... DCSAlias6
    // Each contains 1000 values randomly generated by TRandom::Gaus + 5*nAlias

    TMap* aliasMap = new TMap;
    aliasMap->SetOwner(1);
    
    TRandom random;

    for(int nAlias=0;nAlias<6;nAlias++)
    {
	TObjArray* valueSet = new TObjArray;
	valueSet->SetOwner(1);
	
	TString aliasName="DCSAlias";
	aliasName += nAlias;
	//printf("\n\n alias: %s\n\n",aliasName.Data());
	
	for (int timeStamp=0;timeStamp<1000;timeStamp+=10)
	{
	    AliDCSValue* dcsVal = new AliDCSValue((Float_t) (random.Gaus()+5*nAlias), timeStamp);
	    //printf("%s\n",dcsVal->ToString().Data());
	    valueSet->Add(dcsVal);
	}
	aliasMap->Add(new TObjString(aliasName), valueSet);
    }
    
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
