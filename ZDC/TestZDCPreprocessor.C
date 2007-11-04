/* $Id$ */

// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle

void TestZDCPreprocessor()
{
  // load library
  gSystem->Load("libTestShuttle.so");

  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(7, 0, 1);

  // TODO if needed, change location of OCDB and Reference test folders
  // by default they are set to $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB and TestReference
  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestReference");

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
  //WriteDCSAliasMap();

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
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "PEDESTALS", "LDC0", "ZDCPedestal.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "EMDCALIB",  "LDC0", "ZDCEMDCalib.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "EMDCALIB",  "LDC0", "ZDCEMDEqual.dat");

  // TODO(3)
  //
  // The shuttle can read run type stored in the DAQ logbook.
  // To test it, we must provide the run type manually. They will be retrieved in the preprocessor
  // using GetRunType function.
//  shuttle->SetInputRunType("PEDESTALS");
  shuttle->SetInputRunType("PULSER_RUN");

  // TODO(4)
  //
  // The shuttle can read run parameters stored in the DAQ run logbook.
  // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
  // using GetRunParameter function.
  // In real life the parameters will be retrieved automatically from the run logbook;
  shuttle->AddInputRunParameter("totalEvents", "1000");
  shuttle->AddInputRunParameter("NumberOfGDCs", "1");

  // TODO(5) NEW!
  //
  // This is for preprocessor that require data from HLT.
  // Since HLT may be switched off, the preprocessor should first query the Run logbook where
  // the HLT status is stored. SHUTTLE implements a shortcut function (GetHLTStatus) that returns
  // a bool directly. 1 = HLT ON, 0 = HLT OFF
  //
  Bool_t hltStatus=kFALSE;
  //shuttle->SetInputHLTStatus(hltStatus);

  // TODO(6)
  //
  // The shuttle can query condition parameters valid from the current run from the OCDB
  // To test it, we must first store the object into the OCDB. It will be retrieved in the preprocessor
  // using GetFromOCDB function.
/*
  TObjString obj("This is a condition parameter stored in OCDB");
  AliCDBId id("ZDC/Calib/Data", 0, AliCDBRunRange::Infinity());
  AliCDBMetaData md;
  AliCDBEntry entry(&obj, id, &md);

  shuttle->AddInputCDBEntry(&entry);
*/
  // TODO(6)
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  AliPreprocessor* test = new AliZDCPreprocessor(shuttle);

  // Test the preprocessor
  shuttle->Process();

  // TODO(7)
  // In the preprocessor AliShuttleInterface::Store should be called to put the final
  // data to the CDB. To check if all went fine have a look at the files produced in
  // $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/SHUTTLE/Data
  //
  // Check the file which should have been created
  AliCDBEntry* chkEntry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("ZDC/Calib/Data", 7);
  if (!chkEntry)
  {
    printf("The file is not there. Something went wrong.\n");
    return;
  }

  AliTestDataDCS* output = dynamic_cast<AliTestDataDCS*> (chkEntry->GetObject());
  // If everything went fine, draw the result
  if (output)
    output->Draw();
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
  TString aliasNames[28];

  // ******************************** alignment values
  aliasNames[0] = "ZDC_ZNA_POS.actual.position";
  aliasNames[1] = "ZDC_ZPA_POS.actual.position";
  aliasNames[2] = "ZDC_ZNC_POS.actual.position";
  aliasNames[3] = "ZDC_ZPC_POS.actual.position";
  //
  for(int nAlias=0; nAlias<4; nAlias++)
  {
    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);

    TString aliasName = aliasNames[nAlias];
    printf("\n\n alias: %s\n\n",aliasName.Data());

    Float_t simVal = (Float_t) (random.Rndm()*0.025+random.Rndm()*0.1);
    for(int i=0;i<3;i++)
    {
      int timeStamp1[3] = {0,500,1000};
      AliDCSValue* dcsVal = new AliDCSValue(simVal, timeStamp1[i]);
      printf("%s\n",dcsVal->ToString());
      valueSet->Add(dcsVal);
    }
    aliasMap->Add(new TObjString(aliasName), valueSet);
  }
  // ******************************** HV values
  /*TString ZNAAlias = "ZNA_HV.actual.vMon";
  TString ZPAAlias = "ZPA_HV.actual.vMon";
  TString ZNCAlias = "ZNC_HV.actual.vMon";
  TString ZPCAlias = "ZPC_HV.actual.vMon";
  TString idat[5];
  for(int i=0;i<5;i++)
  {
    idat[i] = i;
    aliasNames[i+3]  = ZNAAlias.Insert(6,idat[i]);
    aliasNames[i+7]  = ZPAAlias.Insert(6,idat[i]);
    aliasNames[i+11] = ZNCAlias.Insert(6,idat[i]);
    aliasNames[i+15] = ZPCAlias.Insert(6,idat[i]);
  }*/
  aliasNames[4]  = "ZDC_ZNA_HV0.actual.vMon";
  aliasNames[5]  = "ZDC_ZNA_HV1.actual.vMon";
  aliasNames[6]  = "ZDC_ZNA_HV2.actual.vMon";
  aliasNames[7]  = "ZDC_ZNA_HV3.actual.vMon";
  aliasNames[8]  = "ZDC_ZNA_HV4.actual.vMon";
  //
  aliasNames[9]   = "ZDC_ZPA_HV0.actual.vMon";
  aliasNames[10]  = "ZDC_ZPA_HV1.actual.vMon";
  aliasNames[11]  = "ZDC_ZPA_HV2.actual.vMon";
  aliasNames[12]  = "ZDC_ZPA_HV3.actual.vMon";
  aliasNames[13]  = "ZDC_ZPA_HV4.actual.vMon";
  //
  aliasNames[14]  = "ZDC_ZNC_HV0.actual.vMon";
  aliasNames[15]  = "ZDC_ZNC_HV1.actual.vMon";
  aliasNames[16]  = "ZDC_ZNC_HV2.actual.vMon";
  aliasNames[17]  = "ZDC_ZNC_HV3.actual.vMon";
  aliasNames[18]  = "ZDC_ZNC_HV4.actual.vMon";
  //
  aliasNames[19]  = "ZDC_ZPC_HV0.actual.vMon";
  aliasNames[20]  = "ZDC_ZPC_HV1.actual.vMon";
  aliasNames[21]  = "ZDC_ZPC_HV2.actual.vMon";
  aliasNames[22]  = "ZDC_ZPC_HV3.actual.vMon";
  aliasNames[23]  = "ZDC_ZPC_HV4.actual.vMon";
  //
  aliasNames[24]  = "ZDC_ZEM_HV0.actual.vMon";
  aliasNames[25]  = "ZDC_ZEM_HV1.actual.vMon";
  //
  aliasNames[26]  = "ZDC_REFA_HV0.actual.vMon";
  aliasNames[27]  = "ZDC_REFC_HV1.actual.vMon";
  //
  for(int nAlias=4;nAlias<28;nAlias++)
  {
    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);
   
    TString aliasName = aliasNames[nAlias];
    printf("\n\n alias: %s\n\n",aliasName.Data());

    for(int timeStamp=0;timeStamp<=1000;timeStamp+=500)
    {
      Float_t simVal = (Float_t) (random.Gaus()*600.+1800.);
      AliDCSValue* dcsVal = new AliDCSValue(simVal, timeStamp);
      printf("%s\n",dcsVal->ToString());
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

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("ZDC/DCS/Data", 0);
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

  AliCDBId id("ZDC/DCS/Data", 0, 0);

  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");

  AliCDBManager::Instance()->Put(dcsAliasMap, id, &metaData);
}
