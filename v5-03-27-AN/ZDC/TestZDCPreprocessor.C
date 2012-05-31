/* $Id$ */

// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle

void TestZDCPreprocessor(const char* runType="PHYSICS")
{
  // load library
  gSystem->Load("libTestShuttle.so");

  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(0, 0, 1);

  // TODO if needed, change location of OCDB and Reference test folders
  // by default they are set to $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB and TestReference
  //AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");
  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/OCDB");
  //AliTestShuttle::SetMainCDB("alien://folder=/alice/data/2009/OCDB/");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestReference");

  printf("\n Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf(" Test Reference storage Uri: %s\n\n", AliShuttleInterface::GetMainRefStorage().Data());

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
  dcsAliasMap->Print("");
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
  //
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "PEDESTALDATA", "LDC", "ZDCPedestal.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "PEDESTALHISTOS", "LDC", "ZDCPedHisto.root");
  //
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "LASERDATA", "LDC", "ZDCLaserCalib.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "LASERHISTOS", "LDC", "ZDCLaserHisto.root");
  //
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "EMDENERGYCALIB", "LDC", "ZDCEnergyCalib.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "EMDTOWERCALIB", "LDC", "ZDCTowerCalib.dat");
  //
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "MBCALIB", "LDC", "ZDCMBCalib.root");
  //
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "MAPPING", "MON", "ZDCChMapping.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "TDCDATA", "MON", "ZDCTDCCalib.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "ZDC", "TDCHISTOS", "MON", "ZDCTDCHisto.root");

  // Todo(3)
  //
  // The shuttle can read run type stored in the DAQ logbook.
  // To test it, we must provide the run type manually. They will be retrieved in the preprocessor
  // using GetRunType function.
  shuttle->SetInputRunType(runType);
  
  // TODO(4)
  //
  // The shuttle can read run parameters stored in the DAQ run logbook.
  // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
  // using GetRunParameter function.
  // In real life the parameters will be retrieved automatically from the run logbook;
  shuttle->AddInputRunParameter("beamType", "A-A");
  shuttle->AddInputRunParameter("beamEnergy", "1380");
  //shuttle->AddInputRunParameter("beamType", "p-p");
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
  shuttle->SetInputHLTStatus(hltStatus);

  // TODO(6)
  //
  // The shuttle can query condition parameters valid from the current run from the OCDB
  // To test it, we must first store the object into the OCDB. It will be retrieved in the preprocessor
  // using GetFromOCDB function.

  /*TObjString obj("This is a condition parameter stored in OCDB");
  AliCDBId id("ZDC/Calib/Data", 0, AliCDBRunRange::Infinity());
  AliCDBMetaData md;
  AliCDBEntry entry(&obj, id, &md);

  shuttle->AddInputCDBEntry(&entry);
  */
  
  // TODO(6)
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  AliPreprocessor* test = new AliZDCPreprocessor(shuttle);
  shuttle->Print();
  
  // Test the preprocessor
  shuttle->Process();
  //printf(" Back to test macro: final checks! \n");

  // TODO(7)
  // In the preprocessor AliShuttleInterface::Store should be called to put the final
  // data to the CDB. To check if all went fine have a look at the files produced in
  // $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/SHUTTLE/Data
  //
  // Check the file which should have been created
  AliCDBEntry* chkEntry0 = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("ZDC/Calib/ChMap", 0);
  TString str(runType);
  AliCDBEntry* chkEntry1;
  if((str.CompareTo("STANDALONE_PEDESTAL")) == 0) chkEntry1 = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("ZDC/Calib/Pedestals", 0);
  else if((str.CompareTo("STANDALONE_LASER")) == 0) chkEntry1 = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("ZDC/Calib/LaserCalib", 0);
  else if((str.CompareTo("CALIBRATION_EMD")) == 0) chkEntry1 = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("ZDC/Calib/EnergyCalib", 0);
  else if((str.CompareTo("PHYSICS")) == 0) chkEntry1 = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("ZDC/Calib/TDCCalib", 0);
  
  
  if(!chkEntry0){
    printf("No file in ZDC/Calib/ChMap\n");
    return;
  }
  if(!chkEntry1){
    if((str.CompareTo("STANDALONE_PEDESTAL")) == 0)  printf("No file in ZDC/Calib/Pedestal\n");
    else if((str.CompareTo("STANDALONE_LASER")) == 0) printf("No file in ZDC/Calib/LaserCalib\n");
    else if((str.CompareTo("CALIBRATION_EMD")) == 0)  printf("No file in ZDC/Calib/EnergyCalib\n");
    else if((str.CompareTo("PHYSICS")) == 0)  printf("No file in ZDC/Calib/TDCCalib\n");
    return;
  }
  

  /*AliTestDataDCS* output = dynamic_cast<AliTestDataDCS*> (chkEntry1->GetObject());
  // If everything went fine, draw the result
  if (output)
    output->Draw();*/
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
  for(int nAlias=0; nAlias<4; nAlias++){
    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);

    TString aliasName = aliasNames[nAlias];
    //printf("\n\n alias: %s\n\n",aliasName.Data());

    Float_t simVal = (Float_t) (random.Rndm()*0.025+random.Rndm()*0.1);
    int timeStamp1[5] = {0,500,1000,1500,2000};
    for(int i=0;i<5;i++){
      AliDCSValue* dcsVal = new AliDCSValue(simVal, timeStamp1[i]);
      //printf("%s\n",dcsVal->ToString());
      valueSet->Add(dcsVal);
    }
    aliasMap->Add(new TObjString(aliasName), valueSet);
  }
  // ******************************** HV values
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
  aliasNames[26]  = "ZDC_REFA_HV.actual.vMon";
  aliasNames[27]  = "ZDC_REFC_HV.actual.vMon";
  //
  for(int nAlias=4;nAlias<28;nAlias++){
//   if(nAlias<14 || nAlias>18){
     TObjArray* valueSet = new TObjArray;
     valueSet->SetOwner(1);
   
     TString aliasName = aliasNames[nAlias];
     //printf("\n\n alias: %s\n\n",aliasName.Data());

     for(int timeStamp=0;timeStamp<=2000;timeStamp+=500){
       Float_t simVal = (Float_t) (random.Gaus()*600.+1800.);
       AliDCSValue* dcsVal = new AliDCSValue(simVal, timeStamp);
       //printf("%s\n",dcsVal->ToString());
       valueSet->Add(dcsVal);
     }
     aliasMap->Add(new TObjString(aliasName), valueSet);
//   }
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

  AliCDBManager *manager = AliCDBManager::Instance();
  AliCDBStorage *sto = manager->GetStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/");
  AliCDBId id("ZDC/DCS/Data",0,999999999);
  AliCDBEntry *entry = sto->Get("ZDC/DCS/Data", 0);
  if(!entry) printf("TestZDCPreprocessor.C -> ERROR! No entry found as DCS Map! \n");
  return dynamic_cast<TMap*> (entry->GetObject());
}

void WriteDCSAliasMap()
{
  // This writes the output from CreateDCSAliasMap to a CDB file

  TMap* dcsAliasMap = CreateDCSAliasMap();

  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Chiara Oppedisano");
	metaData.SetComment("Test object for TestZDCPreprocessor.C");

  AliCDBId id("ZDC/DCS/Data", 0, 999999999);

  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB");

  AliCDBManager::Instance()->Put(dcsAliasMap, id, &metaData);
}
