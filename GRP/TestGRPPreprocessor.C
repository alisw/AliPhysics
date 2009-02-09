
// This macro runs the GRP test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle
  
#include <iostream>
#include <fstream>
using namespace std;

void TestGRPPreprocessor()
{
  // load library
//  gSystem->Load("libSTEER.so");        // needed for AliGRPPreprocessor
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");

  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  Int_t kRun = 7;
  AliTestShuttle* shuttle = new AliTestShuttle(kRun, 1, 10);

  // by default they are set to $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB and TestReference
  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestReference");

  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());


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
  //shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "Period_test000.tag.root", "GDC0", "runTags01.root");
  //shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "Period_test001.tag.root", "GDC1", "runTags02.root");
  //shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "Period_test002.tag.root", "GDC2", "runTags03.root");

  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021390_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root", "GDC0", "$ALICE_ROOT/GRP/ShuttleInput/run000021390_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root");
  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021390_GRP_gdc1_Period_TDSMtest.Seq_0.tag.root", "GDC1", "$ALICE_ROOT/GRP/ShuttleInput/run000021390_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root");
  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021391_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root", "GDC0", "$ALICE_ROOT/GRP/ShuttleInput/run000021391_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root");
  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021391_GRP_gdc1_Period_TDSMtest.Seq_0.tag.root", "GDC1", "$ALICE_ROOT/GRP/ShuttleInput/run000021391_GRP_gdc1_Period_TDSMtest.Seq_0.tag.root");
       
  shuttle->AddInputFile(AliShuttleInterface::kDCS, "GRP", "CTP_runconfig", "", gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/p-p.cfg"));
  shuttle->AddInputFile(AliShuttleInterface::kDCS, "GRP", "CTP_xcounters", "DCS FXS", gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/xcounters.txt"));
  
  Char_t * filename = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/p-p.cfg");
  ifstream is;
  is.open(filename);
  is.seekg(0,ios::end);
  int length = is.tellg();
  const char *buffer = new char[length];
  is.seekg(0,ios::beg);
  is.read(buffer,length);
  is.close();

  shuttle->SetInputTriggerConfiguration(buffer);
  
  // The shuttle can read run type stored in the DAQ logbook.
  // To test it, we must provide the run type manually. They will be retrieved in the preprocessor
  // using GetRunType function.
  shuttle->SetInputRunType("PHYSICS");

  // The shuttle can read run parameters stored in the DAQ run logbook.
  // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
  // using GetRunParameter function.
  shuttle->AddInputRunParameter("time_start", "1233213.22");
  shuttle->AddInputRunParameter("time_end",   "1345645.22");
  shuttle->AddInputRunParameter("beamEnergy", "1400.");
  shuttle->AddInputRunParameter("beamType",    "p-p");
  shuttle->AddInputRunParameter("numberOfDetectors", "5");
  shuttle->AddInputRunParameter("detectorMask", "34555");
  shuttle->AddInputRunParameter("LHCperiod",    "LHC08b");

//  shuttle->AddInputRunParameter("totalEvents", "30000");
//  shuttle->AddInputRunParameter("NumberOfGDCs", "15");

  // TODO(5) NOT NEEDED
  //
  // This is for preprocessor that require data from HLT.
  // Since HLT may be switched off, the preprocessor should first query the Run logbook where
  // the HLT status is stored. SHUTTLE implements a shortcut function (GetHLTStatus) that returns
  // a bool directly. 1 = HLT ON, 0 = HLT OFF
  //

  Bool_t hltStatus = kTRUE;
  shuttle->SetInputHLTStatus(hltStatus);


  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  AliPreprocessor* test = new AliGRPPreprocessor(shuttle);

  // Test the preprocessor
  shuttle->Process();

  printf("\n\n");

  // In the preprocessor AliShuttleInterface::Store should be called to put the final
  // data to the CDB. To check if all went fine have a look at the files produced in
  // $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/SHUTTLE/Data
  //
  // Check the file which should have been created
  AliCDBEntry* chkEntry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get("GRP/GRP/Data", kRun);
  if (!chkEntry) {
    printf("The file is not there. Something went wrong.\n");
    return;
  }
  chkEntry->PrintId();
  chkEntry->GetObject()->Print();
  printf("\n\n");

  AliDCSSensor* sen = (AliDCSSensor*)(((TMap*)chkEntry->GetObject())->GetValue("fP2Pressure"));
  if(sen)
    sen->GetFit()->MakeGraph(1.5, 9.5, 15)->Draw();

//  AliTestDataDCS* output = dynamic_cast<AliTestDataDCS*> (chkEntry->GetObject());
  // If everything went fine, draw the result
//  if (output) {
//    output->Print();
//  }


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
  
  const Int_t fgknDCSDP = 10;
  const char* fgkDCSDataPoints[AliGRPPreprocessor::fgknDCSDP] = {
                   "LHCState",
                   "L3Polarity",
                   "DipolePolarity",
                   "LHCLuminosity",
                   "BeamIntensity",
                   "L3Current",
                   "DipoleCurrent",
                   "CavernTemperature",
                   "CavernAtmosPressure",
                   "SurfaceAtmosPressure"
                 };

  TMap* aliasMap;
  TObjArray* valueSet;
  AliDCSValue* dcsVal;
  
  aliasMap = new TMap;
  aliasMap->SetOwner(1);
  
  // LHCState
  valueSet = new TObjArray;
  valueSet->SetOwner(1);
  dcsVal = new AliDCSValue( 'F', 0 );
  valueSet->Add(dcsVal);
//  aliasMap->Add( new TObjString(fgkDCSDataPoints[0]), valueSet );

  // L3Polarity
  valueSet = new TObjArray;
  valueSet->SetOwner(1);
  dcsVal = new AliDCSValue( kTRUE, 0 );
  valueSet->Add(dcsVal);
  aliasMap->Add( new TObjString(fgkDCSDataPoints[1]), valueSet );
  
  // DipolePolarity
  valueSet = new TObjArray;
  valueSet->SetOwner(1);
  dcsVal = new AliDCSValue( kTRUE, 0 );
  valueSet->Add(dcsVal);
  aliasMap->Add( new TObjString(fgkDCSDataPoints[2]), valueSet );
  
  TRandom random;

  //  for( int nAlias=3; nAlias<fgknDCSDP-1; nAlias++)  {
  for( int nAlias=3; nAlias<fgknDCSDP; nAlias++)  {
    valueSet = new TObjArray;
    valueSet->SetOwner(1);

    for (int timeStamp=0; timeStamp<100; timeStamp++) {
      dcsVal = new AliDCSValue((Float_t) (random.Gaus()+5*nAlias), timeStamp);
      //printf("%s\n",dcsVal->ToString().Data());
      valueSet->Add(dcsVal);
    }
    aliasMap->Add( new TObjString( fgkDCSDataPoints[nAlias]), valueSet );
  }

  return aliasMap;
}

/*
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
*/


void WriteDCSAliasMap()
{
  // This writes the output from CreateDCSAliasMap to a CDB file

  TMap* dcsAliasMap = CreateDCSAliasMap();

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Ernesto Lopez Torres");
  metaData.SetComment("Test object for TestGRPPreprocessor.C");

  AliCDBId id("GRP/Data", 0, 0);

  // look into AliTestShuttle's CDB main folder

  AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Put(dcsAliasMap, id, &metaData);
}
