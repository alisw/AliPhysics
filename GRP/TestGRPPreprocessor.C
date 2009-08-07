
// This macro runs the GRP test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle
  
// Taking as input runtype and errorLevel
// errorLevel used to simulate errors:
// 0 --> no error
// 1 --> DAQ logbook error
// 2 --> DAQ FXS error
// 3 --> DAQ logbook_trigger_config erro
// 4 --> DCS FXS error
// 5 --> DCS DPs error
// 6 --> Missing beamEnergy
// 7 --> null buffer for Trigger Config

// Need to include dummy files in TestShuttle/TestCDB for CTP Configuration and Scalers 
// (see macro $ALICE_ROOT/GRP/MakeCTPDummyEntries.C)

// Modified by C. Zampolli 


#include <iostream>
#include <fstream>
using namespace std;

void TestGRPPreprocessor(const char* runtype="PHYSICS", TString partition="ALICE", TString detector="", TString beamType = "p-p", Int_t errorLevel=0)
{
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");

  Int_t kRun = 7;
  AliTestShuttle* shuttle = new AliTestShuttle(kRun, 1, 10);

  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestReference");

  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

  // setting runtype
  shuttle->SetInputRunType(runtype);

  // simulating DCS DPs
  TMap* dcsAliasMap = CreateDCSAliasMap(errorLevel);
  shuttle->SetDCSInput(dcsAliasMap);

  // simulating input from DAQ FXS
  if (errorLevel != 2){
	  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021390_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root", "GDC0", "$ALICE_ROOT/GRP/ShuttleInput/run000021390_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root");
	  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021390_GRP_gdc1_Period_TDSMtest.Seq_0.tag.root", "GDC1", "$ALICE_ROOT/GRP/ShuttleInput/run000021390_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root");
	  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021391_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root", "GDC0", "$ALICE_ROOT/GRP/ShuttleInput/run000021391_GRP_gdc0_Period_TDSMtest.Seq_0.tag.root");
	  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "GRP", "run000021391_GRP_gdc1_Period_TDSMtest.Seq_0.tag.root", "GDC1", "$ALICE_ROOT/GRP/ShuttleInput/run000021391_GRP_gdc1_Period_TDSMtest.Seq_0.tag.root");
	  
  }

  // simulating input from DCS FXS
  if (errorLevel != 4 && !partition.IsNull() && detector.IsNull()){
	  shuttle->AddInputFile(AliShuttleInterface::kDCS, "GRP", "CTP_xcounters", "", gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/xcounters.txt"));
  }

  Char_t * filename = gSystem->ExpandPathName("$ALICE_ROOT/GRP/CTP/p-p.cfg");
  ifstream is;
  is.open(filename);
  is.seekg(0,ios::end);
  int length = is.tellg();
  const char *buffer = new char[length];
  is.seekg(0,ios::beg);
  is.read(buffer,length);
  is.close();
  const char *emptybuffer = NULL;

  // simulating input from DAQ logbook_trigger_config
  if {
	  (errorLevel != 3 errorLevel != 7 && !partition.IsNull() && detector.IsNull()) shuttle->SetInputTriggerConfiguration(buffer);
  }
  else if (errorLevel == 7) {
	  shuttle->SetInputTriggerConfiguration(emptybuffer);
  }
  
  // simulating input from DAQ logbook
  if (errorLevel != 1){
	shuttle->AddInputRunParameter("DAQ_time_start", "1233213.22");
  }
  if (errorLevel != 6){
	shuttle->AddInputRunParameter("beamEnergy", "1400.");
  }

  shuttle->AddInputRunParameter("DAQ_time_end",   "1345645.22");
  shuttle->AddInputRunParameter("beamType",    beamType);
  shuttle->AddInputRunParameter("numberOfDetectors", "5");
  shuttle->AddInputRunParameter("detectorMask", "34555");
  shuttle->AddInputRunParameter("LHCperiod",    "LHC08b");
  shuttle->AddInputRunParameter("partition",partition);
  shuttle->AddInputRunParameter("detector",detector);

  // simulating HLT
  Bool_t hltStatus = kTRUE;
  shuttle->SetInputHLTStatus(hltStatus);

  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  AliPreprocessor* test = new AliGRPPreprocessor(shuttle);

  // Test the preprocessor
  shuttle->Process();

  printf("\n\n");

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

}

TMap* CreateDCSAliasMap(Int_t errorLevel)
{
  // Creates a DCS structure
  // The structure is the following:
  // TMap (key --> value)
  // <DCSAlias> --> <valueList>
  // <DCSAlias> is a string
  // <valueList> is a TObjArray of AliDCSValue
  // An AliDCSValue consists of timestamp and a value in form of a AliSimpleValue
  
  const Int_t fgknDCSDP = 50;
  const char* fgkDCSDataPoints[AliGRPPreprocessor::fgknDCSDP] = {
                   "LHCState",              // missing in DCS
                   "L3Polarity",
                   "DipolePolarity",
                   "LHCLuminosity",         // missing in DCS
                   "BeamIntensity",         // missing in DCS
                   "L3Current",
                   "DipoleCurrent",
		   "L3_BSF17_H1",
		   "L3_BSF17_H2",
		   "L3_BSF17_H3",
		   "L3_BSF17_Temperature",
		   "L3_BSF4_H1",
		   "L3_BSF4_H2",
		   "L3_BSF4_H3",
		   "L3_BSF4_Temperature",
		   "L3_BKF17_H1",
		   "L3_BKF17_H2",
		   "L3_BKF17_H3",
		   "L3_BKF17_Temperature",
		   "L3_BKF4_H1",
		   "L3_BKF4_H2",
		   "L3_BKF4_H3",
		   "L3_BKF4_Temperature",
		   "L3_BSF13_H1",
		   "L3_BSF13_H2",
		   "L3_BSF13_H3",
		   "L3_BSF13_Temperature",
		   "L3_BSF8_H1",
		   "L3_BSF8_H2",
		   "L3_BSF8_H3",
		   "L3_BSF8_Temperature",
		   "L3_BKF13_H1",
		   "L3_BKF13_H2",
		   "L3_BKF13_H3",
		   "L3_BKF13_Temperature",
		   "L3_BKF8_H1",
		   "L3_BKF8_H2",
		   "L3_BKF8_H3",
		   "L3_BKF8_Temperature",
		   "Dipole_Inside_H1",
		   "Dipole_Inside_H2",
		   "Dipole_Inside_H3",
		   "Dipole_Inside_Temperature",
		   "Dipole_Outside_H1",
		   "Dipole_Outside_H2",
		   "Dipole_Outside_H3",
		   "Dipole_Outside_Temperature",
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
  dcsVal = new AliDCSValue( 'F', 2 );
  valueSet->Add(dcsVal);
  aliasMap->Add( new TObjString(fgkDCSDataPoints[0]), valueSet );

  // L3Polarity
  valueSet = new TObjArray;
  valueSet->SetOwner(1);
  dcsVal = new AliDCSValue( kTRUE, 2 );
  valueSet->Add(dcsVal);
  aliasMap->Add( new TObjString(fgkDCSDataPoints[1]), valueSet );
  
  // DipolePolarity
  valueSet = new TObjArray;
  valueSet->SetOwner(1);
  dcsVal = new AliDCSValue( kTRUE, 2 );
  valueSet->Add(dcsVal);
  aliasMap->Add( new TObjString(fgkDCSDataPoints[2]), valueSet );
  
  TRandom random;

  Int_t maxDPindex = 0;
  if (errorLevel != 5) {
	  maxDPindex = fgknDCSDP;
  }
  else {
	  maxDPindex = 3;  // simulating only a few DP in case errorLevel=5
  }

  for( int nAlias=3; nAlias<maxDPindex; nAlias++)  {
	  if (nAlias>=7 && nAlias < 47) continue; 
    valueSet = new TObjArray;
    valueSet->SetOwner(1);

    for (int timeStamp=0; timeStamp<10; timeStamp++) {
      dcsVal = new AliDCSValue((Float_t) (timeStamp+1+10*nAlias), timeStamp+1);
      valueSet->Add(dcsVal);
    }
    aliasMap->Add( new TObjString( fgkDCSDataPoints[nAlias]), valueSet );
  }

  // Hall Probes
  TString probe1[3] = {"L3_BSF","L3_BKF","Dipole_"};
  TString probe2[6] = {"17_","4_","13_","8_","Inside_","Outside_"};
  TString probe3[4] = {"H1","H2","H3","Temperature"};
  Int_t hp = 0;

  for (Int_t i=0;i<3;i++){
	  for (Int_t j=0;j<6;j++){
		  for (Int_t k=0;k<4;k++){
			  TString dpAlias = probe1[i]+probe2[j]+probe3[k];
			  valueSet = new TObjArray;
			  valueSet->SetOwner(1);
			  for (int timeStamp=0; timeStamp<10; timeStamp++) {
				  dcsVal = new AliDCSValue((Float_t) (timeStamp+1+10*hp), timeStamp+1);
				  valueSet->Add(dcsVal);
				  //cout << " hall probe = " << dpAlias << " with value = " << dcsVal->GetFloat() << endl;
			  }
			  aliasMap->Add( new TObjString(dpAlias), valueSet );
			  hp++;
		  }
	  }
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

  TMap* dcsAliasMap = CreateDCSAliasMap(Int_t errorLevel);

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Ernesto Lopez Torres");
  metaData.SetComment("Test object for TestGRPPreprocessor.C");

  AliCDBId id("GRP/Data", 0, 0);

  // look into AliTestShuttle's CDB main folder

  AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Put(dcsAliasMap, id, &metaData);
}
