/* $Id $ */

// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

void TestTRIGGERPreprocessor()
{
  // load library
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");

  // Defining shuttle instance (run number, start time, end time) 
  AliTestShuttle* shuttle = new AliTestShuttle(7, 0, 1);

  // Setting OCDB
  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestReference");

  // Adding test input files
  shuttle->AddInputFile(AliTestShuttle::kDCS, "TRI", "PITConditions", "", "$ALICE_ROOT/TRIGGER/ShuttleInput/pit_dumpFileSep08.txt");

  // Adding Trigger mask
  shuttle->SetInputTriggerDetectorMask("0100001");
  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

  // Setting RunType
  shuttle->SetInputRunType("PHYSICS");

  // Defining Preprocessor
  AliPreprocessor* test = new AliTRIPreprocessor(shuttle);

  // Test the preprocessor
  shuttle->Process();

}

