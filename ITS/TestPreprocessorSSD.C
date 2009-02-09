/* $Id$ */

// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

void TestPreprocessorSSD()
{
  // load library
  gSystem->Load("libTestShuttle.so");

  // initialize location of CDB
  //  AliCDBManager::Instance()->SetDefaultStorage("local://${ALICE_ROOT}/SHUTTLE/TestShuttle/TestCDB");

  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestReference");

  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime

  AliTestShuttle* shuttle = new AliTestShuttle(7, 0, 1);
  shuttle->SetInputRunType("PEDESTAL_RUN");

  shuttle->AddInputFile(AliTestShuttle::kDAQ, "SSD", "CALIBRATION", "LDC0", "ITSSSDda_LDC0.root");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "SSD", "CALIBRATION", "LDC1", "ITSSSDda_LDC1.root");
  shuttle->AddInputFile(AliTestShuttle::kDAQ, "SSD", "CALIBRATION", "LDC2", "ITSSSDda_LDC2.root");

  // TODO(3)
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  AliPreprocessor *pp = new AliITSPreprocessorSSD(shuttle);

  // Test the preprocessor
  shuttle->Process();

  
  //
  // Check the file which should have been created
  AliCDBManager::Instance()->SetDefaultStorage("local://${ALICE_ROOT}/SHUTTLE/TestShuttle/TestCDB");  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("ITS/Calib/NoiseSSD", 7);
  if (!entry)
  {
    printf("The file is not there. Something went wrong.\n");
    return;
  }
  

}

