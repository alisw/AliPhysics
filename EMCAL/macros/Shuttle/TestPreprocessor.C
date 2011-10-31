/* $Id: TestPreprocessor.C 21848 2007-10-29 18:07:14Z acolla $ */

// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle

void TestPreprocessor(const int physics = 1)
{
  // load library 
  // [compiled with: cd $ALICE_ROOT/SHUTTLE/TestShuttle/; make; cd -]
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");

   // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(7, 0, 1);

  // TODO if needed, change location of OCDB and Reference test folders
  // by default they are set to $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB and TestReference
  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/OCDB");
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB");

  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

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

  // For now, we only check the files from DAQ
  if (physics) {  
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "EMC", "signal", "MON0", "EMCALLED.root");
    printf("EMCALLED.root added\n");
  }
  else {
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "EMC", "pedestals", "MON0", "EMCALPED.root");
    printf("EMCALPED.root added\n");
  }

  // TODO(3)
  //
  // The shuttle can read run type stored in the DAQ logbook.
  // To test it, we must provide the run type manually. They will be retrieved in the preprocessor
  // using GetRunType function.
  if (physics) { 
    shuttle->SetInputRunType("PHYSICS"); 
    printf("RunType PHYSICS\n");
  }
  else { 
    shuttle->SetInputRunType("PEDESTAL"); 
    printf("RunType PEDESTAL\n");
  }

  // TODO(6)
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  AliPreprocessor* test = new AliEMCALPreprocessor(shuttle);
  printf("AliEMCALPreprocessor created\n");

  // Test the preprocessor
  shuttle->Process();
  printf("shuttle->Process() done\n");

  // TODO(7)
  // In the preprocessor AliShuttleInterface::Store should be called to put the final
  // data to the CDB. To check if all went fine have a look at the files produced in
  // $ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB/<detector>/SHUTTLE/Data
  //
  // Check the file which should have been created
  AliCDBEntry* chkEntry;

  if (physics) {
    chkEntry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
      ->Get("EMCAL/Calib/LED", 7);
  }
  else {
    chkEntry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainRefStorage())
      ->Get("EMCAL/Calib/Pedestals", 7);
  }

  if (!chkEntry)
  {
    printf("The file is not there. Something went wrong.\n");
    return;
  }

}
