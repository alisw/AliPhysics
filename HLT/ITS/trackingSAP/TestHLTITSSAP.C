void TestHLTITSSAP() {

  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libhijing");
  gSystem->Load("libTHijing");
  gSystem->Load("libgeant321");

  if (gSystem->Getenv("EVENT"))
   nev = atoi(gSystem->Getenv("EVENT")) ;   
  
  AliSimulation simulator;
  simulator.SetRunGeneration(0);
  simulator.SetRunSimulation(0);
  simulator.SetMakeDigits("");
  simulator.SetMakeSDigits("");
  simulator.SetMakeDigitsFromHits("");
  //simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetSpecificStorage("ITS/Align/Data",
			       "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  
  simulator.SetRunQA(":") ; 

  AliHLTConfigurationHandler::CreateHandler();
  AliHLTConfigurationHandler *handler = AliHLTConfigurationHandler::Instance();

  handler->CreateConfiguration("ITS-SAPtracker1","ITSSAPTracker","DigitClusterFinder ITS-SPD-vertexer",
			       "");


  simulator.SetRunHLT("chains=ITS-SAPtracker1");
  
  TStopwatch timer;
  timer.Start();
  simulator.Run();
  timer.Stop();
  timer.Print();
}
