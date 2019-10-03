void sim(Int_t nevents=100)
{
  
  // User-supplied hadron spectrum
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/EVGEN");
  gROOT->LoadMacro("AliGenPHOSlibPlus.h++") ;
  
  // AliCDBManager::Instance()->SetDefaultStorage("raw://");

  AliSimulation sim;
  sim.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  sim.SetSpecificStorage("GRP/GRP/Data",
			 Form("local://%s",gSystem->pwd()));
  // sim.SetRunNumber(146802);
  
  sim.SetMakeSDigits("PHOS") ;
  sim.SetMakeDigits("PHOS") ;
  
  sim.SetRunQA(":");
  sim.SetRunHLT("");
  
  sim.Run(nevents);
  
}
