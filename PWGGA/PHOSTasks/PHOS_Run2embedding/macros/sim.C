void sim(Int_t nev=50) {

  if (gSystem->Getenv("SIM_EVENTS"))
    nev = atoi(gSystem->Getenv("SIM_EVENTS"));

  printf("GENERATE << %d >> events \n",nev);

  AliSimulation simulator;
  simulator.SetMakeSDigits("PHOS");
  simulator.SetMakeDigits("PHOS");

  simulator.SetDefaultStorage("alien://Folder=/alice/data/2018/OCDB");

  simulator.SetRunHLT("");

  AliPHOSSimParam *simParam =  AliPHOSSimParam::GetInstance() ;
  simParam->SetAPDNoise(0.000001) ;
  simParam->SetCellNonLineairyA(0.001) ;
  simParam->SetCellNonLineairyB(0.2) ;
  simParam->SetCellNonLineairyC(1.031) ; //no NL

  simulator.UseMagFieldFromGRP();
  simulator.SetRunQA(":") ;

  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");        // Pythia 6.2

  simulator.Run(nev);
}
