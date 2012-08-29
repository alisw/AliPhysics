void sim(Int_t nev=50) {

  if (gSystem->Getenv("SIM_EVENTS"))
    nev = atoi(gSystem->Getenv("SIM_EVENTS"));

  printf("GENERATE << %d >> events \n",nev);

  gROOT->LoadMacro("IpPion.C++") ;


  AliSimulation simulator;
  simulator.SetMakeSDigits("PHOS");
  simulator.SetMakeDigits("PHOS");
//
// Ideal OCDB
//  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetDefaultStorage("local://./OCDB");
//  simulator.SetSpecificStorage("GRP/GRP/Data",
//                               Form("local://%s",gSystem->pwd()));

//  simulator.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");

 //simulator.SetSpecificStorage("GRP/Calib/MeanVertexSPD", "alien://folder=/alice/data/2010/OCDB");

  //PHOS bad map from RAW OCDB
  simulator.SetSpecificStorage("PHOS/*/*/","local://./OCDB");
//  simulator.SetSpecificStorage("PHOS/Calib/EmcBadChannels/","local://./OCDB");
//  simulator.SetSpecificStorage("PHOS/Calib/EmcGainPedestals/","local://./OCDB");

  simulator.SetRunHLT("");
//

  simulator.SetSpecificStorage("GRP/GRP/Data", "alien://Folder=/alice/data/2010/OCDB");

  AliPHOSSimParam *simParam =  AliPHOSSimParam::GetInstance() ;
  simParam->SetAPDNoise(0.000001) ;
  simParam->SetCellNonLineairyA(0.001) ;
//  simParam->SetCellNonLineairyA(0.1) ; //Default
  simParam->SetCellNonLineairyB(0.2) ;
//  simParam->SetCellNonLineairyC(0.989) ; //Jan4
//  simParam->SetCellNonLineairyC(0.995) ; //Jan5 - 2GeV
  simParam->SetCellNonLineairyC(1.031) ; //no NL

// Vertex and Mag.field from OCDB

//  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();
  simulator.SetRunQA(":") ;

//
// The rest

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
