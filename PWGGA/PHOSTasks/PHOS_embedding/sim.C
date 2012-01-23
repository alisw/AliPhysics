void sim(Int_t nev=20) {

  if (gSystem->Getenv("SIM_EVENTS"))
    nev = atoi(gSystem->Getenv("SIM_EVENTS"));

  printf("GENERATE << %d >> events \n",nev);


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
