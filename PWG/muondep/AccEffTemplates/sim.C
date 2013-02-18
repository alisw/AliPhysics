void sim(Int_t nev=100) 
{  
  AliSimulation simulator;
  simulator.SetTriggerConfig("MUON");
  simulator.SetRunQA("MUON:ALL");
  simulator.SetRunHLT("");

  simulator.SetMakeSDigits("MUON");
  simulator.SetMakeDigits("MUON");// ITS"); // ITS needed to propagate the simulated vertex
//  simulator.SetMakeDigitsFromHits("ITS"); // ITS needed to propagate the simulated vertex
  
  // raw OCDB
//  simulator.SetDefaultStorage("alien://folder=/alice/data/2011/OCDB?cacheFold=/local/cdb");
  simulator.SetDefaultStorage(VAR_OCDB_PATH);
  
  if ( VAR_OCDB_SNAPSHOT )
  {
    simulator.SetCDBSnapshotMode("OCDB_sim.root");
  }
  // MUON Tracker
  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
  
  // Vertex and Mag.field from OCDB
//  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();
  
  // The rest
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
