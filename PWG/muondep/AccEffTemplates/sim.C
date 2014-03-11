void sim(Int_t nev=100) 
{
  if ( VAR_PURELY_LOCAL) {
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));
  }

  AliSimulation simulator;
  simulator.SetTriggerConfig("MUON");
  simulator.SetRunQA("MUON:ALL");
  simulator.SetRunHLT("");

  simulator.SetMakeSDigits("MUON");
  simulator.SetMakeDigits("MUON");// ITS"); // ITS needed to propagate the simulated vertex
//  simulator.SetMakeDigitsFromHits("ITS"); // ITS needed to propagate the simulated vertex

  simulator.SetDefaultStorage(VAR_OCDB_PATH);
  
  if ( VAR_OCDB_SNAPSHOT )
  {
    simulator.SetCDBSnapshotMode("OCDB_sim.root");
  }
  
  if ( ! VAR_PURELY_LOCAL ) {
    
    // MUON Tracker
    simulator.SetSpecificStorage("MUON/Align/Data",VAR_SIM_ALIGNDATA);
  
    // Mag.field from OCDB
    simulator.UseMagFieldFromGRP();

    //  simulator.UseVertexFromCDB();  
  }
  
  // The rest
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
