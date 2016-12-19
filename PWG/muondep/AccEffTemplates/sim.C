void sim(Int_t nev=VAR_EVENTS_PER_JOB)
{

#if defined(__CINT__)
  gSystem->Load("VAR_LHAPDF");      // Parton density functions
  if ( VAR_NEEDS_PYTHIA6 )
  {
    std::cout << "Setting up Pythia6 required env. variables" << std::endl;
    VAR_PYTHIA6_SETENV
  }
  else  gSystem->Load("libpythia6");     // Pythia 6.2 (for decayer)

  if ( VAR_NEEDS_PYTHIA8 )
  {
    std::cout << "Setting up Pythia8 required libraries and env. variables" << std::endl;
    //    gSystem->Load("libpythia8");
    //    gSystem->Load("libAliPythia8");
    VAR_PYTHIA8_SETENV
  }
#endif

  AliSimulation simulator;
  simulator.SetRunQA("MUON:ALL");
  simulator.SetRunHLT("");

  if ( VAR_OCDB_SNAPSHOT ) 
  {
      simulator.SetCDBSnapshotMode("OCDB_sim.root");
  }
  else
  {
      simulator.SetDefaultStorage(VAR_OCDB_PATH);
  }
  
  if ( VAR_USE_ITS_RECO )
  {
    simulator.SetMakeSDigits("MUON T0 VZERO FMD"); // T0 and VZERO for trigger efficiencies, FMD for diffractive studies
    simulator.SetMakeDigitsFromHits("ITS"); // ITS needed to propagate the simulated vertex
    simulator.SetMakeDigits("MUON T0 VZERO FMD");// ITS"); // ITS needed to propagate the simulated vertex
  }
  else
  {
    simulator.SetTriggerConfig("MUON");
    simulator.SetMakeSDigits("MUON");
    simulator.SetMakeDigits("MUON");
  }

  if ( VAR_PURELY_LOCAL ) {
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));
  }
  else {

    // Mag.field from OCDB
    simulator.UseMagFieldFromGRP();

    // MUON Tracker
    
    if (!VAR_OCDB_SNAPSHOT)
    {
        if ( VAR_USE_RAW_ALIGN )
        {
        simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Full");
        }
        else
        {
        simulator.SetSpecificStorage("MUON/Align/Data",VAR_SIM_ALIGNDATA);
        }
    }

    // ITS
    if ( VAR_USE_ITS_RECO )
    {
      simulator.SetSpecificStorage("ITS/Align/Data",  "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal");
    }
    if ( VAR_USE_ITS_RECO || VAR_USE_MC_VERTEX )
    {
      simulator.UseVertexFromCDB();
    }
  }
  
   // The rest
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
