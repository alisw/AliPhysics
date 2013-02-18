
void rec(int dset=0)
{
  
  AliReconstruction reco;
  reco.SetRunQA("MUON:ALL");
  
  reco.SetCleanESD(kFALSE);
  reco.SetStopOnError(kFALSE);
  
  reco.SetDefaultStorage(VAR_OCDB_PATH);

  if ( dset == 1 )
  {
    reco.SetRunReconstruction("MUON ITS");
  }
  else
  {
    reco.SetRunReconstruction("MUON");
  }
  
  if ( VAR_OCDB_SNAPSHOT )
  {
    reco.SetCDBSnapshotMode("OCDB_rec.root");
  }

  // GRP from local OCDB
  reco.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  
  // MUON Tracker Residual Alignment
  reco.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Residual");

     reco.Run();
}
