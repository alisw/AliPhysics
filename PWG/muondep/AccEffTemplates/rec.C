void rec()
{
  AliReconstruction reco;
  reco.SetRunQA("MUON:ALL");
  
  reco.SetCleanESD(kFALSE);
  reco.SetStopOnError(kFALSE);
  
  reco.SetDefaultStorage(VAR_OCDB_PATH);
  if ( VAR_OCDB_SNAPSHOT ) reco.SetCDBSnapshotMode("OCDB_rec.root");

  if ( VAR_USE_ITS_RECO )
  {
    reco.SetRunReconstruction("VZERO T0 MUON ITS FMD");
  }
  else
  {
    reco.SetRunReconstruction("MUON");
  }


  // GRP from local OCDB
  reco.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));

  if ( ! VAR_USE_RAW_ALIGN )
  {
    // MUON Tracker Residual Alignment
    reco.SetSpecificStorage("MUON/Align/Data",VAR_REC_ALIGNDATA);
  }
  
  if ( VAR_USE_ITS_RECO )
  {
    // ITS
    reco.SetRunPlaneEff(kTRUE);
    reco.SetUseTrackingErrorsForAlignment("ITS");

    reco.SetSpecificStorage("ITS/Align/Data", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
    reco.SetSpecificStorage("ITS/Calib/SPDSparseDead", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");

  }
  else if ( VAR_USE_MC_VERTEX )
  {
    // ITS (use for MC vertex)
    AliITSRecoParam *itsRecoParam = AliITSRecoParam::GetLowFluxParam();
    itsRecoParam->SetVertexerSmearMC();
    itsRecoParam->ReconstructOnlySPD();
    reco.SetRecoParam("ITS",itsRecoParam);
  }

  reco.Run();
}
