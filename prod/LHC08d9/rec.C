void rec() {
  //AliLog::SetGlobalLogLevel(AliLog::kError);

  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();


  AliITSRecoParam * itspar = AliITSRecoParam::GetLowFluxParam();
  itspar->SetStoreLikeSignV0s(kTRUE);
  reco.SetRecoParam("ITS",itspar);
  reco.SetRecoParam("TPC",AliTPCRecoParam::GetLowFluxParam());
  reco.SetRecoParam("TRD",AliTRDrecoParam::GetLowFluxParam());
  reco.SetRecoParam("PHOS",AliPHOSRecoParam::GetDefaultParameters());
  reco.SetRecoParam("MUON",AliMUONRecoParam::GetLowFluxParam());
  reco.SetRecoParam("EMCAL",AliEMCALRecParam::GetLowFluxParam());
  reco.SetRecoParam("GRP",AliGRPRecoParam::GetLowFluxParam());

  reco.SetOption("TRD","sl_tr_1");            // Stream Level for the tracker equal to 1

  // Only in case of Full misalignment
//   AliGRPRecoParam *grpRecoParam = AliGRPRecoParam::GetLowFluxParam();
//   grpRecoParam->SetVertexerTracksConstraintITS(kFALSE);
//   grpRecoParam->SetVertexerTracksConstraintTPC(kFALSE);
//   reco.SetRecoParam("GRP",grpRecoParam);

  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  //reco.SetSpecificStorage("GRP/GRP/Data/","alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  reco.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
  // No write access to the OCDB => local specific storage
//   reco.SetSpecificStorage("GRP/GRP/Data",
//                        Form("local://%s",gSystem->pwd()));

  //-------------------------------------------------------------------------
  // Setting the cuts for the V0 and cascade finding
  // The values of the cuts below are "reasonable" for pp events
  //-------------------------------------------------------------------------

  Double_t v0sels[]={33,    // max allowed chi2
                     0.05,  // min allowed impact parameter for the 1st daughter
                     0.05,  // min allowed impact parameter for the 2nd daughter
                     0.5,   // max allowed DCA between the daughter tracks
                     0.99,  // max allowed cosine of V0's pointing angle
                     0.2,   // min radius of the fiducial volume
                     100    // max radius of the fiducial volume
  };
  AliV0vertexer::SetDefaultCuts(v0sels);

  Double_t xisels[]={33.,   // max allowed chi2 (same as PDC07)
                     0.025, // min allowed V0 impact parameter (PDC07 was 0.05)
                     0.010, // "window" around the Lambda mass (PDC07 was 0.008)
                     0.025, // min allowed bachelor's impact parameter (PDC07 was 0.035)
                     0.2,   // max allowed DCA between the V0 and the bachelor (PDC07 was 0.1)
                     0.998, // max allowed cosine of the cascade pointing angle (PDC07 was 0.9985)
                     0.2,   // min radius of the fiducial volume (PDC07 was 0.9)
                     100    // max radius of the fiducial volume (same as PDC07)
  };
  AliCascadeVertexer::SetDefaultCuts(xisels);

  reco.SetRunQA("ALL:ALL");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}

