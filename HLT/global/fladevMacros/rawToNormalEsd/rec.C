void rec(char* input = "raw.root") {
  AliReconstruction reco;

	reco.SetInput(input);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetFractionFriends(1.0);

  reco.SetRunReconstruction("ITS TPC HLT");

  
//reco.SetDefaultStorage("local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/data/2010/OCDB");
  reco.SetDefaultStorage("local:///home/steffen/ALICE/ocdb10");

  // -------------------------------------------------------

  // ITS (2 objects)
  reco.SetSpecificStorage("ITS/Align/Data",          "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("ITS/Calib/SPDSparseDead", "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");  

  // MUON (1 object)
  reco.SetSpecificStorage("MUON/Align/Data",         "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");

  // TPC (24 objects)
  reco.SetSpecificStorage("TPC/Align/Data",          "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/PadTime0",      "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/ClusterParam",  "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Pedestals",     "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Parameters",    "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/ExB",           "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Mapping",       "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/PadNoise",      "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/PadGainFactor", "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Temperature",   "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/RecoParam",     "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/TimeGain",      "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/AltroConfig",   "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/CE",            "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Pulser",        "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Distortion",    "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Ref",           "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Raw",           "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/QA",            "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/TimeDrift",     "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Goofie",        "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/HighVoltage",   "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/LaserTracks",   "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("TPC/Calib/Correction",    "local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual"); 

  /*
  // TPC
  //  reco.SetSpecificStorage("TPC/Calib/GainFactorDedx", "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/PadTime0",       "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/Pedestals",      "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/Mapping",        "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/PadNoise",       "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/PadGainFactor",  "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/CE",             "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/Pulser",         "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/Distortion",     "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/Ref",            "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/Raw",            "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/QA",             "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/Goofie",         "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/HighVoltage",    "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  //  reco.SetSpecificStorage("TPC/Calib/LaserTracks",    "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");

  // maybe from RAW
  reco.SetSpecificStorage("TPC/Calib/AltroConfig",    "local:///lustre/alice/alien/alice/simulation/2008/v4-15-Release/Ideal/");
  */

  // -------------------------------------------------------
  
  reco.SetRunPlaneEff(kTRUE);
  //  reco.SetRecoParam("ZDC",AliZDCRecoParamPbPb::GetHighFluxParam(2760));

  reco.SetRunQA(":") ;
  
  // -------------------------------------------------------

Bool_t useHLT= kTRUE;

  if (useHLT)  
    reco.SetOption("TPC", "useHLT");
  else
    reco.SetOption("TPC", "useRAW");

  // -------------------------------------------------------

  TStopwatch timer;
  timer.Start();

  reco.Run();

  timer.Stop();
  timer.Print();
}
