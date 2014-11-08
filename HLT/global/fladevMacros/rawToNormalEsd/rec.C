void rec(char* input = "raw.root", char* ocdb="local://$OCDB10") {
  AliReconstruction reco;

	reco.SetInput(input);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetFractionFriends(1.0);

  reco.SetRunReconstruction("ITS TPC HLT");

  
//reco.SetDefaultStorage("local:///home/steffen/gsi/cvmfs/alice.gsi.de/alice/data/2010/OCDB");
  reco.SetDefaultStorage(ocdb);

  // -------------------------------------------------------

  // ITS (2 objects)
  reco.SetSpecificStorage("ITS/Align/Data",          "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("ITS/Calib/SPDSparseDead", "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");  

  // MUON (1 object)
  reco.SetSpecificStorage("MUON/Align/Data",         "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");

  // TPC (24 objects)
  reco.SetSpecificStorage("TPC/Align/Data",          "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/PadTime0",      "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/ClusterParam",  "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Pedestals",     "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Parameters",    "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/ExB",           "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Mapping",       "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/PadNoise",      "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/PadGainFactor", "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Temperature",   "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/RecoParam",     "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/TimeGain",      "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/AltroConfig",   "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/CE",            "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Pulser",        "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Distortion",    "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Ref",           "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Raw",           "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/QA",            "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/TimeDrift",     "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Goofie",        "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/HighVoltage",   "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/LaserTracks",   "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual");
  reco.SetSpecificStorage("TPC/Calib/Correction",    "local:///cvmfs/alice-ocdb.cern.ch/calibration/MC/Residual"); 

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
