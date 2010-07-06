// $Id$

void rec() {

  AliReconstruction reco;

//
// switch off cleanESD, write ESDfriends and Alignment data

  reco.SetCleanESD(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

//
// ITS Efficiency and tracking errors

  reco.SetRunPlaneEff(kTRUE);
  reco.SetUseTrackingErrorsForAlignment("ITS");

//
// Residual OCDB

  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");


//
// Vertex from RAW OCDB

  reco.SetSpecificStorage("GRP/Calib/MeanVertexTPC","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("GRP/Calib/MeanVertex","alien://folder=/alice/data/2010/OCDB");


//
// EMCAL from RAW OCDB

  reco.SetSpecificStorage("EMCAL/Calib/Data","alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("EMCAL/Calib/Pedestals","alien://Folder=/alice/data/2010/OCDB");

//
// PHOS from RAW OCDB

  reco.SetSpecificStorage("PHOS/Calib/EmcBadChannels","alien://Folder=/alice/data/2010/OCDB");

//
// SPD and SDD from RAW OCDB

  reco.SetSpecificStorage("ITS/Calib/SPDDead","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("TRIGGER/SPD/PITConditions","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("ITS/Calib/SPDNoise","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("ITS/Calib/CalibSDD","alien://Folder=/alice/data/2010/OCDB");

//
// TRD from RAW OCDB

  reco.SetSpecificStorage("TRD/Calib/ChamberStatus","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("TRD/Calib/PadStatus","alien://folder=/alice/data/2010/OCDB");

//
// TPC from RAW OCDB

  reco.SetSpecificStorage("TPC/Calib/PadGainFactor","alien://folder=/alice/data/2010/OCDB");

//
// V0 from RAW OCDB

  reco.SetSpecificStorage("VZERO/Trigger/Data","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/RecoParam","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/Data","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/TimeSlewing","alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/TimeDelays","alien://folder=/alice/data/2010/OCDB");

//
// TOF from RAW OCDB

//  reco.SetSpecificStorage("TOF/Calib/Status","alien://folder=/alice/data/2010/OCDB");

//
// QA off

  reco.SetRunQA(":");

//
// The rest
 
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
