// $Id$

void rec() {

  AliReconstruction reco;

//
// switch off cleanESD, write ESDfriends and Alignment data
  
  reco.SetCleanESD(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetRunQA("kFALSE:kFALSE");

//
// ITS Efficiency and tracking errors

  reco.SetRunPlaneEff(kTRUE);
  reco.SetUseTrackingErrorsForAlignment("ITS");

//
// Residual OCDB
  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");


//
// GRP from local OCDB
 reco.SetSpecificStorage("GRP/GRP/Data",
                          Form("local://%s",gSystem->pwd()));

//T0
  reco.SetSpecificStorage("T0/Calib/TimeAdjust", "alien://folder=/alice/data/2010/OCDB"); 

//
// Vertex from RAW OCDB
  reco.SetSpecificStorage("GRP/Calib/MeanVertexTPC",   "alien://folder=/alice/data/2010/OCDB"); 
  reco.SetSpecificStorage("GRP/Calib/MeanVertex",      "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("GRP/Calib/MeanVertexSPD",   "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("GRP/Calib/RecoParam",       "alien://folder=/alice/data/2010/OCDB"); 

// Clock phase from RAW OCDB 
  reco.SetSpecificStorage("GRP/Calib/LHCClockPhase",   "alien://folder=/alice/data/2010/OCDB");

// 
// EMCAL from RAW OCDB
  reco.SetSpecificStorage("EMCAL/Calib/Data",          "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("EMCAL/Calib/Pedestals",     "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("EMCAL/Calib/RecoParam",     "alien://Folder=/alice/data/2010/OCDB"); 

//
// PHOS from RAW OCDB
  reco.SetSpecificStorage("PHOS/Calib/EmcBadChannels", "alien://Folder=/alice/data/2010/OCDB"); 
  reco.SetSpecificStorage("PHOS/Calib/RecoParam",      "alien://Folder=/alice/data/2010/OCDB"); 
  reco.SetSpecificStorage("PHOS/Calib/Mapping",        "alien://Folder=/alice/data/2010/OCDB"); 

//
// ITS
  // Reco Param
  reco.SetSpecificStorage("ITS/Calib/RecoParam",       "alien://folder=/alice/data/2010/OCDB");
  // SPD	
  reco.SetSpecificStorage("ITS/Calib/SPDDead",         "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("TRIGGER/SPD/PITConditions", "alien://folder=/alice/data/2010/OCDB");
  // SDD
  reco.SetSpecificStorage("ITS/Calib/CalibSDD",        "alien://Folder=/alice/data/2010/OCDB"); 
  // SSD
  reco.SetSpecificStorage("ITS/Calib/NoiseSSD",        "alien://Folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("ITS/Calib/BadChannelsSSD",  "alien://Folder=/alice/data/2010/OCDB"); 

//
// TRD from RAW OCDB
  reco.SetSpecificStorage("TRD/Calib/ChamberStatus",   "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("TRD/Calib/PadStatus",       "alien://folder=/alice/data/2010/OCDB");

//
// TPC from RAW OCDB
  reco.SetSpecificStorage("TPC/Calib/PadGainFactor",   "alien://folder=/alice/data/2010/OCDB");

//
// V0 from RAW OCDB
  reco.SetSpecificStorage("VZERO/Trigger/Data",        "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/RecoParam",     "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/Data",          "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/TimeSlewing",   "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("VZERO/Calib/TimeDelays",    "alien://folder=/alice/data/2010/OCDB");

//
// TOF from RAW OCDB
  reco.SetSpecificStorage("TOF/Calib/Status",          "alien://folder=/alice/data/2010/OCDB"); 

//
// FMD from RAW OCDB
  reco.SetSpecificStorage("FMD/Calib/Pedestal",        "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("FMD/Calib/PulseGain",       "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("FMD/Calib/Dead",            "alien://folder=/alice/data/2010/OCDB");
  reco.SetSpecificStorage("FMD/Calib/AltroMap",        "alien://folder=/alice/data/2010/OCDB");

// MUON Trigger LuT:
  reco.SetSpecificStorage("MUON/Calib/TriggerLut",     "alien://folder=/alice/data/2010/OCDB");
// MUON Alignment:
//  reco.SetSpecificStorage("MUON/Align/Data",           "alien://folder=/alice/simulation/2008/v4-15-Release/Full");
  reco.SetSpecificStorage("MUON/Align/Data", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
// MUON  tracking RejectList:
  reco.SetSpecificStorage("MUON/Calib/RejectList",     "alien://folder=/alice/simulation/2008/v4-15-Release/Full");
// MUON reconstruction parameters:
  reco.SetSpecificStorage("MUON/Calib/RecoParam",      "alien://folder=/alice/simulation/2008/v4-15-Release/Full");
// ZDC
  reco.SetSpecificStorage("ZDC/Calib/EnergyCalib",     "alien://folder=/alice/data/2010/OCDB");
//
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
