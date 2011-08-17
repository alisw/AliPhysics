// $Id$

void sim(Int_t nev=20) {
  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");

//
// Ideal OCDB

  simulator.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");

//
// Mean verrtex from RAW OCDB 
 simulator.SetSpecificStorage("GRP/Calib/MeanVertex",          "alien://folder=/alice/data/2010/OCDB");

// Clock phase from RAW OCDB 
 simulator.SetSpecificStorage("GRP/Calib/LHCClockPhase",       "alien://folder=/alice/data/2010/OCDB");

// ITS
//    SDD from RAW OCDB
 simulator.SetSpecificStorage("ITS/Calib/CalibSDD",            "alien://Folder=/alice/data/2010/OCDB");
//    SSD
simulator.SetSpecificStorage("ITS/Calib/NoiseSSD",             "alien://Folder=/alice/data/2010/OCDB");
simulator.SetSpecificStorage("ITS/Calib/BadChannelsSSD",       "alien://Folder=/alice/data/2010/OCDB"); 

//
// EMCAL from RAW OCDB
  simulator.SetSpecificStorage("EMCAL/Calib/Data",             "alien://Folder=/alice/data/2010/OCDB");

//
// TRD from RAW OCDB
  simulator.SetSpecificStorage("TRD/Calib/ChamberStatus",      "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("TRD/Calib/PadStatus",          "alien://folder=/alice/data/2010/OCDB");

//
// V0 from RAW OCDB
  simulator.SetSpecificStorage("VZERO/Trigger/Data",           "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/RecoParam",        "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/Data",             "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/TimeSlewing",      "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/TimeDelays",       "alien://folder=/alice/data/2010/OCDB");

//
// TOF from RAW OCDB
  simulator.SetSpecificStorage("TOF/Calib/Status",             "alien://folder=/alice/data/2010/OCDB");

//
// FMD from RAW OCDB
  simulator.SetSpecificStorage("FMD/Calib/Pedestal",           "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("FMD/Calib/PulseGain",          "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("FMD/Calib/Dead",               "alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("FMD/Calib/AltroMap",           "alien://folder=/alice/data/2010/OCDB");


//
// MUON Trigger (LuT & efficiency)
  simulator.SetSpecificStorage("MUON/Calib/TriggerLut",        "alien://folder=/alice/data/2010/OCDB");
// MUON Trigger Chamber efficiency
  simulator.SetSpecificStorage("MUON/Calib/TriggerEfficiency", "alien://folder=/alice/simulation/2008/v4-15-Release/Full");  

// ZDC
  simulator.SetSpecificStorage("ZDC/Calib/EnergyCalib",        "alien://folder=/alice/data/2010/OCDB");
//
// Read GRP Data from RAW
  simulator.SetSpecificStorage("GRP/GRP/Data",                 "alien://Folder=/alice/data/2010/OCDB");

//
// Vertex and Mag.field from OCDB

  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();

//
// The rest

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
