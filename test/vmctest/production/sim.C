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

 simulator.SetSpecificStorage("GRP/Calib/MeanVertexSPD", "alien://folder=/alice/data/2010/OCDB");

//
// SDD from RAW OCDB

 simulator.SetSpecificStorage("ITS/Calib/CalibSDD","alien://Folder=/alice/data/2010/OCDB");

//
// EMCAL from RAW OCDB

  simulator.SetSpecificStorage("EMCAL/Calib/Data","alien://Folder=/alice/data/2010/OCDB");

//
// TRD from RAW OCDB

  simulator.SetSpecificStorage("TRD/Calib/ChamberStatus","alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("TRD/Calib/PadStatus","alien://folder=/alice/data/2010/OCDB");

//
// V0 from RAW OCDB

  simulator.SetSpecificStorage("VZERO/Trigger/Data","alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/RecoParam","alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/Data","alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/TimeSlewing","alien://folder=/alice/data/2010/OCDB");
  simulator.SetSpecificStorage("VZERO/Calib/TimeDelays","alien://folder=/alice/data/2010/OCDB");

//
// TOF from RAW OCDB

//  simulator.SetSpecificStorage("TOF/Calib/Status","alien://folder=/alice/data/2010/OCDB");

//
// Read GRP Data from RAW

  simulator.SetSpecificStorage("GRP/GRP/Data", "alien://Folder=/alice/data/2010/OCDB");

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
