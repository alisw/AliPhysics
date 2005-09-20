//////////////////////////////////////////////////////////
// Example macro to demonstrate the usage of the IceXtalk
// processor for cross talk correction as a subtask of
// an F2K conversion job. The macro also shows how one
// can interactively invoke one or more subtasks
// (i.e. EvtAna.cxx) to be executed.
// The latter is very convenient in developing/testing
// new reconstruction/analysis algorithms.
//
// To run this macro in batch, just do
//
// root -b -q icextalk.cc
//
// For more details see the docs of class IceXtalk
//
// NvE 20-sep-2005 Utrecht University
//////////////////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");
 gSystem->Load("iceconvert");

 // Interactively compile and load the EvtAna.cxx code
 gROOT->LoadMacro("EvtAna.cxx+");

 // The database loading job (needed for the Xtalk constants)
 IceCal2Root cal("IceCal2Root","Calibration format conversion");
 cal.SetAmacalibFile("amacalib_amanda2_2003.txt");
 cal.ExecuteJob();

 AliObjMatrix* omdb=cal.GetOMdbase();

 // The Xtalk correction processor task
 IceXtalk xtalk("IceXtalk","Cross talk correction");
 xtalk.SetOMdbase(omdb);
 xtalk.SetMinProb(0.5);
 xtalk.SetXtalkPE(1);

 // The event analysis task
 EvtAna evtana("evtana","Event analysis");

 // The F2K event data processing job
 IceF2k q("IceF2k","Processing of the F2K event data");
 q.SetMaxEvents(10);
 q.SetPrintFreq(0);
 q.SetInputFile("real-reco.f2k");

 // Add the Xtalk and EvtAna processors as subtasks to the F2K job
 q.Add(&xtalk);
 q.Add(&evtana);

 // Perform the conversion
 q.ExecuteJob();
}
