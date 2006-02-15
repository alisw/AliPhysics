//////////////////////////////////////////////////////////
// Example macro to demonstrate the usage of the IcePandel
// processor for Pandel fit reconstruction as a subtask
// of an F2K conversion job. The macro also shows how one
// can interactively invoke one or more subtasks
// (i.e. IceXtalk and EvtAna.cxx) to be executed.
// The latter is very convenient in testing/comparing
// new reconstruction/analysis algorithms.
//
// To run this macro in batch, just do
//
// root -b -q icepandel.cc
//
// For more details see the docs of class IcePandel
//
// NvE 09-feb-2006 Utrecht University
//////////////////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");
 gSystem->Load("iceconvert");

 // Interactively compile and load the EvtAna.cxx code
 gROOT->LoadMacro("EvtAna.cxx+");

 // The database conversion job
 IceCal2Root cal("IceCal2Root","Calibration format conversion");
 cal.SetAmacalibFile("amacalib_amanda2_2003.txt");
 cal.SetOutputFile("cal2003.root");
 cal.ExecuteJob();

 // The Pandel fitting procedure
 IcePandel pandel("IcePandel","Pandel fitting");
 pandel.UseTracks("IceDwalk");

 // The direct walk reconstruction task
 IceDwalk dwalk("IceDwalk","Direct walk reconstruction");

 // The calibration processor task
 IceCalibrate calib("IceCalibrate","Signal calibration");
 calib.SetCalibFile("cal2003.root");

 // The Xtalk correction processor task
 IceXtalk xtalk("IceXtalk","Cross talk correction");
 xtalk.SetCalibFile("cal2003.root");
 xtalk.SetMinProb(0.5);
 xtalk.SetXtalkPE(1);

 // The hit cleaning processor task
 IceCleanHits clean("IceCleanHits","Hit cleaning");

 // The event analysis task
 EvtAna evtana("evtana","Event analysis");

 // The F2K event data processing job
 IceF2k q("IceF2k","Processing of the F2K event data");
 q.SetMaxEvents(2);
 q.SetPrintFreq(0);
 q.SetInputFile("real-reco.f2k");
 q.SetOutputFile("real-reco.root");

 // Add the various processors as subtasks to the F2K job
 q.Add(&calib);
 q.Add(&xtalk);
 q.Add(&clean);
 q.Add(&dwalk);
 q.Add(&pandel);
 q.Add(&evtana);

 // Perform the conversion
 q.ExecuteJob();
}
