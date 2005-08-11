//////////////////////////////////////////////////////////
// Example macro to demonstrate the IceCal2Root facility
// and also shows how one can interactively invoke a 
// subtask (i.e. Analyse.cxx) to be executed.
// The latter is very convenient in developing/testing
// new reconstruction/analysis algorithms.
//
// To run this macro in batch, just do
//
// root -b -q icecal2root.cc
//
// For more details see the docs of class IceCal2Root
//
// NvE 09-aug-2005 Utrecht University
//////////////////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");
 gSystem->Load("iceconvert");

 // Interactively compile and load the Analyse.cxx code
 gROOT->LoadMacro("Analyse.cxx+");

 IceCal2Root q("IceCal2Root","Calibration format conversion");

 q.SetAmacalibFile("amacalib_amanda2_2003.txt");
 q.SetOutputFile("cal2003.root");

 // Add Analyse as a sub-task to the IceCal2Root job
 Analyse ana("Analyse","OMDB analyser");
 q.Add(&ana);

 q.ExecuteJob();
}