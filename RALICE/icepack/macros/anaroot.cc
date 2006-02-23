//////////////////////////////////////////////////////////
// Example macro to demonstrate the usage of various
// processors (e.g. IceCalibrate, IceXtalk, IceCleanHits,
// IcePandel) from a single AliJob which reads ROOT
// data structures.
//
// To run this macro in batch, just do
//
// root -b -q anaroot.cc
//
// NvE 22-feb-2006 Utrecht University
//////////////////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");
 gSystem->Load("iceconvert");

 // The input chain
 TChain data("T");
 data.AddFile("run7149-mbias-sldw-nocal.root");

 Int_t nentries=(int)data.GetEntries(); 
 
 // Pointer to the event data structures
 IceEvent* evt=0;
 data.SetBranchAddress("IceEvent",&evt);

 // The database conversion job
 // In case the file cal2003.root is already available,
 // the following 4 statments may be removed.
 IceCal2Root cal("IceCal2Root","Calibration format conversion");
 cal.SetAmacalibFile("amacalib_amanda2_2003.txt");
 cal.SetOutputFile("cal2003.root");
 cal.ExecuteJob();

 // The main data processing job
 AliJob q("AliJob","Processing of the Amanda event ROOT data");

 // The calibration processor task
 IceCalibrate calib("IceCalibrate","Signal calibration");
 calib.SetCalibFile("cal2003.root");

 // The Xtalk correction processor task
 IceXtalk xtalk("IceXtalk","Cross talk correction");
 xtalk.SetCalibFile("cal2003.root");

 // The hit cleaning processor task
 IceCleanHits clean("IceCleanHits","Hit cleaning");

 // The direct walk reconstruction task
 IceDwalk dwalk("IceDwalk","Direct walk reconstruction");

 // The Pandel fitting procedure
 IcePandel pandel("IcePandel","Pandel fitting");
 pandel.UseTracks("IceDwalk");

 // Add the various processors as subtasks to the main job
 q.Add(&calib);
 q.Add(&xtalk);
 q.Add(&clean);
 q.Add(&dwalk);
 q.Add(&pandel);

 q.ListEnvironment();

 // Perform the analysis
 for (Int_t ient=0; ient<nentries; ient++)
 {
  if (evt) evt->Reset();
  data.GetEntry(ient);
  evt->SetOwner();
  q.AddObject(evt);
  q.ExecuteJob(10);
  q.RemoveObject(evt);
  
  // Further inspection of the event data if needed
  evt->Data();
 }
}
