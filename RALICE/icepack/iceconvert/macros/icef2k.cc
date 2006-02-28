/////////////////////////////////////////////////
// Macro to test the IceF2k conversion job class
//
// To run this macro in batch, just do
//
// root -b -q icef2k.cc
//
// For more details and a more user friendly
// output file setting, see the docs of class IceF2k
//
// NvE 11-mar-2005 Utrecht University
/////////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");
 gSystem->Load("iceconvert");

 IceF2k q("IceF2k","Test of IceF2k conversion job");

 // Limit the number of entries for testing
 q.SetMaxEvents(100);

 // Print frequency to produce a short summary print every printfreq events
 q.SetPrintFreq(10);

 // The F2K input filename(s)
 q.AddInputFile("run7825.f2k");

 // Output file for the event structures
 TFile* ofile=new TFile("events.root","RECREATE","F2K data in IceEvent structure");
 q.SetOutputFile(ofile);

 ///////////////////////////////////////////////////////////////////
 // Here the user can specify his/her sub-tasks to be executed
 // on an event-by-event basis after the IceEvent structure
 // has been filled and before the data is written out.
 // Sub-tasks (i.e. a user classes derived from TTask) are entered
 // as follows :
 //
 //    MyXtalk task1("task1","Cross talk correction");
 //    MyClean task2("task2","Hit cleaning");
 //    q.Add(&task1);
 //    q.Add(&task2);
 //
 // The sub-tasks will be executed in the order as they are entered.
 ///////////////////////////////////////////////////////////////////

 // Perform the conversion
 q.ExecuteJob();

 // Select various objects to be added to the output file

 ofile->cd(); // Switch to the output file directory

 AliObjMatrix* omdb=q.GetOMdbase();
 if (omdb) omdb->Write();

 AliDevice* fitdefs=q.GetFitdefs();
 if (fitdefs) fitdefs->Write();

 TDatabasePDG* pdg=q.GetPDG();
 if (pdg) pdg->Write();

 // Flush additional objects to the output file.
 // The output file is not explicitly closed here
 // to allow interactive investigation of the data tree
 // when this macro is run in an interactive ROOT/CINT session.
 ofile->Write();
}
