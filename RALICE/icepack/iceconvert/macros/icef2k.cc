//////////////////////////////////////////////
// Macro to test new IceF2k conversion class
//
// To run this macro in batch, just do
//
// root -b -q 
//
// NvE 11-mar-2005 Utrecht University
//////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");
 gSystem->Load("iceconvert");

 // Output file for the event structures
 TFile* ofile=new TFile("events.root","RECREATE","F2K data in IceEvent structure");
 TTree* otree=new TTree("T","Data of an Amanda run");

 // Limit the number of entries for testing
 Int_t nentries=5;

 // Print frequency to produce a short summary print every printfreq events
 Int_t printfreq=1;

 // Split level for the output structures
 Int_t split=2;

 // Buffer size for the output structures
 Int_t bsize=32000;

 IceF2k q("run7825.f2k",split,bsize);
 q.Loop(otree,nentries,printfreq);

 // Select various objects to be added to the output file

 AliObjMatrix* omdb=q.GetOMdbase();
 if (omdb) omdb->Write();

 AliDevice* fitdefs=q.GetFitdefs();
 if (fitdefs) fitdefs->Write();

 TDatabasePDG* pdg=q.GetPDG();
 if (pdg) pdg->Write();
 
 // Close output file
 ofile->Write();
 ofile->Close();
}
