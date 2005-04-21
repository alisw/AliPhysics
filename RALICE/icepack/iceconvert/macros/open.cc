////////////////////////////////////////////////////////
// Macro to attach and open an IceCube event file
// for interactive investigation.
//
// To run this macro, just do
//
// root> .x open.cc
//
//--- NvE 18-jun-2004 Utrecht University
/////////////////////////////////////////////
{
 gSystem->Load("ralice");
 gSystem->Load("icepack");

 // Access to the input data
 TFile* f=new TFile("events.root");
 TTree* data=(TTree*)f->Get("T");

 // Provide overview of Tree contents
 cout << endl;
 data->Print();

 // Define a pointer for an event
 IceEvent* evt=0;

 // Branch in the tree for the event input
 data->SetBranchAddress("IceEvent",&evt);

 cout << endl;
 cout << " *READ* nentries : " << data->GetEntries() << endl;
 cout << endl;

 cout << " Use data->GetEntry(i) to load the i-th entry." << endl;
 cout << " The event object is called evt " << endl;
 cout << endl;
}
