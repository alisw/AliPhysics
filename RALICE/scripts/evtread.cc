//////////////////////////////////////////////////////////////////////////
// Test of the AliVertex, AliTrack, AliJet and AliVertex functionality.
// Various vertices are read-in from the events.root file created by
// evtwrite.cc.
// These vertices are then analysed in terms of tracks and jets
// using the RALICE facilities.
// In the writing program several identical events were created to enable
// testing of the multiple event structure on the output file.
// 
//--- NvE 28-may-1998 UU-SAP Utrecht
//////////////////////////////////////////////////////////////////////////
{
 gSystem->Load("ralice");

 // Get the input file and a Tree 
 TFile* f=new TFile("events.root");
 TTree* tree=(TTree*)f->Get("T");

 // Provide overview of Tree contents
 cout << endl;
 tree->Print();

 // Define an event (which is also the main Vertex)
 AliEvent* vmain=0;

 // Branch in the tree for the event input
 tree->SetBranchAddress("Events",&vmain);

 Int_t nen=tree->GetEntries();
 cout << endl;
 cout << " *READ* nentries : " << nen << endl; 

 for (Int_t ient=0; ient<nen; ient++)
 {
  cout << endl;
  cout << " === Going for event number : " << (ient+1) << endl;
  cout << endl;

  tree->GetEntry(ient);

  // Print the full event information
  vmain->ListAll();
 }

 // Close input file
 f->Close();
}
