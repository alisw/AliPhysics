#include <iostream.h>


//***************************************************************************

void ReadTree()
{
// This is a basic example on how the tree of reconstructed events
// should be accessed
//Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
	loadlibs();
   }
    
   
// Connect the Root Galice file containing ...
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("tree_reco.root");
   if (!file) file = new TFile("tree_reco.root");
   if (!file) {
      cout << "File tree_reco.root not found\n";
      return;
   }
    
   TFile *galice_file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!galice_file) galice_file = new TFile("galice.root");
   if (!galice_file) {
      cout << "File galice.root not found\n";
      return;
   }
     
// Get AliRun object from file or create it if not on file
   if (!gAlice) {
	gAlice = (AliRun*)galice_file->Get("gAlice");
	if (gAlice) {
         cout << "AliRun object found on file\n";
      } else {
         cout << "AliRun object not found on file !!!\n";
         return;
      }
   }
    
   TTree *tree = (TTree*)file->Get("TreeRecoEvent");
   TBranch *branch = tree->GetBranch("Event");
   AliMUONRecoEvent *event = 0;
   branch->SetAddress(&event);
   Int_t nentries = (Int_t)tree->GetEntries(); // no. of reco events on file

   for (Int_t evNumber=0; evNumber<nentries; evNumber++) {
      tree->GetEntry(evNumber);
      cout << "Event : " << event->GetNoEvent() << " with" << event->GetNoTracks() << " tracks\n";
// print reconstr. event information
//      event->EventInfo();
   }
}
