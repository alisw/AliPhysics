#define RecAna_cxx
#include "RecAna.h"

void RecAna::Loop()
{
  //   In a Root session, you can do:
  //      Root > .L RecAna.C
  //      Root > RecAna t(filename)
  //      Root > t.GetEvent(evt); // Fill t data members with event number evt
  //      Root > t.Show();       // Show values of current event e
  //      Root > t.Show(evt);     // Read and show values of entry evt
  //      Root > t.Loop();       // Loop on all entries
  //
  
  //     This is the loop skeleton
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fTree->SetBranchStatus("*",0);  // disable all branches
  //    fTree->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fTree->GetEntry(i);  // read all branches
  //by  b_branchname->GetEntry(i); //read only this branch
  if (fTree == 0) return;
  
  Int_t nentries = Int_t(fTree->GetEntries());
  
   Int_t nbytes = 0, nb = 0;
   for (Int_t i=0; i<nentries;i++) {
      if (LoadTree(i) < 0) break;
      nb = fTree->GetEntry(i);   nbytes += nb;
   }
}
