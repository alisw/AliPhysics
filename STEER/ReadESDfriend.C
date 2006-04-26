//********************************************************************
//  Example of accessing the information stored in ESD friends
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TFile.h>
  #include <TTree.h>

  #include "AliESD.h"
  #include "AliESDfriend.h"
#endif

void ReadESDfriend(Bool_t readFriend=kTRUE) {
   TFile f("AliESDs.root");
   TTree *esdTree=(TTree*)f.Get("esdTree");
   AliESD *ev=0;
   esdTree->SetBranchAddress("ESD",&ev);

   // Attach the tree with ESD friends
   AliESDfriend *evf=0;
   if (readFriend) {
      esdTree->AddFriend("esdFriendTree","AliESDfriends.root");
      esdTree->SetBranchAddress("ESDfriend",&evf);
   }

   Int_t nev=esdTree->GetEntries();
   for (Int_t i=0; i<nev; i++) {
       cout<<"Event number: "<<i<<endl;
       esdTree->GetEntry(i);
        
       ev->SetESDfriend(evf); //Attach the friend to the ESD

    // Now the attached information can be accessed via pointer to ESD.
    // Example: indices of the TPC clusters associated with the track number 0.
       const AliESDtrack *t=ev->GetTrack(0);
       Int_t idx[AliESDfriendTrack::kMaxTPCcluster], n=t->GetTPCclusters(idx);
       cout<<"Number of TPC clusters: "<<n<<endl;
       cout<<"Index of the 7th TPC cluster: "<<idx[7]<<endl<<endl; 
   }

   delete esdTree;
}
