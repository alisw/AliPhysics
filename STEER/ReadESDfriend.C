//********************************************************************
//  Example of accessing the information stored in ESD friends
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TFile.h>
  #include <TChain.h>

  #include "AliESD.h"
  #include "AliESDfriend.h"
  #include "AliTrackPointArray.h"
#endif

void ReadESDfriend(Bool_t readFriend=kTRUE) {
   Char_t *name[]={
     //Put here the names of the ESD files to be chained
     "AliESDs.root"
   };
   Int_t n=sizeof(name)/sizeof(Char_t *); 
   TChain *esdTree=new TChain("esdTree");
   for (Int_t i=0; i<n; i++) esdTree->AddFile(name[i]);

   AliESD *ev=0;
   esdTree->SetBranchAddress("ESD",&ev);

   // Attach the branch with ESD friends
   AliESDfriend *evf=0;
   if (readFriend) {
      esdTree->SetBranchStatus("ESDfriend*",1);
      esdTree->SetBranchAddress("ESDfriend.",&evf);
   }

   Int_t nev=esdTree->GetEntries();
   for (Int_t i=0; i<nev; i++) {
       esdTree->GetEntry(i);
        
       cout<<endl<<"Event number: "<<i<<endl;
       Int_t n=ev->GetNumberOfTracks();
       cout<<"Number of tracks: "<<n<<endl;

       ev->SetESDfriend(evf); //Attach the friend to the ESD

    // Now the attached information can be accessed via pointer to ESD.
    // Example: indices of the TPC clusters associated with the track number 0.
       if (n > 0) {
          const AliESDtrack *t=ev->GetTrack(0);
          Int_t idx[AliESDfriendTrack::kMaxTPCcluster]; 
          n=t->GetTPCclusters(idx);
          cout<<"Track number 0"<<endl;
          cout<<"   Number of TPC clusters: "<<n<<endl;
          cout<<"   Index of the 7th TPC cluster: "<<idx[7]<<endl;
          UChar_t map=t->GetITSClusterMap();
          cout<<"   ITS cluster map (from SPDs to SSDs): ";
          for (Int_t i=0; i<6; i++) cout<<TESTBIT(map,i)<<' ';
          cout<<endl;

    // Example: track points associated with the track number 0.
          const AliTrackPointArray *pa=t->GetTrackPointArray();
          if (pa != 0) {
             n=pa->GetNPoints();
             const Float_t *x=pa->GetX();
             cout<<"   Number of track points: "<<n<<endl;
             if (n>7)
             cout<<"   X coordinate of the 7th track point: "<<x[7]<<endl;
          }
       }

       delete ev;  ev=0;
       delete evf; evf=0;

   }

   delete esdTree;
}
