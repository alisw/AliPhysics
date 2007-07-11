//********************************************************************
//  Example of accessing the information stored in ESD friends
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//********************************************************************

#if !defined( __CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TFile.h>
  #include <TChain.h>

  #include "AliESDEvent.h"
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

   AliESDEvent *ev= new AliESDEvent();
   ev->ReadFromTree(esdTree);

   // Attach the branch with ESD friends
   AliESDfriend *evf=0;
   if (readFriend) {
      esdTree->SetBranchStatus("ESDfriend*",1);
      evf = (AliESDfriend*)ev->FindListObject("AliESDfriend");
   }

   Int_t nev=esdTree->GetEntries();
   for (Int_t i=0; i<nev; i++) {
       esdTree->GetEntry(i);

      
       cout<<endl<<"Event number: "<<i<<endl;
       Int_t ntr=ev->GetNumberOfTracks();
       cout<<"Number of tracks: "<<ntr<<endl;
       // Now the attached information can be accessed via pointer to ESD.
       // Example: indices of the TPC clusters associated with the track number 0.
       if (ntr > 0) {
	 ev->SetESDfriend(evf); //Attach the friend to the ESD
          const AliESDtrack *t=ev->GetTrack(0);
          Int_t idx[AliESDfriendTrack::kMaxTPCcluster]; 
          n=t->GetTPCclusters(idx);
          cout<<"Track number 0"<<endl;
          cout<<"   Number of TPC clusters: "<<n<<endl;
          cout<<"   Index of the 150th TPC cluster: "<<idx[7]<<endl;
          UChar_t map=t->GetITSClusterMap();
          cout<<"   ITS cluster map (from SPDs to SSDs): ";
          for (Int_t i=0; i<6; i++) printf(" Bit %d: %d\n",i,(map&(1<<i))==(1<<i)) <<' ';
	  
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
       
   }
   delete ev;
   delete esdTree;
}
