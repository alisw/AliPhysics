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
   if(readFriend)esdTree->SetBranchStatus("ESDfriend*",1);
   ev->ReadFromTree(esdTree);

   // Attach the branch with ESD friends
   AliESDfriend *evf=0;
   if (readFriend) {
     evf = (AliESDfriend*)ev->FindListObject("AliESDfriend");
     if(!evf){
       // works for both, we just want to avoid setting the branch adress twice
       // in case of the new ESD
       esdTree->SetBranchAddress("ESDfriend.",&evf); 
     }
   }

   Int_t nev=esdTree->GetEntries();
   for (Int_t i=0; i<nev; i++) {
       esdTree->GetEntry(i);
       if(ev->GetAliESDOld())ev->CopyFromOldESD();
      
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
	 cout<<"Number of friend tracks: "<<evf->GetNumberOfTracks() <<endl;
	 cout<<"Track number 0"<<endl;
	 cout<<"   Number of TPC clusters: "<<n<<endl;
	 cout<<"   Index of the 1st TPC cluster: "  <<idx[1]<<endl;
	 UChar_t map=t->GetITSClusterMap();
	 cout<<"   ITS cluster map (from SPDs to SSDs):\n";
	 for (Int_t k=0; k<6; k++) printf(" Bit %d: %d\n",k,(map&(1<<k))==(1<<k)) <<' ';
	 
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
