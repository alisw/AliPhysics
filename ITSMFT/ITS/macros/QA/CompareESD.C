//***********************************************************
// A simple macro to check an approximate identity of tracks 
// taken from two different ESD files
// (A good test for a multi-threaded tracking)
//***********************************************************   

#if !defined(__CINT__) || defined(__MAKECINT__)
   #include <Riostream.h>
   #include <TFile.h>
   #include <TTree.h>

   #include "AliESDEvent.h"
#endif

void CompareESD() {
  TFile *fold=TFile::Open("AliESDs_old.root");
  TFile *fnew=TFile::Open("AliESDs.root");

  TTree *told=(TTree*)fold->Get("esdTree");
  Int_t nold=told->GetEntries();
  TTree *tnew=(TTree*)fnew->Get("esdTree");
  Int_t nnew=tnew->GetEntries();
  if (nold != nnew) {
    cout<<"Number of events is different !"<<endl;
    return;
  }

  AliESDEvent *eold = new AliESDEvent();
  eold->ReadFromTree(told); 
  AliESDEvent *enew = new AliESDEvent();
  enew->ReadFromTree(tnew); 

  for (Int_t i=0; i<nnew; i++) {
      cout<<"Checking event #"<<i<<endl;
      told->GetEvent(i);
      tnew->GetEvent(i);
      Int_t ntold=eold->GetNumberOfTracks();
      Int_t ntnew=enew->GetNumberOfTracks();
      if (ntold != ntnew) {
	 cout<<"Number of tracks is different ! "<<ntold<<' '<<ntnew<<endl;
         return;
      }
      for (Int_t j=0; j<ntnew; j++) {
	  AliESDtrack *to=eold->GetTrack(j);
	  AliESDtrack *tn=enew->GetTrack(j);
          UInt_t so=to->GetStatus();
          UInt_t sn=tn->GetStatus();
          if (so != sn) {
	     cout<<"Track status is different !"<<endl;
	     return;
          }
          Int_t idxo[20], idxn[20];
          Int_t nco=to->GetITSclusters(idxo);
          Int_t ncn=tn->GetITSclusters(idxn);
          if (nco != ncn) {
	     cout<<"Track #"<<j<<' ';
	     cout<<"Number of clusters is different ! "<<nco<<' '<<ncn<<endl;
	     return;
          }
          for (Int_t k=0; k<ncn; k++) {
	      if (idxo[k] != idxn[k]) {
	         cout<<"Cluster index is different !"<<endl;
	         return;
	      }
	  }
      }
  }

}
