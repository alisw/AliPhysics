void MUONtrackRecoModel() {
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Thu Sep 21 14:53:11 2000 by ROOT version2.25/02)
//   from TTree MUONtrackReco/MUONtrackReco
//   found on file: MUONtrackReco.root
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MUONtrackReco.root");
   if (!f) {
      f = new TFile("MUONtrackReco.root");
   }
   TTree *MUONtrackReco = (TTree*)gDirectory->Get("MUONtrackReco");

//Declaration of leaves types
   Int_t           fEvent;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           Tracks_;
   Int_t           Tracks_fCharge[5];
   Float_t         Tracks_fPxRec[5];
   Float_t         Tracks_fPyRec[5];
   Float_t         Tracks_fPzRec[5];
   Float_t         Tracks_fZRec[5];
   Float_t         Tracks_fZRec1[5];
   Int_t           Tracks_fNHits[5];
   Float_t         Tracks_fChi2[5];
   Float_t         Tracks_fPxGen[5];
   Float_t         Tracks_fPyGen[5];
   Float_t         Tracks_fPzGen[5];
   UInt_t          Tracks_fUniqueID[5];
   UInt_t          Tracks_fBits[5];

//Set branch addresses
 //MUONtrackReco->SetBranchAddress("Header",&Header);
   MUONtrackReco->SetBranchAddress("fEvent",&fEvent);
   MUONtrackReco->SetBranchAddress("fUniqueID",&fUniqueID);
   MUONtrackReco->SetBranchAddress("fBits",&fBits);
   MUONtrackReco->SetBranchAddress("Tracks_",&Tracks_);
   MUONtrackReco->SetBranchAddress("Tracks.fCharge",Tracks_fCharge);
   MUONtrackReco->SetBranchAddress("Tracks.fPxRec",Tracks_fPxRec);
   MUONtrackReco->SetBranchAddress("Tracks.fPyRec",Tracks_fPyRec);
   MUONtrackReco->SetBranchAddress("Tracks.fPzRec",Tracks_fPzRec);
   MUONtrackReco->SetBranchAddress("Tracks.fZRec",Tracks_fZRec);
   MUONtrackReco->SetBranchAddress("Tracks.fZRec1",Tracks_fZRec1);
   MUONtrackReco->SetBranchAddress("Tracks.fNHits",Tracks_fNHits);
   MUONtrackReco->SetBranchAddress("Tracks.fChi2",Tracks_fChi2);
   MUONtrackReco->SetBranchAddress("Tracks.fPxGen",Tracks_fPxGen);
   MUONtrackReco->SetBranchAddress("Tracks.fPyGen",Tracks_fPyGen);
   MUONtrackReco->SetBranchAddress("Tracks.fPzGen",Tracks_fPzGen);
   MUONtrackReco->SetBranchAddress("Tracks.fUniqueID",Tracks_fUniqueID);
   MUONtrackReco->SetBranchAddress("Tracks.fBits",Tracks_fBits);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// MUONtrackReco->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Int_t nentries = MUONtrackReco->GetEntries();

   Int_t nbytes = 0;
//   for (Int_t i=0; i<nentries;i++) {
//      nbytes += MUONtrackReco->GetEntry(i);
//   }
}
