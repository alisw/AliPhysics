//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 20 09:42:44 2009 by ROOT version 5.21/01
// from TTree TrackletsinTRD/TrackletsinTRD
// found on file: TRD.DebugTrackingMultiplicity.root
//////////////////////////////////////////////////////////

#ifndef TrackletsinTRD_h
#define TrackletsinTRD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class TrackletsinTRD {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           standalone;
   Int_t           eventcounter;
   Int_t           layer;
   Float_t         xtracklet;
   Double_t        xtrack;
   Float_t         ytracklet;
   Double_t        ytrack;
   Float_t         ztracklet;
   Double_t        ztrack;
   Int_t           num_tracklets;
   Int_t           dettracklet;

   // List of branches
   TBranch        *b_B0;   //!
   TBranch        *b_B1;   //!
   TBranch        *b_B2;   //!
   TBranch        *b_B3;   //!
   TBranch        *b_B4;   //!
   TBranch        *b_B5;   //!
   TBranch        *b_B6;   //!
   TBranch        *b_B7;   //!
   TBranch        *b_B8;   //!
   TBranch        *b_B9;   //!
   TBranch        *b_B10;  //!

   TrackletsinTRD(TTree *tree=0);
   virtual ~TrackletsinTRD();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TrackletsinTRD_cxx
TrackletsinTRD::TrackletsinTRD(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/alidata80/alice_u/pachmay/AliRoot/v4-16-Rev-01/TRD/qaRec/TRD.DebugTrackingMultiplicity.root");
      if (!f) {
         f = new TFile("/alidata80/alice_u/pachmay/AliRoot/v4-16-Rev-01/TRD/qaRec/TRD.DebugTrackingMultiplicity.root");
      }
      tree = (TTree*)gDirectory->Get("TrackletsinTRD");

   }
   Init(tree);
}

TrackletsinTRD::~TrackletsinTRD()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TrackletsinTRD::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TrackletsinTRD::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TrackletsinTRD::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("standalone", &standalone, &b_B0);
   fChain->SetBranchAddress("eventcounter", &eventcounter, &b_B1);
   fChain->SetBranchAddress("layer", &layer, &b_B2);
   fChain->SetBranchAddress("xtracklet", &xtracklet, &b_B3);
   fChain->SetBranchAddress("xtrack", &xtrack, &b_B4);
   fChain->SetBranchAddress("ytracklet", &ytracklet, &b_B5);
   fChain->SetBranchAddress("ytrack", &ytrack, &b_B6);
   fChain->SetBranchAddress("ztracklet", &ztracklet, &b_B7);
   fChain->SetBranchAddress("ztrack", &ztrack, &b_B8);
   fChain->SetBranchAddress("num_tracklets", &num_tracklets, &b_B9);
   fChain->SetBranchAddress("dettracklet", &dettracklet, &b_B10);
   Notify();
}

Bool_t TrackletsinTRD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TrackletsinTRD::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TrackletsinTRD::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TrackletsinTRD_cxx
