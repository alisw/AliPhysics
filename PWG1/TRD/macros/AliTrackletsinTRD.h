//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 20 09:42:44 2009 by ROOT version 5.21/01
// from TTree AliTrackletsinTRD/AliTrackletsinTRD
// found on file: TRD.DebugTrackingMultiplicity.root
//////////////////////////////////////////////////////////

#ifndef ALITRACKLETSINTRD_h
#define ALITRACKLETSINTRD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class AliTrackletsinTRD {
public :
   TTree           *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           fstandalone;   //!flag is track reconstructed in TPC only or TPC+TRD
   Int_t           feventcounter; //!event number
   Int_t           flayer;        //!layer number
   Float_t         fxtracklet;    //!x position of tracklet
   Double_t        fxtrack;       //!x position of track
   Float_t         fytracklet;    //!y position of tracklet
   Double_t        fytrack;       //!y position of track
   Float_t         fztracklet;    //!z position of tracklet
   Double_t        fztrack;       //!z position of track
   Int_t           fnumtracklets; //!number of tracklets
   Int_t           fdettracklet;  //!detector number

   // List of branches
   TBranch        *fbB0;   //!
   TBranch        *fbB1;   //!
   TBranch        *fbB2;   //!
   TBranch        *fbB3;   //!
   TBranch        *fbB4;   //!
   TBranch        *fbB5;   //!
   TBranch        *fbB6;   //!
   TBranch        *fbB7;   //!
   TBranch        *fbB8;   //!
   TBranch        *fbB9;   //!
   TBranch        *fbB10;  //!

   AliTrackletsinTRD(TTree *tree=0);
   virtual ~AliTrackletsinTRD();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


   AliTrackletsinTRD (AliTrackletsinTRD& p):
       fChain(p.fChain),
       fbB0(p.fbB0),
       fbB1(p.fbB1),
       fbB2(p.fbB2),
       fbB3(p.fbB3),
       fbB4(p.fbB4),
       fbB5(p.fbB5),
       fbB6(p.fbB6),
       fbB7(p.fbB7),
       fbB8(p.fbB8),
       fbB9(p.fbB9),
       fbB10(p.fbB10)
   {
   }

};

  

AliTrackletsinTRD::AliTrackletsinTRD(TTree *tree)
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

AliTrackletsinTRD::~AliTrackletsinTRD()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AliTrackletsinTRD::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AliTrackletsinTRD::LoadTree(Long64_t entry)
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

void AliTrackletsinTRD::Init(TTree *tree)
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

   fChain->SetBranchAddress("standalone", &fstandalone, &fbB0);
   fChain->SetBranchAddress("eventcounter", &feventcounter, &fbB1);
   fChain->SetBranchAddress("layer", &flayer, &fbB2);
   fChain->SetBranchAddress("xtracklet", &fxtracklet, &fbB3);
   fChain->SetBranchAddress("xtrack", &fxtrack, &fbB4);
   fChain->SetBranchAddress("ytracklet", &fytracklet, &fbB5);
   fChain->SetBranchAddress("ytrack", &fytrack, &fbB6);
   fChain->SetBranchAddress("ztracklet", &fztracklet, &fbB7);
   fChain->SetBranchAddress("ztrack", &fztrack, &fbB8);
   fChain->SetBranchAddress("num_tracklets", &fnumtracklets, &fbB9);
   fChain->SetBranchAddress("dettracklet", &fdettracklet, &fbB10);
   Notify();
}

Bool_t AliTrackletsinTRD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AliTrackletsinTRD::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AliTrackletsinTRD::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AliTrackletsinTRD_cxx
