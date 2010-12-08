//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 26 11:48:52 2010 by ROOT version 5.26/00
// from TTree fMyTr/beam tree
// found on file: Memory Directory
//////////////////////////////////////////////////////////

#ifndef AnaEveMyTree_h
#define AnaEveMyTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class AnaEveMyTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           jev;
   Int_t           iev;
   Int_t           ntr;
   Int_t           etr;
   Float_t         z1n;
   Float_t         z2n;
   Float_t         z1p;
   Float_t         z2p;
   Float_t         zvt;
   Float_t         z1s;
   Float_t         z2s;
   Float_t         w1s;
   Float_t         w2s;
   Float_t         t1s;
   Float_t         t2s;
   Float_t         v1s;
   Float_t         v2s;
   Float_t         z1nt[10];
   Float_t         z2nt[10];
   Float_t         z1pt[10];
   Float_t         z2pt[10];
   Float_t         t0amp[24];
   Float_t         t0tim[24];
   Float_t         v0mul[64];

   // List of branches
   TBranch        *b_jev;   //!
   TBranch        *b_iev;   //!
   TBranch        *b_ntr;   //!
   TBranch        *b_etr;   //!
   TBranch        *b_z1n;   //!
   TBranch        *b_z2n;   //!
   TBranch        *b_z1p;   //!
   TBranch        *b_z2p;   //!
   TBranch        *b_zvt;   //!
   TBranch        *b_z1s;   //!
   TBranch        *b_z2s;   //!
   TBranch        *b_w1s;   //!
   TBranch        *b_w2s;   //!
   TBranch        *b_t1s;   //!
   TBranch        *b_t2s;   //!
   TBranch        *b_v1s;   //!
   TBranch        *b_v2s;   //!
   TBranch        *b_z1nt;   //!
   TBranch        *b_z2nt;   //!
   TBranch        *b_z1pt;   //!
   TBranch        *b_z2pt;   //!
   TBranch        *b_t0amp;   //!
   TBranch        *b_t0tim;   //!
   TBranch        *b_v0mul;   //!

   AnaEveMyTree(int irun=137431, int ifile=2, TTree *tree=0);
   virtual ~AnaEveMyTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnaEveMyTree_cxx
AnaEveMyTree::AnaEveMyTree(int irun, int ifile, TTree *tree)
{
// if parameter tree is not specified (or zero), connect the ifile
// used to generate this class and read the Tree.
   char name[80];
   sprintf(name,"%06d/EveMyTree%03d.root",irun,ifile);
   cout << "open file :" << name << endl;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(name);
      if (!f) {
         f = new TFile(name);
      }
      tree = (TTree*)gROOT->FindObject("chist")->FindObject("fMyTr");
   }
   Init(tree);
}

AnaEveMyTree::~AnaEveMyTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaEveMyTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaEveMyTree::LoadTree(Long64_t entry)
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

void AnaEveMyTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("jev", &jev, &b_jev);
   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("ntr", &ntr, &b_ntr);
   fChain->SetBranchAddress("etr", &etr, &b_etr);
   fChain->SetBranchAddress("z1n", &z1n, &b_z1n);
   fChain->SetBranchAddress("z2n", &z2n, &b_z2n);
   fChain->SetBranchAddress("z1p", &z1p, &b_z1p);
   fChain->SetBranchAddress("z2p", &z2p, &b_z2p);
   fChain->SetBranchAddress("zvt", &zvt, &b_zvt);
   fChain->SetBranchAddress("z1s", &z1s, &b_z1s);
   fChain->SetBranchAddress("z2s", &z2s, &b_z2s);
   fChain->SetBranchAddress("w1s", &w1s, &b_w1s);
   fChain->SetBranchAddress("w2s", &w2s, &b_w2s);
   fChain->SetBranchAddress("t1s", &t1s, &b_t1s);
   fChain->SetBranchAddress("t2s", &t2s, &b_t2s);
   fChain->SetBranchAddress("v1s", &v1s, &b_v1s);
   fChain->SetBranchAddress("v2s", &v2s, &b_v2s);
   fChain->SetBranchAddress("z1nt", z1nt, &b_z1nt);
   fChain->SetBranchAddress("z2nt", z2nt, &b_z2nt);
   fChain->SetBranchAddress("z1pt", z1pt, &b_z1pt);
   fChain->SetBranchAddress("z2pt", z2pt, &b_z2pt);
   fChain->SetBranchAddress("t0amp", t0amp, &b_t0amp);
   fChain->SetBranchAddress("t0tim", t0tim, &b_t0tim);
   fChain->SetBranchAddress("v0mul", v0mul, &b_v0mul);
   Notify();
}

Bool_t AnaEveMyTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaEveMyTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaEveMyTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnaEveMyTree_cxx
