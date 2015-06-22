// Tree Analysis Task for Pb-Pb data reduced event trees

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 18 15:43:18 2012 by ROOT version 5.34/02
// from TTree HFEtree/HFE event tree
// found on file: HFEevents.root
//////////////////////////////////////////////////////////

#ifndef IPStudy_h
#define IPStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.
#include "AliHFEreducedEvent.h"
#include "AliHFEreducedTrack.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class IPStudy : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   AliHFEreducedEvent *HFEevent;

   // List of branches
   TBranch        *b_HFEevent;   //!

   IPStudy(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~IPStudy() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   Bool_t passesCuts(AliHFEreducedTrack * hfeTrack);
   Bool_t EventPassesCuts();
   
private :
   TH1D * fTrackDCAerr;
   TH1D * fhMultiplicity;
   TH3D * fTrackDCAerrpT;
   TH3D * fTrackDCAerrpTConv;
   TH3D * fTrackDCApTnoCut;
   TH2D * fTrackDCApTnoCutHadron;
   TH3D * fTrackDCApTCut;
   TH2D * fTrackProdRnoCut;
   TH2D * fTrackProdRnoCutV0;
   TH2D * fTrackProdRwithCut;
   TH2D * fTrackDCAerrpTFake;
   TH2D * fTrackDCAerrpTTrueConv;
   //TH2D * TrackEta;
   TH1D * fConvPhi;
   TH1D * fTOFPhi;
   //TH3D * TPCTOFsigmas;
   TH1D * fEventNr;
   TH2D * fNContribVertex;
   TH2D * fDeuteronTOF;
   TH1D * fBin0Numbers;
   TH2D * fTPCPionV0;
   TH2D * fTPCElectronV0;
   TH2D * fTPCSignal;
   TH2D * fITSSignal;
   TH2D * fTPCPionV0nsigma;
   TH2D * fTPCElectronV0nsigma;
   TH2D * fK0pTDCA;
   TH2D * fK0pTDCARCut;
   TH2D * fTPCElectronV0nsigmaNoCuts;
   
   TH2D * fVertexXY;

   ClassDef(IPStudy,0);
};

#endif

#ifdef IPStudy_cxx
void IPStudy::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   HFEevent = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("HFEevent", &HFEevent, &b_HFEevent);
}

Bool_t IPStudy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef IPStudy_cxx
