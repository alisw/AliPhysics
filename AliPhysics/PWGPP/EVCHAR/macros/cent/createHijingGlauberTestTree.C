// $Id$

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TFolder.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliPDG.h"
#include "THijing/AliGenHijing.h"
#endif

class MyHeader
{
public:
  MyHeader() : fNATT(0), fEATT(0),fJATT(0),fNT(0),fNP(0),
               fN00(0),fN01(0),fN10(0),fN11(0),fBB(0),fRP(0),
               fPSn(0), fPSp(0), fTSn(0), fTSp(0) {;}
  virtual ~MyHeader() {;}
  Int_t       fNATT;
  Double32_t  fEATT;
  Int_t       fJATT;
  Int_t       fNT;
  Int_t       fNP;
  Int_t       fN00;
  Int_t       fN01;
  Int_t       fN10;
  Int_t       fN11;
  Double32_t  fBB;
  Double32_t  fRP;
  Int_t       fPSn;
  Int_t       fPSp;
  Int_t       fTSn;
  Int_t       fTSp;
  ClassDef(MyHeader,1) // My header class
};

class MyResponse
{
public:
  MyResponse() : fEtch0p(0), fEtch1p(0), fEtch2p(0), fEtch3p(0), fEtch4p(0), fEtch5p(0), fEtchrp(0), 
                 fEtch0n(0), fEtch1n(0), fEtch2n(0), fEtch3n(0), fEtch4n(0), fEtch5n(0), fEtchrn(0),
                 fNch0p(0), fNch1p(0), fNch2p(0), fNch3p(0), fNch4p(0), fNch5p(0), fNchrp(0),
                 fNch0n(0), fNch1n(0), fNch2n(0), fNch3n(0), fNch4n(0), fNch5n(0), fNchrn(0) {;}
  virtual ~MyResponse() {;}
  Double32_t  fEtch0p;
  Double32_t  fEtch1p;
  Double32_t  fEtch2p;
  Double32_t  fEtch3p;
  Double32_t  fEtch4p;
  Double32_t  fEtch5p;
  Double32_t  fEtchrp;
  Double32_t  fEtch0n;
  Double32_t  fEtch1n;
  Double32_t  fEtch2n;
  Double32_t  fEtch3n;
  Double32_t  fEtch4n;
  Double32_t  fEtch5n;
  Double32_t  fEtchrn;
  Int_t       fNch0p;
  Int_t       fNch1p;
  Int_t       fNch2p;
  Int_t       fNch3p;
  Int_t       fNch4p;
  Int_t       fNch5p;
  Int_t       fNchrp;
  Int_t       fNch0n;
  Int_t       fNch1n;
  Int_t       fNch2n;
  Int_t       fNch3n;
  Int_t       fNch4n;
  Int_t       fNch5n;
  Int_t       fNchrn;
  ClassDef(MyResponse,1) // My reponse class
};

void createGlauberTree(Int_t nEvents,
                       const char *outFileName = "treeout.root");

//-----------------------------------------------------------------------------------------------------

void createGlauberTree(Int_t nEvents,
                       const char *outFileName) 
{
  AliPDG::AddParticlesToPdgDataBase();
  TDatabasePDG::Instance();

  // Run loader
  TFolder *folder = new TFolder("myfolder","myfolder");
  AliRunLoader* rl = new AliRunLoader(folder);
  rl->MakeHeader();
  rl->MakeStack();
  AliStack* stack = rl->Stack();
  //AliHeader* rheader = rl->GetHeader();

  AliGenHijing *genHi = new AliGenHijing(-1);
  genHi->SetStack(stack);
  genHi->SetEnergyCMS(2760);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A", 208, 82);
  genHi->SetTarget    ("A", 208, 82);
  genHi->SetPtHardMin (2.3);
  genHi->SetImpactParameterRange(0.,30);
  genHi->SetJetQuenching(0); // enable jet quenching
  genHi->SetShadowing(1);    // enable shadowing
  genHi->SetDecaysOff(1);    // neutral pion and heavy particle decays switched off
  genHi->Init();

  MyHeader  *myheader = new MyHeader;
  MyResponse *myresp  = new MyResponse;

  TFile *outFile = TFile::Open(outFileName, "RECREATE");
  outFile->SetCompressionLevel(5);
  TDirectory::TContext context(outFile);

  TTree *tree = new TTree("glaubertree", "Glauber tree");
  tree->Branch("header",&myheader, 32*1024, 99);
  tree->Branch("response",&myresp, 32*1024, 99);

  TNtuple *ntuple = new TNtuple("gnt", "Glauber ntuple", "npart:ncoll:b");

  Double_t etas[] = {-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10};
  TH1D *hNEta = new TH1D("hNeta","",12,etas);
  TH1D *hEtEta = new TH1D("hEteta","",12,etas);

  // create events and fill them
  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent) {

    cout << "Event " << iEvent+1 << "/" << nEvents << endl;;
    stack->Reset();
    hNEta->Reset();
    hEtEta->Reset();
    genHi->Generate();
  
    AliStack *s = genHi->GetStack();
    const TObjArray *parts = s->Particles();
    Int_t nents = parts->GetEntries();
    for (Int_t i = 0; i<nents; ++i) {
      TParticle *p = (TParticle*)parts->At(i);
      //p->Print();
      TParticlePDG *pdg = p->GetPDG(1);
      Int_t c = (Int_t)(TMath::Abs(pdg->Charge()));
      if (c!=0) {
        hNEta->Fill(p->Eta());
        hEtEta->Fill(p->Eta(),p->Pt());
      }
    }

    AliGenHijingEventHeader *h = (AliGenHijingEventHeader*)genHi->CollisionGeometry();
    myheader->fNATT = nents;
    myheader->fEATT = h->TotalEnergy();
    myheader->fJATT = h->HardScatters();
    myheader->fNT   = h->TargetParticipants();
    myheader->fNP   = h->ProjectileParticipants();
    myheader->fN00  = h->NwNw();
    myheader->fN01  = h->NwN();
    myheader->fN10  = h->NNw();
    myheader->fN11  = h->NN();
    myheader->fBB   = h->ImpactParameter();
    myheader->fRP   = h->ReactionPlaneAngle();
    myheader->fPSn  = h->ProjSpectatorsn();
    myheader->fPSp  = h->ProjSpectatorsp();
    myheader->fTSn  = h->TargSpectatorsn();
    myheader->fTSp  = h->TargSpectatorsn();

    myresp->fEtch0p = hEtEta->GetBinContent(hEtEta->FindBin(0.5));
    myresp->fEtch1p = hEtEta->GetBinContent(hEtEta->FindBin(1.5));
    myresp->fEtch2p = hEtEta->GetBinContent(hEtEta->FindBin(2.5));
    myresp->fEtch3p = hEtEta->GetBinContent(hEtEta->FindBin(3.5));
    myresp->fEtch4p = hEtEta->GetBinContent(hEtEta->FindBin(4.5));
    myresp->fEtch5p = hEtEta->GetBinContent(hEtEta->FindBin(5.5));
    myresp->fEtchrp = hEtEta->GetBinContent(hEtEta->FindBin(10.5));
    myresp->fEtch0n = hEtEta->GetBinContent(hEtEta->FindBin(-0.5));
    myresp->fEtch1n = hEtEta->GetBinContent(hEtEta->FindBin(-1.5));
    myresp->fEtch2n = hEtEta->GetBinContent(hEtEta->FindBin(-2.5));
    myresp->fEtch3n = hEtEta->GetBinContent(hEtEta->FindBin(-3.5));
    myresp->fEtch4n = hEtEta->GetBinContent(hEtEta->FindBin(-4.5));
    myresp->fEtch5n = hEtEta->GetBinContent(hEtEta->FindBin(-5.5));
    myresp->fEtchrn = hEtEta->GetBinContent(hEtEta->FindBin(-10.5));
    myresp->fNch0p  = hNEta->GetBinContent(hNEta->FindBin(0.5));
    myresp->fNch1p  = hNEta->GetBinContent(hNEta->FindBin(1.5));
    myresp->fNch2p  = hNEta->GetBinContent(hNEta->FindBin(2.5));
    myresp->fNch3p  = hNEta->GetBinContent(hNEta->FindBin(3.5));
    myresp->fNch4p  = hNEta->GetBinContent(hNEta->FindBin(4.5));
    myresp->fNch5p  = hNEta->GetBinContent(hNEta->FindBin(5.5));
    myresp->fNchrp  = hNEta->GetBinContent(hNEta->FindBin(10.5));
    myresp->fNch0n  = hNEta->GetBinContent(hNEta->FindBin(-0.5));
    myresp->fNch1n  = hNEta->GetBinContent(hNEta->FindBin(-1.5));
    myresp->fNch2n  = hNEta->GetBinContent(hNEta->FindBin(-2.5));
    myresp->fNch3n  = hNEta->GetBinContent(hNEta->FindBin(-3.5));
    myresp->fNch4n  = hNEta->GetBinContent(hNEta->FindBin(-4.5));
    myresp->fNch5n  = hNEta->GetBinContent(hNEta->FindBin(-5.5));
    myresp->fNchrn  = hNEta->GetBinContent(hNEta->FindBin(-10.5));

    tree->Fill();

    if (ntuple) {
      Int_t np = h->TargetParticipants() + h->ProjectileParticipants();
      Int_t nc = h->NwNw() + h->NwN() + h->NNw() + h->NN();
      Double_t b = h->ImpactParameter();
      ntuple->Fill(np,nc,b);
    }

  } // end of event loop

  tree->Write();
  ntuple->Write();
  outFile->Close();
}
