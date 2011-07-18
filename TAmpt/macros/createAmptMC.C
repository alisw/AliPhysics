#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
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
#include "AliGenAmpt.h"
#include "TAmpt.h"
#endif

class MyHeader
{
public:
  MyHeader() : fNATT(0), fEATT(0),fJATT(0),fNT(0),fNP(0),fN0(0),fN01(0),fN10(0),fN11(0),fBB(0),fRP(0) {;}
  virtual ~MyHeader() {;}
  Int_t       fNATT;
  Float_t     fEATT;
  Int_t       fJATT;
  Int_t       fNT;
  Int_t       fNP;
  Int_t       fN0;
  Int_t       fN01;
  Int_t       fN10;
  Int_t       fN11;
  Float_t     fBB;
  Float_t     fRP;
  ClassDef(MyHeader,1) // My header class
};

class MyNuc : public TObject
{
public:
  MyNuc(Int_t pid=0, Int_t st=0, Int_t type=0, Double_t x=0, Double_t y=0, Double_t z=0) :
    TObject(), fPid(pid), fSt(st), fType(type), fX(x), fY(y), fZ(z) {;}
  Double_t    X()    const { return fX; }
  Double_t    Y()    const { return fY; }
  Double_t    Z()    const { return fZ; }
  Int_t       Pid()  const { return fPid; }
  Int_t       St()   const { return fSt; }
  Int_t       Type() const { return fType; }
protected:
  Int_t       fPid;
  Int_t       fSt;
  Int_t       fType;
  Double32_t  fX;
  Double32_t  fY;
  Double32_t  fZ;
  ClassDef(MyNuc,1) // My nucleon class in cartesian coordinates
};

class MyPart : public TObject
{
public:
  MyPart(Int_t pid=0, Int_t st=0, Int_t type=0, Double_t pt=0, Double_t eta=0, Double_t phi=0) :
    TObject(), fPid(pid), fSt(st), fType(type), fPt(pt), fEta(eta), fPhi(phi) {;}
  Double_t    Px()   const { return fPt*TMath::Cos(fPhi);  }
  Double_t    Py()   const { return fPt*TMath::Sin(fPhi);  }
  Double_t    Pz()   const { return fPt*TMath::SinH(fEta); }
  Double_t    Pt()   const { return fPt;  }
  Double_t    Eta()  const { return fEta; }
  Double_t    Phi()  const { return fPhi; }
  Int_t       Pid()  const { return fPid; }
  Int_t       St()   const { return fSt; }
  Int_t       Type() const { return fType; }
protected:
  Int_t       fPid;
  Int_t       fSt;
  Int_t       fType;
  Double32_t  fPt;
  Double32_t  fEta;
  Double32_t  fPhi;
  ClassDef(MyPart,1) // My particle class in cylindrical coordinates
};

void createAmptMC(Int_t nEvents,
                  const char *outFileName = "amptout.root");

void anaAmptMC(Int_t nEvents,
               const char *inFileNames = "/eliza6/alice/loizides/ampt/run_1/ampt*.root");

//-----------------------------------------------------------------------------------------------------

void createAmptMC(Int_t nEvents,
                  const char *outFileName) 
{
  TClass::GetClass("MyNuc")->IgnoreTObjectStreamer();
  TClass::GetClass("MyPart")->IgnoreTObjectStreamer();

  AliPDG::AddParticlesToPdgDataBase();
  TDatabasePDG::Instance();

  // Run loader
  TFolder *folder = new TFolder("myfolder","myfolder");
  AliRunLoader* rl = new AliRunLoader(folder);
  rl->MakeHeader();
  rl->MakeStack();
  AliStack* stack = rl->Stack();
  //AliHeader* rheader = rl->GetHeader();

  AliGenAmpt *genHi = new AliGenAmpt(-1);
  genHi->SetStack(stack);
  genHi->SetEnergyCMS(2760);
  genHi->SetReferenceFrame("CMS");
  genHi->SetProjectile("A", 208, 82);
  genHi->SetTarget    ("A", 208, 82);
  genHi->SetPtHardMin (2);
  genHi->SetImpactParameterRange(0,30);
  genHi->SetJetQuenching(0); // enable jet quenching
  genHi->SetShadowing(1);    // enable shadowing
  genHi->SetDecaysOff(1);    // neutral pion and heavy particle decays switched off
  genHi->SetSpectators(0);   // track spectators 
  genHi->SetIsoft(4);        // 4=string melting, 1=standard AMPT
  genHi->SetXmu(3.2264);     // parton xsection
  genHi->SetNtMax(150);      // time bins
  if (0) { //RHIC settings
    genHi->SetAlpha(0.47140452);
    genHi->SetStringFrag(2.2,0.5);
  }
  genHi->Init();

  TClonesArray *inucs = new TClonesArray("TParticle",10000);
  TClonesArray *parts = new TClonesArray("MyPart",10000);
  TClonesArray *nucs  = new TClonesArray("MyNuc",500);
  MyHeader *myheader = new MyHeader;

  TFile *outFile = TFile::Open(outFileName, "RECREATE");
  outFile->SetCompressionLevel(5);
  TDirectory::TContext context(outFile);

  TTree *tree = new TTree("ampt", "AMPT tree");
  tree->Branch("parts",&parts, 256*1024, 99);
  tree->Branch("nucs",&nucs, 64*1024, 99);
  tree->Branch("info",&myheader, 32*1024, 99);

  // create events and fill them
  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent) {

    cout << "Event " << iEvent+1 << "/" << nEvents << endl;;
    stack->Reset();
    genHi->Generate();
    parts->Clear();
    const TClonesArray *iarr = genHi->GetParticles();
    if (iarr) {
      for(Int_t i=0;i<iarr->GetEntriesFast();++i) {
        TParticle *p = (TParticle*)iarr->At(i);
        new((*parts)[i]) MyPart(p->GetPdgCode(),
                                p->GetStatusCode(),
                                p->GetUniqueID(),
                                p->Pt(),
                                p->Eta(),
                                p->Phi());
      }
    }
    TAmpt *ampt = genHi->Ampt();
    if (ampt) {
      ampt->ImportNucleons(inucs);
      nucs->Clear();
      for(Int_t i=0;i<inucs->GetEntriesFast();++i) {
        TParticle *p = (TParticle*)inucs->At(i);
        new((*nucs)[i]) MyNuc(p->GetPdgCode(),
                              p->GetStatusCode(),
                              p->GetUniqueID(),
                              p->Vx(),
                              p->Vy(),
                              p->Vz());
      }

      myheader->fNATT=ampt->GetNATT();
      myheader->fEATT=ampt->GetEATT();
      myheader->fJATT=ampt->GetJATT();
      myheader->fNT=ampt->GetNT();
      myheader->fNP=ampt->GetNP();
      myheader->fN0=ampt->GetN0();
      myheader->fN01=ampt->GetN01();
      myheader->fN10=ampt->GetN10();
      myheader->fN11=ampt->GetN11();
      myheader->fBB=ampt->GetBB();
      myheader->fRP=ampt->GetPhi();
    }
    tree->Fill();
  } // end of event loop
  tree->Write();
  outFile->Close();
}

//-----------------------------------------------------------------------------------------------------

void anaAmptMC(Int_t nEvents,
               const char *inFileNames)
{

  TChain *c = new TChain("ampt");
  c->Add(inFileNames);


  TClonesArray *parts = 0;
  TClonesArray *nucs  = 0;
  MyHeader *info = 0;

  c->SetBranchAddress("info",&info);
  c->SetBranchAddress("nucs",&nucs);
  c->SetBranchAddress("parts",&parts);

  Int_t nRead = nEvents;
  if (nRead<0)
    nRead = c->GetEntries();
  else if (0 && (nRead>c->GetEntries()))
    nRead = c->GetEntries();

   for (Int_t ev=0;ev<nRead;++ev) {
      c->GetEntry(ev);

      Int_t fAN=0, fBN=0, fAStat[250], fBStat[250];
      Double_t fXA[250], fXB[250], fYA[250], fYB[250];
      for (Int_t k=0;k<nucs->GetEntries();++k) {
        MyNuc *n = (MyNuc*)nucs->At(k);
        if (n->Type()>0) {
          fAStat[fAN] =0;
          fXA[fAN] = n->X();
          fYA[fAN] = n->Y();
          ++fAN;
        } else {
          fBStat[fBN] =0;
          fXB[fBN] = n->X();
          fYB[fBN] = n->Y();
          ++fBN;
        }
      }

      // Glauber calculation
      Double_t fXSect = 65; //mbarn
      Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

      // For each of the A nucleons in nucleus B
      for (Int_t i = 0; i<fBN; i++) {
        for (Int_t j = 0 ; j < fAN ;j++) {
          Double_t dx = fXB[i]-fXA[j];
          Double_t dy = fYB[i]-fYA[j];
          Double_t dij = dx*dx+dy*dy;
          if (dij < d2) {
            fBStat[i]++;
            fAStat[j]++;
          }
        }
      }
      // Calculate npart
      Int_t npart=0;
      for (Int_t i = 0; i<fAN; i++) {
        if (fAStat[i]>0) {
          npart++;
        }
      }
      for (Int_t i = 0; i<fBN; i++) {
        if (fBStat[i]>0) {
          npart++;
        }
      }
      cout << ev << " : "  << info->fBB << " np " << npart << " vs " << info->fNP+info->fNT << endl;
   }
}
