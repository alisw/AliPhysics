// $Id$

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TH2F.h>
#include <TTimeStamp.h>
#include "TreeClasses.h"
#include "EventPoolManager.h"
#endif

void anaCorr(const char *inFileNames,
             Int_t    nEvents=-1,
             Int_t    nmin=10,
             Int_t    nmax=100,
             Double_t vzmin=-8,
             Double_t vzmax=+8,
             Int_t    vzbins=8,
             Int_t    vcmin=10,
             Int_t    psize=100,
             Double_t nmix=10,
             Double_t ptmin=0.4, 
             Double_t ptmax=10, 
             Double_t etamin=-1.4, 
             Double_t etamax=+1.4);

inline Double_t DeltaPhi(const MyPart &t1, const MyPart &t2, 
                         Double_t rangeMin = -TMath::Pi()/2, 
                         Double_t rangeMax = 3*TMath::Pi()/2);
inline Double_t DeltaEta(const MyPart &t1, const MyPart &t2);  
inline Bool_t   InBounds(Double_t val, Double_t min, Double_t max);
inline Double_t InvMass(const MyPart &p1, const MyPart &p2);

class Noti: public TObject
{
public:
  Noti() : fc(0) {;}
  virtual ~Noti() {;}
  Bool_t Notify()         { fc=1; return 1; }
  Bool_t Notified() const { return fc;      }
  void   Reset()          { fc=0;           }
protected:
  Bool_t fc; //=1 when file changed
  ClassDef (Noti,0)
};

//-----------------------------------------------------------------------------------------------------

void anaCorr(const char *inFileNames,
             Int_t    nEvents,
             Int_t    nmin,
             Int_t    nmax,
             Double_t vzmin,
             Double_t vzmax,
             Int_t    vzbins,
             Int_t    vcmin,
             Int_t    psize,
             Double_t nmix,
             Double_t ptmin,
             Double_t ptmax,
             Double_t etamin,
             Double_t etamax)
{
  Noti *nme = new Noti;
  TChain *tree = new TChain("MyTree");
  tree->Add(inFileNames);
  tree->SetNotify(nme);

  Int_t nents = tree->GetEntries();
  cout << "Found " << nents << " entries!" << endl;
  if (nents<=0)
    return;
  if (nEvents>0&&nEvents<nents) { 
    nents=nEvents;
    cout << "Using " << nents << " entries!" << endl;
  }

  if (TClass::GetClass("MyPart"))
    TClass::GetClass("MyPart")->IgnoreTObjectStreamer();
  if (TClass::GetClass("MyTracklet"))
    TClass::GetClass("MyTracklet")->IgnoreTObjectStreamer();

  MyHeader     *header  = 0;
  TClonesArray *tracks  = 0;
  TBranch      *hBranch = 0;
  TBranch      *tBranch = 0;

  Double_t *multbinsa = new Double_t[2]; multbinsa[0] = nmin;  multbinsa[1] = nmax; 
  Double_t *vzbinsa = new Double_t[vzbins+1];
  Double_t vzbinwidth = Double_t(vzmax-vzmin)/vzbins;
  for(Int_t i=0;i<vzbins+1;++i) {
    vzbinsa[i] = vzmin+i*vzbinwidth;
  }
  TH1F* hMultBins = new TH1F("hMultBins", "Event multiplicity binning", 
			     1, multbinsa);
  TH1F* hZvtxBins = new TH1F("hZvtxBins", "Event Z-vertex binning", 
			     vzbins, vzbinsa);


  EventPoolManager *pm = new EventPoolManager;
  pm->InitEventPools(psize, 1, multbinsa, vzbins, vzbinsa);

  TH2F *hSig = new TH2F("hSig",";#Delta#eta;#Delta#phi",60,-3,3,64,-TMath::Pi()/2,3*TMath::Pi()/2);
  TH2F *hBkg = new TH2F("hBgk",";#Delta#eta;#Delta#phi",60,-3,3,64,-TMath::Pi()/2,3*TMath::Pi()/2);

  Double_t dpmax = 0.03;
  Double_t demax = 0.01;

  for (Int_t i=0;i<nents;++i) {
    Int_t li = tree->LoadTree(i);
    if (nme->Notified()) {
      hBranch = tree->GetBranch("header");
      hBranch->SetAddress(&header);
      tBranch = tree->GetBranch("parts");
      tBranch->SetAddress(&tracks);
      nme->Reset();
    }

    hBranch->GetEntry(li);
    if (header->fIsPileupSPD)
      continue;
    if ((header->fVz<vzmin) ||
        (header->fVz>vzmax) )
      continue;
    if ((header->fVc<vcmin) ||
        (header->fVc>nmax) )
      continue;

    tBranch->GetEntry(li);
    Int_t ntracks = tracks->GetEntries();
    TClonesArray *seltracks = new TClonesArray("MyPart",ntracks);
    Int_t nseltracks = 0;
    for (Int_t t1=0;t1<ntracks;++t1) {
      MyPart *part1 = (MyPart*)tracks->At(t1);
      if (!InBounds(part1->fPt,ptmin,ptmax))
        continue;
      if (!InBounds(part1->fEta,etamin,etamax))
        continue;
      if (TMath::Abs(part1->fD)>3)
        continue;
      if (TMath::Abs(part1->fZ)>3)
        continue;
      new((*seltracks)[nseltracks]) MyPart(*part1);
      ++nseltracks;
    }

    Int_t mBin = hMultBins->FindBin(nseltracks) - 1;
    Int_t zBin = hZvtxBins->FindBin(header->fVz) - 1;
    GenericEventPool *pool = pm->GetEventPool(mBin, zBin);
    if (!pool) {
      continue;
    }

    if (pool->IsReady()) {
      cout << "Working on event " << i << " from " << header->fRun << endl;

      for (Int_t t1=0;t1<nseltracks;++t1) {
        MyPart *part1 = (MyPart*)seltracks->At(t1);
        for (Int_t t2=t1+1;t2<nseltracks;++t2) {
          MyPart *part2 = (MyPart*)seltracks->At(t2);

          Double_t dphi = DeltaPhi(*part1,*part2);
          Double_t deta = DeltaEta(*part1,*part2);
          Double_t dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
          if(dr>1)
            hSig->Fill(deta,dphi);

          Int_t mtracks = nmix;
          for (Int_t t3=0;t3<mtracks;) {
            MyPart *part2 = (MyPart*)pool->GetRandomTrack();
            Double_t dphi = DeltaPhi(*part1,*part2);
            Double_t deta = DeltaEta(*part1,*part2);
            Double_t dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
            if(dr>1)
              hBkg->Fill(deta,dphi);
            ++t3;
          }
        }
      }
    }
    pool->UpdatePool(i,header,seltracks);
    if (!pool->IsReady())
      pool->PrintInfo();
  }

  TCanvas *c1 = new TCanvas("cSig");
  hSig->Draw("surf1");
  TCanvas *c2 = new TCanvas("cBkg");
  hBkg->Draw("surf1");
  TH2F *hMix1=(TH2F*)hSig->Clone("mix1");
  hMix1->Divide(hSig,hBkg,nmix,1);
  TCanvas *c3 = new TCanvas("cMix1");
  hMix1->Draw("surf1");
  TH2F *hMix2=(TH2F*)hSig->Clone("mix2");
  hMix2->Divide(hSig,hBkg,hBkg->Integral(),hSig->Integral());
  TCanvas *c4 = new TCanvas("cMix2");
  hMix2->Draw("surf1");
  TTimeStamp t;
  TString fname(Form("histout-%d.root",(Int_t)t.GetSec()));
  TFile outfile(fname,"recreate");
  hSig->Write();
  hBkg->Write();
  hMix1->Write();
  hMix2->Write();
}

//-----------------------------------------------------------------------------------------------------

Bool_t InBounds(Double_t val, Double_t min, Double_t max)
{
  if (val<min)
    return 0;
  if (val>max)
    return 0;
  return 1;
}

Double_t DeltaPhi(const MyPart &t1, const MyPart &t2,
                  Double_t rangeMin, Double_t rangeMax)
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();
  Double_t phia = t1.fPhi;  
  Double_t phib = t2.fPhi;  
  
  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;
  
  return dphi;
}

Double_t DeltaEta(const MyPart &t1, const MyPart &t2)
{
  return t1.fEta - t2.fEta;
}

Double_t InvMass(const MyPart &p1, const MyPart &p2)
{
  Double_t px1 = p1.Px();
  Double_t py1 = p1.Py();
  Double_t pz1 = p1.Pz();
  Double_t px2 = p2.Px();
  Double_t py2 = p2.Py();
  Double_t pz2 = p2.Pz();

  Double_t pm1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
  Double_t pm2 = TMath::Sqrt(px2*px2+py2*py2+pz1*pz2);
  Double_t p12 = px1*px2+py1*py2+pz1*pz2;
  
  Double_t m = TMath::Sqrt(pm1*pm2-p12);
  return m;
}

