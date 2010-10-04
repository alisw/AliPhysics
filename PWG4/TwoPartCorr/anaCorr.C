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
             Double_t vzmin=-10,
             Double_t vzmax=+10,
             Int_t    vzbins=40,
             Int_t    vcmin=10,
             Int_t    vcmax=500,
             Int_t    psize=5,
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

//-----------------------------------------------------------------------------------------------------

void anaCorr(const char *inFileNames,
             Int_t    nEvents,
             Double_t vzmin,
             Double_t vzmax,
             Int_t    vzbins,
             Int_t    vcmin,
             Int_t    vcmax,
             Int_t    psize,
             Double_t ptmin,
             Double_t ptmax,
             Double_t etamin,
             Double_t etamax)
{
  Bool_t doExclCut = 0;
  Double_t dpmax = 0.03;
  Double_t demax = 0.01;

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

  Int_t multbinsa[] = {1,2,3,5,10,20,30,40,50,60,70,80,90,100,120,150,200,500};
  Int_t multbins=17;
  Double_t *multbinsa2 = new Double_t[multbins];
  for(Int_t i=0;i<multbins;++i) {
    multbinsa2[i] = multbinsa[i];
  }

  TH1F* hMultBins = new TH1F("hMultBins", "Event multiplicity binning", 
			     multbins-1, multbinsa2);
  Double_t *vzbinsa = new Double_t[vzbins+1];
  Double_t vzbinwidth = Double_t(vzmax-vzmin)/vzbins;
  for(Int_t i=0;i<vzbins+1;++i) {
    vzbinsa[i] = vzmin+i*vzbinwidth;
  }
  TH1F* hZvtxBins = new TH1F("hZvtxBins", "Event Z-vertex binning", 
			     vzbins, vzbinsa);

  EventPoolManager *pm = new EventPoolManager;
  pm->InitEventPools(psize, multbins, multbinsa, vzbins, vzbinsa);


  TH2F **hSig = new TH2F*[multbins];
  TH2F **hBkg = new TH2F*[multbins];
  TH1F **hMul = new TH1F*[multbins];
  for (Int_t i = 0; i<multbins; ++i) {
    hSig[i] = new TH2F(Form("hSig%d",i),";#Delta#eta;#Delta#phi",60,-3,3,64,-TMath::Pi()/2,3*TMath::Pi()/2);
    hBkg[i] = new TH2F(Form("hBgk%d",i),"#;Delta#eta;#Delta#phi",60,-3,3,64,-TMath::Pi()/2,3*TMath::Pi()/2);
    hMul[i] = new TH1F(Form("hMul%d",i),"#;sel tracks",500,0,500);
  }

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
        (header->fVc>vcmax) )
      continue;

    tBranch->GetEntry(li);
    Int_t ntracks = tracks->GetEntries();
    TClonesArray *seltracks = new TClonesArray("MyPart",ntracks);
    Int_t nseltracks = 0;
    for (Int_t t=0;t<ntracks;++t) {
      MyPart *part = (MyPart*)tracks->At(t);
      if (part->IsITSRefit() && !part->IsTPCRefit())
        continue;
      if (!InBounds(part->fPt,ptmin,ptmax))
        continue;
      if (!InBounds(part->fEta,etamin,etamax))
        continue;
      if (TMath::Abs(part->fD)>0.3)
        continue;
      if (TMath::Abs(part->fZ)>0.3)
        continue;
      new((*seltracks)[nseltracks]) MyPart(*part);
      ++nseltracks;
    }

    Int_t mBin = hMultBins->FindBin(nseltracks) - 1;
    Int_t zBin = hZvtxBins->FindBin(header->fVz) - 1;
    GenericEventPool *pool = pm->GetEventPool(mBin, zBin);
    if (!pool) {
      continue;
    }

    if (pool->IsReady()) {
      if (i%100==0)
        cout << "Working on event " << i << " from " << header->fRun << endl;

      hMul[mBin]->Fill(nseltracks);

      Int_t sigpairs = 0; // count pairs
      for (Int_t t1=0;t1<nseltracks;++t1) {
        MyPart *part1 = (MyPart*)seltracks->At(t1);
        for (Int_t t2=t1+1;t2<nseltracks;++t2) {
          MyPart *part2 = (MyPart*)seltracks->At(t2);
          Double_t dphi = DeltaPhi(*part1,*part2);
          Double_t deta = DeltaEta(*part1,*part2);
          Double_t dr = 1e10;
          if (doExclCut)
            dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
          if(dr>1) {
            ++sigpairs;
          }
        }
      }
      Double_t weight = 1./sigpairs;
      for (Int_t t1=0;t1<nseltracks;++t1) {
        MyPart *part1 = (MyPart*)seltracks->At(t1);
        for (Int_t t2=t1+1;t2<nseltracks;++t2) {
          MyPart *part2 = (MyPart*)seltracks->At(t2);
          Double_t dphi = DeltaPhi(*part1,*part2);
          Double_t deta = DeltaEta(*part1,*part2);
          Double_t dr = 1e10;
          if (doExclCut)
            dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
          if(dr>1) {
            ++sigpairs;
            hSig[mBin]->Fill(deta,dphi,weight);
          }
        }
      }

      TObjArray *mseltracks = pool->GetRandomEvent();
      Int_t nmseltracks = mseltracks->GetEntries();
      Int_t mixpairs = 0; // count mix pairs
      for (Int_t t1=0;t1<nseltracks;++t1) {
        MyPart *part1 = (MyPart*)seltracks->At(t1);
        for (Int_t t2=1;t2<nmseltracks;++t2) {
          MyPart *part2 = (MyPart*)mseltracks->At(t2);
          Double_t dphi = DeltaPhi(*part1,*part2);
          Double_t deta = DeltaEta(*part1,*part2);
          Double_t dr = 1e10;
          if (doExclCut)
            dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
          if(dr>1) {
            ++mixpairs;
          }
        }
      }
      Double_t weight2 = 1./mixpairs;
      for (Int_t t1=0;t1<nseltracks;++t1) {
        MyPart *part1 = (MyPart*)seltracks->At(t1);
        for (Int_t t2=1;t2<nmseltracks;++t2) {
          MyPart *part2 = (MyPart*)mseltracks->At(t2);
          Double_t dphi = DeltaPhi(*part1,*part2);
          Double_t deta = DeltaEta(*part1,*part2);
          Double_t dr = 1e10;
          if (doExclCut)
            dr = dphi*dphi/(dpmax*dpmax) + deta*deta/(demax*demax);
          if(dr>1) {
            ++sigpairs;
            hBkg[mBin]->Fill(deta,dphi,weight2);
          }
        }
      }
    }
    pool->UpdatePool(i,header,seltracks);
    if (!pool->IsReady())
      pool->PrintInfo();
  }

  TTimeStamp t;
  TString fname(Form("histout-%d.root",(Int_t)t.GetSec()));
  TFile outfile(fname,"recreate");
  for (Int_t i = 0; i<multbins; ++i) {
    hSig[i]->Write();
    hBkg[i]->Write();
    hMul[i]->Write();
    TH2F *hMix=(TH2F*)hSig[i]->Clone(Form("hMix%d",i));
    hMix->Divide(hSig[i],hBkg[i],1.,1.);
    hMix->Write();
    TH2F *hMix2=(TH2F*)hSig[i]->Clone(Form("hMix2%d",i));
    for(Int_t x=1;x<=hMix->GetNbinsX();++x) {
      for(Int_t y=1;y<=hMix->GetNbinsY();++y) {
        Int_t bin = hMix->GetBin(x,y);
        Double_t val = hMix->GetBinContent(bin);
        val = hMul[i]->GetMean()*(val-1);
        hMix2->SetBinContent(bin,val);
      }
    }
    hMix2->Write();
  }
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

