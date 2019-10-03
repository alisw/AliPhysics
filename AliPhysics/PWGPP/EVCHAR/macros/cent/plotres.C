#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TParticle.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TColor.h>
#include <TLegend.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TTimeStamp.h>
#endif

TObjArray *gList = 0;
TTree     *gGTree = 0;

// from PPR
const Int_t nclasses = 6;
Double_t bmin[nclasses] = {0,3,6,9,12,15};
Double_t bmax[nclasses] = {3,6,9,12,15,18};
Double_t fxs[nclasses];

// analysis cuts
const int nclassesan = 11;
Double_t bminan[nclassesan];
Double_t bmaxan[nclassesan];
Double_t fxsan[nclassesan];
Double_t fxsan1[nclassesan] = {5,5,10,10,10,10,10,10,10,10,10};
Double_t npminan[nclassesan];
Double_t npmaxan[nclassesan];
Double_t npfxsan[nclassesan];
const char *lan[nclassesan] = 
  {"0-5","5-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100"};
Int_t colorcl[nclassesan] = 
  {kYellow-9,kYellow,kOrange-4,kOrange+6,kOrange+8,kRed,kRed+1,kRed+2,kMagenta+3,kBlue+3,kBlue+4};
Double_t npmean[nclassesan];
Double_t nprms[nclassesan];

TCanvas *Canvas(const char *name, const char *title=0, Int_t ww=600, Int_t wh=400);
TH1 *Hist(const char *hname, const char *name, const char *title, const char *xtitle, 
          const char *ytitle, Bool_t stats=0, Int_t lc=0, Int_t ls=0);
TH1 *Hist(TH1 *h, const char *name, const char *title, const char *xtitle, 
          const char *ytitle, Bool_t stats=0, Int_t lc=0, Int_t ls=0);
TH1 *Hist(const char *name, const char *title, const char *xtitle, const char *ytitle, 
          Bool_t stats=0, Int_t lc=1, Int_t ls=1);
TObjArray *Draw(const char *expr, const char *sel=0, const char *opt=0, Int_t type=1,
                const char *scl=0, Double_t *smin=0, Double_t *smax=0, TTree *t=0);
void Classes(TH1 *h, Double_t *resmin, Double_t *resmax, Double_t *fxs, Int_t pos=1, Int_t verbose=0);
TCanvas *FitNpartDists(const char *name, TObjArray *arr, Int_t verbose=0);
void Store(TObject *o, const char *name=0, const char *fname="gres");

// macro starts here

void glauber() 
{
  Bool_t pImDist               = 0;
  Bool_t pNpDist               = 0;
  Bool_t pNpDistSelWithImp     = 1;
  Bool_t pNpDistSelWithImpFits = 0;
  Bool_t pMidRecResStudy       = 1;

  Double_t fwdres=0.50;
  Double_t midres=0.25;

  gStyle->SetOptFit(1);
  gROOT->ForceStyle();
  TH1::SetDefaultSumw2(1);

  if (gList) 
    delete gList;
  gList = new TObjArray;
  gList->SetOwner(1);

  TFile *f = TFile::Open("hj-unquenched.root");
  TTree *t = (TTree*)f->Get("glaubertree");
  if (!t) {
    cerr << " not find glaubertree" <<endl;
    return;
  }
  if (gGTree)
    delete gGTree;
  gGTree = t;

  if (1) {
    TNtuple *nt = new TNtuple("nt","nt","g1:g2:g3");
    nt->SetDirectory(0);
    Int_t nents = t->GetEntries();
    for (Int_t i=0;i<nents;++i) {
      nt->Fill(gRandom->Gaus(),gRandom->Gaus(),gRandom->Gaus());
    }
    t->AddFriend(nt,"nt");
  }

  t->SetAlias("Etmidn","response.fEtch0n");
  t->SetAlias("Etmidp","response.fEtch0p");
  t->SetAlias("Etfwdn","response.fEtch3n+response.fEtch4n");
  t->SetAlias("Etfwdp","response.fEtch3p+response.fEtch4p");
  t->SetAlias("Nmidn","response.fNch0n");
  t->SetAlias("Nmidp","response.fNch0p");
  t->SetAlias("Nfwdn","response.fNch3n+response.fNch4n");
  t->SetAlias("Nfwdp","response.fNch3p+response.fNch4p");
  t->SetAlias("Nmid","(Nmidn+Nmidp)/2.");

  t->SetAlias("Etfwdnres",Form("Etfwdn*(1+%f*nt.g1)",fwdres));
  t->SetAlias("Etfwdpres",Form("Etfwdp*(1+%f*nt.g2)",fwdres));
  t->SetAlias("Nmidrec",Form("Nmid*(1+%f*nt.g3)",midres));
  t->SetAlias("npart","header.fNT+header.fNP");
  t->SetAlias("ncoll","header.fN00+header.fN01+header.fN10+header.fN11");
  t->SetAlias("bb","header.fBB");

  t->SetAlias("tresp","1+npart*0.");
  //t->SetAlias("tresp","1-exp(-npart/2.)");
  t->SetAlias("trig","1+Nmid*0.");
  //t->SetAlias("trig","rndm<tresp&&rndm<tresp&&Etfwdnres>2.5&&Etfwdpres>2.5");

  TString name;
  if (1) {
    name="ImpactDistFine";
    t->Draw("bb>>htemp(3000,0,30)","1","goff");
    TH1F *hbb=(TH1F*)Hist(name,"","impact parameter [fm]","counts per bin");

    Double_t totxs = 0;
    for (Int_t i=0;i<nclasses;++i) {
      fxs[i] = 100. * hbb->Integral(hbb->FindBin(bmin[i]), hbb->FindBin(bmax[i])) / hbb->Integral();
      totxs += fxs[i];
      printf("Class PPR %d: %.1f - %.1f -> %.1f\n",i+1,bmin[i],bmax[i],fxs[i]);
    }
    cout << "Total from PPR: " << totxs << endl;

    Classes(hbb, bminan, bmaxan, fxsan, 1, 0);
    totxs = 0;
    for (Int_t i=0;i<nclassesan;++i) {
      printf("Centrality Class %d: %.2ffm - %.2ffm -> %.2f (%s)\n",i+1,bminan[i],bmaxan[i],fxsan[i],lan[i]);
      totxs+=fxsan[i];
    }
    cout << "Total: " << totxs << endl;

    if (pImDist) {
      name="ImpactDist";
      TCanvas *c = Canvas(name);
      t->Draw("bb>>htemp(440,0,22)");
      Hist(name,"","impact parameter [fm]","counts per bin");
      TObjArray *arr=Draw("bb>>htemp(440,0,22)",0,"hist");
      TLegend *l = dynamic_cast<TLegend*>(arr->At(0));
      if (l) {
        l->SetX1(0.15); l->SetX2(0.35); l->Draw();
      }
      Store(c);
    }
  }
  if (1) {
    name="NpartDistFine";
    t->Draw("npart>>htemp(440,0,440)","1","goff");
    TH1F *hnp=(TH1F*)Hist(name,"","number participants","counts per bin");

    Classes(hnp, npminan, npmaxan, npfxsan, -1, 0);
    Double_t totxs = 0;
    for (Int_t i=0;i<nclassesan;++i) {
      printf("Npart Class %d: %.1f - %.1f -> %.1f (%s)\n",i+1,npminan[i],npmaxan[i],npfxsan[i],lan[i]);
      totxs+=npfxsan[i];
    }
    cout << "Total: " << totxs << endl;
    if (pNpDist) {
      name="NpartDist";
      TCanvas *c1=Canvas(name);
      c1->SetLogy(1);
      t->Draw("npart>>htemp(440,0,440)","1","goff");
      TH1 *h=Hist(name,"","number participants","counts per bin");
      h->SetMinimum(1);
      h->Draw();
      /*TObjArray *arr=*/Draw("npart>>htemp(440,0,440)",0,"hist",-2);
      Store(c1);
    }
    if (pNpDistSelWithImp) {
      name="NpartDistsWithImpact";
      TCanvas *c1=Canvas(name);
      c1->SetLogy(1);
      t->Draw("npart>>htemp(440,0,440)","1","goff");
      TH1 *h=Hist(name,"","number participants","counts per bin");
      h->SetMinimum(1);
      h->Draw();
      TObjArray *arr=Draw("npart>>htemp(440,0,440)",0,"hist",-1);
      TGraph *ge1 = new TGraph(nclassesan);
      TGraph *ge2 = new TGraph(nclassesan);
      ge1->SetMarkerSize(1.2);
      ge1->SetMarkerStyle(20);
      ge2->SetMarkerSize(1.2);
      ge2->SetMarkerStyle(20);
      for (Int_t i=1;i<arr->GetEntries();++i) {
        TH1F *h = (TH1F*)arr->At(i);
        Int_t N=nclassesan-i;
        Double_t mean  = h->GetMean();
        Double_t width = h->GetRMS();
        npmean[N] = mean;
        nprms[N]  = width;
        ge1->SetPoint(N,N,mean);
        ge2->SetPoint(N,N,width);
      }
      Store(c1);
      Store(ge1,"gNpartMean");
      Store(ge2,"gNpartRms");
      if (pNpDistSelWithImpFits) {
        TCanvas *c2 = FitNpartDists("NpartDistsWithImpactFits",arr,1);
        Store(c2);
      }
    }
  }
  if (pMidRecResStudy) {
    Int_t resint = 100*midres;
    name=Form("MidRecDistribution_res%d",resint);
    t->Draw("Nmidrec>>htemp(3000,0,3000)","1","goff");
    TH1 *hNm = Hist(name,"","Nch in -0.5<#eta<0.5","counts per bin");
    Double_t resmin[nclassesan];
    Double_t resmax[nclassesan];
    Double_t fxs[nclassesan];
    Classes(hNm,resmin,resmax,fxs,-1,1);
    if (1) {
      TCanvas *c = Canvas(name);
      hNm->SetMinimum(1);
      hNm->Draw();
      Draw("Nmidrec>>htemp(3000,0,3000)",0,"hist",-3,"Nmidrec",resmin,resmax);
      Store(c);
    } else {
      delete hNm;
    }

    name=Form("NpartDistsWithTracks_res%d",resint);
    TCanvas *c1 = Canvas(name);
    t->Draw("npart>>htemp(440,0,440)","1","goff");
    TH1 *hNpart = Hist(name,"","number participants","counts per bin");
    hNpart->SetMinimum(1);
    hNpart->Draw();
    TObjArray *arr=Draw("npart>>htemp(440,0,440)",0,"hist",-3,"Nmidrec",resmin,resmax);
    Store(c1);
    if (1) {
      TCanvas *c2 = FitNpartDists(Form("NpartDistsWithTracksFits_res%d",resint),arr,1);
      Store(c2);
    }
    TGraph *ge1 = new TGraph(nclassesan);
    TGraph *ge2 = new TGraph(nclassesan);
    ge1->SetMarkerSize(1.2);
    ge1->SetMarkerStyle(20);
    ge2->SetMarkerSize(1.2);
    ge2->SetMarkerStyle(20);
    for (Int_t i=1;i<arr->GetEntries();++i) {
      Int_t N=nclassesan-i;
      TH1F *h = (TH1F*)arr->At(i);
      Double_t mean  = h->GetMean();
      Double_t rms = h->GetRMS();
      ge1->SetPoint(N,N,mean/npmean[N]-1);
      ge2->SetPoint(N,N,rms/mean*npmean[N]/nprms[N]-1);
      cout << i << " " << mean << " " << rms << " " << npmean[N] << " " << nprms[N] << endl;
    }
    Store(ge1,Form("gMidrecMean_res%d",resint));
    Store(ge2,Form("gMidrecRms_res%d",resint));
  }
  if (0) {
    name="FwdSumCorr";
    Canvas(name);
    t->Draw("Etfwdn:Etfwdp");
    Hist(name,"","sum p_{t} [GeV] in 3<#eta<5","sum p_{t} [GeV] in -3<#eta<-5");
  }
  if (0) {
    name="FwdNvsNpart";
    Canvas(name);
    t->Draw("Etfwdn:npart","1","prof");
    Hist(name,"","Npart","sum p_{t} [GeV] in -3<#eta<-5");
  }
  if (0) {
    name="FwdPvsNpart";
    Canvas(name);
    t->Draw("Etfwdp:npart","1","prof");
    Hist(name,"","Npart","sum p_{t} [GeV] in 3<#eta<5");
  }
  if (0) {
    name="FwdSumvsNpart";
    Canvas(name);
    t->Draw("Etfwdn+Etfwdp:npart","1","prof");
    Hist(name,"","Npart","sum p_{t} [GeV] in 3<|#eta|<5");
  }
  if (0) {
    name="FwdSumDistribution";
    Canvas(name);
    t->Draw("Etfwdn+Etfwdp","1","");
    TH1 *h1=Hist(name,"","sum p_{t} [GeV] in 3<|#eta|<5","counts per bin");
    name="FwdSumTriggered";
    t->Draw("Etfwdn+Etfwdp","trig>0","same");
    TH1 *h2=Hist(name,"","sum p_{t} [GeV] in 3<|#eta|<5","counts per bin",0,2,1);
    Double_t percent = h2->Integral()/h1->Integral()*100.;
    cout << "Recorded " << percent << "% of tot. cross section" << endl;
  }
  if (0) {
    name="Mid2Distribution";
    Canvas(name);
    t->Draw("Nmidn+Nmidp","1","");
    Hist(name,"","Nch in -1<#eta<1","counts per bin");
    name="Mid2Triggered";
    t->Draw("Nmidn+Nmidp","trig>0","same");
    Hist(name,"","Nch in -1<#eta<1","counts per bin",0,2,1);
  }
  if (0) {
    name="MidDistribution";
    Canvas(name);
    t->Draw("Nmid","1","");
    Hist(name,"","Nch in -0.5<#eta<0.5","counts per bin");
    name="MidTriggered";
    t->Draw("Nmid","trig>0","same");
    Hist(name,"","Nch in -0.5<#eta<0.5","counts per bin",0,2,1);
  }
}

//--------------------------------------------------------------------------------------------------

TCanvas *Canvas(const char *name, const char *title, Int_t ww, Int_t wh)
{
  if (!name)
    return 0;
  TString hname(Form("c%s",name));
  if (!title)
    title = name;
  TCanvas *c = new TCanvas(hname,title,ww,wh);
  return c;
}

//--------------------------------------------------------------------------------------------------

TH1 *Hist(const char *name, const char *title, const char *xtitle, const char *ytitle, 
          Bool_t stats, Int_t lc, Int_t ls)
{
  TH1 *h=dynamic_cast<TH1*>(gROOT->FindObject("htemp"));
  return Hist(h,name,title,xtitle,ytitle,stats,lc,ls);
}

//--------------------------------------------------------------------------------------------------

TH1 *Hist(const char *hname, const char *name, const char *title, const char *xtitle, 
          const char *ytitle, Bool_t stats, Int_t lc, Int_t ls)
{
  TH1 *h=dynamic_cast<TH1*>(gROOT->FindObject(hname));
  return Hist(h,name,title,xtitle,ytitle,stats,lc,ls);
}

//--------------------------------------------------------------------------------------------------

TH1 *Hist(TH1 *h, const char *name, const char *title, const char *xtitle, 
          const char *ytitle, Bool_t stats, Int_t lc, Int_t ls)
{
  if (!h && !name) 
    return 0;
  TString hname(Form("h%s",name));
  h->SetName(hname);
  if (!title)
    title = name;
  h->SetTitle(title);
  if (xtitle)
    h->SetXTitle(xtitle);
  if (ytitle)
    h->SetYTitle(ytitle);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->SetStats(stats);
  h->SetLineColor(lc);
  h->SetLineStyle(ls);
  gList->Add(h);
  return h;
}

//--------------------------------------------------------------------------------------------------

TObjArray *Draw(const char *expr, const char *sel, const char *opt, Int_t type,
                const char *scl, Double_t *smin, Double_t *smax, TTree *t)
{
  if (!t) 
    t = gGTree;
  
  TObjArray *oarr = new TObjArray;
  oarr->SetOwner();

  TLegend *l = new TLegend(0.2,0.4,0.4,0.9);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  oarr->Add(l);

  TString tsel("1");
  if (sel)
    tsel=sel;
  TString doopt("same");
  if (opt)
    doopt=Form("%s,same",opt);

  Int_t beg=0;
  Int_t end=nclassesan;
  if (type<0) {
    beg=nclassesan-1;
    end=-1;
  }
  Double_t *cmin = bminan;
  Double_t *cmax = bmaxan;
  const char *varname="bb";
  if (TMath::Abs(type)==2) {
    cmin = npminan;
    cmax = npmaxan;
    varname="npart";
  } else if (TMath::Abs(type)>2) {
    cmin    = smin;
    cmax    = smax;
    varname = scl;
  }

  while(beg!=end) {
    Int_t i=beg;
    TString dosel(Form("(%s>%.2f&&%s<%.2f)&&(%s)",varname,cmin[i],varname,cmax[i],tsel.Data()));
    t->Draw(expr,dosel,doopt);
    TH1 *h=dynamic_cast<TH1*>(gROOT->FindObject("htemp"));
    if (h) {
      h->SetName(Form("cl_%s_%s",expr,dosel.Data()));
      h->SetLineColor(colorcl[i]);
      h->SetMarkerColor(colorcl[i]);
      h->SetFillColor(colorcl[i]);
      h->SetFillStyle(1000);
      l->AddEntry(h,Form("%s%%",lan[i]),"f");
      oarr->Add(h);
    } else {
      cerr << "Could not obtain htemp for: " << expr << " " << dosel << " " << doopt << endl;
    }
    if(type<0)
      --beg;
    else
      ++beg;
  }
  return oarr;
}

//--------------------------------------------------------------------------------------------------

void Classes(TH1 *h, Double_t *resmin, Double_t *resmax, Double_t *fxs, Int_t pos, Int_t verbose)
{
  Int_t bfrom = 0;
  Int_t cbin  = h->GetNbinsX();
  if (pos<0) {
    pos   = -1;
    bfrom = h->GetNbinsX()+1;
    cbin  = 1;
  } else {
    pos = 1;
  }

  Double_t totxs=0;
  for (Int_t i=0,lbin=bfrom,bin=lbin;i<nclassesan;++i) {
    lbin = bin+pos;
    Int_t lxs = 0;
    Double_t norm = h->Integral();
    while (1) {
      bin += pos;
      lxs += h->GetBinContent(bin);
      Double_t pxs = lxs/norm*100;
      Double_t tdiff = (lxs+h->GetBinContent(bin+pos))/norm*100;
      //cout << pos << " " << bin << " " << pxs << " " << tdiff << " " << lbin << endl;
      if ((pxs>1&&TMath::Abs(pxs-fxsan1[i])<=TMath::Abs(tdiff-fxsan1[i]))||(bin==cbin)) {
        if (pos>0) {
          resmin[i] = h->GetBinLowEdge(lbin);
          resmax[i] = h->GetBinLowEdge(bin+1);
        } else {
          resmin[i] = h->GetBinLowEdge(bin);
          resmax[i] = h->GetBinLowEdge(lbin+1);
        }
        fxs[i] = pxs;
        if (verbose)
          printf("Class %d: %.1f - %.1f -> %.1f (%s)\n",i+1,resmin[i],resmax[i],pxs,lan[i]);
        totxs += pxs;
        break;
      }
    }
  }
  if (verbose)
    cout << "Total: " << totxs << endl;
}

//--------------------------------------------------------------------------------------------------

TCanvas *FitNpartDists(const char *name, TObjArray *arr, Int_t verbose)
{
  TCanvas *c = Canvas(name,0,1200,900);
  c->Divide(4,3);
  TGraphErrors *ge = new TGraphErrors(nclassesan);
  for (Int_t i=1;i<arr->GetEntries();++i) {
    c->cd(i);
    Int_t N=nclassesan-i;
    TH1F *h = (TH1F*)arr->At(i);
    Double_t mean  = h->GetMean();
    Double_t width = h->GetRMS();
    TH1 *h1 = h->DrawCopy("hist");
    h1->SetName(lan[N]);
    h1->SetTitle(Form("%s%% most central",lan[N]));
    h1->SetStats(1);
    h1->SetXTitle("Npart");
    h1->SetYTitle("counts per bin");
    TF1 *fit = new TF1(Form("fit%d",i),"gaus(0)",0,440);
    fit->SetParameters(1,mean,width);
    fit->SetLineWidth(3);
    h1->Fit(fit,"QM0");
    fit->Draw("same");
    ge->SetPoint(N,mean,width);
    //ge->SetPointError(N,0,width);
    if (verbose) 
      cout << i << " hist: " << mean << " " << width 
           << " fit:  " << fit->GetParameter(1) << " " << fit->GetParameter(2) << endl;
  }
  c->cd(12);
  TH2F *h2f = new TH2F(Form("%s_summary",name),";#LTNpart#GT;RMS",1,0,440,1,0,49);
  h2f->SetStats(0);
  h2f->Draw();
  ge->SetMarkerSize(1.2);
  ge->SetMarkerStyle(20);
  ge->Draw("P");
  return c;
}

//--------------------------------------------------------------------------------------------------

void Store(TObject *o, const char *name, const char *fname)
{
  if (!o || !fname)
    return;

  TTimeStamp t;
  TString filename(Form("%s_%d.root",fname, t.GetDate()));
  TFile *f = new TFile(filename,"update");
  if (!name)
    name = o->GetName();
  if (o) {
    f->Delete(Form("%s;1",name));
    o->Write(name);
  }

  f->Close();
  delete f;
}

void plotRes(const char *fname)
{
  TFile *f = new TFile(fname);
  
  TGraph *gm0 = (TGraph*)f->Get("gMidrecMean_res0");
  TGraph *gm5 = (TGraph*)f->Get("gMidrecMean_res5");
  TGraph *gm10 = (TGraph*)f->Get("gMidrecMean_res10");
  TGraph *gm15 = (TGraph*)f->Get("gMidrecMean_res14");
  TGraph *gm25 = (TGraph*)f->Get("gMidrecMean_res25");

  gm0->SetMarkerColor(1);
  gm5->SetMarkerColor(2);
  gm10->SetMarkerColor(3);
  gm15->SetMarkerColor(4);
  gm25->SetMarkerColor(5);

  TH2F *h2f = new TH2F("h2f",";classes (0=most central);#LTNpart#GT/#LTNpart^{true}#GT-1",1,-0.5,10.5,1,-0.25,0.25);
  h2f->SetStats(0);
  h2f->GetXaxis()->SetTitleOffset(1.1);
  h2f->GetYaxis()->SetTitleOffset(1.2);
  h2f->Draw();
  gm0->Draw("P");
  gm5->Draw("P");
  gm10->Draw("P");
  gm15->Draw("P");
  gm25->Draw("P");
  TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
  leg->AddEntry(gm0,"0% smearing","p");
  leg->AddEntry(gm5,"5% smearing","p");
  leg->AddEntry(gm10,"10% smearing","p");
  leg->AddEntry(gm15,"15% smearing","p");
  leg->AddEntry(gm25,"25% smearing","p");
  leg->Draw();
}
