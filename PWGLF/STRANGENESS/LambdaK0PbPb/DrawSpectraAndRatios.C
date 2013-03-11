#if !defined(__CINT__) || defined(__MAKECINT__)
   #include <TMath.h>
   #include <TROOT.h>
   #include <Riostream.h>
   #include <TCanvas.h>
   #include <TColor.h>
   #include <TLatex.h>
   #include <TLegend.h>
   #include <TLegendEntry.h>

   #include <TStyle.h>
   #include <TString.h>
   #include <TASImage.h>

   #include <TFile.h>
   #include <TList.h>
   #include <TH1F.h>
   #include <TH1D.h>
   #include <TF2.h>
   #include <TFitResult.h>
   #include <TFitResultPtr.h>
   #include <TH2F.h>
   #include <TH3F.h>
#endif

extern TStyle *gStyle;

static Double_t xBins[]={
   0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
   1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
   2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,
   4.5,5.0,5.5,6.5,8.0,10.0,12.0
};
const Int_t nBins=sizeof(xBins)/sizeof(Double_t) - 1; //37

//*** The systematic uncertainties for combining
// cos(PA)
// DCA between V0 daughters
// TPC crossed pad rows
// DCA daughters <-> PV
// c*tau
static 
Double_t sysEffK0s[nBins]={//Efficiency, combined over cuts mentioned above
  0.0,
  0.05,0.05,0.04,0.04,  //Dominated by cos(PA)
  0.04,0.04,0.04,0.04,0.04,0.04,
  0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,
  0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,
  0.04,0.04,0.04,0.04,0.04
};
/*
static Double_t sysSigK0s[nBins]={//Signal extraction
  0.00728589, 0.00728589, 0.00728539, 0.0073469, 0.00737846, 0.00741705,
  0.00750887, 0.00753641, 0.00769012, 0.00789154, 0.00796624, 0.00822856,
  0.00838203, 0.00864603, 0.00906498, 0.00923208, 0.00931179, 0.0100081,
  0.0100768, 0.0105292, 0.0112067, 0.0122143, 0.0130162, 0.0139799, 0.0155215,
  0.0158903, 0.0168517, 0.0189316, 0.0188103, 0.020233, 0.0223704, 0.0260327,
  0.0260327, 0.0260327, 0.0260327, 0.0260327, 0.0260327
};
*/
static Double_t sysSigK0s[nBins]={//Signal extraction
  0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
  0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
  0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
  0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
  0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
  0.03, 0.03, 0.03, 0.03, 0.03
};

static 
Double_t sysEffLam[nBins]={//Efficiency, combined over cuts mentioned above
  0.0,0.0,0.0,
  0.20,0.12,0.05, //Dominated by cos(PA)
  0.06,0.06,0.06,0.06,0.06,  //Dominated by c*tau
  0.05,0.05, //Dominated by cos(PA)
  0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
  0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
  0.05,0.05,0.05
};
static Double_t sysSigLam[nBins]={//Signal extraction
  0.0447214, 0.0447214, 0.0447214, 0.0447214, 0.0447214,
  0.0447214, 0.0447214, 0.0447214, 0.0447214, 0.0447214,
  0.0447214, 0.0447214, 0.0447214, 0.0447214, 0.0447214,
  0.0447214, 0.0447214, 0.0447214, 0.0447214, 0.0447214,
  0.0447214, 0.0447214, 0.0447214, 0.0447214, 0.0447214,
  0.0447214, 0.0447214, 0.0447214, 0.0447214, 0.0465862,
  0.0472600, 0.0488576, 0.0526797, 0.0575355, 0.0651647,
  0.0721110, 0.072111
};

const Double_t sysPID=0.02; //PID
const Double_t sysArm=0.01; //Armenteros cut
const Double_t sysFD=0.05;  //Feed down

const Double_t sysRatio=0.05;//Efficiency systematics for the L/K ratio

Double_t fd(Double_t x) {
  //Effective FD correction
  return 0.1619 + 0.05295*x - 0.01749*x*x + 0.001425*x*x*x - 3.446e-05*x*x*x*x;
}

void FeedDown(TH1 *spe) {
  for (Int_t i=1; i<=spe->GetNbinsX(); i++) {
      Double_t pt=spe->GetBinCenter(i);
      Double_t c=spe->GetBinContent(i);
      c -= (c*fd(pt));
      spe->SetBinContent(i,c);
  }
}

TH1 *MapHisto(const TH1 *h) {
  const Double_t eps=0.0001;
  TString name("m");
  name = name + h->GetName(); 
  TH1F *mh=new TH1F(name.Data(),h->GetTitle(),nBins,xBins);
  mh->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  mh->GetYaxis()->SetTitle("1/N_{ev}d^{2}N/dp_{t}/dy (GeV/c)^{-1}");

  Double_t xh=h->GetBinCenter(1), xmh=0.;

  Int_t n=1;
  for (; n<=nBins; n++) {
    xmh=mh->GetBinCenter(n);
    if (TMath::Abs(xh-xmh)<eps) break; 
  }

  Int_t iii=h->GetNbinsX();
  for (Int_t i=1; i<=iii; i++) {
    Int_t ni1=n+i-1;

    if (ni1 > nBins) {
      cerr<<"Input number of bins is larger than the output number of bins!\n";
      delete mh;
      return 0;
    }

    xh = h->GetBinCenter(i);
    xmh=mh->GetBinCenter(ni1);
    if (TMath::Abs(xh-xmh)>eps) {
      cerr<<"Wrong binning !\n";
      delete mh;
      return 0;
    }

    Double_t c=h->GetBinContent(i);
    Double_t e=h->GetBinError(i);
    mh->SetBinContent(ni1,c);
    mh->SetBinError(ni1,e);
  }
  
  return mh;
}

Bool_t 
GetHistos(const Char_t *rName[], const Char_t *eName[], TH1 *&raw, TH1 *&eff) {

  /*TFile *fr=*/TFile::Open(rName[0]);
  TList *lst=(TList*)gFile->Get("c1DataYields");
  //TList *lst=(TList*)gFile->Get("cLK0Spectra");

  raw=(TH1F*)lst->FindObject(rName[1]);
  if (!raw) {
     cerr<<"No raw yield !"<<eName[0]<<' '<<eName[1]<<endl; 
     return kFALSE;
  }

  /*TFile *fe=*/TFile::Open(eName[0]);
  eff=(TH1F*)gFile->Get(eName[1]);
  //eff=(TH1F*)lst->FindObject(eName[1]);
  if (!eff) {
     cerr<<"No efficiency ! "<<eName[0]<<' '<<eName[1]<<endl; 
     return kFALSE;
  }
  return kTRUE;
}

void SetAttributes(TH1 *h,const Char_t *tit,Int_t col,Int_t mar,Float_t siz,
Float_t min=1e-7, Float_t max=1000., Int_t range=nBins) {
  h->SetTitle(tit);
  h->SetLineColor(col); 
  h->SetMarkerColor(col);
  h->SetMarkerStyle(mar);
  h->SetMarkerSize(siz);
  h->SetMaximum(max);
  h->SetMinimum(min);
  h->GetXaxis()->SetRange(1,range);
}

void 
DrawHisto(const TH1 *h, const Option_t *option, Double_t *sysEff, 
	  Double_t *sysSig, Double_t scale=1) {
  TH1F *hh=new TH1F(*((TH1F*)h));
  Int_t nb=hh->GetNbinsX();

  for (Int_t i=1; i<=nb; i++) {
      Double_t c=hh->GetBinContent(i);
      Double_t e=hh->GetBinError(i);
      Int_t j=i-1;
      e = sysEff[j]*sysEff[j] + sysSig[j]*sysSig[j];

      if (sysEff==sysEffLam) {// for Lambda
	 e += sysFD*sysFD;
         if (i<13) e += sysPID*sysPID;
      } else {// for K0s
         e += sysArm*sysArm;
      }

      e=c*TMath::Sqrt(e); 

      hh->SetBinError(i,e);
  }
  hh->SetFillColor(17);
  //TString opt("E5"); opt+=option;
  TString opt("E2"); opt+=option;
  /*
  TFile *f=TFile::Open("systematics.root","update");
  hh->Write();
  f->Close();
  */
  hh->Scale(scale);
  hh->SetMinimum(h->GetMinimum());
  hh->SetMaximum(h->GetMaximum());
  hh->Draw(opt.Data());

  TH1F *ch=new TH1F(*((TH1F*)h));
  /*  
  TFile *f=TFile::Open("k0s_lambda_5cm_with_0510.root","update");
  ch->Write();
  f->Close();
  */
  ch->Scale(scale);
  ch->SetMinimum(h->GetMinimum());
  ch->SetMaximum(h->GetMaximum());
  ch->Draw("e x0 same");
  //opt="e x0";  opt+=option;
  //ch->Draw(opt.Data());
}

void DrawRatio(TH1 *h, const Option_t *option) {
  TH1F *hh=new TH1F(*((TH1F*)h));
  Int_t nb=hh->GetNbinsX();

  for (Int_t i=1; i<=nb; i++) {
      Double_t c=hh->GetBinContent(i);
      Double_t e=hh->GetBinError(i);
      Int_t j=i-1;
      e = sysSigK0s[j]*sysSigK0s[j] + sysSigLam[j]*sysSigLam[j];

      e += sysRatio*sysRatio;

      e += sysArm*sysArm;
      e += sysFD*sysFD;
      if (i<13) e += sysPID*sysPID;
 
      e=c*TMath::Sqrt(e); 

      hh->SetBinError(i,e);
  }
  hh->SetFillColor(17);
  //TString opt("E5"); opt+=option;
  TString opt("E2"); opt+=option;
  hh->Draw(opt.Data());
  h->Draw("e x0 same");
  //opt="e x0";  opt+=option;
  //h->Draw(opt.Data());
}

void DrawALICELogo(Float_t x1, Float_t y1, Float_t x2, Float_t y2)
{
// Correct for aspect ratio of figure plus aspect ratio of pad.
// Coordinates are NDC!

  x2 = x1 + (y2 - y1)*0.891*gPad->GetCanvas()->GetWindowHeight()*gPad->GetHNDC() / (gPad->GetWNDC() * gPad->GetCanvas()->GetWindowWidth());
  
  TPad *myPadLogo = new TPad("myPadLogo","Pad for ALICE Logo", x1, y1, x2, y2);
  myPadLogo->SetLeftMargin(0);
  myPadLogo->SetTopMargin(0);
  myPadLogo->SetRightMargin(0);
  myPadLogo->SetBottomMargin(0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = 
  new TASImage("alice_logo_preliminary.eps");
  myAliceLogo->Draw("same");
}

void DrawFit(const Char_t *nam[], const Float_t *fac, Int_t n){
  for (Int_t i=0; i<n; i++) {
      const Char_t *name=nam[i];
      Float_t factor=fac[i];
      TF1 *f=(TF1*)gROOT->FindObject(name);
      TH1 *h=f->GetHistogram();
      h->Scale(factor);
      h->SetLineColor(1);
      h->SetLineWidth(1);
      h->SetLineStyle(1);
      h->Draw("Lsame");
  }
}

void DrawSpectraAndRatios() {

  const Int_t nCent=7;

  const Char_t *title[nCent]={
    "0-5 %",
    "5-10 %",
    "10-20 %",
    "20-40 %",
    "40-60 %",
    "60-80 %",
    "80-90 %"
  };
  const Int_t   colour[nCent]={2,   635, 797, 419, 4 , 6,  1  };
  const Int_t   marker[nCent]={22,  29, 34,  21,  23, 33, 20 };
  const Float_t masize[nCent]={1.3, 1.6, 1.3, 1.2, 1.4, 1.8,  1.3};
  const Float_t factor[nCent]={1.0, 1.0, 1.1, 1.5, 3.0,7.5,15.0}; //scale for drawing
  const Float_t factor2[nCent]={1.0, 0.9, 0.9, 1.1, 2.0,5.5,15.0}; //scale for lin lambda drawing
  const Float_t factor1[nCent]={1.0, 1/2., 1/4., 1/8., 1/16., 1/32., 1/64.}; //scale for log drawing
  
  const Char_t *rNameL[2*nCent]={ // file name, histo name
    "raw.root", "YieldLambda_0005", 
    "raw.root", "YieldLambda_0510", 
    "raw.root", "YieldLambda_1020", 
    "raw.root", "YieldLambda_2040", 
    "raw.root", "YieldLambda_4060", 
    "raw.root", "YieldLambda_6080", 
    "raw.root", "YieldLambda_8090" 
  };
  const Char_t *eNameL[2*nCent]={ // file name, histo name
    "eff.root", "eff_Lambda_comb_0005",
    "eff.root", "eff_Lambda_comb_0510",
    "eff.root", "eff_Lambda_comb_1020",
    "eff.root", "eff_Lambda_comb_2040",
    "eff.root", "eff_Lambda_comb_4060",
    "eff.root", "eff_Lambda_comb_6080",
    "eff.root", "eff_Lambda_comb_8090"
  };

  const Char_t *rNameK[2*nCent]={ // file name, histo name
    "raw.root", "YieldK0Short_0005", 
    "raw.root", "YieldK0Short_0510", 
    "raw.root", "YieldK0Short_1020", 
    "raw.root", "YieldK0Short_2040", 
    "raw.root", "YieldK0Short_4060", 
    "raw.root", "YieldK0Short_6080", 
    "raw.root", "YieldK0Short_8090" 
  };
  const Char_t *eNameK[2*nCent]={ // file name, histo name
    "eff.root", "eff_K0s_comb_0005",
    "eff.root", "eff_K0s_comb_0510",
    "eff.root", "eff_K0s_comb_1020",
    "eff.root", "eff_K0s_comb_2040",
    "eff.root", "eff_K0s_comb_4060",
    "eff.root", "eff_K0s_comb_6080",
    "eff.root", "eff_K0s_comb_8090"
  };

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFillColor(0);

  TH1 *raw=0;
  TH1 *eff=0;
  TString option(""), ratio("ratio");

  TCanvas *c1=new TCanvas; c1->SetLogy();
  c1->SetLeftMargin(0.13); c1->SetBottomMargin(0.13);
  TCanvas *c1lin=new TCanvas;
  c1lin->SetLeftMargin(0.13); c1lin->SetBottomMargin(0.13);

  TCanvas *c2=new TCanvas; c2->SetLogy();
  c2->SetLeftMargin(0.13); c2->SetBottomMargin(0.13);
  TCanvas *c2lin=new TCanvas;
  c2lin->SetLeftMargin(0.13); c2lin->SetBottomMargin(0.13);

  TCanvas *c3=new TCanvas;
  c3->SetLeftMargin(0.13); c3->SetBottomMargin(0.13);

  TH1 *lkRatio[nCent]={0};

  for (Int_t cent=0; cent<nCent; cent++) {
      const Char_t *tit=title[cent];
      Int_t col=colour[cent];
      Int_t mar=marker[cent];
      Float_t siz=masize[cent];
 
      // Lambda
      if (!GetHistos(rNameL+2*cent, eNameL+2*cent, raw, eff)) return;
      TH1 *rawHl=MapHisto(raw);
      TH1 *effHl=MapHisto(eff);

      //Feed down
      FeedDown(rawHl);

      rawHl->Divide(effHl);
      SetAttributes(rawHl,tit,col,mar,siz);
      c1->cd();
      DrawHisto(rawHl, option.Data(), sysEffLam, sysSigLam, factor1[cent]);

      TH1 *linHl=(TH1*)rawHl->Clone();
      SetAttributes(linHl,tit,col,mar,siz,0.,20.,32); 
      c1lin->cd();
      DrawHisto(linHl, option.Data(), sysEffLam, sysSigLam, factor2[cent]);

      // K0s
      if (!GetHistos(rNameK+2*cent, eNameK+2*cent, raw, eff)) return;
      TH1 *rawHk=MapHisto(raw);
      TH1 *effHk=MapHisto(eff);
      rawHk->Divide(effHk);
      SetAttributes(rawHk,tit,col,mar,siz,1e-7);
      c2->cd();
      DrawHisto(rawHk, option.Data(), sysEffK0s, sysSigK0s, factor1[cent]);

      TH1 *linHk=(TH1*)rawHk->Clone();
      SetAttributes(linHk,tit,col,mar,siz,0.,120.,32); 
      c2lin->cd();
      DrawHisto(linHk, option.Data(), sysEffK0s, sysSigK0s, factor[cent]);

      // Lambda/K0s
      TH1 *rawHlk=(TH1*)rawHl->Clone();
      lkRatio[cent]=rawHlk;
      TString name=ratio+rawHlk->GetName();
      rawHlk->SetName(name.Data());      
      rawHlk->SetMaximum(1.7);      
      rawHlk->Divide(rawHk);
      rawHlk->GetYaxis()->SetTitle("#Lambda/K^{0}_{S}");
      c3->cd();
      //if (cent!=1)
      DrawRatio(rawHlk,option.Data());

      option+="same";
  }

  for (Int_t cent=0; cent<nCent; cent++) {
    //if (cent != 1) 
    lkRatio[cent]->Draw("same");
  }


  TLegend *leg=c1->BuildLegend(0.68,0.46,0.88,0.82,"Centrality:");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  TLegendEntry *entry=leg->AddEntry("NULL","systematic uncertainty","lpf");
  Int_t ci = TColor::GetColor("#cccccc");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(10);
  entry->SetMarkerColor(ci);

  c1->cd(); 
  TLatex *   tex = new TLatex(0.5,0.65,"#Lambda");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.11);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.27,0.83,"Pb-Pb at #sqrt{s_{NN}}=2.76 TeV, |y|<0.5");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();
  //Float_t offx=0.15, offy=0.16, sizx=0.22, sizy=0.22;
  //DrawALICELogo(offx,offy,offx+sizx,offy+sizy);

   leg=c1lin->BuildLegend(0.69,0.43,0.88,0.80,"Centrality:");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);

   entry=leg->AddEntry("NULL","systematic uncertainty","lpf");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(10);
   entry->SetMarkerColor(ci);
   entry=leg->AddEntry("NULL","BGBW fit","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);

   c1lin->cd();
   tex=new TLatex(0.27,0.83,"Pb-Pb at #sqrt{s_{NN}}=2.76 TeV, |y|<0.5");
   tex->SetNDC();
   tex->Draw();
      tex = new TLatex(0.5,0.65,"#Lambda");  
   tex->SetNDC();
   tex->SetTextSize(0.11);
   tex->Draw();
   {
     TFile::Open(" BWFitResults_Lambda0_stat.root");
     const Char_t *name[nCent]={
       "BWFit_0005",
       "BWFit_0510","BWFit_1020","BWFit_2040",
       "BWFit_4060","BWFit_6080","BWFit_8090"
     };
     DrawFit(name, factor2, nCent);
   }
   //Float_t offx1=0.70, offy1=0.18;
   //DrawALICELogo(offx1,offy1,offx1+sizx,offy1+sizy);
       



  leg=c2->BuildLegend(0.68,0.46,0.88,0.82,"Centrality:");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  entry=leg->AddEntry("NULL","systematic uncertainty","lpf");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(10);
  entry->SetMarkerColor(ci);

  c2->cd(); 
  tex = new TLatex(0.5,0.65,"K^{0}_{S}");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.089);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.27,0.83,"Pb-Pb at #sqrt{s_{NN}}=2.76 TeV, |y|<0.5");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();
  //DrawALICELogo(offx,offy,offx+sizx,offy+sizy);

   leg=c2lin->BuildLegend(0.69,0.43,0.88,0.80,"Centrality:");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);

   entry=leg->AddEntry("NULL","systematic uncertainty","lpf");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(10);
   entry->SetMarkerColor(ci);
   entry=leg->AddEntry("NULL","BGBW fit","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);

   c2lin->cd();
      tex = new TLatex(0.27,0.83,"Pb-Pb at #sqrt{s_{NN}}=2.76 TeV, |y|<0.5");
   tex->SetNDC();
   tex->Draw();
      tex = new TLatex(0.5,0.65,"K^{0}_{S}");  
   tex->SetNDC();
   tex->SetTextSize(0.089);
   tex->Draw();
   {
     TFile::Open("BWFitResults_K0_stat.root");
     const Char_t *name[nCent]={
       "BWFit_0005",
       "BWFit_0510","BWFit_1020","BWFit_2040",
       "BWFit_4060","BWFit_6080","BWFit_8090"
     };
     DrawFit(name, factor, nCent);
    }
   //DrawALICELogo(offx1,offy1,offx1+sizx,offy1+sizy);
       
   
   //leg=c3->BuildLegend(0.74,0.62,0.88,0.88,"Centrality:");
  leg=c3->BuildLegend(0.55,0.55,0.88,0.88,"Centrality:");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  entry=leg->AddEntry("NULL","systematic uncertainty","lpf");
  entry->SetLineColor(ci);
  entry->SetLineStyle(1);
  entry->SetLineWidth(10);
  entry->SetMarkerColor(ci);

  return;
}

