#if !defined(__CINT__) || defined(__MAKECINT__)
   #include <TMath.h>
   #include <TROOT.h>
   #include <Riostream.h>
   #include <TCanvas.h>
   #include <TColor.h>
   #include <TLatex.h>
   #include <TLegend.h>

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
   //4.5,5.0,5.5,6.5,8.0,12.0
};
static Int_t nBins=sizeof(xBins)/sizeof(Double_t) - 1;

static Bool_t gFlag=kFALSE;

//*** The systematic uncertainties for combining
// cos(PA)
// DCA between V0 daughters
// TPC crossed/findable
static Double_t sysErrK0s[]={ 
  0.0,
  0.05,0.05,0.04,0.03,  //Dominated by cos(PA)
  0.04,0.05,0.06,0.07,0.08,0.09, //Dominated by TPC crossed/findable
  0.10,0.10,0.11,0.12,0.12,0.13,0.13,0.13,0.14,0.14,0.14,
  0.14,0.13,0.13,0.12,0.11,0.11,0.11,0.10,0.09,0.08,
  0.08,0.05,0.04,0.03,
  0.3,0.3,0.3 // "trailers"
};
const Double_t sysPID=0.02; //PID
const Double_t sysSig=0.03; //Signal extraction
const Double_t sysArm=0.01; //Armenteros cut


TH1 *MapHisto(const TH1 *h) {
  const Double_t eps=0.0001;
  TString name("m");
  name = name + h->GetName(); 
  TH1F *mh=new TH1F(name.Data(),h->GetTitle(),nBins,xBins);
  mh->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  mh->GetYaxis()->SetTitle("1/N_{ev}dN/dp_{T}/dy (GeV/c)^{-1}");

  Double_t xh=h->GetBinCenter(1), xmh=0.;

  Int_t n=1;
  for (; n<=nBins; n++) {
    xmh=mh->GetBinCenter(n);
    if (TMath::Abs(xh-xmh)<eps) break; 
  }

  Int_t iii=h->GetNbinsX();
  if (gFlag) iii--;
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
  raw=(TH1F*)lst->FindObject(rName[1]);
  if (!raw) {
     cerr<<"No raw yield !"<<eName[0]<<' '<<eName[1]<<endl; 
     return kFALSE;
  }
  /*TFile *fe=*/TFile::Open(eName[0]);
  eff=(TH1F*)gFile->Get(eName[1]);
  if (!eff) {
     cerr<<"No efficiency ! "<<eName[0]<<' '<<eName[1]<<endl; 
     return kFALSE;
  }
  return kTRUE;
}

void SetAttributes(TH1 *h,const Char_t *tit,Int_t col,Int_t mar,Float_t siz) {
  h->SetTitle(tit);
  h->SetLineColor(col); 
  h->SetMarkerColor(col);
  h->SetMarkerStyle(mar);
  h->SetMarkerSize(siz);
  h->SetMaximum(1000);
  h->SetMinimum(1e-5);
}

void DrawHisto(TH1 *h, const Option_t *option, Double_t *sysErr) {
  TH1F *hh=new TH1F(*((TH1F*)h));
  Int_t nb=hh->GetNbinsX();
  for (Int_t i=1; i<=nb; i++) {
      Double_t c=hh->GetBinContent(i);
      Double_t e=hh->GetBinError(i);
      e = c*TMath::Sqrt(sysErr[i]*sysErr[i] + 
                        sysPID*sysPID + 
                        sysSig*sysSig +
                        sysArm*sysArm);
      hh->SetBinError(i,e);
  }
  hh->SetFillColor(17);
  TString opt("E5"); opt+=option;
  //TString opt("E2"); opt+=option;
  hh->Draw(opt.Data());
  h->Draw("e x0 same");
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


void DrawSpectraAndRatios() {

  const Int_t nCent=5;

  const Char_t *title[nCent]={
    "00-05 %",
    "20-40 %",
    "40-60 %",
    "60-80 %",
    "80-90 %"
  };
  const Int_t   colour[nCent]={2  , 419, 4  , 6 , 1  };
  const Int_t   marker[nCent]={22 , 21 , 23 , 33, 20 };
  const Float_t masize[nCent]={1.6, 1.3, 1.6, 2 , 1.3};
  
  const Char_t *rNameL[2*nCent]={ // file name, histo name
    "luke/DataYields_5sigCubic_070512.root", "YieldLambda0005", 
    "luke/DataYields_200312.root", "YieldLambda2040", 
    "luke/DataYields_5sigCubic_070512.root", "YieldLambda4060", 
    "luke/DataYields_200312.root", "YieldLambda6080", 
    "luke/DataYields_200312.root", "YieldLambda8090" 
  };
  const Char_t *eNameL[2*nCent]={ // file name, histo name
    "marian/combined/EFF_Lambda_PbPb_00_05.root", "Efficiency",
    "marian/combined/EFF_Lambda_PbPb_20_40.root", "Efficiency",
    "marian/combined/EFF_Lambda_PbPb_40_60.root", "Efficiency",
    "marian/combined/EFF_Lambda_PbPb_60_80.root", "Efficiency",
    "marian/combined/EFF_Lambda_PbPb_80_90.root", "Efficiency"
  };
  /*
  const Char_t *rNameL[2*nCent]={ // file name, histo name
    "luke/DataYields_5sigCubic_070512.root", "YieldAntiLambda0005", 
    "luke/DataYields_200312.root", "YieldAntiLambda2040", 
    "luke/DataYields_5sigCubic_070512.root", "YieldAntiLambda4060", 
    "luke/DataYields_200312.root", "YieldAntiLambda6080", 
    "luke/DataYields_200312.root", "YieldAntiLambda8090" 
  };
  const Char_t *eNameL[2*nCent]={ // file name, histo name
    "marian/combined/EFF_AntiLambda_PbPb_00_05.root", "Efficiency",
    "marian/combined/EFF_AntiLambda_PbPb_20_40.root", "Efficiency",
    "marian/combined/EFF_AntiLambda_PbPb_40_60.root", "Efficiency",
    "marian/combined/EFF_AntiLambda_PbPb_60_80.root", "Efficiency",
    "marian/combined/EFF_AntiLambda_PbPb_80_90.root", "Efficiency"
  };
  */
  const Char_t *rNameK[2*nCent]={ // file name, histo name
    "luke/DataYields_5sigCubic_070512.root", "YieldK0Short0005", 
    "luke/DataYields_5sigCubic_070512.root", "YieldK0Short2040", 
    "luke/DataYields_5sigCubic_070512.root", "YieldK0Short4060", 
    "luke/DataYields_5sigCubic_070512.root", "YieldK0Short6080", 
    "luke/DataYields_5sigCubic_070512.root", "YieldK0Short8090" 
  };
  const Char_t *eNameK[2*nCent]={ // file name, histo name
    "marian/combined/EFF_K0s_PbPb_00_05.root", "Efficiency",
    "marian/combined/EFF_K0s_PbPb_20_40.root", "Efficiency",
    "marian/combined/EFF_K0s_PbPb_40_60.root", "Efficiency",
    "marian/combined/EFF_K0s_PbPb_60_80.root", "Efficiency",
    "marian/combined/EFF_K0s_PbPb_80_90.root", "Efficiency"
  };

//   const Char_t *fdName[nCent]={
//     "fd_0005", "fd_2040", "fd_4060", "fd_6080", "fd_8090"
//   };
//   TFile *fdFile=TFile::Open("lambda_fd.root");  

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFillColor(0);

  TH1 *raw=0;
  TH1 *eff=0;
  TString option(""), ratio("ratio");

  TCanvas *c1=new TCanvas; c1->SetLogy();
  c1->SetLeftMargin(0.13); c1->SetBottomMargin(0.13);

  TCanvas *c2=new TCanvas; c2->SetLogy();
  c2->SetLeftMargin(0.13); c2->SetBottomMargin(0.13);

  TCanvas *c3=new TCanvas;
  //TCanvas *c4=new TCanvas;


  for (Int_t cent=0; cent<nCent; cent++) {
      // Lambda
      if (!GetHistos(rNameL+2*cent, eNameL+2*cent, raw, eff)) return;
      TH1 *rawHl=MapHisto(raw);
      TH1 *effHl=MapHisto(eff);
      rawHl->Divide(effHl);
      SetAttributes(rawHl,title[cent],colour[cent],marker[cent],masize[cent]);
      c1->cd();
      //rawHl->Draw(option.Data());
      DrawHisto(rawHl, option.Data(), sysErrK0s);

      // K0s
    if (cent==nCent-1) gFlag=kTRUE;
      if (!GetHistos(rNameK+2*cent, eNameK+2*cent, raw, eff)) return;
      TH1 *rawHk=MapHisto(raw);
      TH1 *effHk=MapHisto(eff);
      rawHk->Divide(effHk);
      SetAttributes(rawHk,title[cent],colour[cent],marker[cent],masize[cent]);
      c2->cd();
      //rawHk->Draw(option.Data());
      DrawHisto(rawHk, option.Data(), sysErrK0s);

      // Lambda/K0s
      TH1 *rawHlk=(TH1*)rawHl->Clone();
      //FD
      rawHlk->Scale(0.8);
//       fdFile->cd();
//       TH1 *fd=(TH1*)gFile->Get(fdName[cent]);
//       TH1 *fdM=MapHisto(fd);
//       SetAttributes(fdM,title[cent],colour[cent],marker[cent],masize[cent]);
//       c4->cd();
//       fdM->Draw(option.Data());
//       rawHlk->Add(fdM,-1);
      //
      TString name=ratio+rawHlk->GetName();
      rawHlk->SetName(name.Data());      
      rawHlk->Divide(rawHk);
      c3->cd();
      rawHlk->Draw(option.Data());

      option+="same";
  }

  Float_t offx=0.15, offy=0.16, sizx=0.22, sizy=0.22;
  TLegend *leg=c1->BuildLegend(0.68,0.52,0.88,0.88,"Centrality:");
  leg->SetFillColor(0);
  c1->cd(); 
  TLatex *   tex = new TLatex(5.5,19.,"#Lambda");
  tex->SetTextFont(42);
  tex->SetTextSize(0.11);
  tex->SetLineWidth(2);
  tex->Draw();
  DrawALICELogo(offx,offy,offx+sizx,offy+sizy);

  //c2->BuildLegend(0.74,0.62,0.88,0.88);
  leg=c2->BuildLegend(0.68,0.52,0.88,0.88,"Centrality:");
  leg->SetFillColor(0);
  c2->cd(); 
  tex = new TLatex(5.5,19.,"K^{0}_{S}");
  tex->SetTextFont(42);
  tex->SetTextSize(0.089);
  tex->SetLineWidth(2);
  tex->Draw();
  DrawALICELogo(offx,offy,offx+sizx,offy+sizy);

  c3->BuildLegend(0.74,0.62,0.88,0.88);

  return;
}

