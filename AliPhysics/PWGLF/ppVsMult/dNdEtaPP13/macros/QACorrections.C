// Plots the multiplicity dependence of the corrections used in the dndeta analysis
#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TFile.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TLegend.h"

#endif

void PrettifyAndAddEntry(TLegend * l, TH1*h, Int_t ibin) ;

void QACorrections (const char * resultsFile, Int_t nbin = 12){
  TFile * fres = TFile::Open(resultsFile);
  TObjArray * objArr = (TObjArray*) fres->Get("TObjArrayAux");

  TProfile ** profAlpha = new TProfile*[nbin];
  TProfile ** profBeta  = new TProfile*[nbin];
  Int_t firstbin = 0;
  for(Int_t ibin = firstbin; ibin < nbin; ibin++){
    
    TH2F * hAlpha = (TH2F*) objArr->FindObject(Form("bin%2.2d_Alpha" ,ibin));
    TH2F * hBeta  = (TH2F*) objArr->FindObject(Form("bin%d_mc_h1mBetaMC",ibin));
    std::cout << hBeta << std::endl;
    profAlpha[ibin] = hAlpha->ProfileX(Form("alpha_%d_px", ibin), 1, -1, "i");
    profBeta [ibin] = hBeta ->ProfileX(Form("beta_%d_px" , ibin), 1, -1, "i");
  }
  

  TCanvas * c1 = new TCanvas("cAlphaBetaVsCentr", "cAlphaBetaVsCentr", 1200, 500);
  c1->Divide(3,1);
  TLegend * l = new TLegend(0.1,0.1,0.8,0.8);
  for(Int_t ibin = firstbin; ibin < nbin; ibin++){
    TString drawOpt = ibin?"same":"";
    c1->cd(1);
    profAlpha[ibin]->Draw(drawOpt);
    PrettifyAndAddEntry(l, profAlpha[ibin], ibin);
    c1->cd(2);
    profBeta [ibin]->Draw(drawOpt);    
    PrettifyAndAddEntry(0, profBeta[ibin], ibin);
  }
  c1->cd(3);
  l->Draw();
  
}

void PrettifyAndAddEntry(TLegend * l, TH1*h, Int_t ibin) {
  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2,
                              kBlack-2, kRed-2 , kBlue-2, kGreen-2, kMagenta-2, kOrange+3,kCyan-2,kYellow-2,
  };
  const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar,
                              kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};


  if(l) {
    l->AddEntry(h, h->GetName(), "lp");    
  }
  h->SetMarkerStyle(markers[ibin]);
  h->SetMarkerColor(colors[ibin]);
  h->SetLineColor(colors[ibin]);
  
}
