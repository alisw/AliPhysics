/* $Id$ */

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"

void PPRstyle();

Double_t PeakPosition(Double_t pt);
Double_t PeakWidth(Double_t pt);
Double_t CB(Double_t * x, Double_t * par);
Double_t CB2(Double_t * x, Double_t * par);
Double_t CBs(Double_t * x, Double_t * par);
Double_t BG1(Double_t * x, Double_t * par);
Double_t BG2(Double_t * x, Double_t * par);


// const Int_t nPadX = 3, nPadY = 2;
// const Int_t nPtBins=6 ;
// const Double_t ptBinEdges[21]={1., 2., 3., 4., 5., 7., 10.} ;

const Int_t nPadX = 5, nPadY = 4;
Int_t nPtBins=0;
Double_t ptBinEdges[1000] = {0};
double GetPtBin(int bin){
  if( bin==0 )
    return 1.;

  // return GetPtBin(bin-1) * 1.1;

  // if ( bin % 2 )
  //   return GetPtBin(bin-1) + 0.4;
  // else
  //   return GetPtBin(bin-1) + 0.2;

  double previousBin = GetPtBin(bin-1);
  double linInc = 0.2;
  double threshold = 5.;
  double logFact = 1 + linInc/threshold;
  if ( previousBin < threshold ) // linear
    return previousBin + linInc;
  else { // logarithmic
    return previousBin * logFact;
  }
}
void MakePtBins() {
  int bin = -1;
  do {
    ++bin;
    ptBinEdges[bin] = GetPtBin(bin);
  } while(ptBinEdges[bin] < 40.);
  nPtBins = bin -2;

  printf("Making Pt Bins:\n");
  for(int b=0; b < nPtBins+1; ++b)
    printf("%.1f, ", ptBinEdges[b]);
  printf("\n N. Bins: %d\n", nPtBins);


  // for(int bin = 0; bin <= nPtBins; ++bin){
  //   ptBinEdges[bin] = GetPtBin(bin);
  //   printf("%.1f, ", ptBinEdges[bin]);
  // }
  // printf("\n");
}


const int kNCentralityBins = 3;

const Double_t kMean=0.136 ; //Approximate peak position to facilitate error estimate

const char format[] = ".pdf"; // say if you want .pdf

TH2F* FindPi0(TList *histoList, bool mix, const char* pid, int centrality);

//-----------------------------------------------------------------------------
void MakeMmixPi0(const TString filename,
		 const TString listPath = "PHOSPi0Flow/PHOSPi0FlowCoutput1", // lego train
		 const Int_t centrality=0,
		 const char* pid="CPV",
		 const char* trigger="kCentral",
		 const char* saveToDir="")
{
  MakePtBins();
  Printf("\nMakeMmixPi0(%s, %s, %i, %s, %s)", filename.Data(), listPath.Data(), centrality, pid, saveToDir);

  if( TString(saveToDir).Length() )
    gSystem->mkdir(saveToDir, true);

  //Fit Real/Mixed ratio, normalize Mixed and subtract it from Real

  TFile * file = new TFile(filename) ;
  TList *histoList = (TList*)file->Get(listPath);

  char key[125] ;

  TH1F * hev = (TH1F*)histoList->FindObject("hTotSelEvents") ;
  TH2F * hCentrality  = (TH2F*)histoList->FindObject("hCenPHOSCells") ;
  TH1D * hCentralityX = hCentrality->ProjectionX();

  printf("TotSelEvents (4): %.0f \n",hev->GetBinContent(4)) ;
  printf("Centrality:   %.0f \n",hCentralityX->Integral()) ;

  TH2F *hPi0 = FindPi0(histoList, false, pid, centrality);
  TH2F *hPi0Mix = FindPi0(histoList, true, pid, centrality);
  if( !hPi0 || !hPi0Mix || hPi0->GetEntries() < 10000) {
    Printf(Form("no histogram(0x%p, 0x%p) or to low number of entries(%.1f) in hPi0, skipping this cent", hPi0, hPi0Mix, hPi0->GetEntries()));
    file->Close();
    delete file;
    return;
  }


  TFile* outFile = new TFile(Form("%sPi0_FitResult.root", saveToDir, centrality),"update");

  PPRstyle();
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.01);
  //gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.08);

  //Fit real only
  //Linear Bg
  char kkey[55];
  sprintf(kkey, Form("%s_cent%d_%s",pid, centrality, trigger)) ;
  char key2[155];
  sprintf(key,"Mix%s",kkey) ;
  sprintf(key2,"%s_mr1",key) ;
  TH1D * mr1 = new TH1D(key2,"Mass",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_sr1",key) ;
  TH1D * sr1 = new TH1D(key2,"Width",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_ar1",key) ;
  TH1D * ar1 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_br1",key) ;
  TH1D * br1 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr1",key) ;
  TH1D * nr1 = new TH1D(key2,"Raw yield",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr1int",key) ;
  TH1D * nr1int = new TH1D(key2,"Raw yield, integrated",nPtBins,ptBinEdges) ;

  //Quadratic Bg
  sprintf(key2,"%s_mr2",key) ;
  TH1D * mr2 = new TH1D(key2,"Mass",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_sr2",key) ;
  TH1D * sr2 = new TH1D(key2,"Width",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_ar2",key) ;
  TH1D * ar2 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_br2",key) ;
  TH1D * br2 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_cr2",key) ;
  TH1D * cr2 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr2",key) ;
  TH1D * nr2 = new TH1D(key2,"Raw yield",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr2int",key) ;
  TH1D * nr2int = new TH1D(key2,"Raw yield, integrated",nPtBins,ptBinEdges) ;

  // Functions:
  TF1 * funcRatioFit1 = new TF1("funcRatioFit1",CB,0.,1.,6) ;
  funcRatioFit1->SetParName(0,"A") ;
  funcRatioFit1->SetParName(1,"m_{0}") ;
  funcRatioFit1->SetParName(2,"#sigma") ;
  funcRatioFit1->SetParName(3,"a_{0}") ;
  funcRatioFit1->SetParName(4,"a_{1}") ;
  funcRatioFit1->SetParName(5,"a_{2}") ;
  funcRatioFit1->SetLineWidth(2) ;
  funcRatioFit1->SetLineColor(2) ;
  TF1 * fgs = new TF1("gs",CBs,0.,1.,4) ;
  fgs->SetParName(0,"A") ;
  fgs->SetParName(1,"m_{0}") ;
  fgs->SetParName(2,"#sigma") ;
  fgs->SetParName(3,"B") ;
  fgs->SetLineColor(2) ;
  fgs->SetLineWidth(1) ;
  TF1 * funcRatioFit2 = new TF1("funcRatioFit2",CB2,0.,1.,7) ;
  funcRatioFit2->SetParName(0,"A") ;
  funcRatioFit2->SetParName(1,"m_{0}") ;
  funcRatioFit2->SetParName(2,"#sigma") ;
  funcRatioFit2->SetParName(3,"a_{0}") ;
  funcRatioFit2->SetParName(4,"a_{1}") ;
  funcRatioFit2->SetParName(5,"a_{2}") ;
  funcRatioFit2->SetParName(6,"a_{3}") ;
  funcRatioFit2->SetLineWidth(2) ;
  funcRatioFit2->SetLineColor(4) ;
  funcRatioFit2->SetLineStyle(2) ;
  TF1 * fbg1 = new TF1("bg1",BG1,0.,1.,3) ;
  TF1 * fbg2 = new TF1("bg2",BG2,0.,1.,4) ;


  // Canvases
  TCanvas * rawCanvas = new TCanvas("rawCanvas","rawCanvas",10,10,1200,800);
  rawCanvas->Divide(nPadX, nPadY);
  TCanvas * ratioCanvas = new TCanvas("mggFit1","mggFit1",10,10,1200,800) ;
  ratioCanvas->Divide(nPadX,nPadY) ;
  TCanvas * sbsCanvas = new TCanvas("mggFit1_Signal","mggFitCB",10,10,1200,800) ;
  sbsCanvas->Divide(nPadX,nPadY) ;


  for(Int_t ptBin=1; ptBin<=nPtBins; ptBin++){
    ratioCanvas->cd(ptBin) ;
    TAxis * pta=hPi0->GetYaxis() ;
    Int_t imin=pta->FindBin(ptBinEdges[ptBin-1]+0.0001);
    Int_t imax=pta->FindBin(ptBinEdges[ptBin]-0.0001) ;
    Double_t pt=(ptBinEdges[ptBin]+ptBinEdges[ptBin-1])/2. ;

    TH1D * hPi0Proj = hPi0->ProjectionX(Form("hPi0_ptBin%d",ptBin), imin, imax) ;
    TH1D * hPi0ProjMix= hPi0Mix->ProjectionX(Form("hPi0Mix_ptBin%d",ptBin), imin, imax) ;

    hPi0Proj->SetTitle(Form("M_{#gamma#gamma}, PID=%s, %.1f<p_{T}<%.1f GeV/c",pid,ptBinEdges[ptBin-1],ptBinEdges[ptBin]));
    hPi0Proj->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
    hPi0ProjMix->SetTitle(Form("M_{#gamma#gamma}^{Mix}, PID=%s, %.1f<p_{T}<%.1f GeV/c",pid,ptBinEdges[ptBin-1],ptBinEdges[ptBin]));
    hPi0ProjMix->SetXTitle("M_{#gamma#gamma}^{Mix} (GeV/c^{2})");

    Printf("\n\t%.1f<pt<%.1f GeV/c, entries: %.0f",ptBinEdges[ptBin-1],ptBinEdges[ptBin],  hPi0Proj->GetEntries());
    if( hPi0Proj->GetEntries() < 100. ) {
      Printf("skipping this bin as n. entries is to low");
      continue;
    }

    hPi0Proj ->Rebin(4) ;
    hPi0ProjMix->Rebin(4) ;
    hPi0Proj->SetName(Form("%s_rebin", hPi0Proj->GetName()));
    hPi0ProjMix->SetName(Form("%s_rebin", hPi0ProjMix->GetName()));

    // Error Fix
    for(Int_t ib=1; ib<=hPi0Proj->GetNbinsX();ib++){if(hPi0Proj ->GetBinContent(ib)==0)hPi0Proj ->SetBinError(ib,1.);}
    for(Int_t ib=1; ib<=hPi0Proj->GetNbinsX();ib++){if(hPi0ProjMix->GetBinContent(ib)==0)hPi0ProjMix->SetBinError(ib,1.);}

    // Signal-Mix Ratio
    TH1D * hPi0Ratio = (TH1D*)hPi0Proj->Clone( Form("hPi0Ratio_ptBin%d",ptBin) ) ;
    hPi0Ratio->SetTitle(Form("#frac{M_{#gamma#gamma}}{M_{#gamma#gamma}^{Mix}}, %.1f<p_{T}<%.1f GeV/c", ptBinEdges[ptBin-1], ptBinEdges[ptBin]));
    hPi0Ratio->Divide(hPi0ProjMix) ;
    hPi0Ratio->SetMarkerStyle(20) ;
    hPi0Ratio->SetMarkerSize(0.7) ;

    funcRatioFit1->SetParameters(0.001,0.136,0.0055,0.0002,-0.002,0.0) ;
    funcRatioFit1->SetParLimits(0,0.000,1.000) ;
    funcRatioFit1->SetParLimits(1,0.120,0.145) ;
    funcRatioFit1->SetParLimits(2,0.004,0.012) ;

    Double_t rangeMin=0.05 ;
    Double_t rangeMax=0.3 ;
    if(centrality==0) rangeMax=0.4 ;
    if(ptBin==1){
      rangeMin=0.06 ;
      rangeMax=0.25 ;
    }
    int error = hPi0Ratio->Fit(funcRatioFit1,"Q" ,"",rangeMin,rangeMax) ;
    if( error ) {
      Printf("fit (ratio) error: %d ", error);
      continue;
    }
    error = hPi0Ratio->Fit(funcRatioFit1,"MQ","",rangeMin,rangeMax) ;
    if( error % 4000) {
      Printf("fit (ratio more) error: %d ", error);
      continue;
    }


    ar1->SetBinContent(ptBin,funcRatioFit1->GetParameter(3)) ;
    ar1->SetBinError  (ptBin,funcRatioFit1->GetParError(3)) ;
    br1->SetBinContent(ptBin,funcRatioFit1->GetParameter(4)) ;
    br1->SetBinError  (ptBin,funcRatioFit1->GetParError(4)) ;

    funcRatioFit2->SetParameters(funcRatioFit1->GetParameters()) ;
    funcRatioFit2->SetParLimits(0,0.000,1.000) ;
    funcRatioFit2->SetParLimits(1,0.110,0.145) ;
    funcRatioFit2->SetParLimits(2,0.004,0.012) ;

    error = hPi0Ratio->Fit(funcRatioFit2,"+NQ","",rangeMin,rangeMax) ;
    if( error )  {
      Printf("fit (funcRatioFit2) error: %d ", error);
      continue;
    }
    error =hPi0Ratio->Fit(funcRatioFit2,"+MQ","",rangeMin,rangeMax) ;
    if( error % 4000) {
      Printf("fit (funcRatioFit2 more) error: %d ", error);
      continue;
    }

    ar2->SetBinContent(ptBin,funcRatioFit2->GetParameter(3)) ;
    ar2->SetBinError  (ptBin,funcRatioFit2->GetParError(3)) ;
    br2->SetBinContent(ptBin,funcRatioFit2->GetParameter(4)) ;
    br2->SetBinError  (ptBin,funcRatioFit2->GetParError(4)) ;
    cr2->SetBinContent(ptBin,funcRatioFit2->GetParameter(5)) ;
    cr2->SetBinError  (ptBin,funcRatioFit2->GetParError(5)) ;
    hPi0Ratio->GetXaxis()->SetRangeUser(0.05,0.30) ;
    hPi0Ratio->DrawCopy() ;
    hPi0Ratio->Write(0,TObject::kOverwrite) ;
    ratioCanvas->Update() ;

    fbg1->SetParameters(funcRatioFit1->GetParameter(3),funcRatioFit1->GetParameter(4),funcRatioFit1->GetParameter(5));
    fbg2->SetParameters(funcRatioFit2->GetParameter(3),funcRatioFit2->GetParameter(4),funcRatioFit2->GetParameter(5),
			funcRatioFit2->GetParameter(6));

    Double_t intRangeMin = PeakPosition(pt)-3.*PeakWidth(pt) ;
    Double_t intRangeMax = PeakPosition(pt)+3.*PeakWidth(pt) ;
    Int_t    intBinMin   = hPi0Proj->GetXaxis()->FindBin(intRangeMin) ;
    Int_t    intBinMax   = hPi0Proj->GetXaxis()->FindBin(intRangeMax) ;
    Double_t errStat     = hPi0ProjMix->Integral(intBinMin,intBinMax);

    rawCanvas->cd(ptBin);
    TH1D * hPi0MixScaled = (TH1D*)hPi0ProjMix ->Clone(Form("%s_clone2", hPi0Proj->GetName())) ;
    hPi0MixScaled->Scale(fbg1->Eval(0.1349));
    hPi0Proj->SetLineColor(kBlack);
    hPi0Proj->SetAxisRange(0.0, 0.3);
    hPi0Proj->DrawCopy();
    hPi0Proj->Write(0,TObject::kOverwrite) ;
    hPi0MixScaled->SetLineColor(kRed);
    hPi0MixScaled->SetTitle(Form("M_{#gamma#gamma}^{Mix,scaled}, PID=%s, %.1f<p_{T}<%.1f GeV/c",pid,ptBinEdges[ptBin-1],ptBinEdges[ptBin]));
    hPi0MixScaled->DrawCopy("same");
    rawCanvas->Update();
    hPi0MixScaled->Write(0,TObject::kOverwrite) ;

    // Linear background subtraction
    TH1D * hPi0MixScaledPol1 = (TH1D*)hPi0ProjMix->Clone(Form("hPi0MixScaledPol1_ptBin%d", ptBin)) ;
    hPi0MixScaledPol1 ->Multiply(fbg1) ;
    TH1D * hPi0BSPol1 = (TH1D*)hPi0Proj->Clone(Form("hPi0BSPol1_ptBin%d", ptBin)) ;
    hPi0BSPol1->Add(hPi0MixScaledPol1 ,-1.) ;

    // Quadratic background
    TH1D * hPi0MixScaledPol2 = (TH1D*)hPi0ProjMix->Clone(Form("hPi0MixScaledPol2_ptBin%d", ptBin)) ;
    hPi0MixScaledPol2->Multiply(fbg2) ;
    TH1D * hPi0BSPol2     = (TH1D*)hPi0Proj    ->Clone(Form("hPi0BSPol2_ptBin%d", ptBin)) ;
    hPi0BSPol2 ->Add(hPi0MixScaledPol2,-1.) ;

    sbsCanvas->cd(ptBin) ;

    Int_t binPi0 = hPi0BSPol1->FindBin(funcRatioFit1->GetParameter(1));
    Int_t nWidPi0 = 2 * (Int_t) (funcRatioFit1->GetParameter(2)/hPi0BSPol1->GetBinWidth(1));
    fgs->SetParameters(hPi0BSPol1->Integral(binPi0-nWidPi0,binPi0+nWidPi0)/5., funcRatioFit1->GetParameter(1), funcRatioFit1->GetParameter(2)) ;
    fgs->SetParLimits(0,0.,1.e+5) ;
    fgs->SetParLimits(1,0.110,0.145) ;
    fgs->SetParLimits(2,0.004,0.02) ;
    error = hPi0BSPol1->Fit(fgs,"Q","",rangeMin,rangeMax) ;
    if( error ) {
      Printf("fit (hPi0BSPol1 fgs) error: %d ", error);
      continue;
    }
    hPi0BSPol1->SetMaximum(hPi0BSPol2->GetMaximum()*1.5) ;
    hPi0BSPol1->SetMinimum(hPi0BSPol2->GetMinimum()*1.1) ;
    hPi0BSPol1->SetMarkerStyle(20) ;
    hPi0BSPol1->SetMarkerSize(0.7) ;
    mr1->SetBinContent(ptBin,fgs->GetParameter(1)) ;
    mr1->SetBinError  (ptBin,fgs->GetParError(1) ) ;
    sr1->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
    sr1->SetBinError  (ptBin,fgs->GetParError(2) ) ;

    Double_t y=fgs->GetParameter(0)/hPi0BSPol1->GetXaxis()->GetBinWidth(1) ;
    nr1->SetBinContent(ptBin,y) ;
    Double_t ey=fgs->GetParError(0)/hPi0BSPol1->GetXaxis()->GetBinWidth(1) ;
    nr1->SetBinError(ptBin,ey) ;

    Double_t npiInt = hPi0BSPol1->Integral(intBinMin,intBinMax) ;
    Double_t norm   = fbg1->GetParameter(0) ;
    Double_t normErr= fbg1->GetParError(0) ;
    if(npiInt>0.){
      nr1int->SetBinContent(ptBin,npiInt) ;
      nr1int->SetBinError(ptBin,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
    }
    hPi0BSPol2->GetXaxis()->SetRangeUser(0.05,0.3) ;
    hPi0BSPol2->SetMaximum(hPi0BSPol2->GetMaximum()*1.5) ;
    hPi0BSPol2->SetMinimum(hPi0BSPol2->GetMinimum()*1.1) ;
    hPi0BSPol2->SetMarkerStyle(20) ;
    hPi0BSPol2->SetMarkerSize(0.7) ;
    hPi0BSPol2->Fit(fgs,"Q","",rangeMin,rangeMax) ;
    if( error ) {
      Printf("fit (hPi0BSPol2 fgs) error: %d ", error);
      continue;
    }
    mr2->SetBinContent(ptBin,fgs->GetParameter(1)) ;
    mr2->SetBinError  (ptBin,fgs->GetParError(1)) ;
    sr2->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
    sr2->SetBinError  (ptBin,fgs->GetParError(2)) ;
    y=fgs->GetParameter(0)/hPi0BSPol1->GetXaxis()->GetBinWidth(1) ;
    nr2->SetBinContent(ptBin,y) ;
    ey= fgs->GetParError(0)/hPi0BSPol1->GetXaxis()->GetBinWidth(1) ;
    nr2->SetBinError(ptBin,ey) ;
    Double_t npiInt2=hPi0BSPol2->Integral(intBinMin,intBinMax) ;
    norm=fbg2->GetParameter(0) ;
    normErr=fbg2->GetParError(0) ;
    if(npiInt2>0.){
      nr2int->SetBinContent(ptBin,npiInt2) ;
      nr2int->SetBinError(ptBin,TMath::Sqrt(npiInt2 + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
    }
    hPi0BSPol2->SetTitle(Form("M_{#gamma#gamma}^{BS_{2}}, PID=%s, %.1f<p_{T}<%.1f GeV/c",pid,ptBinEdges[ptBin-1],ptBinEdges[ptBin]));
    hPi0BSPol2->DrawCopy() ;
    sbsCanvas->Update() ;
    hPi0BSPol2->Write(0,TObject::kOverwrite) ;
  }

  char name[55] ;
  sprintf(name,"%sPi0_ratio_%s%s", saveToDir, kkey, format) ;
  ratioCanvas->Print(name) ;
  sprintf(name,"%sPi0_signal_%s%s", saveToDir, kkey, format) ;
  sbsCanvas->Print(name) ;
  sprintf(name,"%sPi0_raw_%s%s", saveToDir, kkey, format) ;
  rawCanvas->Print(name);

  //Normalize by the number of events
  Int_t cMin=0, cMax=100;
  if      (centrality == 0) {
    cMin=1;
    cMax=10;
  }
  else if (centrality == 1) {
    cMin=11;
    cMax=40;
  }
  else if (centrality == 2) {
    cMin=41;
    cMax=80;
  }
  else if (centrality == -1) {
    cMin=1;
    cMax=80;
  }
  else {
    Printf("ERROR: Centrality: %d not defined !!! ERROR", centrality);
    return;
  }

  Double_t nevents = hCentralityX->Integral(cMin,cMax);
  if ( nevents > 0.9 ) {
    nr1   ->Scale(1./nevents) ;
    nr1int->Scale(1./nevents) ;
    nr2   ->Scale(1./nevents) ;
    nr2int->Scale(1./nevents) ;
  } else {
    Printf("WARNING: non positive nEvents in centrality range, cMin:%d, cMax:%d, nEvents:%f", cMin, cMax, nevents );

  }

  mr1->Write(0,TObject::kOverwrite) ;
  sr1->Write(0,TObject::kOverwrite) ;
  ar1->Write(0,TObject::kOverwrite) ;
  br1->Write(0,TObject::kOverwrite) ;
  nr1->Write(0,TObject::kOverwrite) ;
  nr1int->Write(0,TObject::kOverwrite) ;
  ar2->Write(0,TObject::kOverwrite) ;
  br2->Write(0,TObject::kOverwrite) ;
  cr2->Write(0,TObject::kOverwrite) ;
  mr2->Write(0,TObject::kOverwrite) ;
  sr2->Write(0,TObject::kOverwrite) ;
  nr2->Write(0,TObject::kOverwrite) ;
  nr2int->Write(0,TObject::kOverwrite) ;


  outFile->Close();
  delete outFile;

  file->Close();
  delete file;

  delete ratioCanvas;
  delete sbsCanvas;
  delete rawCanvas;
}

//-----------------------------------------------------------------------------
Double_t PeakPosition(Double_t pt){
  //Fit to the measured peak position
  return 4.99292e-003*exp(-pt*9.32300e-001)+1.34944e-001 ;
}
//-----------------------------------------------------------------------------
Double_t PeakWidth(Double_t pt){
  //fit to the measured peak width
  Double_t a=0.0068 ;
  Double_t b=0.0025 ;
  Double_t c=0.000319 ;
  return TMath::Sqrt(a*a+b*b/pt/pt+c*c*pt*pt) ;
}

//-----------------------------------------------------------------------------
Double_t CB(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean) ;
}

//-----------------------------------------------------------------------------
Double_t CB2(Double_t * x, Double_t * par){
  //Another parameterization of Real/Mixed ratio7TeV/README
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean)+par[5]*(x[0]-kMean)*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
Double_t CBs(Double_t * x, Double_t * par){
  //Parameterizatin of signal
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3] ;
}
//-----------------------------------------------------------------------------
Double_t BG1(Double_t * x, Double_t * par){
  //Normalizatino of Mixed
  return par[0]+par[1]*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
Double_t BG2(Double_t * x, Double_t * par){
  //Another normalization of  Mixed
  return par[0]+par[1]*(x[0]-kMean)+par[2]*(x[0]-kMean)*(x[0]-kMean) ;
}


//-----------------------------------------------------------------------------
void PPRstyle()
{

  //////////////////////////////////////////////////////////////////////
  //
  // ROOT style macro for the TRD TDR
  //
  //////////////////////////////////////////////////////////////////////

  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(-1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(10);

  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(-1);
  gStyle->SetFrameLineWidth(1.2);
  gStyle->SetFrameLineColor(1);

  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(10);
  gStyle->SetPadBorderSize(1);
  gStyle->SetPadBorderMode(-1);

  gStyle->SetStatColor(10);
  gStyle->SetTitleColor(kBlack,"X");
  gStyle->SetTitleColor(kBlack,"Y");

  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetLabelSize(0.04,"Z");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelFont(42,"Z");
  gStyle->SetStatFont(42);

  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.4,"Y");

  gStyle->SetFillColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);

}

// For different tasks (e.g. triggers)
void MakeMmixPi0()
{
  // Take care when uncommenting, seems to leak memory consumption.

  // // All Centrality
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow/PHOSPi0FlowCoutput1", -1, "CPV", "out/kCentral/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_SemiCentral/PHOSPi0Flow_SemiCentralCoutput1", -1, "CPV", "out/kSemiCentral/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_MB/PHOSPi0Flow_MBCoutput1", -1, "CPV", "out/kMB/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_PHOSPb/PHOSPi0Flow_PHOSPbCoutput1", -1, "CPV", "out/kPHOSPb/");

  // 0-10%
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow/PHOSPi0FlowCoutput1", 0, "CPV", "out/kCentral/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_MB/PHOSPi0Flow_MBCoutput1", 0, "CPV", "out/kMB/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_PHOSPb/PHOSPi0Flow_PHOSPbCoutput1", 0, "CPV", "out/kPHOSPb/");

  // // 10-40%
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_SemiCentral/PHOSPi0Flow_SemiCentralCoutput1", 1, "CPV", "out/kSemiCentral/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_MB/PHOSPi0Flow_MBCoutput1", 1, "CPV", "out/kMB/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_PHOSPb/PHOSPi0Flow_PHOSPbCoutput1", 1, "CPV", "out/kPHOSPb/");

  // // 40-80%
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_SemiCentral/PHOSPi0Flow_SemiCentralCoutput1", 2, "CPV", "out/kSemiCentral/");
  // MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_MB/PHOSPi0Flow_MBCoutput1", 2, "CPV", "out/kMB/");
  MakeMmixPi0("AnalysisResults.root", "PHOSPi0Flow_PHOSPb/PHOSPi0Flow_PHOSPbCoutput1", 2, "CPV", "out/kPHOSPb/");
}



TH2F* FindPi0(TList *histoList, bool mix, const char* pid, int centrality)
{
  // consider mixed or not
  TString mixStr("");
  if( mix )
    mixStr = "Mi";

  // If centrality is some integer in range [0, kNCentralityBins] return Pi0 plot for that cent.
  if( centrality >= 0 && centrality < kNCentralityBins ) {
    TString name = Form("h%sPi0%s_cen%d", mixStr.Data(), pid, centrality);
    TH2F* hist = (TH2F*) histoList->FindObject(name.Data());
    hist->Sumw2();
    return hist;
  }
  // If centrality is '-1' Merge [0, kNCentralityBins)
  else if ( centrality == -1 ) {
    TString name = Form("h%sPi0%s_cen%d", mixStr.Data(), pid, 0);
    TH2F* histMerge = (TH2F*) histoList->FindObject(name.Data()) -> Clone(Form("h%sPi0%s_merged", mixStr.Data(), pid));
    histMerge->Sumw2();

    for( int cent = 1; cent < kNCentralityBins; ++cent ) {
      name = Form("h%sPi0%s_cen%d", mixStr.Data(), pid, cent);
      TH2F* other = (TH2F*) histoList->FindObject(name.Data());
      other->Sumw2();
      histMerge->Add( other );
    }
    return histMerge;
  }
  else { // case not defined
    Printf("ERROR: Centrality must be in range: [%d,%d]", -1, kNCentralityBins - 1 );
    return 0;
  }
}
