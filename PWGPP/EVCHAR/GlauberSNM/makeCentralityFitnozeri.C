#if !defined(__CINT__) || defined(__MAKECINT__)

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include "AliCentralityGlauberFit.h"

#endif

void makeCentralityFitnozeri(const int nRun=195483, const char *system = "ZNA",
	int ntrials=1, int Rebin=1, int Nevt=1.e6)
{
 //load libraries
  gSystem->SetBuildDir(".");
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("AliCentralityGlauberFit.cxx++g");

  const char *finnameGlau ="GlauberMC_pPb_ntuple_sigma70_mind4_r662_a546_Rpro6.root";
  char histname[8];
  double chi2min=0.;
	
  if((strncmp(system,"ZNA",3))==0){
    printf("\n  Glauber fit on ZNA spectrum\n\n");
    sprintf(histname,"hZNA");
    chi2min = 200.;
  }
  else if((strncmp (system,"ZPA",3)) == 0){
    printf("\n  Glauber fit on ZPA spectrum\n\n");
    sprintf(histname,"hZPA");
    chi2min=5000.;
  }
  
  TString finname = Form("centrHistos%d-dNdeta.root",nRun);
  printf(" Opening file %s\n",finname.Data());
  //
  TString foutname = Form("%s_fit_%d.root",system,nRun);
  TString foutnameGlau = Form("%s_ntuple_%d.root",system,nRun);

  
  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit(finnameGlau);
  mPM->SetInputFile(finname);        
  mPM->SetInputNtuple(finnameGlau);     
  mPM->SetOutputFile(foutname);  
  mPM->SetOutputNtuple(foutnameGlau);
  mPM->AddHisto(histname);

  mPM->SetRebin(Rebin);
  mPM->SetNevents(Nevt);
  mPM->UseChi2(kTRUE);     // If TRUE minimize Chi2
  mPM->SetNTrials(ntrials);
  mPM->SetChi2Min(chi2min);
  //
  // COSY ORIGINAL!!!!
  //mPM->SetNParam(51.5, 469.2, 8.762);
  // COSY-like
//  mPM->SetNParam(51.5, 620., 10.);
// change 2
  mPM->SetNParam(50., 230., 4.4, 0.48);
  //mPM->SetNParam(40., 230., 4.5, 0.45);
  //
  // ALICE ex-novo
  //mPM->SetNParam(78., 700., 7.);
  //mPM->SetNParam(61., 470., 7.);

  if(strncmp(system,"ZNA",3) == 0) {
    printf(" Setting parameters for ZNA Glauber fit\n\n");
    mPM->SetIsZN();
    mPM->SetRangeToFit(0.8, 90.5);   
    mPM->SetRangeToScale(0);  
    // original
    if(nRun==195483) mPM->SetGlauberParam(1,0.0,0.0, 1,0.957,1., 1,0.26,0.30, 1,0.65,0.65, 1,0.585,0.585, 1, 0.25,0.3);
    //to study systematics
    //if(nRun==195483) mPM->SetGlauberParam(1,0.0,0.0, 1,0.956,1., 1,0.29,0.30, 1,1.80,0.65, 1,0.585,0.585, 1, 0.25,0.3);
    else             mPM->SetGlauberParam(1,0.0,0.0, 1,0.957,1., 1,0.25,0.30, 1,0.65,0.65, 1,0.585,0.585, 1, 0.25,0.3);
  }
  else if (strncmp(system,"ZPA",3) == 0) {
    mPM->SetIsZP();
    mPM->SetRangeToFit(1., 30.);  
    mPM->SetRangeToScale(0);  
    // original
    //if(nRun==195483) mPM->SetGlauberParam(1,0.0,0.0, 1,0.60,0.8, 1,0.6,0.65, 1,0.65,0.65, 1,0.585,0.585, 1, 0.25,0.4);
    if(nRun==195483) mPM->SetGlauberParam(1,0.0,0.0, 1,0.40,0.8, 1,0.4,0.65, 1,1.80,0.65, 1,0.585,0.585, 1, 0.25,0.4);
    else             mPM->SetGlauberParam(1,0.0,0.0, 1,0.60,0.8, 1,0.5,0.65, 1,0.65,0.65, 1,0.585,0.585, 1, 0.25,0.4);
  }
  mPM->MakeFits();  

  char hnam[10], hnamg[20];
  double xt=0;
  if(strncmp(system,"ZNA",3) == 0){
    sprintf(hnam,"hZNA");
    sprintf(hnamg,"hZNA_GLAU");
    xt=20.;
  }
  else if(strncmp(system,"ZPA",3) == 0) {
    sprintf(hnam,"hZPA");
    sprintf(hnamg,"hZPA_GLAU");
    xt=5.;
  }

  TFile *f = TFile::Open(foutname);
  TH1 * hd = dynamic_cast<TH1*> (f->Get((hnam)));
  TH1 * hg = dynamic_cast<TH1*> (f->Get((hnamg)));
  hg->SetLineColor(kPink-2);
  hg->SetLineWidth(2);
  hd->SetMarkerStyle(20);
  //hd->SetMarkerSize(1.2);
  hd->SetMarkerColor(kBlue+3);
  hd->SetLineWidth(2);
  //hd->SetMinimum(10.);
  
  TCanvas *g = new TCanvas("g","g",0,0,700,700);
  g->cd();
  gPad->SetLogy(1);
  hg->Draw("E");
  hd->Draw("PEsame");
  hd->SetXTitle("E_{ZNA} (TeV)");
  hg->Draw("Esame");
  
  printf("\n Entries with zero energy %1.0f  -> %f of the total (%d)\n\n",
  	hg->GetBinContent(hg->FindBin(0.)), hg->GetBinContent(hg->FindBin(0.))/hg->Integral(),
	hg->Integral());

  /*double zncut[7]  = {0., 16.3423, 38.2513, 52.7698, 63.7601, 70.2207, 75.0126};
  TH1F *hist[7];
  for(int ih=0; ih<7; ih++){
    char hnam[24];
    sprintf(hnam,"hZNA%d",ih);
    hist[ih] = new TH1F(hnam,hnam,hd->GetNbinsX(),0.,142.5);
  }
  int index=0;
  for(int ic=0; ic<7; ic++){
    for(int ib=1; ib<hd->GetNbinsX(); ib++){
      if(hd->GetBinLowEdge(ib)>=zncut[ic] && hd->GetBinLowEdge(ib+1)<zncut[ic+1]){
        index = ic;
        hist[index]->SetBinContent(ib, hd->GetBinContent(ib));
        if(hd->GetBinError(ib)>0) hist[index]->SetBinError(ib, 1./(hd->GetBinError(ib)*hd->GetBinError(ib)));
      }
    }
  }
  for(int ih=0; ih<7; ih++){
    if(ih%2==0) hist[ih]->SetFillColor(kAzure+6);
    hist[ih]->Draw("hist SAME");
  }
  hg->Draw("Esame");*/
  
  TLatex text0;
  text0.SetTextSize(0.04);
  text0.SetTextColor(kBlue+3);
  char ch[60];
  if(strncmp(system,"ZNA",3) == 0) sprintf(ch,"<E_{ZNA}> DATA = %1.2f TeV ",hd->GetMean());
  else  sprintf(ch,"<E_{ZPA}> DATA = %1.2f TeV ",hd->GetMean());
  text0.DrawLatex(xt,20.,ch);
  char chd[60];
  if(strncmp(system,"ZNA",3) == 0) sprintf(chd,"<E_{ZNA}> Glauber = %1.2f TeV",hg->GetMean());
  else  sprintf(chd,"<E_{ZPA}> Glauber = %1.2f TeV",hg->GetMean());
  text0.SetTextColor(kPink-2);
  text0.DrawLatex(xt,10.,chd);

  char ct[60];
  sprintf(ct,"RUN %d",nRun);
  text0.SetTextColor(kTeal+2);
  text0.DrawLatex(xt+80.,1.e4,ct);
  
  char psn[30];
  if(strncmp(system,"ZNA",3) == 0) sprintf(psn,"ZNA_fit%d.gif",nRun);
  else  sprintf(psn,"ZPA_fit%d.gif",nRun);
  g->Print(psn);
	
	printf("\n Everything is OK: ntuple and fit results saved!!! \n\n");
}

