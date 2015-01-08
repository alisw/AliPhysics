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

void makeCentralityFit(const char * run="188359", const char *system = "ZNA", int Rebin=1, int Nevt=1e6)
{
 //load libraries
  gSystem->SetBuildDir("/tmp/");
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("AliCentralityGlauberFit.cxx+");

  const char *finnameGlau ="GlauberMC_pPb_ntuple_sigma70_mind4_r662_a546_Rpro6.root";
  char histname[8];
  
  if((strncmp(system,"ZNA",3))==0){
    printf("\n  Glauber fit on ZNA spectrum\n\n");
    sprintf(histname,"hZNA");
  }
  else if((strncmp (system,"ZNs",3)) == 0){
    printf("\n  Glauber fit on ZNA spectrum subtracting ZNC contribution (~SD)\n\n");
    sprintf(histname,"hZPA");
  }
  else if((strncmp (system,"ZPA",3)) == 0){
    printf("\n  Glauber fit on ZPA spectrum\n\n");
    sprintf(histname,"hZPA");
  }
  
  TString finname = Form("Histos%s.root",run);
  TString foutname = Form("%s_fit_%s.root",system,run);
  TString foutnameGlau = Form("%s_ntuple_%s.root",system,run);

  
  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit(finnameGlau);
  mPM->SetInputFile(finname);        
  mPM->SetInputNtuple(finnameGlau);     
  mPM->SetOutputFile(foutname);  
  mPM->SetOutputNtuple(foutnameGlau);
  mPM->AddHisto(histname);

  mPM->SetRebin(Rebin);
  mPM->SetNevents(Nevt);
  mPM->UseChi2(kTRUE);     // If TRUE minimize Chi2

  if (strncmp(system,"ZNA",3) == 0) {
    mPM->SetIsZN();
    mPM->SetRangeToFit(1., 300.);   
    mPM->SetRangeToScale(0);  
    mPM->SetGlauberParam(1,0.,1., 1,0.96,1., 2,0.22,0.25, 1,0.65,0.7, 1,0.56,0.585);
  }
  else if (strncmp(system,"ZPA",3) == 0) {
    mPM->SetIsZP();
    mPM->SetRangeToFit(1., 25.);  
    mPM->SetRangeToScale(1.);  
    mPM->SetGlauberParam(1,0.0,1., 1,0.6,1, 1,0.25,0.3, 1,0.65,0.8, 1,0.58,0.59);
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
  hg->Draw("E");
  hd->Draw("PEsame");
  hd->SetXTitle("E_{ZNA} (TeV)");
  hg->Draw("Esame");
  gPad->SetLogy(1);
  
  TLatex text0;
  text0.SetTextSize(0.04);
  text0.SetTextColor(kBlue+3);
  char ch[60];
  sprintf(ch,"<E_{ZNA}> DATA = %1.2f TeV ",hd->GetMean());
  text0.DrawLatex(xt,10.,ch);
  char chd[60];
  sprintf(chd,"<E_{ZNA}> Glauber = %1.2f TeV",hg->GetMean());
  text0.SetTextColor(kPink-2);
  text0.DrawLatex(xt,5.,chd);

  char ct[60];
  sprintf(ct,"RUN %s",run);
  text0.SetTextColor(kAzure);
  text0.DrawLatex(xt,1.,ct);

}

