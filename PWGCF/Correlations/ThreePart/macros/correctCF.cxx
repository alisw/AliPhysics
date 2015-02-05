#include <iostream>
#include "TFile.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStyle.h"
#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TNamed.h"
#include "TList.h"
#include "TDirectory.h"
#include <X3DDefs.h>
#include "TPaveText.h"
#include "TMinuit.h"
#include "TLegend.h"

using namespace std;
//global pointers for the histograms:

TH2D* gSameEvent =0x0;
TH2D* gMETA = 0x0;
TH2D* gMETrigger = 0x0;

const double gkEtaFitRange = 1.5;

//functions and a class I need.
TDirectory* resultsdirectory(TDirectory* motherdir, const char* name){
  motherdir->cd();
  TDirectory * outdir = motherdir->GetDirectory(name);
  if(outdir) motherdir->rmdir(name);
  outdir = motherdir->mkdir(name);
  return outdir;
}

//class for containing bin directories
class BinDirs : public TObject{
public:
  BinDirs(TFile* file, const char* bin, bool empty = false){
    if(file){
      Bin = TString(bin);
      file->cd();
      Samedir = file->GetDirectory(bin);
      if(!Samedir)Samedir =  file->mkdir(bin);
      Samedir->cd();
      if(!Samedir->GetDirectory("bin_stats")) Samedir->mkdir("bin_stats");
      if(!Samedir->GetDirectory("divided")) Samedir->mkdir("divided");
      if(empty) ::resultsdirectory(Samedir,"bin_stats");
      if(empty) ::resultsdirectory(Samedir,"divided");
      TDirectory * METAparentdir = file->GetDirectory("META");
      if(!METAparentdir)METAparentdir = file->mkdir("META");
      METAparentdir->cd();
      METAdir = METAparentdir->GetDirectory(bin);
      if(!METAdir)METAdir= METAparentdir->mkdir(bin);
      METAdir->cd();
      if(!METAdir->GetDirectory("bin_stats")) METAdir->mkdir("bin_stats");
      if(!METAdir->GetDirectory("divided")) METAdir->mkdir("divided");
      if(empty) ::resultsdirectory(METAdir,"bin_stats");
      if(empty) ::resultsdirectory(METAdir,"divided");
      TDirectory * METriggerparentdir = file->GetDirectory("METrigger");
      if(!METriggerparentdir)METriggerparentdir = file->mkdir("METrigger");
      METriggerparentdir->cd();
      METriggerdir = METriggerparentdir->GetDirectory(bin);
      if(!METriggerdir)METriggerdir= METriggerparentdir->mkdir(bin);
      METriggerdir->cd();
      if(!METriggerdir->GetDirectory("bin_stats")) METriggerdir->mkdir("bin_stats");
      if(!METriggerdir->GetDirectory("divided")) METriggerdir->mkdir("divided");
      if(empty) ::resultsdirectory(METriggerdir,"bin_stats");
      if(empty) ::resultsdirectory(METriggerdir,"divided");
      ready = true;
    }
    else{
      Samedir = NULL;
      METAdir = NULL;
      METriggerdir = NULL;
      Bin = TString("");
      ready = false;
    }
  }
  BinDirs(TDirectory* samed, TDirectory * METAd, TDirectory * METriggerd, bool empty = false){
    if(samed&&METAd&&METriggerd){
      Bin = TString(samed->GetName());
      samed->cd();
      Samedir = samed;
      if(!Samedir->GetDirectory("divided")) Samedir->mkdir("divided");
      if(empty) ::resultsdirectory(Samedir,"divided");
      METAd->cd();
      METAdir = METAd;
      if(!METAdir->GetDirectory("divided")) METAdir->mkdir("divided");
      if(empty) ::resultsdirectory(METAdir,"divided");
      METriggerd->cd();
      METriggerdir = METriggerd;
      if(!METriggerdir->GetDirectory("divided")) METriggerdir->mkdir("divided");
      if(empty) ::resultsdirectory(METriggerdir,"divided");
      ready = true;
    }
    else{
      Samedir = NULL;
      METAdir = NULL;
      METriggerdir = NULL;
      Bin = TString("");
      ready = false;
    }
  }  
//   ~BinDirs(){
//     delete Samedir;
//     delete METAdir;
//     delete METriggerdir;
//   }
  bool CompareTo(const char* name){return (Bin.CompareTo(name)==0);}
  TDirectory * Same(){return Samedir;}
  TDirectory * Samediv(){return Samedir->GetDirectory("divided");}
  TDirectory * SameDir(const char* name){return dynamic_cast<TDirectory*>(Samedir->Get(name));}
  TDirectory * META(){return METAdir;}
  TDirectory * METAdiv(){return METAdir->GetDirectory("divided");}
  TDirectory * METADir(const char* name){return dynamic_cast<TDirectory*>(METAdir->Get(name));}
  TDirectory * METrigger(){return METriggerdir;}
  TDirectory * METriggerdiv(){return METriggerdir->GetDirectory("divided");}
  TDirectory * METriggerDir(const char* name){return dynamic_cast<TDirectory*>(METriggerdir->Get(name));}
  BinDirs resultsdirectory(const char* dir){
    TDirectory * samedirdir = ::resultsdirectory(Samedir,dir);
    TDirectory * METAdirdir = ::resultsdirectory(METAdir,dir); 
    TDirectory * METriggerdirdir = ::resultsdirectory(METriggerdir,dir);
    return BinDirs(samedirdir,METAdirdir,METriggerdirdir);}
  const char* path(){return Samedir->GetPath();}
  const char* pathMETA(){return METAdir->GetPath();}
  const char* pathMETrigger(){return METriggerdir->GetPath();}
private:
  TDirectory * Samedir;
  TDirectory * METAdir;
  TDirectory * METriggerdir;
  TString    Bin;
  bool ready;
  ClassDef(BinDirs, 1);

};

  Double_t CGausPol0(Double_t * x, Double_t * par){
    //Parameterizatin of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1] ;
    Double_t s=par[2] ;
    Double_t dx=(x[0]-c)/s ;
    return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3] ;
  }

  Double_t CGausAPol1(Double_t * x, Double_t * par){
    //Parameterizatin of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1] ;
    Double_t s=par[2] ;
    Double_t dx=(x[0]-c)/s ;
    return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3]-abs(dx)*par[4] ;
  }
  Double_t CGausAPol2(Double_t * x, Double_t * par){
    //Parameterizatin of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1] ;
    Double_t s=par[2] ;
    Double_t dx=(x[0]-c)/s ;
    return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3]-abs(dx)*par[4] + dx*dx*par[5] ;
  }
TH2D * GetHistFrom3d(TH3D* hist,const char* which,const char* name)
{//Function to create the 2d histograms from the fully corrected 3d histogram.
  //hist has x=dEta_12, y=dPhi_1 and z=dPhi_2
  TString key;
  TH2D* hist2d;
  key = "DPHIDPHI";
  if(key.CompareTo(which)==0){
    hist2d = new TH2D(Form("DPHIDPHI%s",name),"#Delta#Phi for Associated 1 vs #Delta#Phi for Associated 2",50,-0.5*TMath::Pi(),1.5*TMath::Pi(),50 ,-0.5*TMath::Pi(),1.5*TMath::Pi());
    hist2d->Sumw2();
  }
  
  for(int x=1;x<hist->GetNbinsX();x++){
    for(int y=1;y<hist->GetNbinsY();y++){
      for(int z=1;z<hist->GetNbinsZ();z++){
	key = "DPHIDPHI";
	cout << "x="<<x<<", y="<<y<<", z="<<z<<" has "<< hist->GetBinContent(x,y,z) <<endl;
	if(key.CompareTo(which)==0){
	  hist2d->Fill(y,z,hist->GetBinContent(x,y,z));
	}
      }
    }
  }
  return hist2d;

}

Double_t Diagonal(TH2D* hist, Int_t xmin, Int_t xmax, Int_t ymin, Int_t ymax){
  Double_t out =0.0;
  if(xmin>xmax) return 0.0;
  if(ymin>ymax) return 0.0;
  int x = xmin;
  int y = ymin;
  while(x<=xmax&&y<=ymax){
    out += hist->GetBinContent(x,y);
    x++;y++;
  }
  return out;
}

Double_t gradient(TH2D* hist, Int_t minbin, Int_t maxbin){
  //calculates the gradient between the bins in the region and the 4 adjacent bins.
  Double_t out=0.0;
  for(int x=minbin+1;x<maxbin;x++){
    for(int y= minbin+1;y<maxbin;y++){
      out += TMath::Abs(hist->GetBinContent(x,y)-0.25*(hist->GetBinContent(x-1,y)+hist->GetBinContent(x+1,y)+hist->GetBinContent(x,y-1)+hist->GetBinContent(x,y+1)));
    }
  }
  return out;
}

TCanvas* Makecanvas(TH2D* hist, const char* name, Bool_t Stats)
{
  TCanvas * Canvas = new TCanvas(name);
  Canvas->Divide(2,2);
  if(!Stats)hist->SetStats(0);
  Canvas->cd(1);
  hist->Draw("surf2");
  Canvas->cd(2);
  hist->Draw("surf3");
  Canvas->cd(3);
  Int_t binpihn = hist->GetYaxis()->FindBin(0+0.2);
  Int_t binpiln = hist->GetYaxis()->FindBin(0-0.2);
  TH1D* projY = hist->ProjectionX(Form("%s%s",hist->GetName(),"_nearside"),binpiln,binpihn);
  if(!Stats) projY->SetStats(0);
//   projY->GetYaxis()->SetRangeUser(0., 1.1*projY->GetBinContent(projY->GetMaximumBin()));
  projY->SetTitle("Integral over the near side peak with #Delta#eta_{12} = 0#pm 0.2");
  projY->GetYaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  projY->SetTitleSize(0.04,"xy");
  projY->SetTitleOffset(1.05,"xy");
  projY->Draw("E");
  Canvas->cd(4);
  Int_t binpih = hist->GetYaxis()->FindBin(TMath::Pi()+0.2);
  Int_t binpil = hist->GetYaxis()->FindBin(TMath::Pi()-0.2);
  TH1D* projX = hist->ProjectionX(Form("%s%s",hist->GetName(),"_px"),binpil,binpih);
  if(!Stats) projX->SetStats(0);
  projX->SetTitle("Integral over a slice of the away side peak around with #Delta#eta_{12} = #pi#pm 0.2");
  projX->GetYaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  projX->SetTitleSize(0.04,"xy");
  projX->SetTitleOffset(1.05,"xy");
  projX->Draw("E");
  return Canvas;
}

TCanvas* Makecanvas(TH2D* histtopl, TH2D* histtopr, TH2D* histbotl, TH2D* histbotr, const char* name, Bool_t Stats)
{
  TCanvas * Canvas = new TCanvas(name);
  Canvas->Divide(2,2);
  Canvas->cd(1);
  if(!Stats)histtopl->SetStats(0);
  histtopl->Draw("surf3");
  Canvas->cd(2);
  if(!Stats)histtopr->SetStats(0);
  histtopr->Draw("surf3");
  Canvas->cd(3);
  if(!Stats)histbotl->SetStats(0);
  histbotl->Draw("surf3");
  Canvas->cd(4);
  if(!Stats)histbotr->SetStats(0);
  histbotr->Draw("surf3");
  return Canvas;
}

TStringToken GetHistTokens(TDirectory * dir){
  TList * histlist = dir->GetListOfKeys();
  TString * histlistS = new TString("");
  for(int i=0; i<histlist->GetEntries();i++)
  {
    TKey* key = dynamic_cast<TKey*>(histlist->At(i));
    if((!dynamic_cast<TH1D*>(key->ReadObj()))&&(!dynamic_cast<TH2D*>(key->ReadObj()))&&(!dynamic_cast<TH3D*>(key->ReadObj())))continue;//exclude non histograms
    if(dynamic_cast<TH3D*>(key->ReadObj()))continue;//exclude 3d histograms
    if(TString(key->GetName()).Contains("scaled"))continue;
    histlistS->Append(Form("%s ",histlist->At(i)->GetName()));
    }
  TStringToken histtokens(histlistS->Data()," ");
  return histtokens;
}

void fillwithbinnr(TH2D* fillfrom, TH1D* fillto, int binnr){
  for(int x=1;x<fillto->GetNbinsX();x++){
    fillto->SetBinContent(x,fillfrom->GetBinContent(x,binnr));
    fillto->SetBinError(x,fillfrom->GetBinError(x,binnr));
  }
}

void RemovePlateau(Double_t plateauheight, TH2D * hist){
  Double_t loccontent=0.0;
  for(int x = hist->GetXaxis()->FindBin(-1.6)+1;x<=hist->GetXaxis()->FindBin(1.6)-1;x++){
    for(int y = 0; y<hist->GetNbinsY()+1;y++){
      loccontent = hist->GetBinContent(x,y)-plateauheight;
      hist->SetBinContent(x,y,loccontent);
    }
  }
}

void extractbinyield(TDirectory* dir, TDirectory* yielddir, int type =0){
  if(!dir||!yielddir)return;
//   dir->pwd();
//   yielddir->GetPath();
  TDirectory * bindir = resultsdirectory(yielddir,"bins");

  TH2D* dphideta12ss;
  if(dynamic_cast<TH2D*>(dir->Get("DPhi_1_DEta_12_SameSide"))) dphideta12ss = dynamic_cast<TH2D*>(dir->Get("DPhi_1_DEta_12_SameSide")->Clone("DPhi_1_DEta_12_SameSidec"));
  else return;
  yielddir->cd();
  dphideta12ss->Write("DPhi_1_DEta_12_SameSide");

  TH1D* deta12ss = dphideta12ss->ProjectionX("DEta12");
  TH1D* deta12ssdraw = dphideta12ss->ProjectionX("DEta12d");
  deta12ss->Reset();
  TString title = TString("#Delta#eta_{12} distribution in bin");
  //flat background and a gaussian:
  int flatcolor = 2;
  const char* flatlable = "flat bg + gaussian";
  TH1D* yieldfg = dphideta12ss->ProjectionY("dphiyieldfg");
  yieldfg->Reset();
  yieldfg->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a flat background and a gaussian.");
  yieldfg->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  yieldfg->SetLineColor(flatcolor);
  TH1D* backgroundfg = dynamic_cast<TH1D*>(yieldfg->Clone("dphibackgroundfg"));
  backgroundfg->SetTitle("Height of the background as a function of #Delta#Phi#.");
  backgroundfg->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  backgroundfg->SetLineColor(flatcolor);
  TH1D* widthfg = dynamic_cast<TH1D*>(yieldfg->Clone("dphiwidthfg")); 
  widthfg->SetTitle("Width of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  widthfg->GetYaxis()->SetTitle("width of the peak (rad)");
  widthfg->SetLineColor(flatcolor);
  TH1D* peakposfg = dynamic_cast<TH1D*>(yieldfg->Clone("dphiposfg")); 
  peakposfg->SetTitle("Position of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  peakposfg->SetYTitle("Position of the peak (rad)");
  peakposfg->SetLineColor(flatcolor);
  TH1D* chisqfg = dynamic_cast<TH1D*>(yieldfg->Clone("dphichisqfg"));
  chisqfg->SetTitle("#Chi^2 of the fit as a function of #Delta#Phi");
  chisqfg->SetYTitle("#Chi^2 ");
  chisqfg->SetLineColor(flatcolor);
  TH1D* probfg = dynamic_cast<TH1D*>(yieldfg->Clone("dphiprobfg"));
  probfg->SetTitle("Probability of the fit as a function of #Delta#Phi");
  probfg->SetYTitle("Probability ");
  probfg->SetLineColor(flatcolor);
  //flat + gaussian:
  TF1 * fgs = new TF1("fgs",CGausPol0,-gkEtaFitRange,gkEtaFitRange,4) ;
  fgs->SetParNames("peakhight", "peakpos", "peakwidth", "B") ;
  fgs->SetLineColor(flatcolor);
  fgs->SetLineWidth(1);
  fgs->SetParameters(0.0,0.0,0.1,0.0);
  fgs->SetParLimits(0,0.0,0.3);
  fgs->SetParLimits(1,-0.1,0.1);
  fgs->SetParLimits(2,0.1,0.8);  
  //extract near and away width:
  fillwithbinnr(dphideta12ss,deta12ss,dphideta12ss->FindBin(0.0));
  TFitResultPtr fitresultfgn = deta12ss->Fit(fgs,"SQ","",-gkEtaFitRange,gkEtaFitRange);
  if(int(fitresultfgn)!=4000)//if "error", try again with more:
    fitresultfgn = deta12ss->Fit(fgs,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
  double widthnearfg = fgs->GetParameter(2);  
  deta12ss->Reset();
  fillwithbinnr(dphideta12ss,deta12ss,dphideta12ss->FindBin(TMath::Pi()));
  TFitResultPtr fitresultfga = deta12ss->Fit(fgs,"SQ","",-gkEtaFitRange,gkEtaFitRange);
  if(int(fitresultfga)!=4000)//if "error", try again with more:
    fitresultfga = deta12ss->Fit(fgs,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
  double widthawayfg = fgs->GetParameter(2);
  deta12ss->Reset();      
  //pol1 background and a gaussian:
  int p1color = 3;
  const char* p1lable = "pol1 bg + gaussian";
  TH1D* yieldp1g = dphideta12ss->ProjectionY("dphiyieldp1g");
  yieldp1g->Reset();
  yieldp1g->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a pol1 background and a gaussian.");
  yieldp1g->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  yieldp1g->SetLineColor(p1color);
  TH1D* backgroundp1g = dynamic_cast<TH1D*>(yieldfg->Clone("dphibackgroundp1g"));
  backgroundp1g->SetTitle("Height of the background as a function of #Delta#Phi#.");
  backgroundp1g->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  backgroundp1g->SetLineColor(p1color);
  TH1D* backgroundlp1g = dynamic_cast<TH1D*>(yieldfg->Clone("dphibackgroundlp1g"));
  backgroundlp1g->SetTitle("Slope of the background as a function of #Delta#Phi#.");
  backgroundlp1g->GetYaxis()->SetTitle("??#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  backgroundlp1g->SetLineColor(p1color);
  TH1D* widthp1g = dynamic_cast<TH1D*>(yieldfg->Clone("dphiwidthp1g")); 
  widthp1g->SetTitle("Width of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  widthp1g->GetYaxis()->SetTitle("width of the peak (rad)");
  widthp1g->SetLineColor(p1color);
  TH1D* peakposp1g = dynamic_cast<TH1D*>(yieldfg->Clone("dphiposp1g")); 
  peakposp1g->SetTitle("Position of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  peakposp1g->SetYTitle("Position of the peak (rad)");
  peakposp1g->SetLineColor(p1color);  
  TH1D* chisqp1g = dynamic_cast<TH1D*>(yieldfg->Clone("dphichisqp1g"));
  chisqp1g->SetTitle("#Chi^2 of the fit as a function of #Delta#Phi");
  chisqp1g->SetYTitle("#Chi^2 ");
  chisqp1g->SetLineColor(p1color);  
  TH1D* probp1g = dynamic_cast<TH1D*>(yieldfg->Clone("dphiprobp1g"));
  probp1g->SetTitle("Probability of the fit as a function of #Delta#Phi");
  probp1g->SetYTitle("Probability ");
  probp1g->SetLineColor(p1color);
  //pol1 + gaussian:
  TF1 * p1gs = new TF1("p1gs",CGausAPol1,-gkEtaFitRange,gkEtaFitRange,5) ;
  p1gs->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C") ;
  p1gs->SetLineColor(p1color);
  p1gs->SetLineWidth(1);
  p1gs->SetParameters(0.0,0.0,0.1,0.0,0.0);
  p1gs->SetParLimits(0,0.0,0.3);
  p1gs->SetParLimits(1,-0.1,0.1);
  p1gs->SetParLimits(2,0.1,0.8);  
  p1gs->SetParLimits(4,0.0,1.0);
  //extract near and away width:
  fillwithbinnr(dphideta12ss,deta12ss,dphideta12ss->FindBin(0.0));
  TFitResultPtr fitresultp1gn = deta12ss->Fit(p1gs,"SQ","",-gkEtaFitRange,gkEtaFitRange);
  if(int(fitresultp1gn)!=4000)//if "error", try again with more:
    fitresultp1gn = deta12ss->Fit(p1gs,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
  double widthnearp1g = p1gs->GetParameter(2);  
  deta12ss->Reset();
  fillwithbinnr(dphideta12ss,deta12ss,dphideta12ss->FindBin(TMath::Pi()));
  TFitResultPtr fitresultp1ga = deta12ss->Fit(p1gs,"SQ","",-1.5,1.5);
  if(int(fitresultp1ga)!=4000)//if "error", try again with more:
    fitresultp1ga = deta12ss->Fit(p1gs,"MSQ","",-1.5,1.5);
  double widthawayp1g = p1gs->GetParameter(2);
  deta12ss->Reset();      
  //pol2 background and a gaussian:
  int p2color = 4;
  const char* p2lable = "pol2 bg + gaussian";
  TH1D* yieldp2g = dphideta12ss->ProjectionY("dphiyieldp2g");
  yieldp2g->Reset();
  yieldp2g->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a pol1 background and a gaussian.");
  yieldp2g->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  yieldp2g->SetLineColor(p2color);
  TH1D* backgroundp2g = dynamic_cast<TH1D*>(yieldfg->Clone("dphibackgroundp2g"));
  backgroundp2g->SetTitle("Height of the background as a function of #Delta#Phi#.");
  backgroundp2g->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  backgroundp2g->SetLineColor(p2color);
  TH1D* backgroundlp2g = dynamic_cast<TH1D*>(yieldfg->Clone("dphibackgroundlp2g"));
  backgroundlp2g->SetTitle("Slope of the background as a function of #Delta#Phi#.");
  backgroundlp2g->GetYaxis()->SetTitle("??#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  backgroundlp2g->SetLineColor(p2color);  
  TH1D* backgroundcp2g = dynamic_cast<TH1D*>(yieldfg->Clone("dphibackgroundcp2g"));
  backgroundcp2g->SetTitle("Curvature of the background as a function of #Delta#Phi#.");
  backgroundcp2g->GetYaxis()->SetTitle("??#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  backgroundcp2g->SetLineColor(p2color);  
  TH1D* widthp2g = dynamic_cast<TH1D*>(yieldfg->Clone("dphiwidthp2g")); 
  widthp2g->SetTitle("Width of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  widthp2g->GetYaxis()->SetTitle("width of the peak (rad)");
  widthp2g->SetLineColor(p2color);
  TH1D* peakposp2g = dynamic_cast<TH1D*>(yieldfg->Clone("dphiposp2g")); 
  peakposp2g->SetTitle("Position of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  peakposp2g->SetYTitle("Position of the peak (rad)");
  peakposp2g->SetLineColor(p2color);
  TH1D* chisqp2g = dynamic_cast<TH1D*>(yieldfg->Clone("dphichisqp2g"));
  chisqp2g->SetTitle("#Chi^2 of the fit as a function of #Delta#Phi");
  chisqp2g->SetYTitle("#Chi^2 ");
  chisqp2g->SetLineColor(p2color);
  TH1D* probp2g = dynamic_cast<TH1D*>(yieldfg->Clone("dphiprobp2g"));
  probp2g->SetTitle("Probability of the fit as a function of #Delta#Phi");
  probp2g->SetYTitle("Probability ");
  probp2g->SetLineColor(p2color);
  //pol2 + gaussian:
  TF1 * p2gs = new TF1("p2gs",CGausAPol2,-gkEtaFitRange,gkEtaFitRange,6) ;
  p2gs->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C", "D") ;
  p2gs->SetLineColor(p2color);
  p2gs->SetLineWidth(1);
  p2gs->SetParameters(0.0,0.0,0.1,0.0,0.0);
  p2gs->SetParLimits(0,0.0,0.3);
  p2gs->SetParLimits(1,-0.1,0.1);
  p2gs->SetParLimits(2,0.1,0.8);  
  p2gs->SetParLimits(4,0.0,1.0);
//   p2gs->SetParLimits(5,0.0,1.0);
  //extract near and away width:
  fillwithbinnr(dphideta12ss,deta12ss,dphideta12ss->FindBin(0.0));
  TFitResultPtr fitresultp2gn = deta12ss->Fit(p2gs,"SQ","",-gkEtaFitRange,gkEtaFitRange);
  if(int(fitresultp2gn)!=4000)//if "error", try again with more:
    fitresultp2gn = deta12ss->Fit(p2gs,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
  double widthnearp2g = p2gs->GetParameter(2);  
  deta12ss->Reset();
  fillwithbinnr(dphideta12ss,deta12ss,dphideta12ss->FindBin(TMath::Pi()));
  TFitResultPtr fitresultp2ga = deta12ss->Fit(p2gs,"SQ","",-1.5,1.5);
  if(int(fitresultp2ga)!=4000)//if "error", try again with more:
    fitresultp2ga = deta12ss->Fit(p2gs,"MSQ","",-1.5,1.5);
  double widthawayp2g = p2gs->GetParameter(2);
  deta12ss->Reset();   

  double deltaphi = yieldfg->GetXaxis()->GetBinCenter(2)-yieldfg->GetXaxis()->GetBinCenter(1);
  bindir->cd();
  for(int dphi=1;dphi<=dphideta12ss->GetNbinsY();dphi++){
    deta12ss->SetTitle(Form("%s %4.2f #pi < #Delta#Phi < %4.2f #pi",title.Data(),dphideta12ss->GetYaxis()->GetBinLowEdge(dphi)/TMath::Pi(),dphideta12ss->GetYaxis()->GetBinUpEdge(dphi)/TMath::Pi()));
    fillwithbinnr(dphideta12ss,deta12ss,dphi);
    deta12ss->SetStats(false);
    
    if(dphi<TMath::Pi()/2.0) fgs->FixParameter(2,widthnearfg);
    if(dphi>=TMath::Pi()/2.0) fgs->FixParameter(2,widthawayfg);
    fgs->FixParameter(1,0.0);
    TFitResultPtr fitresult = deta12ss->Fit(fgs,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    if(int(fitresult)!=4000)//if "error", try again with more:
      fitresult = deta12ss->Fit(fgs,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
    deta12ss->Write(Form("%sflatg_%i",deta12ss->GetName(),dphi));
//     fgs->Draw("same");
    if(dphi<TMath::Pi()/2.0) p1gs->FixParameter(2,widthnearp1g);
    if(dphi>=TMath::Pi()/2.0) p1gs->FixParameter(2,widthawayp1g);
    p1gs->FixParameter(1,0.0);
    TFitResultPtr fitresultp1 = deta12ss->Fit(p1gs,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    if(int(fitresultp1)!=4000)//if "error", try again with more:
      fitresultp1 = deta12ss->Fit(p1gs,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
//     p1gs->Draw("same");
    deta12ss->Write(Form("%sp1g_%i",deta12ss->GetName(),dphi));
    if(dphi<TMath::Pi()/2.0) p2gs->FixParameter(2,widthnearp2g);
    if(dphi>=TMath::Pi()/2.0) p2gs->FixParameter(2,widthawayp2g);
    p2gs->FixParameter(1,0.0);
    TFitResultPtr fitresultp2 = deta12ss->Fit(p2gs,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    if(int(fitresultp2)!=4000)//if "error", try again with more:
      fitresultp2 = deta12ss->Fit(p2gs,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
//     p2gs->Draw("same");    
    deta12ss->Write(Form("%sp2g_%i",deta12ss->GetName(),dphi));
    
    TCanvas * bincanvas = new TCanvas(Form("%sCanvas_%i",deta12ss->GetName(),dphi));
    deta12ssdraw->SetTitle(Form("%s %4.2f #pi < #Delta#Phi < %4.2f #pi",title.Data(),dphideta12ss->GetYaxis()->GetBinLowEdge(dphi)/TMath::Pi(),dphideta12ss->GetYaxis()->GetBinUpEdge(dphi)/TMath::Pi()));
    deta12ssdraw->Reset();
    fillwithbinnr(dphideta12ss,deta12ssdraw,dphi);
    deta12ssdraw->SetStats(false);    
    deta12ssdraw->Draw("ESAME");

    TF1 * plfgs = new TF1("plfgs",CGausPol0,-gkEtaFitRange,gkEtaFitRange,4) ;
    plfgs->SetParNames("peakhight", "peakpos", "peakwidth", "B") ;
    plfgs->SetLineColor(flatcolor);
    plfgs->SetLineWidth(1);
    plfgs->SetParameters(fgs->GetParameter(0),fgs->GetParameter(1),fgs->GetParameter(2),fgs->GetParameter(3));
    plfgs->Draw("LSAME");
    TF1 * plp1gs = new TF1("plp1gs",CGausAPol1,-gkEtaFitRange,gkEtaFitRange,5) ;
    plp1gs->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C") ;
    plp1gs->SetLineColor(p1color);
    plp1gs->SetLineWidth(1);
    plp1gs->SetParameters(p1gs->GetParameter(0),p1gs->GetParameter(1),p1gs->GetParameter(2),p1gs->GetParameter(3),p1gs->GetParameter(4));
    plp1gs->Draw("LSAME");
    TF1 * plp2gs = new TF1("plp2gs",CGausAPol2,-gkEtaFitRange,gkEtaFitRange,6) ;
    plp2gs->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C","D") ;
    plp2gs->SetLineColor(p2color);
    plp2gs->SetLineWidth(1);
    plp2gs->SetParameters(p2gs->GetParameter(0),p2gs->GetParameter(1),p2gs->GetParameter(2),p2gs->GetParameter(3),p2gs->GetParameter(4),p2gs->GetParameter(5));
    plp2gs->Draw("LSAME");    
//     double diffmaxmin =  deta12ssdraw->GetBinContent(deta12ssdraw->FindBin(0.0)) - deta12ssdraw->GetBinContent(deta12ssdraw->FindBin(-1.0));
//     deta12ssdraw->SetMinimum(deta12ssdraw->GetBinContent(deta12ssdraw->GetXaxis()->FindBin(-1.5)) - 1.5*deta12ssdraw->GetBinError(deta12ssdraw->GetXaxis()->FindBin(-1.5)));//deta12ssdraw->GetBinContent(deta12ssdraw->GetMinimumBin())-10*diffmaxmin);
//     deta12ssdraw->SetMaximum(deta12ssdraw->GetBinContent(deta12ssdraw->GetXaxis()->FindBin(0.0)) + 1.5*deta12ssdraw->GetBinError(deta12ssdraw->GetXaxis()->FindBin(0.0)));
//    deta12ssdraw->SetMaximum();//deta12ssdraw->GetBinContent(deta12ssdraw->GetMaximumBin())-10*diffmaxmin);

   TLegend * leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->SetHeader("Different fits:");
   leg->AddEntry(plfgs,flatlable);
   leg->AddEntry(plp1gs,p1lable);
   leg->AddEntry(plp2gs,p2lable);
   leg->Draw("SAME");    
    bincanvas->Update();
    bincanvas->Write();
    delete bincanvas;
//     deta12ss->Write(Form("%s_%i",deta12ss->GetName(),dphi));
//     if(!(int(fitresult)%4000)){
    yieldfg->SetBinContent(dphi,fgs->GetParameter(0)/deltaphi);
    yieldfg->SetBinError(dphi,fgs->GetParError(0)/deltaphi);
    widthfg->SetBinContent(dphi,fgs->GetParameter(2)/deltaphi);
    widthfg->SetBinError(dphi,fgs->GetParError(2)/deltaphi);
    peakposfg->SetBinContent(dphi,fgs->GetParameter(1)/deltaphi);
    peakposfg->SetBinError(dphi,fgs->GetParError(1)/deltaphi);
    backgroundfg->SetBinContent(dphi,fgs->GetParameter(3)/deltaphi);
    backgroundfg->SetBinError(dphi,fgs->GetParError(3)/deltaphi);
    chisqfg->SetBinContent(dphi,fgs->GetChisquare());
    probfg->SetBinContent(dphi,fgs->GetProb());
//     }
//     if(!int(fitresultp1)%4000){
    yieldp1g->SetBinContent(dphi,p1gs->GetParameter(0)/deltaphi);
    yieldp1g->SetBinError(dphi,p1gs->GetParError(0)/deltaphi);
    widthp1g->SetBinContent(dphi,p1gs->GetParameter(2)/deltaphi);
    widthp1g->SetBinError(dphi,p1gs->GetParError(2)/deltaphi);
    peakposp1g->SetBinContent(dphi,p1gs->GetParameter(1)/deltaphi);
    peakposp1g->SetBinError(dphi,p1gs->GetParError(1)/deltaphi);
    backgroundp1g->SetBinContent(dphi,p1gs->GetParameter(3)/deltaphi);
    backgroundp1g->SetBinError(dphi,p1gs->GetParError(3)/deltaphi);
    backgroundlp1g->SetBinContent(dphi,p1gs->GetParameter(4)/deltaphi);
    backgroundlp1g->SetBinError(dphi,p1gs->GetParError(4)/deltaphi);
    chisqp1g->SetBinContent(dphi,p1gs->GetChisquare());
    probp1g->SetBinContent(dphi,p1gs->GetProb());
//     }    
//     if(!int(fitresultp2)%4000){
    yieldp2g->SetBinContent(dphi,p2gs->GetParameter(0)/deltaphi);
    yieldp2g->SetBinError(dphi,p2gs->GetParError(0)/deltaphi);
    widthp2g->SetBinContent(dphi,p2gs->GetParameter(2)/deltaphi);
    widthp2g->SetBinError(dphi,p2gs->GetParError(2)/deltaphi);
    peakposp2g->SetBinContent(dphi,p2gs->GetParameter(1)/deltaphi);
    peakposp2g->SetBinError(dphi,p2gs->GetParError(1)/deltaphi);
    backgroundp2g->SetBinContent(dphi,p2gs->GetParameter(3)/deltaphi);
    backgroundp2g->SetBinError(dphi,p2gs->GetParError(3)/deltaphi);
    backgroundlp2g->SetBinContent(dphi,p2gs->GetParameter(4)/deltaphi);
    backgroundlp2g->SetBinError(dphi,p2gs->GetParError(4)/deltaphi);
    backgroundcp2g->SetBinContent(dphi,p2gs->GetParameter(5)/deltaphi);
    backgroundcp2g->SetBinError(dphi,p2gs->GetParError(5)/deltaphi);
    chisqp2g->SetBinContent(dphi,p2gs->GetChisquare());
    probp2g->SetBinContent(dphi,p2gs->GetProb());
//     }        
    deta12ss->Reset();
  }
  
  yielddir->cd();
  yieldfg->Write();
  backgroundfg->Write();
  widthfg->Write();
  peakposfg->Write();
  chisqfg->Write();
  probfg->Write();
  yieldp1g->Write();
  backgroundp1g->Write();
  widthp1g->Write();
  peakposp1g->Write();
  backgroundlp1g->Write();
  chisqp1g->Write();
  probp1g->Write();
  yieldp2g->Write();
  backgroundp2g->Write();
  widthp2g->Write();
  peakposp2g->Write();
  backgroundlp2g->Write();
  chisqp2g->Write();
  probp2g->Write();
  
  TCanvas * totcanvas = new TCanvas("YieldCanvas");
  totcanvas->cd();
  TPad * histpad = new TPad("histpad","yields in different fits",0.05,0.3,0.95,0.95);
  histpad->SetFillColor(16);
  histpad->cd();
  TLegend * histleg = new TLegend(0.62,0.7,0.9,0.9);
  histleg->SetHeader("Different fits:");
  histleg->AddEntry(yieldfg,flatlable);
  histleg->AddEntry(yieldp1g,p1lable);
  histleg->AddEntry(yieldp2g,p2lable);
  yieldfg->SetStats(false);
  yieldfg->Draw("E");
  yieldp1g->SetStats(false);
  yieldp1g->Draw("ESAME");
  yieldp2g->SetStats(false);
  yieldp2g->Draw("ESAME");
  histleg->Draw("SAME");    
  totcanvas->cd();
  TPad* chipad = new TPad("chisqpad","#Chi^2 in different fits.",0.05,0.05,0.95,0.3);
  chipad->SetFillColor(16);
  chipad->cd();
  gPad->SetLogy();
  chisqfg->SetStats(false);
  chisqfg->Draw("E");
  chisqp1g->SetStats(false);
  chisqp1g->Draw("ESAME");
  chisqp2g->SetStats(false);
  chisqp2g->Draw("ESAME");
  totcanvas->cd();
  histpad->Draw();
  chipad->Draw();
  totcanvas->Update();
  totcanvas->Write();
  
//   delete yielddir;
//   delete dphideta12ss; delete fgs;
//   delete p1gs; delete p2gs;
}

TList * GetMZDirectories(BinDirs* same){
  TList * keys = same->Same()->GetListOfKeys();
  TList * directories = new TList();
  for(int i=0; i<keys->GetEntries();i++){if(TString(keys->At(i)->GetName()).BeginsWith("BinM(")&&TString(keys->At(i)->GetName()).Contains("Z(")) directories->Add(keys->At(i) );}
  return directories;
}

void canvasmaker(const char* histname, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs * Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs * Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs * Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  if(TString(histname).CompareTo("DPhi_1_DPHI_2_far")==0){
    TH2D* hist;TH2D* histnear;TH2D* histmid; TH2D* histfar;
    TCanvas* canvas;
    hist = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2"));
    histnear = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2_near"));
    histmid = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2_mid"));
    histfar = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2_far"));
    if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);All->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2"));
    histnear = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2_near"));
    histmid = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2_mid"));
    histfar = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2_far"));
    if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin1->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2"));
    histnear = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_near"));
    histmid = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_mid"));
    histfar = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_far"));
    if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin2->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2"));
    histnear = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_near"));
    histmid = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_mid"));
    histfar = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_far"));
    if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin3->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2"));
    histnear = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_near"));
    histmid = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_mid"));
    histfar = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_far"));
    if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin4->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2"));
    histnear = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_near"));
    histmid = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_mid"));
    histfar = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_far"));
    if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin5->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
  }
  else if(TString(histname).Contains("DPhi_1_DEta")||TString(histname).Contains("DPhi_DEta")){
    TH2D* hist;
    TCanvas* canvas;
    hist = dynamic_cast<TH2D*>(All->Samediv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(All->METAdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(All->METriggerdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin1->Samediv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin1->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin1->METAdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin1->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin1->METriggerdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin1->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin2->Samediv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin2->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin2->METAdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin2->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin2->METriggerdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin2->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin3->Samediv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin3->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin3->METAdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin3->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin3->METriggerdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin3->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin4->Samediv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin4->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin4->METAdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin4->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin4->METriggerdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin4->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin5->Samediv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin5->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin5->METAdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin5->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
    hist = dynamic_cast<TH2D*>(Bin5->METriggerdiv()->Get(histname));
    if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin5->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
  }
  return;
}

void savedircontent(BinDirs* same, double METAScale, double METriggerScale, const char* iteration){
  TPaveText * whatinthisbin;
  TCanvas * whatinthisbincanvas = new TCanvas("Thisversion");
  whatinthisbincanvas->cd();
  whatinthisbin = new TPaveText(.05,.1,.95,.8);
  whatinthisbin->AddText("The corrected same event histogram");
  whatinthisbin->AddText(Form("- %4.2f times META     ",METAScale));
  whatinthisbin->AddText(Form("- %4.2f times METrigger",METriggerScale));
  whatinthisbin->Draw();
  same->SameDir(iteration)->cd();
  whatinthisbincanvas->Write();
  delete whatinthisbincanvas;
  delete whatinthisbin;
}

void CollectHistbs(TH1D* histo, TList * directories, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs * Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs * Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs * Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH1D* hist 			= dynamic_cast<TH1D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH1D* histMETA 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH1D* histMETrigger 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  //Bin 1:
  TH1D* histbin1 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH1D* histMETAbin1 	  	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH1D* histMETriggerbin1 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  //Bin 2:
  TH1D* histbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
  TH1D* histMETAbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
  TH1D* histMETriggerbin2 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));
  //Bin 3:
  TH1D* histbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
  TH1D* histMETAbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
  TH1D* histMETriggerbin3 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  //Bin 4:
  TH1D* histbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
  TH1D* histMETAbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
  TH1D* histMETriggerbin4 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  //Bin 5:
  TH1D* histbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
  TH1D* histMETAbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
  TH1D* histMETriggerbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  //Doubles to hold the values until they are put into the hists.
  TH1D* histtmp;TH1D* histMETAtmp;TH1D* histMETriggertmp;
  //loop over the bins:
  for(int i=0;i<directories->GetEntries();i++){
    int Mbin = 0;
    //find the Multiplicity bin we are in
    for(int j = 1;j<=5;j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
    
    //test if the relevant histogram exists in this bin:
    if(!dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/bin_stats",directories->At(i)->GetName()))->Get(histo->GetName())))
    {// 	cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
      continue;}
    //extract histogram in relevant bin:
    histtmp 		= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/bin_stats",directories->At(i)->GetName()))->Get(histo->GetName()));
    histMETAtmp 	= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/bin_stats",directories->At(i)->GetName()))->Get(histo->GetName()));
    histMETriggertmp 	= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/bin_stats",directories->At(i)->GetName()))->Get(histo->GetName()));
    hist->Add(histtmp); histMETA->Add(histMETAtmp); histMETrigger->Add(histMETriggertmp);
    if(Mbin == 1){	histbin1->Add(hist);histMETAbin1->Add(histMETAtmp); histMETriggerbin1->Add(histMETriggertmp);}
    else if(Mbin == 2){	histbin2->Add(hist);histMETAbin2->Add(histMETAtmp); histMETriggerbin2->Add(histMETriggertmp);}	      
    else if(Mbin == 3){	histbin3->Add(hist);histMETAbin3->Add(histMETAtmp); histMETriggerbin3->Add(histMETriggertmp);}
    else if(Mbin == 4){	histbin4->Add(hist);histMETAbin4->Add(histMETAtmp); histMETriggerbin4->Add(histMETriggertmp);}	      
    else if(Mbin == 5){ histbin5->Add(hist);histMETAbin5->Add(histMETAtmp); histMETriggerbin5->Add(histMETriggertmp);}
  }
  //save the histograms in the relevant directories:
  All->Same()->GetDirectory("bin_stats")->cd();
  hist->Write(histo->GetName());
  delete hist;
  All->META()->GetDirectory("bin_stats")->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  All->METrigger()->GetDirectory("bin_stats")->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
  Bin1->Same()->GetDirectory("bin_stats")->cd();
  histbin1->Write(histo->GetName());
  delete histbin1;
  Bin1->META()->GetDirectory("bin_stats")->cd();
  histMETAbin1->Write(histo->GetName());
  delete histMETAbin1;
  Bin1->METrigger()->GetDirectory("bin_stats")->cd();
  histMETriggerbin1->Write(histo->GetName());
  delete histMETriggerbin1;  
  Bin2->Same()->GetDirectory("bin_stats")->cd();
  histbin2->Write(histo->GetName());
  delete histbin2;
  Bin2->META()->GetDirectory("bin_stats")->cd();
  histMETAbin2->Write(histo->GetName());
  delete histMETAbin2;
  Bin2->METrigger()->GetDirectory("bin_stats")->cd();
  histMETriggerbin2->Write(histo->GetName());
  delete histMETriggerbin2;  
  Bin3->Same()->GetDirectory("bin_stats")->cd();
  histbin3->Write(histo->GetName());
  delete histbin3;
  Bin3->META()->GetDirectory("bin_stats")->cd();
  histMETAbin3->Write(histo->GetName());
  delete histMETAbin3;
  Bin3->METrigger()->GetDirectory("bin_stats")->cd();
  histMETriggerbin3->Write(histo->GetName());
  delete histMETriggerbin3;  
  Bin4->Same()->GetDirectory("bin_stats")->cd();
  histbin4->Write(histo->GetName());
  delete histbin4;
  Bin4->META()->GetDirectory("bin_stats")->cd();
  histMETAbin4->Write(histo->GetName());
  delete histMETAbin4;
  Bin4->METrigger()->GetDirectory("bin_stats")->cd();
  histMETriggerbin4->Write(histo->GetName());
  delete histMETriggerbin4;  
  Bin5->Same()->GetDirectory("bin_stats")->cd();
  histbin5->Write(histo->GetName());
  delete histbin5;
  Bin5->META()->GetDirectory("bin_stats")->cd();
  histMETAbin5->Write(histo->GetName());
  delete histMETAbin5;
  Bin5->METrigger()->GetDirectory("bin_stats")->cd();
  histMETriggerbin5->Write(histo->GetName());
  delete histMETriggerbin5;  
}

void CollectHist(TH1D* histo, TList * directories, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs * Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs * Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs * Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH1D* hist 			= dynamic_cast<TH1D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH1D* histMETA 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH1D* histMETrigger 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  //Bin 1:
  TH1D* histbin1 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH1D* histMETAbin1 	  	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH1D* histMETriggerbin1 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  //Bin 2:
  TH1D* histbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
  TH1D* histMETAbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
  TH1D* histMETriggerbin2 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));
  //Bin 3:
  TH1D* histbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
  TH1D* histMETAbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
  TH1D* histMETriggerbin3 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  //Bin 4:
  TH1D* histbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
  TH1D* histMETAbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
  TH1D* histMETriggerbin4 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  //Bin 5:
  TH1D* histbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
  TH1D* histMETAbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
  TH1D* histMETriggerbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  //Doubles to hold the values until they are put into the hists.
  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;
  Double_t BinContent = 0.0;Double_t BinError   = 0.0; Double_t BinContentMETA = 0.0;Double_t BinErrorMETA   = 0.0;
  Double_t BinContentMETrigger = 0.0;Double_t BinErrorMEtrigger   = 0.0; Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  //loop over the bins:
  for(int x=0;x<=histo->GetNbinsX()+1;x++){
    //loop over all Multiplicity-Vertex bins:
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
      //find the Multiplicity bin we are in
      for(int j = 1;j<=5;j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName())))
      {
// 	cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
	continue;}
      //extract bin content and error in the relevant bin:
      bincontl 			= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x);
      binerrorl 		= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x);	
      bincontlMETA 		= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x);
      binerrorlMETA 		= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x);	      
      bincontlMEtrigger 	= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x);
      binerrorlMETrigger 	= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x);      
      if(bincontl>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContent += bincontl/(binerrorl*binerrorl);
	BinError   += 1.0/(binerrorl*binerrorl);
	if(Mbin == 1){		BinContentbin1 += bincontl/(binerrorl*binerrorl); BinErrorbin1   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 2){	BinContentbin2 += bincontl/(binerrorl*binerrorl); BinErrorbin2   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 3){	BinContentbin3 += bincontl/(binerrorl*binerrorl); BinErrorbin3   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 4){	BinContentbin4 += bincontl/(binerrorl*binerrorl); BinErrorbin4   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 5){ 	BinContentbin5 += bincontl/(binerrorl*binerrorl); BinErrorbin5   += 1.0/(binerrorl*binerrorl);}
      }
      else{
	BinError += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorbin1   += 1.0;
	else if(Mbin == 2)	BinErrorbin2   += 1.0;
	else if(Mbin == 3)	BinErrorbin3   += 1.0;
	else if(Mbin == 4)	BinErrorbin4   += 1.0;	      
	else if(Mbin == 5)	BinErrorbin5   += 1.0;
      }
      if(bincontlMETA>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETA += bincontlMETA/(binerrorlMETA*binerrorlMETA);
	BinErrorMETA   += 1.0/(binerrorlMETA*binerrorlMETA);
	if(Mbin == 1){		BinContentMETAbin1 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin1   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 2){	BinContentMETAbin2 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin2   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 3){	BinContentMETAbin3 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin3   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 4){	BinContentMETAbin4 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin4   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 5){	BinContentMETAbin5 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin5   += 1.0/(binerrorlMETA*binerrorlMETA);}
      }
      else{
	BinErrorMETA += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMETAbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMETAbin2   += 1.0;	      
	else if(Mbin == 3)	BinErrorMETAbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMETAbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMETAbin5   += 1.0;
      }
      if(bincontlMEtrigger>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETrigger += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger);
	BinErrorMEtrigger   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);
	if(Mbin == 1){		BinContentMETriggerbin1 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin1   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 2){	BinContentMETriggerbin2 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin2   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 3){	BinContentMETriggerbin3 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin3   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 4){	BinContentMETriggerbin4 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin4   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 5){	BinContentMETriggerbin5 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin5   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
      }
      else{
	BinErrorMEtrigger += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMEtriggerbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMEtriggerbin2   += 1.0;
	else if(Mbin == 3)	BinErrorMEtriggerbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMEtriggerbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMEtriggerbin5   += 1.0;
      }
      //reset the local bins to zero just to be shure:
      bincontl  = 0.0;		binerrorl = 0.0;
      bincontlMETA =0.0;	binerrorlMETA=0.0;
      bincontlMEtrigger=0.0;	binerrorlMETrigger=0.0;
    }//end loop over M-V bins
    //normalize the bin and the error for same:
    if(BinError>1.0e-10){		BinContent = BinContent/BinError;						BinError = 1.0/TMath::Sqrt(BinError);				}
    else{				BinContent=0.0;									BinError=0.0;							}
    if(BinErrorbin1>1.0e-10){		BinContentbin1 = BinContentbin1/BinErrorbin1;					BinErrorbin1 = 1.0/TMath::Sqrt(BinErrorbin1);			}
    else{				BinContentbin1=0.0;								BinErrorbin1=0.0;						}
    if(BinErrorbin2>1.0e-10){		BinContentbin2 = BinContentbin2/BinErrorbin2;					BinErrorbin2 = 1.0/TMath::Sqrt(BinErrorbin2);			}
    else{				BinContentbin2=0.0;								BinErrorbin2=0.0;						}
    if(BinErrorbin3>1.0e-10){		BinContentbin3 = BinContentbin3/BinErrorbin3;					BinErrorbin3 = 1.0/TMath::Sqrt(BinErrorbin3);			}
    else{				BinContentbin3=0.0;								BinErrorbin3=0.0;						}
    if(BinErrorbin4>1.0e-10){		BinContentbin4 = BinContentbin4/BinErrorbin4;					BinErrorbin4 = 1.0/TMath::Sqrt(BinErrorbin4);			}
    else{				BinContentbin4=0.0;								BinErrorbin4=0.0;						}
    if(BinErrorbin5>1.0e-10){		BinContentbin5 = BinContentbin5/BinErrorbin5;					BinErrorbin5 = 1.0/TMath::Sqrt(BinErrorbin5);			}
    else{				BinContentbin5=0.0;								BinErrorbin5=0.0;						}
    //normalize the bin and the error for META:    
    if(BinErrorMETA>1.0e-10){		BinContentMETA = BinContentMETA/BinErrorMETA;					BinErrorMETA = 1.0/TMath::Sqrt(BinErrorMETA);			}
    else{				BinContentMETA=0.0;								BinErrorMETA=0.0;						}
    if(BinErrorMETAbin1>1.0e-10){	BinContentMETAbin1 = BinContentMETAbin1/BinErrorMETAbin1;			BinErrorMETAbin1 = 1.0/TMath::Sqrt(BinErrorMETAbin1);		}
    else{				BinContentMETAbin1=0.0;								BinErrorMETAbin1=0.0;						}
    if(BinErrorMETAbin2>1.0e-10){	BinContentMETAbin2 = BinContentMETAbin2/BinErrorMETAbin2;			BinErrorMETAbin2 = 1.0/TMath::Sqrt(BinErrorMETAbin2);		}
    else{				BinContentMETAbin2=0.0;								BinErrorMETAbin2=0.0;						}
    if(BinErrorMETAbin3>1.0e-10){	BinContentMETAbin3 = BinContentMETAbin3/BinErrorMETAbin3;			BinErrorMETAbin3 = 1.0/TMath::Sqrt(BinErrorMETAbin3);		}
    else{				BinContentMETAbin3=0.0;								BinErrorMETAbin3=0.0;						}
    if(BinErrorMETAbin4>1.0e-10){	BinContentMETAbin4 = BinContentMETAbin4/BinErrorMETAbin4;			BinErrorMETAbin4 = 1.0/TMath::Sqrt(BinErrorMETAbin4);		}
    else{				BinContentMETAbin4=0.0;								BinErrorMETAbin4=0.0;						}
    if(BinErrorMETAbin5>1.0e-10){	BinContentMETAbin5 = BinContentMETAbin5/BinErrorMETAbin5;			BinErrorMETAbin5 = 1.0/TMath::Sqrt(BinErrorMETAbin5);		}
    else{				BinContentMETAbin5=0.0;								BinErrorMETAbin5=0.0;						}
    //normalize the bin and the error for METrigger:    
    if(BinErrorMEtrigger>1.0e-10){	BinContentMETrigger = BinContentMETrigger/BinErrorMEtrigger;			BinErrorMEtrigger = 1.0/TMath::Sqrt(BinErrorMEtrigger);		}
    else{				BinContentMETrigger=0.0;							BinErrorMEtrigger=0.0;						}
    if(BinErrorMEtriggerbin1>1.0e-10){	BinContentMETriggerbin1 = BinContentMETriggerbin1/BinErrorMEtriggerbin1;	BinErrorMEtriggerbin1 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin1);	}
    else{				BinContentMETriggerbin1=0.0;							BinErrorMEtriggerbin1=0.0;					}
    if(BinErrorMEtriggerbin2>1.0e-10){	BinContentMETriggerbin2 = BinContentMETriggerbin2/BinErrorMEtriggerbin2;	BinErrorMEtriggerbin2 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin2);	}
    else{				BinContentMETriggerbin2=0.0;							BinErrorMEtriggerbin2=0.0;					}
    if(BinErrorMEtriggerbin3>1.0e-10){	BinContentMETriggerbin3 = BinContentMETriggerbin3/BinErrorMEtriggerbin3;	BinErrorMEtriggerbin3 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin3);	}
    else{				BinContentMETriggerbin3=0.0;							BinErrorMEtriggerbin3=0.0;					}
    if(BinErrorMEtriggerbin4>1.0e-10){	BinContentMETriggerbin4 = BinContentMETriggerbin4/BinErrorMEtriggerbin4;	BinErrorMEtriggerbin4 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin4);	}
    else{				BinContentMETriggerbin4=0.0;							BinErrorMEtriggerbin4=0.0;					}
    if(!BinErrorMEtriggerbin5>1.0e-10){	BinContentMETriggerbin5 = BinContentMETriggerbin5/BinErrorMEtriggerbin5;	BinErrorMEtriggerbin5 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin5);	}
    else{				BinContentMETriggerbin5=0.0;							BinErrorMEtriggerbin5=0.0;					}
    //Set the bin content and error in every histogram.
    if(BinContent>1.0e-10){
      hist->SetBinContent(x,BinContent);
      hist->SetBinError(x,BinError);
      histbin1->SetBinContent(x,BinContentbin1);
      histbin1->SetBinError(x,BinErrorbin1);
      histbin2->SetBinContent(x,BinContentbin2);
      histbin2->SetBinError(x,BinErrorbin2);
      histbin3->SetBinContent(x,BinContentbin3);
      histbin3->SetBinError(x,BinErrorbin3);
      histbin4->SetBinContent(x,BinContentbin4);
      histbin4->SetBinError(x,BinErrorbin4);
      histbin5->SetBinContent(x,BinContentbin5);
      histbin5->SetBinError(x,BinErrorbin5);
      histMETA->SetBinContent(x,BinContentMETA);
      histMETA->SetBinError(x,BinErrorMETA);
      histMETAbin1->SetBinContent(x,BinContentMETAbin1);
      histMETAbin1->SetBinError(x,BinErrorMETAbin1);
      histMETAbin2->SetBinContent(x,BinContentMETAbin2);
      histMETAbin2->SetBinError(x,BinErrorMETAbin2);
      histMETAbin3->SetBinContent(x,BinContentMETAbin3);
      histMETAbin3->SetBinError(x,BinErrorMETAbin3);
      histMETAbin4->SetBinContent(x,BinContentMETAbin4);
      histMETAbin4->SetBinError(x,BinErrorMETAbin4);
      histMETAbin5->SetBinContent(x,BinContentMETAbin5);
      histMETAbin5->SetBinError(x,BinErrorMETAbin5);
      histMETrigger->SetBinContent(x,BinContentMETrigger);
      histMETrigger->SetBinError(x,BinErrorMEtrigger);
      histMETriggerbin1->SetBinContent(x,BinContentMETriggerbin1);
      histMETriggerbin1->SetBinError(x,BinErrorMEtriggerbin1);
      histMETriggerbin2->SetBinContent(x,BinContentMETriggerbin2);
      histMETriggerbin2->SetBinError(x,BinErrorMEtriggerbin2);
      histMETriggerbin3->SetBinContent(x,BinContentMETriggerbin3);
      histMETriggerbin3->SetBinError(x,BinErrorMEtriggerbin3);
      histMETriggerbin4->SetBinContent(x,BinContentMETriggerbin4);
      histMETriggerbin4->SetBinError(x,BinErrorMEtriggerbin4);
      histMETriggerbin5->SetBinContent(x,BinContentMETriggerbin5);
      histMETriggerbin5->SetBinError(x,BinErrorMEtriggerbin5);
    }
    //Reset all:
    BinContent = 0.0;BinError = 0.0;BinContentMETA = 0.0;BinErrorMETA = 0.0;BinContentMETrigger = 0.0;BinErrorMEtrigger = 0.0;
    BinContentbin1 = 0.0;BinErrorbin1 = 0.0;BinContentMETAbin1 = 0.0;BinErrorMETAbin1 = 0.0;BinContentMETriggerbin1 = 0.0;BinErrorMEtriggerbin1 = 0.0;
    BinContentbin2 = 0.0;BinErrorbin2 = 0.0;BinContentMETAbin2 = 0.0;BinErrorMETAbin2 = 0.0;BinContentMETriggerbin2 = 0.0;BinErrorMEtriggerbin2 = 0.0;
    BinContentbin3 = 0.0;BinErrorbin3 = 0.0;BinContentMETAbin3 = 0.0;BinErrorMETAbin3 = 0.0;BinContentMETriggerbin3 = 0.0;BinErrorMEtriggerbin3 = 0.0;
    BinContentbin4 = 0.0;BinErrorbin4 = 0.0;BinContentMETAbin4 = 0.0;BinErrorMETAbin4 = 0.0;BinContentMETriggerbin4 = 0.0;BinErrorMEtriggerbin4 = 0.0;
    BinContentbin5 = 0.0;BinErrorbin5 = 0.0;BinContentMETAbin5 = 0.0;BinErrorMETAbin5 = 0.0;BinContentMETriggerbin5 = 0.0;BinErrorMEtriggerbin5 = 0.0;	
  }//end loop over bins.
  //save the histograms in the relevant directories:
  All->Samediv()->cd();
  hist->Write(histo->GetName());
  delete hist;
  All->METAdiv()->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  All->METriggerdiv()->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
  Bin1->Samediv()->cd();
  histbin1->Write(histo->GetName());
  delete histbin1;
  Bin1->METAdiv()->cd();
  histMETAbin1->Write(histo->GetName());
  delete histMETAbin1;
  Bin1->METriggerdiv()->cd();
  histMETriggerbin1->Write(histo->GetName());
  delete histMETriggerbin1;  
  Bin2->Samediv()->cd();
  histbin2->Write(histo->GetName());
  delete histbin2;
  Bin2->METAdiv()->cd();
  histMETAbin2->Write(histo->GetName());
  delete histMETAbin2;
  Bin2->METriggerdiv()->cd();
  histMETriggerbin2->Write(histo->GetName());
  delete histMETriggerbin2;  
  Bin3->Samediv()->cd();
  histbin3->Write(histo->GetName());
  delete histbin3;
  Bin3->METAdiv()->cd();
  histMETAbin3->Write(histo->GetName());
  delete histMETAbin3;
  Bin3->METriggerdiv()->cd();
  histMETriggerbin3->Write(histo->GetName());
  delete histMETriggerbin3;  
  Bin4->Samediv()->cd();
  histbin4->Write(histo->GetName());
  delete histbin4;
  Bin4->METAdiv()->cd();
  histMETAbin4->Write(histo->GetName());
  delete histMETAbin4;
  Bin4->METriggerdiv()->cd();
  histMETriggerbin4->Write(histo->GetName());
  delete histMETriggerbin4;  
  Bin5->Samediv()->cd();
  histbin5->Write(histo->GetName());
  delete histbin5;
  Bin5->METAdiv()->cd();
  histMETAbin5->Write(histo->GetName());
  delete histMETAbin5;
  Bin5->METriggerdiv()->cd();
  histMETriggerbin5->Write(histo->GetName());
  delete histMETriggerbin5;  
}
void CollectHist(TH2D* histo, TList * directories, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs * Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs * Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs * Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));  
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH2D* hist 			= dynamic_cast<TH2D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH2D* histMETA 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH2D* histMETrigger 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  //Bin 1:
  TH2D* histbin1 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH2D* histMETAbin1 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH2D* histMETriggerbin1 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  //Bin 2:
  TH2D* histbin2		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
  TH2D* histMETAbin2		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
  TH2D* histMETriggerbin2 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));
  //Bin 3:
  TH2D* histbin3		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
  TH2D* histMETAbin3		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
  TH2D* histMETriggerbin3 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  //Bin 4:
  TH2D* histbin4		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
  TH2D* histMETAbin4		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
  TH2D* histMETriggerbin4 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  //Bin 5:
  TH2D* histbin5		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
  TH2D* histMETAbin5		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
  TH2D* histMETriggerbin5	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  //Doubles to hold the values until they are put into the hists.
  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;
  Double_t BinContent = 0.0;Double_t BinError   = 0.0; Double_t BinContentMETA = 0.0;Double_t BinErrorMETA   = 0.0;
  Double_t BinContentMETrigger = 0.0;Double_t BinErrorMEtrigger   = 0.0; Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  //loop over the bins:
  for(int x=0;x<=histo->GetNbinsX()+1;x++){for(int y=0;y<=histo->GetNbinsY();y++){
    //loop over all Multiplicity-Vertex bins:
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
      //find the Multiplicity bin we are in
      for(int j = 1;j<=5;j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName())))
      {
// 	cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
	continue;}
      //extract bin content and error in the relevant bin:
      bincontl 			= dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
      binerrorl 		= dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y);	
      bincontlMETA 		= dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
      binerrorlMETA 		= dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y);	      
      bincontlMEtrigger 	= dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
      binerrorlMETrigger 	= dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y); 
      
      if(bincontl>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContent += bincontl/(binerrorl*binerrorl);
	BinError   += 1.0/(binerrorl*binerrorl);
	if(Mbin == 1){		BinContentbin1 += bincontl/(binerrorl*binerrorl); BinErrorbin1   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 2){	BinContentbin2 += bincontl/(binerrorl*binerrorl); BinErrorbin2   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 3){	BinContentbin3 += bincontl/(binerrorl*binerrorl); BinErrorbin3   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 4){	BinContentbin4 += bincontl/(binerrorl*binerrorl); BinErrorbin4   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 5){ 	BinContentbin5 += bincontl/(binerrorl*binerrorl); BinErrorbin5   += 1.0/(binerrorl*binerrorl);}
      }
      else{
	BinError += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorbin1   += 1.0;
	else if(Mbin == 2)	BinErrorbin2   += 1.0;
	else if(Mbin == 3)	BinErrorbin3   += 1.0;
	else if(Mbin == 4)	BinErrorbin4   += 1.0;	      
	else if(Mbin == 5)	BinErrorbin5   += 1.0;
      }
            if(bincontlMETA>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETA += bincontlMETA/(binerrorlMETA*binerrorlMETA);
	BinErrorMETA   += 1.0/(binerrorlMETA*binerrorlMETA);
	if(Mbin == 1){		BinContentMETAbin1 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin1   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 2){	BinContentMETAbin2 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin2   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 3){	BinContentMETAbin3 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin3   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 4){	BinContentMETAbin4 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin4   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 5){	BinContentMETAbin5 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin5   += 1.0/(binerrorlMETA*binerrorlMETA);}
      }
      else{
	BinErrorMETA += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMETAbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMETAbin2   += 1.0;	      
	else if(Mbin == 3)	BinErrorMETAbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMETAbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMETAbin5   += 1.0;
      }
      if(bincontlMEtrigger>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETrigger += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger);
	BinErrorMEtrigger   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);
	if(Mbin == 1){		BinContentMETriggerbin1 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin1   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 2){	BinContentMETriggerbin2 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin2   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 3){	BinContentMETriggerbin3 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin3   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 4){	BinContentMETriggerbin4 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin4   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 5){	BinContentMETriggerbin5 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin5   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
      }
      else{
	BinErrorMEtrigger += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMEtriggerbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMEtriggerbin2   += 1.0;
	else if(Mbin == 3)	BinErrorMEtriggerbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMEtriggerbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMEtriggerbin5   += 1.0;
      }
      //reset the local bins to zero just to be shure:
      bincontl  = 0.0;		binerrorl = 0.0;
      bincontlMETA =0.0;	binerrorlMETA=0.0;
      bincontlMEtrigger=0.0;	binerrorlMETrigger=0.0;
    }//end loop over M-V bins.
    //normalize the bin and the error for same:
    if(BinError>1.0e-10){		BinContent = BinContent/BinError;						BinError = 1.0/TMath::Sqrt(BinError);				}
    else{				BinContent=0.0;									BinError=0.0;							}
    if(BinErrorbin1>1.0e-10){		BinContentbin1 = BinContentbin1/BinErrorbin1;					BinErrorbin1 = 1.0/TMath::Sqrt(BinErrorbin1);			}
    else{				BinContentbin1=0.0;								BinErrorbin1=0.0;						}
    if(BinErrorbin2>1.0e-10){		BinContentbin2 = BinContentbin2/BinErrorbin2;					BinErrorbin2 = 1.0/TMath::Sqrt(BinErrorbin2);			}
    else{				BinContentbin2=0.0;								BinErrorbin2=0.0;						}
    if(BinErrorbin3>1.0e-10){		BinContentbin3 = BinContentbin3/BinErrorbin3;					BinErrorbin3 = 1.0/TMath::Sqrt(BinErrorbin3);			}
    else{				BinContentbin3=0.0;								BinErrorbin3=0.0;						}
    if(BinErrorbin4>1.0e-10){		BinContentbin4 = BinContentbin4/BinErrorbin4;					BinErrorbin4 = 1.0/TMath::Sqrt(BinErrorbin4);			}
    else{				BinContentbin4=0.0;								BinErrorbin4=0.0;						}
    if(BinErrorbin5>1.0e-10){		BinContentbin5 = BinContentbin5/BinErrorbin5;					BinErrorbin5 = 1.0/TMath::Sqrt(BinErrorbin5);			}
    else{				BinContentbin5=0.0;								BinErrorbin5=0.0;						}
    //normalize the bin and the error for META:    
    if(BinErrorMETA>1.0e-10){		BinContentMETA = BinContentMETA/BinErrorMETA;					BinErrorMETA = 1.0/TMath::Sqrt(BinErrorMETA);			}
    else{				BinContentMETA=0.0;								BinErrorMETA=0.0;						}
    if(BinErrorMETAbin1>1.0e-10){	BinContentMETAbin1 = BinContentMETAbin1/BinErrorMETAbin1;			BinErrorMETAbin1 = 1.0/TMath::Sqrt(BinErrorMETAbin1);		}
    else{				BinContentMETAbin1=0.0;								BinErrorMETAbin1=0.0;						}
    if(BinErrorMETAbin2>1.0e-10){	BinContentMETAbin2 = BinContentMETAbin2/BinErrorMETAbin2;			BinErrorMETAbin2 = 1.0/TMath::Sqrt(BinErrorMETAbin2);		}
    else{				BinContentMETAbin2=0.0;								BinErrorMETAbin2=0.0;						}
    if(BinErrorMETAbin3>1.0e-10){	BinContentMETAbin3 = BinContentMETAbin3/BinErrorMETAbin3;			BinErrorMETAbin3 = 1.0/TMath::Sqrt(BinErrorMETAbin3);		}
    else{				BinContentMETAbin3=0.0;								BinErrorMETAbin3=0.0;						}
    if(BinErrorMETAbin4>1.0e-10){	BinContentMETAbin4 = BinContentMETAbin4/BinErrorMETAbin4;			BinErrorMETAbin4 = 1.0/TMath::Sqrt(BinErrorMETAbin4);		}
    else{				BinContentMETAbin4=0.0;								BinErrorMETAbin4=0.0;						}
    if(BinErrorMETAbin5>1.0e-10){	BinContentMETAbin5 = BinContentMETAbin5/BinErrorMETAbin5;			BinErrorMETAbin5 = 1.0/TMath::Sqrt(BinErrorMETAbin5);		}
    else{				BinContentMETAbin5=0.0;								BinErrorMETAbin5=0.0;						}
    //normalize the bin and the error for METrigger:    
    if(BinErrorMEtrigger>1.0e-10){	BinContentMETrigger = BinContentMETrigger/BinErrorMEtrigger;			BinErrorMEtrigger = 1.0/TMath::Sqrt(BinErrorMEtrigger);		}
    else{				BinContentMETrigger=0.0;							BinErrorMEtrigger=0.0;						}
    if(BinErrorMEtriggerbin1>1.0e-10){	BinContentMETriggerbin1 = BinContentMETriggerbin1/BinErrorMEtriggerbin1;	BinErrorMEtriggerbin1 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin1);	}
    else{				BinContentMETriggerbin1=0.0;							BinErrorMEtriggerbin1=0.0;					}
    if(BinErrorMEtriggerbin2>1.0e-10){	BinContentMETriggerbin2 = BinContentMETriggerbin2/BinErrorMEtriggerbin2;	BinErrorMEtriggerbin2 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin2);	}
    else{				BinContentMETriggerbin2=0.0;							BinErrorMEtriggerbin2=0.0;					}
    if(BinErrorMEtriggerbin3>1.0e-10){	BinContentMETriggerbin3 = BinContentMETriggerbin3/BinErrorMEtriggerbin3;	BinErrorMEtriggerbin3 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin3);	}
    else{				BinContentMETriggerbin3=0.0;							BinErrorMEtriggerbin3=0.0;					}
    if(BinErrorMEtriggerbin4>1.0e-10){	BinContentMETriggerbin4 = BinContentMETriggerbin4/BinErrorMEtriggerbin4;	BinErrorMEtriggerbin4 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin4);	}
    else{				BinContentMETriggerbin4=0.0;							BinErrorMEtriggerbin4=0.0;					}
    if(!BinErrorMEtriggerbin5>1.0e-10){	BinContentMETriggerbin5 = BinContentMETriggerbin5/BinErrorMEtriggerbin5;	BinErrorMEtriggerbin5 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin5);	}
    else{				BinContentMETriggerbin5=0.0;							BinErrorMEtriggerbin5=0.0;					}
    //Set the bin content and error in every histogram.
    if(BinContent>1.0e-10){
      hist->SetBinContent(x,y,BinContent);
      hist->SetBinError(x,y,BinError);
      histbin1->SetBinContent(x,y,BinContentbin1);
      histbin1->SetBinError(x,y,BinErrorbin1);
      histbin2->SetBinContent(x,y,BinContentbin2);
      histbin2->SetBinError(x,y,BinErrorbin2);
      histbin3->SetBinContent(x,y,BinContentbin3);
      histbin3->SetBinError(x,y,BinErrorbin3);
      histbin4->SetBinContent(x,y,BinContentbin4);
      histbin4->SetBinError(x,y,BinErrorbin4);
      histbin5->SetBinContent(x,y,BinContentbin5);
      histbin5->SetBinError(x,y,BinErrorbin5);
      histMETA->SetBinContent(x,y,BinContentMETA);
      histMETA->SetBinError(x,y,BinErrorMETA);
      histMETAbin1->SetBinContent(x,y,BinContentMETAbin1);
      histMETAbin1->SetBinError(x,y,BinErrorMETAbin1);
      histMETAbin2->SetBinContent(x,y,BinContentMETAbin2);
      histMETAbin2->SetBinError(x,y,BinErrorMETAbin2);
      histMETAbin3->SetBinContent(x,y,BinContentMETAbin3);
      histMETAbin3->SetBinError(x,y,BinErrorMETAbin3);
      histMETAbin4->SetBinContent(x,y,BinContentMETAbin4);
      histMETAbin4->SetBinError(x,y,BinErrorMETAbin4);
      histMETAbin5->SetBinContent(x,y,BinContentMETAbin5);
      histMETAbin5->SetBinError(x,y,BinErrorMETAbin5);
      histMETrigger->SetBinContent(x,y,BinContentMETrigger);
      histMETrigger->SetBinError(x,y,BinErrorMEtrigger);
      histMETriggerbin1->SetBinContent(x,y,BinContentMETriggerbin1);
      histMETriggerbin1->SetBinError(x,y,BinErrorMEtriggerbin1);
      histMETriggerbin2->SetBinContent(x,y,BinContentMETriggerbin2);
      histMETriggerbin2->SetBinError(x,y,BinErrorMEtriggerbin2);
      histMETriggerbin3->SetBinContent(x,y,BinContentMETriggerbin3);
      histMETriggerbin3->SetBinError(x,y,BinErrorMEtriggerbin3);
      histMETriggerbin4->SetBinContent(x,y,BinContentMETriggerbin4);
      histMETriggerbin4->SetBinError(x,y,BinErrorMEtriggerbin4);
      histMETriggerbin5->SetBinContent(x,y,BinContentMETriggerbin5);
      histMETriggerbin5->SetBinError(x,y,BinErrorMEtriggerbin5);
    }
    //Reset all:
    BinContent = 0.0;BinError = 0.0;BinContentMETA = 0.0;BinErrorMETA = 0.0;BinContentMETrigger = 0.0;BinErrorMEtrigger = 0.0;
    BinContentbin1 = 0.0;BinErrorbin1 = 0.0;BinContentMETAbin1 = 0.0;BinErrorMETAbin1 = 0.0;BinContentMETriggerbin1 = 0.0;BinErrorMEtriggerbin1 = 0.0;
    BinContentbin2 = 0.0;BinErrorbin2 = 0.0;BinContentMETAbin2 = 0.0;BinErrorMETAbin2 = 0.0;BinContentMETriggerbin2 = 0.0;BinErrorMEtriggerbin2 = 0.0;
    BinContentbin3 = 0.0;BinErrorbin3 = 0.0;BinContentMETAbin3 = 0.0;BinErrorMETAbin3 = 0.0;BinContentMETriggerbin3 = 0.0;BinErrorMEtriggerbin3 = 0.0;
    BinContentbin4 = 0.0;BinErrorbin4 = 0.0;BinContentMETAbin4 = 0.0;BinErrorMETAbin4 = 0.0;BinContentMETriggerbin4 = 0.0;BinErrorMEtriggerbin4 = 0.0;
    BinContentbin5 = 0.0;BinErrorbin5 = 0.0;BinContentMETAbin5 = 0.0;BinErrorMETAbin5 = 0.0;BinContentMETriggerbin5 = 0.0;BinErrorMEtriggerbin5 = 0.0;	
  }}//end binloop
  //save the histograms in the relevant directories:
  All->Samediv()->cd();
  hist->Write(histo->GetName());
  delete hist;
  All->METAdiv()->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  All->METriggerdiv()->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
  Bin1->Samediv()->cd();
  histbin1->Write(histo->GetName());
  delete histbin1;
  Bin1->METAdiv()->cd();
  histMETAbin1->Write(histo->GetName());
  delete histMETAbin1;
  Bin1->METriggerdiv()->cd();
  histMETriggerbin1->Write(histo->GetName());
  delete histMETriggerbin1;  
  Bin2->Samediv()->cd();
  histbin2->Write(histo->GetName());
  delete histbin2;
  Bin2->METAdiv()->cd();
  histMETAbin2->Write(histo->GetName());
  delete histMETAbin2;
  Bin2->METriggerdiv()->cd();
  histMETriggerbin2->Write(histo->GetName());
  delete histMETriggerbin2;  
  Bin3->Samediv()->cd();
  histbin3->Write(histo->GetName());
  delete histbin3;
  Bin3->METAdiv()->cd();
  histMETAbin3->Write(histo->GetName());
  delete histMETAbin3;
  Bin3->METriggerdiv()->cd();
  histMETriggerbin3->Write(histo->GetName());
  delete histMETriggerbin3;  
  Bin4->Samediv()->cd();
  histbin4->Write(histo->GetName());
  delete histbin4;
  Bin4->METAdiv()->cd();
  histMETAbin4->Write(histo->GetName());
  delete histMETAbin4;
  Bin4->METriggerdiv()->cd();
  histMETriggerbin4->Write(histo->GetName());
  delete histMETriggerbin4;  
  Bin5->Samediv()->cd();
  histbin5->Write(histo->GetName());
  delete histbin5;
  Bin5->METAdiv()->cd();
  histMETAbin5->Write(histo->GetName());
  delete histMETAbin5;
  Bin5->METriggerdiv()->cd();
  histMETriggerbin5->Write(histo->GetName());
  delete histMETriggerbin5;  
}
void CollectHist(TH3D* histo, TList * directories, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs * Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs * Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs * Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH3D* hist 			= dynamic_cast<TH3D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH3D* histMETA 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH3D* histMETrigger 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  //Bin 1:
  TH3D* histbin1 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH3D* histMETAbin1 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH3D* histMETriggerbin1 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  //Bin 2:
  TH3D* histbin2		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
  TH3D* histMETAbin2		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
  TH3D* histMETriggerbin2 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));
  //Bin 3:
  TH3D* histbin3		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
  TH3D* histMETAbin3		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
  TH3D* histMETriggerbin3 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  //Bin 4:
  TH3D* histbin4		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
  TH3D* histMETAbin4		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
  TH3D* histMETriggerbin4 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  //Bin 5:
  TH3D* histbin5		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
  TH3D* histMETAbin5		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
  TH3D* histMETriggerbin5	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  //Doubles to hold the values until they are put into the hists.
  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;
  Double_t BinContent = 0.0;Double_t BinError   = 0.0; Double_t BinContentMETA = 0.0;Double_t BinErrorMETA   = 0.0;
  Double_t BinContentMETrigger = 0.0;Double_t BinErrorMEtrigger   = 0.0; Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  //loop over the bins:'
  for(int x=0;x<=histo->GetNbinsX()+1;x++){for(int y=0;y<=histo->GetNbinsY();y++){for(int z=0;z<=histo->GetNbinsZ();z++){
    //loop over all Multiplicity-Vertex bins:
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
      //find the Multiplicity bin we are in
      for(int j = 1;j<=5;j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH3D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName())))
      {
// 	cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
	continue;}
      //extract bin content and error in the relevant bin:
      bincontl 			= dynamic_cast<TH3D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y,z);
      binerrorl 		= dynamic_cast<TH3D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y,z);	
      bincontlMETA 		= dynamic_cast<TH3D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y,z);
      binerrorlMETA 		= dynamic_cast<TH3D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y,z);	      
      bincontlMEtrigger 	= dynamic_cast<TH3D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y,z);
      binerrorlMETrigger 	= dynamic_cast<TH3D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y,z); 
      
      if(bincontl>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContent += bincontl/(binerrorl*binerrorl);
	BinError   += 1.0/(binerrorl*binerrorl);
	if(Mbin == 1){		BinContentbin1 += bincontl/(binerrorl*binerrorl); BinErrorbin1   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 2){	BinContentbin2 += bincontl/(binerrorl*binerrorl); BinErrorbin2   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 3){	BinContentbin3 += bincontl/(binerrorl*binerrorl); BinErrorbin3   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 4){	BinContentbin4 += bincontl/(binerrorl*binerrorl); BinErrorbin4   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 5){ 	BinContentbin5 += bincontl/(binerrorl*binerrorl); BinErrorbin5   += 1.0/(binerrorl*binerrorl);}
      }
      else{
	BinError += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorbin1   += 1.0;
	else if(Mbin == 2)	BinErrorbin2   += 1.0;
	else if(Mbin == 3)	BinErrorbin3   += 1.0;
	else if(Mbin == 4)	BinErrorbin4   += 1.0;	      
	else if(Mbin == 5)	BinErrorbin5   += 1.0;
      }
      if(bincontlMETA>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETA += bincontlMETA/(binerrorlMETA*binerrorlMETA);
	BinErrorMETA   += 1.0/(binerrorlMETA*binerrorlMETA);
	if(Mbin == 1){		BinContentMETAbin1 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin1   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 2){	BinContentMETAbin2 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin2   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 3){	BinContentMETAbin3 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin3   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 4){	BinContentMETAbin4 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin4   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 5){	BinContentMETAbin5 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin5   += 1.0/(binerrorlMETA*binerrorlMETA);}
      }
      else{
	BinErrorMETA += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMETAbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMETAbin2   += 1.0;	      
	else if(Mbin == 3)	BinErrorMETAbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMETAbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMETAbin5   += 1.0;
      }
      if(bincontlMEtrigger>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETrigger += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger);
	BinErrorMEtrigger   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);
	if(Mbin == 1){		BinContentMETriggerbin1 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin1   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 2){	BinContentMETriggerbin2 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin2   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 3){	BinContentMETriggerbin3 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin3   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 4){	BinContentMETriggerbin4 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin4   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 5){	BinContentMETriggerbin5 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin5   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
      }
      else{
	BinErrorMEtrigger += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMEtriggerbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMEtriggerbin2   += 1.0;
	else if(Mbin == 3)	BinErrorMEtriggerbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMEtriggerbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMEtriggerbin5   += 1.0;
      }
      //reset the local bins to zero just to be shure:
      bincontl  = 0.0;		binerrorl = 0.0;
      bincontlMETA =0.0;	binerrorlMETA=0.0;
      bincontlMEtrigger=0.0;	binerrorlMETrigger=0.0;
    }//end loop over M-V bins.
    //normalize the bin and the error for same:
    if(BinError>1.0e-10){		BinContent = BinContent/BinError;						BinError = 1.0/TMath::Sqrt(BinError);				}
    else{				BinContent=0.0;									BinError=0.0;							}
    if(BinErrorbin1>1.0e-10){		BinContentbin1 = BinContentbin1/BinErrorbin1;					BinErrorbin1 = 1.0/TMath::Sqrt(BinErrorbin1);			}
    else{				BinContentbin1=0.0;								BinErrorbin1=0.0;						}
    if(BinErrorbin2>1.0e-10){		BinContentbin2 = BinContentbin2/BinErrorbin2;					BinErrorbin2 = 1.0/TMath::Sqrt(BinErrorbin2);			}
    else{				BinContentbin2=0.0;								BinErrorbin2=0.0;						}
    if(BinErrorbin3>1.0e-10){		BinContentbin3 = BinContentbin3/BinErrorbin3;					BinErrorbin3 = 1.0/TMath::Sqrt(BinErrorbin3);			}
    else{				BinContentbin3=0.0;								BinErrorbin3=0.0;						}
    if(BinErrorbin4>1.0e-10){		BinContentbin4 = BinContentbin4/BinErrorbin4;					BinErrorbin4 = 1.0/TMath::Sqrt(BinErrorbin4);			}
    else{				BinContentbin4=0.0;								BinErrorbin4=0.0;						}
    if(BinErrorbin5>1.0e-10){		BinContentbin5 = BinContentbin5/BinErrorbin5;					BinErrorbin5 = 1.0/TMath::Sqrt(BinErrorbin5);			}
    else{				BinContentbin5=0.0;								BinErrorbin5=0.0;						}
    //normalize the bin and the error for META:    
    if(BinErrorMETA>1.0e-10){		BinContentMETA = BinContentMETA/BinErrorMETA;					BinErrorMETA = 1.0/TMath::Sqrt(BinErrorMETA);			}
    else{				BinContentMETA=0.0;								BinErrorMETA=0.0;						}
    if(BinErrorMETAbin1>1.0e-10){	BinContentMETAbin1 = BinContentMETAbin1/BinErrorMETAbin1;			BinErrorMETAbin1 = 1.0/TMath::Sqrt(BinErrorMETAbin1);		}
    else{				BinContentMETAbin1=0.0;								BinErrorMETAbin1=0.0;						}
    if(BinErrorMETAbin2>1.0e-10){	BinContentMETAbin2 = BinContentMETAbin2/BinErrorMETAbin2;			BinErrorMETAbin2 = 1.0/TMath::Sqrt(BinErrorMETAbin2);		}
    else{				BinContentMETAbin2=0.0;								BinErrorMETAbin2=0.0;						}
    if(BinErrorMETAbin3>1.0e-10){	BinContentMETAbin3 = BinContentMETAbin3/BinErrorMETAbin3;			BinErrorMETAbin3 = 1.0/TMath::Sqrt(BinErrorMETAbin3);		}
    else{				BinContentMETAbin3=0.0;								BinErrorMETAbin3=0.0;						}
    if(BinErrorMETAbin4>1.0e-10){	BinContentMETAbin4 = BinContentMETAbin4/BinErrorMETAbin4;			BinErrorMETAbin4 = 1.0/TMath::Sqrt(BinErrorMETAbin4);		}
    else{				BinContentMETAbin4=0.0;								BinErrorMETAbin4=0.0;						}
    if(BinErrorMETAbin5>1.0e-10){	BinContentMETAbin5 = BinContentMETAbin5/BinErrorMETAbin5;			BinErrorMETAbin5 = 1.0/TMath::Sqrt(BinErrorMETAbin5);		}
    else{				BinContentMETAbin5=0.0;								BinErrorMETAbin5=0.0;						}
    //normalize the bin and the error for METrigger:    
    if(BinErrorMEtrigger>1.0e-10){	BinContentMETrigger = BinContentMETrigger/BinErrorMEtrigger;			BinErrorMEtrigger = 1.0/TMath::Sqrt(BinErrorMEtrigger);		}
    else{				BinContentMETrigger=0.0;							BinErrorMEtrigger=0.0;						}
    if(BinErrorMEtriggerbin1>1.0e-10){	BinContentMETriggerbin1 = BinContentMETriggerbin1/BinErrorMEtriggerbin1;	BinErrorMEtriggerbin1 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin1);	}
    else{				BinContentMETriggerbin1=0.0;							BinErrorMEtriggerbin1=0.0;					}
    if(BinErrorMEtriggerbin2>1.0e-10){	BinContentMETriggerbin2 = BinContentMETriggerbin2/BinErrorMEtriggerbin2;	BinErrorMEtriggerbin2 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin2);	}
    else{				BinContentMETriggerbin2=0.0;							BinErrorMEtriggerbin2=0.0;					}
    if(BinErrorMEtriggerbin3>1.0e-10){	BinContentMETriggerbin3 = BinContentMETriggerbin3/BinErrorMEtriggerbin3;	BinErrorMEtriggerbin3 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin3);	}
    else{				BinContentMETriggerbin3=0.0;							BinErrorMEtriggerbin3=0.0;					}
    if(BinErrorMEtriggerbin4>1.0e-10){	BinContentMETriggerbin4 = BinContentMETriggerbin4/BinErrorMEtriggerbin4;	BinErrorMEtriggerbin4 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin4);	}
    else{				BinContentMETriggerbin4=0.0;							BinErrorMEtriggerbin4=0.0;					}
    if(!BinErrorMEtriggerbin5>1.0e-10){	BinContentMETriggerbin5 = BinContentMETriggerbin5/BinErrorMEtriggerbin5;	BinErrorMEtriggerbin5 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin5);	}
    else{				BinContentMETriggerbin5=0.0;							BinErrorMEtriggerbin5=0.0;					}
    //Set the bin content and error in every histogram.
    if(BinContent>1.0e-10){
      hist->SetBinContent(x,y,z,BinContent);
      hist->SetBinError(x,y,z,BinError);
      histbin1->SetBinContent(x,y,z,BinContentbin1);
      histbin1->SetBinError(x,y,z,BinErrorbin1);
      histbin2->SetBinContent(x,y,z,BinContentbin2);
      histbin2->SetBinError(x,y,z,BinErrorbin2);
      histbin3->SetBinContent(x,y,z,BinContentbin3);
      histbin3->SetBinError(x,y,z,BinErrorbin3);
      histbin4->SetBinContent(x,y,z,BinContentbin4);
      histbin4->SetBinError(x,y,z,BinErrorbin4);
      histbin5->SetBinContent(x,y,z,BinContentbin5);
      histbin5->SetBinError(x,y,z,BinErrorbin5);
      histMETA->SetBinContent(x,y,z,BinContentMETA);
      histMETA->SetBinError(x,y,z,BinErrorMETA);
      histMETAbin1->SetBinContent(x,y,z,BinContentMETAbin1);
      histMETAbin1->SetBinError(x,y,z,BinErrorMETAbin1);
      histMETAbin2->SetBinContent(x,y,z,BinContentMETAbin2);
      histMETAbin2->SetBinError(x,y,z,BinErrorMETAbin2);
      histMETAbin3->SetBinContent(x,y,z,BinContentMETAbin3);
      histMETAbin3->SetBinError(x,y,z,BinErrorMETAbin3);
      histMETAbin4->SetBinContent(x,y,z,BinContentMETAbin4);
      histMETAbin4->SetBinError(x,y,z,BinErrorMETAbin4);
      histMETAbin5->SetBinContent(x,y,z,BinContentMETAbin5);
      histMETAbin5->SetBinError(x,y,z,BinErrorMETAbin5);
      histMETrigger->SetBinContent(x,y,z,BinContentMETrigger);
      histMETrigger->SetBinError(x,y,z,BinErrorMEtrigger);
      histMETriggerbin1->SetBinContent(x,y,z,BinContentMETriggerbin1);
      histMETriggerbin1->SetBinError(x,y,z,BinErrorMEtriggerbin1);
      histMETriggerbin2->SetBinContent(x,y,z,BinContentMETriggerbin2);
      histMETriggerbin2->SetBinError(x,y,z,BinErrorMEtriggerbin2);
      histMETriggerbin3->SetBinContent(x,y,z,BinContentMETriggerbin3);
      histMETriggerbin3->SetBinError(x,y,z,BinErrorMEtriggerbin3);
      histMETriggerbin4->SetBinContent(x,y,z,BinContentMETriggerbin4);
      histMETriggerbin4->SetBinError(x,y,z,BinErrorMEtriggerbin4);
      histMETriggerbin5->SetBinContent(x,y,z,BinContentMETriggerbin5);
      histMETriggerbin5->SetBinError(x,y,z,BinErrorMEtriggerbin5);
    }
    //Reset all:
    BinContent = 0.0;BinError = 0.0;BinContentMETA = 0.0;BinErrorMETA = 0.0;BinContentMETrigger = 0.0;BinErrorMEtrigger = 0.0;
    BinContentbin1 = 0.0;BinErrorbin1 = 0.0;BinContentMETAbin1 = 0.0;BinErrorMETAbin1 = 0.0;BinContentMETriggerbin1 = 0.0;BinErrorMEtriggerbin1 = 0.0;
    BinContentbin2 = 0.0;BinErrorbin2 = 0.0;BinContentMETAbin2 = 0.0;BinErrorMETAbin2 = 0.0;BinContentMETriggerbin2 = 0.0;BinErrorMEtriggerbin2 = 0.0;
    BinContentbin3 = 0.0;BinErrorbin3 = 0.0;BinContentMETAbin3 = 0.0;BinErrorMETAbin3 = 0.0;BinContentMETriggerbin3 = 0.0;BinErrorMEtriggerbin3 = 0.0;
    BinContentbin4 = 0.0;BinErrorbin4 = 0.0;BinContentMETAbin4 = 0.0;BinErrorMETAbin4 = 0.0;BinContentMETriggerbin4 = 0.0;BinErrorMEtriggerbin4 = 0.0;
    BinContentbin5 = 0.0;BinErrorbin5 = 0.0;BinContentMETAbin5 = 0.0;BinErrorMETAbin5 = 0.0;BinContentMETriggerbin5 = 0.0;BinErrorMEtriggerbin5 = 0.0;	
  }}}//end binloop
  //save the histograms in the relevant directories:
  //save the histograms in the relevant directories:
  All->Samediv()->cd();
  hist->Write(histo->GetName());
  delete hist;
  All->METAdiv()->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  All->METriggerdiv()->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
  Bin1->Samediv()->cd();
  histbin1->Write(histo->GetName());
  delete histbin1;
  Bin1->METAdiv()->cd();
  histMETAbin1->Write(histo->GetName());
  delete histMETAbin1;
  Bin1->METriggerdiv()->cd();
  histMETriggerbin1->Write(histo->GetName());
  delete histMETriggerbin1;  
  Bin2->Samediv()->cd();
  histbin2->Write(histo->GetName());
  delete histbin2;
  Bin2->METAdiv()->cd();
  histMETAbin2->Write(histo->GetName());
  delete histMETAbin2;
  Bin2->METriggerdiv()->cd();
  histMETriggerbin2->Write(histo->GetName());
  delete histMETriggerbin2;  
  Bin3->Samediv()->cd();
  histbin3->Write(histo->GetName());
  delete histbin3;
  Bin3->METAdiv()->cd();
  histMETAbin3->Write(histo->GetName());
  delete histMETAbin3;
  Bin3->METriggerdiv()->cd();
  histMETriggerbin3->Write(histo->GetName());
  delete histMETriggerbin3;  
  Bin4->Samediv()->cd();
  histbin4->Write(histo->GetName());
  delete histbin4;
  Bin4->METAdiv()->cd();
  histMETAbin4->Write(histo->GetName());
  delete histMETAbin4;
  Bin4->METriggerdiv()->cd();
  histMETriggerbin4->Write(histo->GetName());
  delete histMETriggerbin4;  
  Bin5->Samediv()->cd();
  histbin5->Write(histo->GetName());
  delete histbin5;
  Bin5->METAdiv()->cd();
  histMETAbin5->Write(histo->GetName());
  delete histMETAbin5;
  Bin5->METriggerdiv()->cd();
  histMETriggerbin5->Write(histo->GetName());
  delete histMETriggerbin5;  
}
void CollectHist(const char* histname,TList * directories, TObjArray* multdirlist){
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(multdirlist->At(0))->Same()->GetDirectory(Form("%s/divided",directories->At(1)->GetName()))->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  TH2D* hist2d = dynamic_cast<TH2D*>(hist);
  TH3D* hist3d = dynamic_cast<TH3D*>(hist);
  if(hist1d){cout << "TH1D " <<histname<<endl;CollectHist(hist1d,directories,multdirlist);}
  if(hist2d){cout << "TH2D " <<histname<<endl;CollectHist(hist2d,directories,multdirlist);}
  if(hist3d){cout << "TH3D " <<histname<<endl;CollectHist(hist3d,directories,multdirlist);}
  canvasmaker(histname,multdirlist);
}
void CollectHistbinstats(const char* histname,TList * directories, TObjArray* multdirlist){
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(multdirlist->At(0))->Same()->GetDirectory(Form("%s/bin_stats",directories->At(1)->GetName()))->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  if(hist1d){cout << "TH1D " <<histname<<endl;CollectHistbs(hist1d,directories,multdirlist);}
  canvasmaker(histname,multdirlist);
}

///////////////////////////////////////


//functions to correct with ME types.
  
//function for minimizing the histogram edge:
Double_t Chi2(Double_t scaleMETA, Double_t scaleMETrigger){
  Double_t grad = 0.0;
  TH2D* testhist = dynamic_cast<TH2D*>(gSameEvent->Clone("testhist"));
  testhist->Add(gMETA,-1.0*scaleMETA);
  testhist->Add(gMETrigger,-1.0*scaleMETrigger);
  for(int x=1;x<testhist->GetNbinsX();x++){
    Double_t gradloc = testhist->GetBinContent(x,1) - testhist->GetBinContent(x+1,1);
    grad += gradloc*gradloc;
    gradloc = testhist->GetBinContent(x,testhist->GetNbinsY()) - testhist->GetBinContent(x+1,testhist->GetNbinsY());
    grad += gradloc*gradloc;
  }
  return grad;
}
//____________________________________________________________________________________
void FcnFitScale(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t) {
  //
  // Chi2 function
  // par[0] - scale factor for the META 	background histogram
  // par[0] - scale factor for the METrigger 	background histogram
  //
  Double_t chi2 = Chi2(par[0],par[1]);
  f = chi2;
} 

  
void GetScalingFactors(BinDirs * BinDir,const char* in, const char* out,Double_t * METAscale, Double_t * METriggerscale,bool ispp, bool isPbPb, int iter)
{
  //Function to find the relevant scaling factors.
  //Clean and create the output directory
  resultsdirectory(BinDir->Same(),out);
  if(!BinDir->SameDir(in)) return;
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * SameHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METAHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_META")));
  tmp = BinDir->METriggerdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METRiggerHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_METrigger")));
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DEta_12");if(!tmp){return;}
  TH2D * SameHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DEta_12");if(!tmp){return;}  
  TH2D * METAHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_META")));
   tmp = BinDir->METriggerdiv()->Get("DPhi_1_DEta_12");if(!tmp){return;}
  TH2D * METRiggerHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_METrigger")));  
  delete tmp;

  gSameEvent 	= SameHistPhiPhi;
  gMETA 	= METAHistPhiPhi;
  gMETrigger 	= METRiggerHistPhiPhi;
  
  
  TMinuit* minuit=new TMinuit(1);
  minuit->SetFCN(FcnFitScale);

  Double_t arglist[2];
  Int_t ierflg=0;
  arglist[0] = 1;
  minuit->mnexcm("SET ERR", arglist, 1, ierflg); 
  // Set starting values and step sizes for parameters
  Double_t vstart[2];
  Double_t step[2];
  vstart[0] = *METAscale;
  step[0] = 0.01; 
  vstart[1] = *METriggerscale;
  step[1] = 0.01;  
  
  minuit->mnparm(0, "METAscale", vstart[0], step[0], 0.0, 1000000.0, ierflg);
  minuit->mnparm(1, "METriggerscale", vstart[1], step[1], 0.0, 1000000.0, ierflg);
  
  
  arglist[0] = 500;
  arglist[1] = 1.;

  minuit->SetMaxIterations(10000);
  minuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
  minuit->mnexcm("MIGRAD", arglist ,2,ierflg); 
  
  
  if(ispp){
    if(iter ==1){
      Double_t ScalingMETA =SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
      ScalingMETA +=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
      ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      Double_t Divideby = METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
      Divideby +=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
      Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      if(Divideby >0.0) *METAscale = ScalingMETA/Divideby;    
      else *METAscale = 0.0;
      if(*METAscale<0.0)*METAscale = 0.0;
      
      Double_t ScalingMETrigger = SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(0.0),SameHistPhiEta->GetXaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
      ScalingMETrigger -=0.5* SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 + 0.2));
      ScalingMETrigger -=0.5* SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 + 0.2));
      Double_t DividebyMETrigger = METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(0.0),METRiggerHistPhiEta->GetXaxis()->FindBin(0.0),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
      DividebyMETrigger -=0.5* METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(-1.0),METRiggerHistPhiEta->GetXaxis()->FindBin(-1.0),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
      DividebyMETrigger -=0.5* METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(1.0),METRiggerHistPhiEta->GetXaxis()->FindBin(1.0),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
      if(DividebyMETrigger>0.0) *METriggerscale = ScalingMETrigger/DividebyMETrigger;
      else *METriggerscale = 0.0;
      if(*METriggerscale<0.0)*METriggerscale = 0.0;
    }
  }
  if(isPbPb){
    if(iter ==1||iter==2){
      Double_t ScalingMETA =SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
      ScalingMETA +=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
      ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      Double_t Divideby = METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
      Divideby +=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
      Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
      if(Divideby >0.0) *METAscale = ScalingMETA/Divideby;    
      else *METAscale = 0.0;
      if(*METAscale<0.0)*METAscale = 0.0;

      Double_t ScalingMETrigger =SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetNbinsX(),1);
      ScalingMETrigger +=SameHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetNbinsX());
      ScalingMETrigger -=SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-1.0),SameHistPhiPhi->GetXaxis()->FindBin(4.0));
      ScalingMETrigger -=SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(4.0),SameHistPhiPhi->GetXaxis()->FindBin(-1.0));
      Double_t DividebyMETrigger = METRiggerHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetNbinsX()); 
      DividebyMETrigger += METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetNbinsX(),1);
      DividebyMETrigger -=METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-1.0),SameHistPhiPhi->GetXaxis()->FindBin(4.0));
      DividebyMETrigger -=METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(4.0),SameHistPhiPhi->GetXaxis()->FindBin(-1.0));
      if(DividebyMETrigger>0.0) *METriggerscale = ScalingMETrigger/DividebyMETrigger;
      else *METriggerscale = 0.0;
      if(*METriggerscale<0.0)*METriggerscale = 0.0;
    }    
  }
  Double_t temp =0.0;
  Double_t tempe = 0.0;
  minuit->GetParameter(0, temp, tempe); 
  *METAscale = temp;
  minuit->GetParameter(1, temp, tempe); 
  *METriggerscale = temp;
  cout << *METAscale <<" " <<*METriggerscale<<endl;


  delete SameHistPhiPhi;delete METAHistPhiPhi;delete METRiggerHistPhiPhi;delete SameHistPhiEta;delete METAHistPhiEta;delete METRiggerHistPhiEta;

//   
//   if(isPbPb){
//     Double_t ScalingMETA =SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
//     ScalingMETA +=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
//     ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
//     ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
//     Double_t Divideby = METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
//     Divideby +=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
//     Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
//     Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
//     if(Divideby >0.0) *METAscale = ScalingMETA/Divideby;    
//     else *METAscale = 0.0;
//     if(*METAscale<0.0)*METAscale = 0.0;
// 
//     
//     Double_t ScalingMETrigger =SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetNbinsX(),1);// Diagonal(SameHistPhiPhi,SameHistPhiPhi->GetXaxis()->FindBin(1.2),SameHistPhiPhi->GetXaxis()->FindBin(1.4),SameHistPhiPhi->GetXaxis()->FindBin(1.2),SameHistPhiPhi->GetXaxis()->FindBin(1.4))/(1.0*(SameHistPhiPhi->GetXaxis()->FindBin(1.2)-SameHistPhiPhi->GetXaxis()->FindBin(1.4)));//SameHistPhiPhi->Integral(SameHistPhiPhi->GetXaxis()->FindBin(-0.1+TMath::Pi()/2.0),SameHistPhiPhi->GetXaxis()->FindBin(0.1+TMath::Pi()/2.0),SameHistPhiPhi->GetXaxis()->FindBin(-0.1+TMath::Pi()/2.0),SameHistPhiPhi->GetXaxis()->FindBin(0.1+TMath::Pi()/2.0));
//     ScalingMETrigger +=SameHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetNbinsX());
// //     cout << ScalingMETrigger<<" ";
//     ScalingMETrigger -=SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-1.0),SameHistPhiPhi->GetXaxis()->FindBin(4.0));//Diagonal(SameHistPhiPhi,SameHistPhiPhi->GetXaxis()->FindBin(2.0),SameHistPhiPhi->GetXaxis()->FindBin(2.0)+1,SameHistPhiPhi->GetXaxis()->FindBin(1.0),SameHistPhiPhi->GetXaxis()->FindBin(1.0)+1)/(2.0);
//     ScalingMETrigger -=SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(4.0),SameHistPhiPhi->GetXaxis()->FindBin(-1.0));//Diagonal(SameHistPhiPhi,SameHistPhiPhi->GetXaxis()->FindBin(1.0),SameHistPhiPhi->GetXaxis()->FindBin(1.0)+1,SameHistPhiPhi->GetXaxis()->FindBin(2.0),SameHistPhiPhi->GetXaxis()->FindBin(2.0)+1)/(2.0);
// //     cout << ScalingMETrigger<<endl;
//     Double_t DividebyMETrigger = METRiggerHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetNbinsX()); //Diagonal(METRiggerHistPhiPhi,METRiggerHistPhiPhi->GetXaxis()->FindBin(1.2),METRiggerHistPhiPhi->GetXaxis()->FindBin(1.4),METRiggerHistPhiPhi->GetXaxis()->FindBin(1.2),METRiggerHistPhiPhi->GetXaxis()->FindBin(1.4))/(1.0*(SameHistPhiPhi->GetXaxis()->FindBin(1.2)-SameHistPhiPhi->GetXaxis()->FindBin(1.4)));
//     DividebyMETrigger += METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetNbinsX(),1);
// //     cout << DividebyMETrigger<<" ";
//     DividebyMETrigger -=METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-1.0),SameHistPhiPhi->GetXaxis()->FindBin(4.0));//Diagonal(METRiggerHistPhiPhi,METRiggerHistPhiPhi->GetXaxis()->FindBin(2.0),METRiggerHistPhiPhi->GetXaxis()->FindBin(2.0)+1,METRiggerHistPhiPhi->GetXaxis()->FindBin(1.0),METRiggerHistPhiPhi->GetXaxis()->FindBin(1.0)+1)/2.0;
//     DividebyMETrigger -=METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(4.0),SameHistPhiPhi->GetXaxis()->FindBin(-1.0));//Diagonal(METRiggerHistPhiPhi,METRiggerHistPhiPhi->GetXaxis()->FindBin(1.0),METRiggerHistPhiPhi->GetXaxis()->FindBin(1.0)+1,METRiggerHistPhiPhi->GetXaxis()->FindBin(2.0),METRiggerHistPhiPhi->GetXaxis()->FindBin(2.0)+1)/2.0;
// //     cout << DividebyMETrigger<<endl;
// 
//     if(DividebyMETrigger>0.0) *METriggerscale = ScalingMETrigger/DividebyMETrigger;
//     else *METriggerscale = 0.0;
//     if(*METriggerscale<0.0)*METriggerscale = 0.0;
//     
//     if(first==2){
//       ScalingMETrigger =SameHistPhiPhi->GetBinContent(2,2);
//       ScalingMETrigger -=0.5*SameHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetXaxis()->FindBin(-0.3));
//       ScalingMETrigger -=0.5*SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-0.3),1);
// //       cout << ScalingMETrigger<<endl;
//       DividebyMETrigger = METRiggerHistPhiPhi->GetBinContent(2,2); 
//       DividebyMETrigger -=0.5*METRiggerHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetXaxis()->FindBin(-0.3));
//       DividebyMETrigger -=0.5*METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-0.3),1);
// //       cout << DividebyMETrigger<<endl;
//       if(DividebyMETrigger>0.0) *METriggerscale = ScalingMETrigger/DividebyMETrigger;
//       else *METriggerscale = 0.0;
// //       if(*METriggerscale<0.0)*METriggerscale = 0.0;
// //       cout <<*METriggerscale<<endl;
//     }
// //     Double_t ScalingMETA =SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(0.0),5);
// //     ScalingMETA +=SameHistPhiPhi->GetBinContent(5,SameHistPhiPhi->GetYaxis()->FindBin(0.0));
// //     ScalingMETA -=SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(4.0),5);
// //     ScalingMETA -=SameHistPhiPhi->GetBinContent(5,SameHistPhiPhi->GetYaxis()->FindBin(4.0));
// //     Double_t Divideby = METAHistPhiPhi->GetBinContent(METAHistPhiPhi->GetXaxis()->FindBin(0.0),5);
// //     Divideby +=METAHistPhiPhi->GetBinContent(5,METAHistPhiPhi->GetYaxis()->FindBin(0.0));
// //     Divideby -=METAHistPhiPhi->GetBinContent(METAHistPhiPhi->GetXaxis()->FindBin(4.0),5);
// //     Divideby -=METAHistPhiPhi->GetBinContent(5,METAHistPhiPhi->GetYaxis()->FindBin(4.0));
// //     if(Divideby >0.0) *METAscale = ScalingMETA/Divideby;    
// //     else *METAscale = 0.0;
// //     
// //     
// //     Double_t ScalingMETrigger = SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(0.0),SameHistPhiEta->GetXaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/3.0),SameHistPhiEta->GetYaxis()->FindBin(2.0*TMath::Pi()/3.0));
// //     ScalingMETrigger -=0.5* SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-0.75),SameHistPhiEta->GetXaxis()->FindBin(-0.75),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/3.0),SameHistPhiEta->GetYaxis()->FindBin(2.0*TMath::Pi()/3.0));
// //     ScalingMETrigger -=0.5* SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(0.75),SameHistPhiEta->GetXaxis()->FindBin(0.75),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/3.0),SameHistPhiEta->GetYaxis()->FindBin(2.0*TMath::Pi()/3.0));
// //     Double_t DividebyMETrigger = METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(0.0),METRiggerHistPhiEta->GetXaxis()->FindBin(0.0),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/3.0),METRiggerHistPhiEta->GetYaxis()->FindBin(2.0*TMath::Pi()/3.0));
// //     DividebyMETrigger -=0.5* METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(-0.75),METRiggerHistPhiEta->GetXaxis()->FindBin(-0.75),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/3.0),METRiggerHistPhiEta->GetYaxis()->FindBin(2.0*TMath::Pi()/3.0));
// //     DividebyMETrigger -=0.5* METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(0.75),METRiggerHistPhiEta->GetXaxis()->FindBin(0.75),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/3.0),METRiggerHistPhiEta->GetYaxis()->FindBin(2.0*TMath::Pi()/3.0));
// //     if(DividebyMETrigger>0.0) *METriggerscale = ScalingMETrigger/DividebyMETrigger;
// //     else *METriggerscale = 0.0;
//   }
  
  
}


void Correct(const char* histname, BinDirs * BinDir,const char* in,const char* out, Double_t  METAscale, Double_t  METriggerscale)
{
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get(histname);if(!tmp){return;}
  TH2D * SameHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname)));
  tmp = BinDir->METAdiv()->Get(histname);if(!tmp){return;}
  TH2D * METAHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname)));
  tmp = BinDir->METriggerdiv()->Get(histname);if(!tmp){return;}
  TH2D * METriggerHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname)));  
  delete tmp;
  //Scale the Mixed hists with the correct scales:
  METAHist->Scale(METAscale);
  METriggerHist->Scale(METriggerscale);
 //perform the removal:  
  SameHist->Add(METAHist,-1.0);SameHist->Add(METriggerHist,-1.0);
//   RemovePlateau(SameHist->Integral(SameHist->GetXaxis()->FindBin(-1.6),SameHist->GetXaxis()->FindBin(1.6),SameHist->GetYaxis()->FindBin(TMath::Pi()*0.5),SameHist->GetYaxis()->FindBin(TMath::Pi()*0.5))/(SameHist->GetXaxis()->FindBin(1.6)-SameHist->GetXaxis()->FindBin(-1.6)),SameHist);
  TCanvas * canvas = Makecanvas(SameHist,Form("%sCanvas",histname),false);
  BinDir->SameDir(out)->cd();
  SameHist->Write(histname);
  canvas->Write();
  delete canvas;delete SameHist;delete METAHist; delete METriggerHist;
}
void Correct(const char* histname1,const char* histname2,const char* histname3,const char* histname4, BinDirs * BinDir,const char* in,const char* out, Double_t  METAscale, Double_t  METriggerscale)
{
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get(histname1);if(!tmp){return;}
  TH2D * SameHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname1)));
  tmp = BinDir->METAdiv()->Get(histname1);if(!tmp){return;}
  TH2D * METAHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname1)));
  tmp = BinDir->METriggerdiv()->Get(histname1);if(!tmp){return;}
  TH2D * METriggerHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname1)));
  tmp = BinDir->SameDir(in)->Get(histname2);if(!tmp){return;}
  TH2D * SameHist2 = dynamic_cast<TH2D*>(BinDir->SameDir(in)->Get(histname2)->Clone(Form("%sSame",histname2)));
  tmp = BinDir->METAdiv()->Get(histname2);if(!tmp){return;}
  TH2D * METAHist2 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname2)));
  tmp = BinDir->METriggerdiv()->Get(histname2);if(!tmp){return;}
  TH2D * METriggerHist2 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname2)));
  tmp = BinDir->SameDir(in)->Get(histname3);if(!tmp){return;}
  TH2D * SameHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname3)));
  tmp = BinDir->METAdiv()->Get(histname3);if(!tmp){return;}
  TH2D * METAHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname3)));
  tmp = BinDir->METriggerdiv()->Get(histname3);if(!tmp){return;}
  TH2D * METriggerHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname3)));  
  tmp = BinDir->SameDir(in)->Get(histname4);if(!tmp){return;}
  TH2D * SameHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname4)));
  tmp = BinDir->METAdiv()->Get(histname4);if(!tmp){return;}
  TH2D * METAHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname4)));
  tmp = BinDir->METriggerdiv()->Get(histname4);if(!tmp){return;}
  TH2D * METriggerHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname4)));
  delete tmp;
  //Scale the Mixed hists with the correct scales:
  METAHist1->Scale(METAscale);
  METAHist2->Scale(METAscale);
  METAHist3->Scale(METAscale);
  METAHist4->Scale(METAscale);
  METriggerHist1->Scale(METriggerscale);
  METriggerHist2->Scale(METriggerscale);
  METriggerHist3->Scale(METriggerscale);
  METriggerHist4->Scale(METriggerscale);
  //perform the removal:  
  SameHist1->Add(METAHist1,-1.0);SameHist1->Add(METriggerHist1,-1.0);
  SameHist2->Add(METAHist2,-1.0);SameHist2->Add(METriggerHist2,-1.0);
  SameHist3->Add(METAHist3,-1.0);SameHist3->Add(METriggerHist3,-1.0);
  SameHist4->Add(METAHist4,-1.0);SameHist4->Add(METriggerHist4,-1.0);
  TCanvas * canvas = Makecanvas(SameHist1,SameHist2,SameHist3,SameHist4,"DPHIDPHI",false);
  BinDir->SameDir(out)->cd();
  SameHist1->Write(histname1);
  SameHist2->Write(histname2);
  SameHist3->Write(histname3);
  SameHist4->Write(histname4);
  canvas->Write();
  delete SameHist1;delete SameHist2;delete SameHist3;delete SameHist4;
  delete METAHist1;delete METAHist2;delete METAHist3;delete METAHist4;
  delete METriggerHist1;delete METriggerHist2;delete METriggerHist3;delete METriggerHist4;
  delete canvas;
}

void Correctpp(BinDirs* BinDir)
{
  //If the bin is empty, skip the entire bin:
  if(!dynamic_cast<TH1D*>(BinDir->Samediv()->Get("number_of_triggers"))){cout << "skipping empty bin "<<BinDir->path()<<endl;return;}
  if(dynamic_cast<TH1D*>(BinDir->Samediv()->Get("number_of_triggers"))->GetEffectiveEntries()<1){cout << "skipping empty bin "<<BinDir->path()<<endl;return;}
  //create the scaling doubles.
  Double_t METAScale 		= 0.0;
  Double_t METriggerScale 	= 0.0;
  //First iteration: 
  GetScalingFactors(BinDir,"divided","iteration1",&METAScale,&METriggerScale,true,false,1);
  cout << BinDir->path()<<" META:"<<METAScale<<" METriggerScale:"<<METriggerScale<<endl;
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,METriggerScale,"iteration1");
  Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","iteration1",METAScale,METriggerScale);
  
}

void CorrectPbPb(BinDirs* BinDir)
{
  //If the bin is empty, skip the entire bin:
  if(!dynamic_cast<TH1D*>(BinDir->Samediv()->Get("number_of_triggers"))){cout << "skipping empty bin "<<BinDir->path()<<endl;return;}
  //create the scaling doubles.
  Double_t METAScale 		= 0.0;
  Double_t METriggerScale 	= 0.0;
  //First iteration: 
  GetScalingFactors(BinDir,"divided","iteration1",&METAScale,&METriggerScale,false,true,1);
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,METriggerScale,"iteration1");
  Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","iteration1",METAScale,METriggerScale);  
//   GetScalingFactors(BinDir,"iteration1","iteration2",&METAScale,&METriggerScale,false,true,2);
//   //Create a canvas that tells what is substracted.
//   savedircontent(BinDir,0.0,METriggerScale,"iteration2");
//   Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"iteration1","iteration2",0.0,METriggerScale);
//   Correct("DPhi_1_DEta_12",BinDir,"iteration1","iteration2",0.0,METriggerScale);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"iteration1","iteration2",0.0,METriggerScale);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"iteration1","iteration2",0.0,METriggerScale);
//   Correct("DPhi_1_DEta_12_SameSide",BinDir,"iteration1","iteration2",0.0,METriggerScale);  
//   GetScalingFactors(BinDir,"iteration2","iteration3",&METAScale,&METriggerScale,false,true,2);
//   //Create a canvas that tells what is substracted.
//   savedircontent(BinDir,METAScale,METriggerScale,"iteration3");
//   Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"iteration2","iteration3",METAScale,METriggerScale);
//   Correct("DPhi_1_DEta_12",BinDir,"iteration2","iteration3",METAScale,METriggerScale);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"iteration2","iteration3",METAScale,METriggerScale);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"iteration2","iteration3",METAScale,METriggerScale);
//   Correct("DPhi_1_DEta_12_SameSide",BinDir,"iteration2","iteration3",METAScale,METriggerScale);  

}   
////////////////////////////////////



void CollectMVbins(){
  //take the histograms from the bins and average over them.
  TFile * outfile = new TFile("results.root","UPDATE");
  TObjArray * multdirlist = new TObjArray(6);
  BinDirs * divsame = new BinDirs(outfile->GetDirectory("/"),outfile->GetDirectory("META"),outfile->GetDirectory("METrigger"),true);
  multdirlist->Add(divsame);
  //List of directories for multiplicity bins:
  TList * directories = GetMZDirectories(divsame);
  //Go through the list and find the multiplicity/centrality binning in order to sum over them and add them to the TObjArray:
  TString  s = TString("");
  for(int i=0; i<directories->GetEntries();i++){
    if(TString(directories->At(i)->GetName()).Contains("Z(")){
      TString * bin = new TString(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName());
      if(bin->CompareTo(s.Data())!=0&&bin->Contains("BinM")){
	BinDirs * dir = new BinDirs(outfile,bin->Data(),true);
	multdirlist->Add(dir);
	s.Clear();
	s.Append(bin->Data());
      }
    }
  }
  
   //Get Tokens for all histograms:
   TStringToken histtokensbinstats = GetHistTokens(divsame->Same()->GetDirectory(Form("%s/bin_stats",directories->At(1)->GetName())));
   while(histtokensbinstats.NextToken()){CollectHistbinstats(histtokensbinstats.Data(),directories,multdirlist);}
   TStringToken histtokens = GetHistTokens(divsame->Same()->GetDirectory(Form("%s/divided",directories->At(1)->GetName())));
   while(histtokens.NextToken()){CollectHist(histtokens.Data(),directories,multdirlist);}
  outfile->Close();
  delete outfile;
}

void correct(const char* options)
{
  bool ispp = false;
  bool isPbPb = false;
  TString delimiter(" ");
  TStringToken token(options, delimiter);
  while (token.NextToken()) {
    TString argument=token;
    if (argument.CompareTo("-h")==0 ||argument.CompareTo("--help")==0) {
      cout<<Form(" options:"
		    "\n\t  pp   - uses pp fitting ranges to substract the correlated background."
		    "\n\t  PbPb     - uses PbPb fitting ranges to substract the correlated background."
		    )<<endl;
      return;
      }
    if(argument.CompareTo("pp")==0){
      ispp = true;
      isPbPb = false;
      
    }
    if(argument.CompareTo("PbPb")==0){
      isPbPb = true;
      ispp = false;
   }
  }
  //open the file in UPDATE mode:
  TFile * rfile = TFile::Open("results.root","UPDATE");
  //Find a list over all directories beginning with BinM and create a BinDirs object for each:
  TList * dirs = rfile->GetListOfKeys();
  TObjArray * BinDirAr = new TObjArray();
  for(int i = 0;i<dirs->GetEntries();i++)
  {
    if(TString(dirs->At(i)->GetName()).BeginsWith("BinM")){
      BinDirs * dir = new BinDirs(rfile,dirs->At(i)->GetName(),false);
      BinDirAr->Add(dir);
    }
    else if(TString(dirs->At(i)->GetName()).CompareTo("divided")==0){
      BinDirs * dir = new BinDirs(rfile->GetDirectory(""),rfile->GetDirectory("META"),rfile->GetDirectory("METrigger"),false);
      BinDirAr->Add(dir);
    }
  }
  for(int i=0;i<BinDirAr->GetEntriesFast();i++){
    if(ispp)	Correctpp(	dynamic_cast<BinDirs*>(BinDirAr->At(i)));
    if(isPbPb)	CorrectPbPb(	dynamic_cast<BinDirs*>(BinDirAr->At(i)));
  }
  rfile->Close();
}

void yield(){
  TFile * rfile = TFile::Open("results.root","UPDATE");
  TObjArray * dirarray = new TObjArray();
  TObjArray * yieldarray = new TObjArray();
//   TObjArray * yieldarrayap1 = new TObjArray();

  TList * dirs = rfile->GetListOfKeys();

  
  for(int i = 0;i<dirs->GetEntries();i++){
    if(TString(dirs->At(i)->GetName()).BeginsWith("BinM")){
      unsigned int j=1;
      cout << dirs->At(i)->GetName()<<endl;
      while(j<10){
	TDirectory * tmp1 = dynamic_cast<TDirectory*>(rfile->Get(Form("%s/iteration%u",dirs->At(i)->GetName(),j)));
	if(!tmp1){j = 11;continue;}
	TDirectory * tmp2 = dynamic_cast<TDirectory*>(rfile->Get(Form("%s/iteration%u",dirs->At(i)->GetName(),j+1)));
	if(tmp1&&!tmp2){
	  dirarray->Add(tmp1);
	  tmp1->pwd();
	  yieldarray->Add(resultsdirectory(rfile->GetDirectory(Form("%s",dirs->At(i)->GetName())),"yield"));
// 	  yieldarrayap1->Add(resultsdirectory(rfile->GetDirectory(Form("%s",dirs->At(i)->GetName())),"yieldap1"));
	  j = 11;
	}
	if(tmp1&&tmp2){j +=1;}
      }
    }
  }
  unsigned int k = 1;
  while(k<10){
    TDirectory * tmp1 = dynamic_cast<TDirectory*>(rfile->Get(Form("iteration%u",k)));
    if(!tmp1){k = 11;continue;}
    TDirectory * tmp2 = dynamic_cast<TDirectory*>(rfile->Get(Form("iteration%u",k+1)));
    if(tmp1&&!tmp2){
      dirarray->Add(rfile->GetDirectory(Form("iteration%u",k)));
      yieldarray->Add(resultsdirectory(rfile->GetDirectory(""),"yield"));
//       yieldarrayap1->Add(resultsdirectory(rfile->GetDirectory(""),"yieldap1"));
      k = 11;
    }
    if(tmp1&&tmp2){k +=1;}
  }
  TDirectory * dir;
  TDirectory * yielddir;
  for(int i=0;i<dirarray->GetEntriesFast();i++){
    dir = dynamic_cast<TDirectory*>(dirarray->At(i));
    yielddir = dynamic_cast<TDirectory*>(yieldarray->At(i));
    extractbinyield(dir,yielddir);
  }
  
}
void draw(){
  TFile * rfile = TFile::Open("results.root","UPDATE");

  BinDirs * dir = new BinDirs(rfile, "BinM(0.00)->(44.44)");
  dir->resultsdirectory("hello");

  
//   TFile * rfile = TFile::Open("results.root","READ");
//   TList * dirs = rfile->GetListOfKeys();
//   TString * dirss = new TString("");
//   for(int i = 0;i<dirs->GetEntries();i++){if(TString(dirs->At(i)->GetName()).BeginsWith("BinM("))dirss->Append(Form("%s ",dirs->At(i)->GetName()));}
//   TStringToken dirtokens(dirss->Data()," ");
//   TDirectory * dir;
// //   TDirectory * ydir;
//   TDirectory * ndir;
// //   TH1D* plothist;
//   TH1D* dirhist;
// //   TH1D* savehist = NULL;
// //   TH1D* numberoftracks = new TH1D("ntracks","Average number of associated particles in event category.",1,0.0,1.0);
//   Int_t colorin = 1;
//   Int_t symbolint =20;
//   TCanvas * averagecanvas = new TCanvas("average");
//   averagecanvas->cd();
// 
// //   TCanvas * yieldcanvas = new TCanvas("yield");
// //   yieldcanvas->cd();
//   while(dirtokens.NextToken()){
//     dir = rfile->GetDirectory(dirtokens.Data());
//     if(dirtokens.Contains("-10.00)->(-5.00)")){symbolint = 20;colorin+=1;if(colorin==5)colorin+=1;}
// //     ydir = dir->GetDirectory("yield");
// //     if(ydir){
// //       plothist =dynamic_cast<TH1D*>(ydir->Get("dphiyields"));
// //       if(plothist){
// // 	plothist->SetBit(TH1::kIsAverage);
// // 	if(!savehist){
// // 	  savehist = dynamic_cast<TH1D*>(plothist->Clone("yieldcol"));
// // 	}
// // 	else{
// // 	  savehist->Add(plothist);
// // 	  }
// // // 	plothist->Draw("sameE");
// //       }
// //     }
//     ndir = dir->GetDirectory("bin_stats");
//     if(ndir){
//       dirhist=dynamic_cast<TH1D*>(ndir->Get("number_of_associated"));
//       if(dirhist){
// 	double n = dirhist->GetEffectiveEntries();
// 	if(n!=0){
// 	  Double_t average = 0.0;
// 	  for(int i=1;i<=dirhist->GetNbinsX();i++){
// 	    average += dirhist->GetBinCenter(i)*dirhist->GetBinContent(i)/n;
// 	  }
// 	  TH1D* numberoftracks = new TH1D("ntracks","Average number of associated particles in event category.",1,0.0,1.0);
// 	  numberoftracks->SetStats(0);
// 	  numberoftracks->SetBinContent(1,average);
// 	  numberoftracks->SetMarkerStyle(symbolint);
// 	  numberoftracks->SetMarkerColor(colorin);
// 	  numberoftracks->Draw("sameP");
// // 	  cout << dirtokens.Data() <<" " <<average<<endl;
// 	  
// 	}
//       }
//       
//     }
//     symbolint+=1;
//     if(symbolint==24)symbolint =29;
//   }
// //   cout << savehist<<endl;
// //   yieldcanvas->cd();
// //   savehist->Draw("E");
}