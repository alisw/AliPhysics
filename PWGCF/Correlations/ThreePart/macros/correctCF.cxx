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
#include "TPaveLabel.h"
#include "TObjString.h"
#include "TPaletteAxis.h"
#include "TParameter.h"

using namespace std;
//global pointers for the histograms:

TH2D* gSameEvent =0x0;
TH2D* gMETA = 0x0;
TH2D* gMETrigger = 0x0;

TH2D* gSameEventETA =0x0;
TH2D* gMETAETA = 0x0;
TH2D* gMETriggerETA = 0x0;

bool gislowpTbin = false;

double gkEtaFitRange = 1.7;
double gEtasigRange = 1.5;

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
  BinDirs(TDirectory* dir, const char* bin, bool empty = false){
    if(dir){
      Bin = TString(bin);
      dir->cd();
      Samedir = dir->GetDirectory(bin);
      if(!Samedir)Samedir =  dir->mkdir(bin);
      Samedir->cd();
      if(!Samedir->GetDirectory("bin_stats")) Samedir->mkdir("bin_stats");
      if(!Samedir->GetDirectory("divided")) Samedir->mkdir("divided");
      if(empty) ::resultsdirectory(Samedir,"bin_stats");
      if(empty) ::resultsdirectory(Samedir,"divided");
      TDirectory * METAparentdir = dir->GetDirectory("META");
      if(!METAparentdir)METAparentdir = dir->mkdir("META");
      METAparentdir->cd();
      METAdir = METAparentdir->GetDirectory(bin);
      if(!METAdir)METAdir= METAparentdir->mkdir(bin);
      METAdir->cd();
      if(!METAdir->GetDirectory("bin_stats")) METAdir->mkdir("bin_stats");
      if(!METAdir->GetDirectory("divided")) METAdir->mkdir("divided");
      if(empty) ::resultsdirectory(METAdir,"bin_stats");
      if(empty) ::resultsdirectory(METAdir,"divided");
      TDirectory * METriggerparentdir = dir->GetDirectory("METrigger");
      if(!METriggerparentdir)METriggerparentdir = dir->mkdir("METrigger");
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
    if(samed&&METriggerd&&METAd){
      Bin = TString(samed->GetName());
      samed->cd();
      Samedir = samed;
      if(!Samedir->GetDirectory("bin_stats")) Samedir->mkdir("bin_stats");      
      else if(empty){::resultsdirectory(Samedir,"bin_stats");}
      if(!Samedir->GetDirectory("divided")){Samedir->mkdir("divided");}
      else if(empty){::resultsdirectory(Samedir,"divided");}
      METriggerd->cd();
      METriggerdir = METriggerd;
      if(!METriggerdir->GetDirectory("bin_stats")) METriggerdir->mkdir("bin_stats");
      else if(empty) ::resultsdirectory(METriggerdir,"bin_stats");
      if(!METriggerdir->GetDirectory("divided")) METriggerdir->mkdir("divided");
      else if(empty) ::resultsdirectory(METriggerdir,"divided");
      METAd->cd();
      METAdir = METAd;
      if(!METAdir->GetDirectory("bin_stats")) METAdir->mkdir("bin_stats");      
      else if(empty) ::resultsdirectory(METAdir,"bin_stats");
      if(!METAdir->GetDirectory("divided")) METAdir->mkdir("divided");
      else if(empty) ::resultsdirectory(METAdir,"divided");
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

class ResHists : public TObject{
public:
  ResHists(TString str, TFile * file,bool collected=false){
    TDirectory * dir = file->GetDirectory(str.Data());
    if(collected) dir = dir->GetDirectory("Collected");
    if(!dir){yieldbc = NULL;RMS=NULL;minTpT=0;maxTpT=0;minApT=0;maxApT=0;}
    else{
      TObjArray * getlims = TString(str).Tokenize("_");
      int minTpTi = dynamic_cast<TObjString*>(getlims->At(1))->GetString().Atoi();
      minTpT = minTpTi;
      maxTpT = dynamic_cast<TObjString*>(getlims->At(2))->GetString().Atoi();
      int minApTi = dynamic_cast<TObjString*>(getlims->At(3))->GetString().Atoi();
      if(minApTi == 0) minApT = 0.5;
      else minApT=minApTi; 
      maxApT = dynamic_cast<TObjString*>(dynamic_cast<TObjString*>(getlims->At(4))->GetString().Tokenize("C")->At(0))->GetString().Atoi();
      TDirectory* yielddir = dir->GetDirectory("yield");
      TDirectory* unbinned = yielddir->GetDirectory("original_binning");
      TDirectory* G0bin = unbinned->GetDirectory("GP0");
      yieldbc = dynamic_cast<TH1D*>(G0bin->Get("dphiyieldbc"));
      TH1D* yieldfitc  = dynamic_cast<TH1D*>(G0bin->Get("dphiyield"));
      TH1D* chisq  = dynamic_cast<TH1D*>(G0bin->Get("dphichisq"));
      yieldfitc->Multiply(chisq);
      yieldfit = yieldfitc;
      RMS     = dynamic_cast<TH1D*>(G0bin->Get("dphiRMS"));
      TDirectory* itdir = dir->GetDirectory("iteration1");
      dphideta12 = dynamic_cast<TH2D*>(itdir->Get("DPhi_1_DEta_12_SameSide"));
    }
  }
  TH1D * GetYield(){return yieldbc;}
  TH1D * GetYieldfit(){return yieldfit;}
  TH2D * Get2Dhist(){return dphideta12;}
  Double_t GetMinTpT(){return minTpT;}
  Double_t GetMaxTpT(){return maxTpT;}
  Double_t GetMinApT(){return minApT;}
  Double_t GetMaxApT(){return maxApT;}

private:
  TH1D * yieldbc;
  TH1D * RMS;
  TH1D * yieldfit;
  TH1D * RMSchisq;
  TH2D * dphideta12;
  Double_t minTpT;
  Double_t maxTpT;
  Double_t minApT;
  Double_t maxApT; 
  ClassDef(ResHists, 1);
};

  Double_t CGausPol0(Double_t * x, Double_t * par){
    //Parameterization of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1];
    Double_t s=TMath::Abs(par[2]);
    Double_t dx=(x[0]-c);
    return par[0]*exp(-dx*dx/(2.*s*s))/TMath::Sqrt(TMath::Pi())/s+par[3];
  }
  Double_t CBKGPol0(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    //reject if inside +- gEtasigRange:
    Double_t c=par[1];
    Double_t dx=(x[0]-c);
    
    if(dx<gEtasigRange&&dx>-gEtasigRange){TF1::RejectPoint();return 0;}
    return par[0] ;
  }
  Double_t CBKGPol0p(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    if(TMath::Abs(x[0])>2.0)return 0.0;
    return par[0] ;
  }  
  Double_t CGausAPol1(Double_t * x, Double_t * par){
    //Parameterization of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1];
    Double_t s=TMath::Abs(par[2]);
    Double_t dx=(x[0]-c);
    return par[0]*exp(-dx*dx/(2.*s*s))/TMath::Sqrt(TMath::Pi())/s+par[3]-TMath::Abs(dx)*par[4];
  }
  Double_t CBKGAPol1(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    Double_t c=par[1];
    Double_t dx=(x[0]-c);
    if(dx<gEtasigRange&&dx>-gEtasigRange){TF1::RejectPoint();return 0;}
    return par[0]-TMath::Abs(dx)*par[2];
  }
  Double_t CBKGAPol1p(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    Double_t c=par[1];
    Double_t dx=(x[0]-c);
    return par[0]-TMath::Abs(dx)*par[2];
  }  
  Double_t CGausAPol2(Double_t * x, Double_t * par){
    //Parameterizatin of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1];
    Double_t s=TMath::Abs(par[2]);
    Double_t dx=(x[0]-c);
    return par[0]*exp(-dx*dx/(2.*s*s))/TMath::Sqrt(TMath::Pi())/s+par[3]-TMath::Abs(dx)*par[4] + dx*dx*par[5];
  }
  Double_t CBKGAPol2(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    Double_t c=par[1];
    Double_t dx=(x[0]-c);
    if(dx<gEtasigRange&&dx>-gEtasigRange){TF1::RejectPoint();return 0;}
    return par[0]-TMath::Abs(dx)*par[2] + dx*dx*par[3];
  }  
  Double_t CBKGAPol2p(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    Double_t c=par[1];
    Double_t dx=(x[0]-c);
    return par[0]-TMath::Abs(dx)*par[2] + dx*dx*par[3];
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

void AddSigBins(TH2D* hist, TH2D* METAhist, TH2D* METriggerhist){
  //Calculates hist -METAhist -METriggerhist for nonzero bins.
  if(!hist||!METAhist||!METriggerhist)return;
  for(int i = 0; i<=hist->GetNbinsX();i++){
    for(int j=0; j<=hist->GetNbinsY();j++){
      if(hist->GetBinContent(i,j)>1.0E-10){
	Double_t content = hist->GetBinContent(i,j);
	Double_t error = hist->GetBinError(i,j);
	Double_t METAcontent = METAhist->GetBinContent(i,j);
	Double_t METAerror = METAhist->GetBinError(i,j);
	Double_t METriggercontent = METriggerhist->GetBinContent(i,j);
	Double_t METriggererror = METriggerhist->GetBinError(i,j);
	content -= METAcontent + METriggercontent;
	hist->SetBinContent(i,j,content);
	hist->SetBinError(i,j,TMath::Sqrt(error*error + METAerror*METAerror+METriggererror*METriggererror));
      }
    }
  }
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

void removeconstant(TH1D * hist, Double_t plateau, Double_t erroronit){
  for(int x =0;x<=hist->GetNbinsX();x++){
    if(TMath::Abs(hist->GetBinContent(x))>1.0E-10){
      hist->SetBinContent(x,hist->GetBinContent(x) + plateau);
      hist->SetBinError(x,TMath::Sqrt(hist->GetBinError(x)*hist->GetBinError(x) + erroronit*erroronit));
    }
  }
}

void fitwith(TDirectory * dir, const char* type, TH2D * histo,Double_t etalimit){
  //Fit the hist with the method specified in type.
  //Types: fitfunc/rebin:
  //Fitfunc: fgauspol1 , fgauspol2 , fgauspol0
  //rebin: adds together n bins for /n.
  TObjArray * types = TString(type).Tokenize("/");
  int rebin = dynamic_cast<TObjString*>(types->At(1))->GetString().Atoi();
  dir->cd();
  TH2D* hist = dynamic_cast<TH2D*>(histo->Clone("DPhi_1_DEta_12_SameSidecc"));
  if(rebin>1){hist->RebinX(3);hist->RebinY(2);}
  hist->Write("DPhi_1_DEta_12_SameSide");

  TH1D* deta12ss = hist->ProjectionX("DEta12");
  TH1D* deta12ssdraw = hist->ProjectionX("DEta12d");
  TH1D* deta12ssbinc = hist->ProjectionX("DEta12d");

  deta12ss->Reset();deta12ssdraw->Reset();deta12ssbinc->Reset();
  
  TF1* fitsig;
  TF1* fitbg;
  TF1* rembg;
  
  TDirectory* typedir;
  TDirectory* bindir;
  
  TString title = TString("#Delta#eta_{12} distribution in bin");
  int color;
//   const char* lable;
  TH1D* yield = hist->ProjectionY("dphiyield");
  yield->Reset();
  yield->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  TH1D* yieldbc = hist->ProjectionY("dphiyieldbc");
  yieldbc->Reset();
  yield->SetTitle("Yield from bin counting vs #Delta#Phi_1");
  yieldbc->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  
  TH1D* background = dynamic_cast<TH1D*>(yield->Clone("dphibackground"));
  background->SetTitle("Height of the background as a function of #Delta#Phi#.");
  background->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  TH1D* width = dynamic_cast<TH1D*>(yield->Clone("dphiwidth")); 
  width->SetTitle("Width of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  width->GetYaxis()->SetTitle("width of the peak (rad)");
  TH1D* hRMS = dynamic_cast<TH1D*>(yield->Clone("dphiRMS")); 
  width->SetTitle("RMS of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  width->GetYaxis()->SetTitle("RMS of the peak (rad)");

  TH1D* peakpos = dynamic_cast<TH1D*>(yield->Clone("dphipos")); 
  peakpos->SetTitle("Position of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  peakpos->SetYTitle("Position of the peak (rad)");
  TH1D* chisq = dynamic_cast<TH1D*>(yield->Clone("dphichisq"));
  chisq->SetTitle("#Chi^2/NDF of the fit as a function of #Delta#Phi");
  chisq->SetYTitle("#Chi^2/NDF ");
  TH1D* prob = dynamic_cast<TH1D*>(yield->Clone("dphiprob"));
  prob->SetTitle("Probability of the fit as a function of #Delta#Phi");
  prob->SetYTitle("Probability ");

  TH1D * etasigrange = dynamic_cast<TH1D*>(yield->Clone("detasigrange"));
  etasigrange->SetTitle("etarange bin by bin");
  etasigrange->SetYTitle("eta range ");
  
  Double_t AvWidth = 1.0;
  if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol0")==0){
    typedir = resultsdirectory(dir,"GP0");    
    bindir = resultsdirectory(typedir,"bins");
    //flat background and a gaussian:
    color = 2;
//     lable = "flat bg + gaussian";
    yield->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a flat background and a gaussian.");
    yield->SetLineColor(color);
    background->SetLineColor(color);
    width->SetLineColor(color);
    peakpos->SetLineColor(color);  
    chisq->SetLineColor(color);  
    prob->SetLineColor(color);    
    //Initialize the fit function:
    fitsig = new TF1("fgs",CGausPol0,-gkEtaFitRange,gkEtaFitRange,4) ;
    fitsig->SetParNames("peakhight", "peakpos", "peakwidth", "B") ;
    fitsig->SetLineColor(color);
    fitsig->SetLineWidth(1);
    fitsig->SetParameters(0.0,0.0,0.1,0.0);
    fitsig->SetParLimits(0,0.0,0.3);
    fitsig->SetParLimits(1,-0.1,0.1);
    fitsig->SetParLimits(2,0.01,0.8);  
    //Initialize the fit function for BKG:
    fitbg = new TF1("fgs",CBKGPol0,-gkEtaFitRange,gkEtaFitRange,2) ;
    fitbg->SetParNames("B", "peakpos" ) ;
    fitbg->SetLineColor(color);
    fitbg->SetLineWidth(1);
    fitbg->SetParameters(0.0,0.0);
//     fitbg->SetParLimits(0,0.0,1.0E10);
    fitbg->SetParLimits(1,-0.1,0.1);    
    //Initialize the fit function for BKG:
    rembg = new TF1("fgs",CBKGPol0p,-gkEtaFitRange,gkEtaFitRange,1) ;
    rembg->SetParNames("B", "peakpos" ) ;
    rembg->SetLineColor(color);
    rembg->SetLineWidth(1);
    rembg->SetParameters(0.0,0.0);    
    
    //find the width at 0:
    int bin1 = yield->GetXaxis()->FindBin(-0.001);
    int bin2 = yield->GetXaxis()->FindBin(0.001);
    
    fillwithbinnr(hist,deta12ss,bin1);
    fitsig->FixParameter(1,0.0);
    TFitResultPtr fitresultwidth1 = deta12ss->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    Double_t width1 = fitsig->GetParameter(2);
    fillwithbinnr(hist,deta12ss,bin2);
    TFitResultPtr fitresultwidth2 = deta12ss->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    Double_t width2 = fitsig->GetParameter(2);  
    deta12ss->Reset();
    fitsig->ReleaseParameter(1);
    
    fitsig->FixParameter(2,width1+width2/2.0);
    AvWidth = 3.0*fitsig->GetParameter(2);
    cout << AvWidth<<endl;
    if(3.0*fitsig->GetParameter(2)<etalimit){AvWidth = 3.0*fitsig->GetParameter(2);}//set the range to 3*/sigma
    else {AvWidth = etalimit;}
    cout << AvWidth<<endl;
    fitsig->ReleaseParameter(2);
    fitsig->SetParLimits(2,0.01,0.8);  

  }
  if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol1")==0){
    typedir = resultsdirectory(dir,"GP1");
    bindir = resultsdirectory(typedir,"bins");
    //pol1 background and a gaussian:
    color = 3;
//     lable = "pol1 bg + gaussian";
    yield->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a pol1 background and a gaussian.");
    yield->SetLineColor(color);
    background->SetLineColor(color);
    width->SetLineColor(color);
    peakpos->SetLineColor(color);  
    chisq->SetLineColor(color);  
    prob->SetLineColor(color);
    //Initialize the fit function:
    fitsig = new TF1("p1gs",CGausAPol1,-gkEtaFitRange,gkEtaFitRange,5) ;
    fitsig->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C") ;
    fitsig->SetLineColor(color);
    fitsig->SetLineWidth(1);
    fitsig->SetParameters(0.0,0.0,0.1,0.0,0.0);
    fitsig->SetParLimits(0,0.0,0.3);
    fitsig->SetParLimits(1,-0.1,0.1);
    fitsig->SetParLimits(2,0.1,0.8);  
    fitsig->SetParLimits(4,0.0,1.0);
    //Initialize the fit function for BKG:
    fitbg = new TF1("p1gs",CBKGAPol1,-gkEtaFitRange,gkEtaFitRange,3) ;
    fitbg->SetParNames("B","peakpos", "C") ;
    fitbg->SetLineColor(color);
    fitbg->SetLineWidth(1);
    fitbg->SetParameters(0.0,0.0,0.0);
//     fitbg->SetParLimits(0,0.0,1.0E10);
    fitbg->SetParLimits(1,-0.1,0.1);
    fitbg->SetParLimits(2,0.0,1.0);
    //Initialize the fit function for BKG:
    rembg = new TF1("p1gs",CBKGAPol1p,-gkEtaFitRange,gkEtaFitRange,3) ;
    rembg->SetParNames("B","peakpos", "C") ;
    rembg->SetLineColor(color);
    rembg->SetLineWidth(1);
    rembg->SetParameters(0.0,0.0,0.0);
//     fitbg->SetParLimits(0,0.0,1.0E10);
    rembg->SetParLimits(1,-0.1,0.1);
    rembg->SetParLimits(2,0.0,1.0);
  }
  if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol2")==0){
    typedir = resultsdirectory(dir,"GP2");  
    bindir  = resultsdirectory(typedir,"bins");    
    //pol2 background and a gaussian:
    color = 4;
//     lable = "pol2 bg + gaussian";
    yield->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a pol2 background and a gaussian.");
    yield->SetLineColor(color);
    background->SetLineColor(color);
    width->SetLineColor(color);
    peakpos->SetLineColor(color);  
    chisq->SetLineColor(color);  
    prob->SetLineColor(color);
    //Initialize the fit function:
    fitsig= new TF1("p2gs",CGausAPol2,-gkEtaFitRange,gkEtaFitRange,6) ;
    fitsig->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C", "D") ;
    fitsig->SetLineColor(color);
    fitsig->SetLineWidth(1);
    fitsig->SetParameters(0.0,0.0,0.1,0.0,0.0,0.0);
    fitsig->SetParLimits(0,0.0,0.3);
    fitsig->SetParLimits(1,-0.1,0.1);
    fitsig->SetParLimits(2,0.1,0.8);  
    fitsig->SetParLimits(4,0.0,1.0);
    //Initialize the fit function:
    fitbg= new TF1("p2gs",CBKGAPol2,-gkEtaFitRange,gkEtaFitRange,4) ;
    fitbg->SetParNames("B", "peakpos",  "C", "D") ;
    fitbg->SetLineColor(color);
    fitbg->SetLineWidth(1);
    fitbg->SetParameters(0.0,0.0,0.0,0.0);
//     fitbg->SetParLimits(0,0.0,1.0E10);
    fitbg->SetParLimits(1,-0.1,0.1);
    fitbg->SetParLimits(2,0.0,1.0);
    //Initialize the fit function:
    rembg= new TF1("p2gs",CBKGAPol2p,-gkEtaFitRange,gkEtaFitRange,4) ;
    rembg->SetParNames("B", "peakpos",  "C", "D") ;
    rembg->SetLineColor(color);
    rembg->SetLineWidth(1);
    rembg->SetParameters(0.0,0.0,0.0,0.0);
//     fitbg->SetParLimits(0,0.0,1.0E10);
    rembg->SetParLimits(1,-0.1,0.1);
    rembg->SetParLimits(2,0.0,1.0);    
  }
  double deltaphi = (yield->GetXaxis()->GetBinCenter(2)-yield->GetXaxis()->GetBinCenter(1));//Width of the bin in phi
  bindir->cd();
  for(int dphi=1;dphi<=hist->GetNbinsY();dphi++){
    gEtasigRange = AvWidth;
    deta12ss->SetTitle(Form("%s %4.2f #pi < #Delta#Phi < %4.2f #pi",title.Data(),hist->GetYaxis()->GetBinLowEdge(dphi)/TMath::Pi(),hist->GetYaxis()->GetBinUpEdge(dphi)/TMath::Pi()));
    fillwithbinnr(hist,deta12ss,dphi);
    deta12ss->SetStats(false);
    fitsig->FixParameter(1,0.0);
    TFitResultPtr fitresult = deta12ss->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    if(int(fitresult)!=4000)//if "error", try again with more:
      fitresult = deta12ss->Fit(fitsig,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
    deta12ss->Write(Form("%sflatg_%i_first",deta12ss->GetName(),dphi));
    fitbg->FixParameter(1,0.0);
//     if(3.0*fitsig->GetParameter(2)<1.0) cout << fitsig->GetParameter(2);
    
    
//     if(3.0*fitsig->GetParameter(2)<etalimit)gEtasigRange = 3.0*fitsig->GetParameter(2);//set the range to 3*/sigma
//     else gEtasigRange = etalimit;
    if(fitsig->GetParameter(0)<1.0E-6) gEtasigRange = 0.3;//if there is no peak, most of it.

    TFitResultPtr bkgresult = deta12ss->Fit(fitbg,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    if(int(bkgresult)!=4000)//if "error", try again with more:
      bkgresult = deta12ss->Fit(fitbg,"MSQ","",-gkEtaFitRange,gkEtaFitRange);    
    deta12ss->Write(Form("%sflatg_%i_bkgonly",deta12ss->GetName(),dphi));
    fitsig->FixParameter(3,0.0);//fitbg->GetParameter(0));
    rembg->SetParameter(0,fitbg->GetParameter(0));

    if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol1")==0){fitsig->FixParameter(4,fitbg->GetParameter(2));rembg->SetParameter(2,fitbg->GetParameter(2));}
    if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol2")==0){fitsig->FixParameter(5,fitbg->GetParameter(3));rembg->SetParameter(2,fitbg->GetParameter(3));}
    
    fillwithbinnr(hist,deta12ssbinc,dphi);
    Double_t ErrBg = 0.0;
    if(fitbg->GetNDF()!=0)ErrBg= fitbg->GetParError(0)*fitbg->GetChisquare()/fitbg->GetNDF();
    removeconstant(deta12ssbinc,-1.0*fitbg->GetParameter(0),ErrBg);
    
    Double_t binerr;
    etasigrange->SetBinContent(dphi,gEtasigRange);
    Double_t binc = deta12ssbinc->IntegralAndError(deta12ssbinc->FindBin(-gEtasigRange),deta12ssbinc->FindBin(gEtasigRange),binerr);
    Double_t rmsv = deta12ssbinc->GetRMS();
    yieldbc->SetBinContent(dphi,binc/deltaphi);
    yieldbc->SetBinError(dphi,binerr/deltaphi);
    hRMS->SetBinContent(dphi,rmsv/deltaphi);
    deta12ssbinc->Write(Form("%sminbgbin_%i",deta12ss->GetName(),dphi));
    
    TFitResultPtr fitresult2 = deta12ssbinc->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
    if(int(fitresult2)!=4000)//if "error", try again with more:
      fitresult2 = deta12ssbinc->Fit(fitsig,"MSQ","",-gkEtaFitRange,gkEtaFitRange);    
    deta12ssbinc->Write(Form("%sflatg_%i_full",deta12ss->GetName(),dphi));
    fitsig->ReleaseParameter(1);fitsig->ReleaseParameter(3);
    if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol1")==0){fitsig->ReleaseParameter(4);}
    if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol2")==0){fitsig->ReleaseParameter(5);}    

    TCanvas * bincanvas = new TCanvas(Form("%sCanvas_%i",deta12ss->GetName(),dphi));
    deta12ssdraw->SetTitle(Form("%s %4.2f #pi < #Delta#Phi < %4.2f #pi",title.Data(),hist->GetYaxis()->GetBinLowEdge(dphi)/TMath::Pi(),hist->GetYaxis()->GetBinUpEdge(dphi)/TMath::Pi()));
    deta12ssdraw->Reset();
    fillwithbinnr(hist,deta12ssdraw,dphi);
    deta12ssdraw->SetStats(false);    
    removeconstant(deta12ssdraw,-1.0*fitbg->GetParameter(0),ErrBg);
    deta12ssdraw->Draw("ESAME");
    fitsig->Draw("LSAME");
    bincanvas->Update();
    bincanvas->Write();
    delete bincanvas;

    yield->SetBinContent(dphi,fitsig->GetParameter(0)/deltaphi);
    yield->SetBinError(dphi,fitsig->GetParError(0)/deltaphi);
    width->SetBinContent(dphi,fitsig->GetParameter(2)/deltaphi);
    width->SetBinError(dphi,fitsig->GetParError(2)/deltaphi);
    peakpos->SetBinContent(dphi,fitsig->GetParameter(1)/deltaphi);
    peakpos->SetBinError(dphi,fitsig->GetParError(1)/deltaphi);
    background->SetBinContent(dphi,fitbg->GetParameter(0)/deltaphi);
    background->SetBinError(dphi,fitbg->GetParError(0)/deltaphi);
    if(fitsig->GetNDF()>=1.0)chisq->SetBinContent(dphi,fitsig->GetChisquare()/fitsig->GetNDF());
    prob->SetBinContent(dphi,fitsig->GetProb());
    deta12ss->Reset();
    
    
//     fitresult.~TFitResultPtr();fitresult2.~TFitResultPtr();bkgresult.~TFitResultPtr();
  }
  typedir->cd();
  yieldbc->Write();
  hRMS->Write();
  yield->Write();
  background->Write();
  width->Write();
  peakpos->Write();
  chisq->Write();
  prob->Write();
  etasigrange->Write();
  
  
  delete types;delete yield; delete background;delete width; delete peakpos; delete chisq; delete prob;delete etasigrange;
}



void extractbinyield(TDirectory* dir, TDirectory* yielddir, Double_t etalimit){
  if(!dir||!yielddir)return;
  TDirectory * bindir = resultsdirectory(yielddir,"original_binning");
  TDirectory * binrebindir = resultsdirectory(yielddir,"Rebinned_3_to_1");
  gkEtaFitRange = 1.4;
  // if(TString(dir->GetPath()).Contains("4_8_0_1"))gkEtaFitRange = 1.5;
  // if(TString(dir->GetPath()).Contains("4_8_1_2"))gkEtaFitRange = 1.5;
  // if(TString(dir->GetPath()).Contains("4_8_2_4"))gkEtaFitRange = 1.3;
  // if(TString(dir->GetPath()).Contains("8_16_0_1"))gkEtaFitRange = 1.5;
  // if(TString(dir->GetPath()).Contains("8_16_1_2"))gkEtaFitRange = 1.5;
  // if(TString(dir->GetPath()).Contains("8_16_2_4"))gkEtaFitRange = 1.2;
  // if(TString(dir->GetPath()).Contains("8_16_4_8"))gkEtaFitRange = 1.0;

  TH2D* dphideta12ss;
  if(dynamic_cast<TH2D*>(dir->Get("DPhi_1_DEta_12_SameSide"))) dphideta12ss = dynamic_cast<TH2D*>(dir->Get("DPhi_1_DEta_12_SameSide")->Clone("DPhi_1_DEta_12_SameSidec"));
  else return;
  yielddir->cd();
  fitwith(bindir,"fgauspol0/1",dphideta12ss,etalimit);
  fitwith(binrebindir,"fgauspol0/3",dphideta12ss,etalimit);
//   fitwith(bindir,"fgauspol1/1",dphideta12ss);
//   fitwith(binrebindir,"fgauspol1/3",dphideta12ss);
  //  gkEtaFitRange = 1.6;
//   fitwith(bindir,"fgauspol2/1",dphideta12ss);
//   fitwith(binrebindir,"fgauspol2/3",dphideta12ss);
  
}



TList * GetMZDirectories(TDirectory* same){
  TList * keys = same->GetListOfKeys();
  TList * directories = new TList();
  for(int i=0; i<keys->GetEntries();i++){
    if(TString(keys->At(i)->GetName()).BeginsWith("BinM(")&&TString(keys->At(i)->GetName()).Contains("Z(")){
      BinDirs* tmp = new BinDirs(same->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("META")->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("METrigger")->GetDirectory(keys->At(i)->GetName()));
      directories->Add(tmp);
    }
    
  }
  return directories;
}

TList * GetMZDirectories(BinDirs* same){
  TList * keys = same->Same()->GetListOfKeys();
  TList * directories = new TList();
  for(int i=0; i<keys->GetEntries();i++){if(TString(keys->At(i)->GetName()).BeginsWith("BinM(")&&TString(keys->At(i)->GetName()).Contains("Z(")) directories->Add(same->Same()->GetDirectory(keys->At(i)->GetName()) );}
  return directories;
}

void canvasmaker(const char* histname, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs* Bin2 = NULL;
  if(multdirlist->GetEntries()>2)Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs* Bin3 = NULL;
  if(multdirlist->GetEntries()>3)Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs* Bin4 = NULL;
  if(multdirlist->GetEntries()>4)Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs* Bin5 = NULL;
  if(multdirlist->GetEntries()>5)Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  BinDirs * Bin6 = NULL;
  if(multdirlist->GetEntries()>6) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  BinDirs * Bin7 = NULL;  
  if(multdirlist->GetEntries()>7) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(7));
  
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
    if(Bin2){
      hist = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin2->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    }
    if(Bin3){
      hist = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin3->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    }
    if(Bin4){
      hist = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin4->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    }
    if(Bin5){
      hist = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin5->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    }
    if(Bin6){
      hist = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin6->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    }
    if(Bin7){
      hist = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin7->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    }
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
    if(Bin2){
      hist = dynamic_cast<TH2D*>(Bin2->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin2->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin2->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin2->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin2->METriggerdiv()->Get(histname));
    }
    if(Bin3){
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin2->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin3->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin3->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin3->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin3->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin3->METriggerdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin3->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    }
    if(Bin4){
      hist = dynamic_cast<TH2D*>(Bin4->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin4->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin4->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin4->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin4->METriggerdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin4->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    }
    if(Bin5){
      hist = dynamic_cast<TH2D*>(Bin5->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin5->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin5->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin5->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin5->METriggerdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin5->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    }
    if(Bin6){
      hist = dynamic_cast<TH2D*>(Bin6->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin6->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin6->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin6->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin6->METriggerdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin6->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    }
    if(Bin7){
      hist = dynamic_cast<TH2D*>(Bin7->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin7->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin7->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin7->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      hist = dynamic_cast<TH2D*>(Bin7->METriggerdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);Bin7->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    }
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
  BinDirs * Bin2 =NULL;
  if(multdirlist->GetEntries()>2) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs* Bin3 = NULL;
  if(multdirlist->GetEntries()>3) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs* Bin4 = NULL;
  if(multdirlist->GetEntries()>4) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs* Bin5 = NULL;
  if(multdirlist->GetEntries()>5) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  BinDirs * Bin6 =NULL; 
  if(multdirlist->GetEntries()>6) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  BinDirs * Bin7=NULL;  
  if(multdirlist->GetEntries()>7) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(7));
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
  TH1D* histbin2;TH1D* histMETAbin2;TH1D* histMETriggerbin2;
  if(Bin2){
    histbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histMETAbin2	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETriggerbin2 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));    
  }
  //Bin 3:
  TH1D* histbin3;TH1D* histMETAbin3;TH1D* histMETriggerbin3;
  if(Bin3){
    histbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histMETAbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETriggerbin3 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  }
  //Bin 4:
  TH1D* histbin4;TH1D* histMETAbin4;TH1D* histMETriggerbin4;
  if(Bin4){
    histbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histMETAbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETriggerbin4 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  }
  //Bin 5:
   TH1D* histbin5;TH1D* histMETAbin5;TH1D* histMETriggerbin5;
  if(Bin4){
    histbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histMETAbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETriggerbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  }
  //Bin 6:
  TH1D* histbin6;TH1D* histMETAbin6;TH1D* histMETriggerbin6;
  if(Bin6){
    histbin6		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6"		,histo->GetName())));
    histMETAbin6		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6META"	,histo->GetName())));
    histMETriggerbin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6METrigger"	,histo->GetName())));
  }
  //Bin 7:
  TH1D* histbin7;TH1D* histMETAbin7;TH1D* histMETriggerbin7;
  if(Bin7){
    histbin7		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7"		,histo->GetName())));
    histMETAbin7		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7META"	,histo->GetName())));
    histMETriggerbin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7METrigger"	,histo->GetName())));
  }  
  //Doubles to hold the values until they are put into the hists.
  TH1D* histtmp;TH1D* histMETAtmp;TH1D* histMETriggertmp;
  //loop over the bins:
  for(int i=0;i<directories->GetEntries();i++){
    int Mbin = 0;
//     if((TString(directories->At(i)->GetName()).Contains("Z(-10.00)"))&&(TString(directories->At(i)->GetName()).Contains("Z(5.00)")))continue;
    //find the Multiplicity bin we are in
    for(int j = 1;j<multdirlist->GetEntries();j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
    
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
    else if(Mbin == 2&&Bin2){ histbin2->Add(hist);histMETAbin2->Add(histMETAtmp); histMETriggerbin2->Add(histMETriggertmp);}	      
    else if(Mbin == 3&&Bin3){ histbin3->Add(hist);histMETAbin3->Add(histMETAtmp); histMETriggerbin3->Add(histMETriggertmp);}
    else if(Mbin == 4&&Bin4){ histbin4->Add(hist);histMETAbin4->Add(histMETAtmp); histMETriggerbin4->Add(histMETriggertmp);}	      
    else if(Mbin == 5&&Bin5){ histbin5->Add(hist);histMETAbin5->Add(histMETAtmp); histMETriggerbin5->Add(histMETriggertmp);}
    else if(Mbin == 6&&Bin6){ histbin6->Add(hist);histMETAbin6->Add(histMETAtmp); histMETriggerbin6->Add(histMETriggertmp);}
    else if(Mbin == 7&&Bin7){ histbin7->Add(hist);histMETAbin7->Add(histMETAtmp); histMETriggerbin7->Add(histMETriggertmp);}

    
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
  if(Bin2){
    Bin2->Same()->GetDirectory("bin_stats")->cd();
    histbin2->Write(histo->GetName());
    delete histbin2;
    Bin2->META()->GetDirectory("bin_stats")->cd();
    histMETAbin2->Write(histo->GetName());
    delete histMETAbin2;
    Bin2->METrigger()->GetDirectory("bin_stats")->cd();
    histMETriggerbin2->Write(histo->GetName());
    delete histMETriggerbin2;  
  }
  if(Bin3){
    Bin3->Same()->GetDirectory("bin_stats")->cd();
    histbin3->Write(histo->GetName());
    delete histbin3;
    Bin3->META()->GetDirectory("bin_stats")->cd();
    histMETAbin3->Write(histo->GetName());
    delete histMETAbin3;
    Bin3->METrigger()->GetDirectory("bin_stats")->cd();
    histMETriggerbin3->Write(histo->GetName());
    delete histMETriggerbin3;  
  }
  if(Bin4){
    Bin4->Same()->GetDirectory("bin_stats")->cd();
    histbin4->Write(histo->GetName());
    delete histbin4;
    Bin4->META()->GetDirectory("bin_stats")->cd();
    histMETAbin4->Write(histo->GetName());
    delete histMETAbin4;
    Bin4->METrigger()->GetDirectory("bin_stats")->cd();
    histMETriggerbin4->Write(histo->GetName());
    delete histMETriggerbin4;  
  }
  if(Bin4){
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
  if(Bin6){
    Bin6->Same()->GetDirectory("bin_stats")->cd();
    histbin6->Write(histo->GetName());
    delete histbin6;
    Bin6->META()->GetDirectory("bin_stats")->cd();
    histMETAbin6->Write(histo->GetName());
    delete histMETAbin6;
    Bin6->METrigger()->GetDirectory("bin_stats")->cd();
    histMETriggerbin6->Write(histo->GetName());
    delete histMETriggerbin6;  
  }
  if(Bin7){
    Bin7->Same()->GetDirectory("bin_stats")->cd();
    histbin7->Write(histo->GetName());
    delete histbin7;
    Bin7->META()->GetDirectory("bin_stats")->cd();
    histMETAbin7->Write(histo->GetName());
    delete histMETAbin7;
    Bin7->METrigger()->GetDirectory("bin_stats")->cd();
    histMETriggerbin7->Write(histo->GetName());
    delete histMETriggerbin7;  
  }
}

void CollectHist(TH1D* histo, TList * directories, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 =NULL;
  if(multdirlist->GetEntries()>2) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs* Bin3 = NULL;
  if(multdirlist->GetEntries()>3) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs* Bin4 = NULL;
  if(multdirlist->GetEntries()>4) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs* Bin5 = NULL;
  if(multdirlist->GetEntries()>5) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  BinDirs * Bin6 =NULL; 
  if(multdirlist->GetEntries()>6) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  BinDirs * Bin7=NULL;  
  if(multdirlist->GetEntries()>7) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(7));
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
  TH1D* histbin2;TH1D* histMETAbin2;TH1D* histMETriggerbin2;
  if(Bin2){
    histbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histMETAbin2	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETriggerbin2 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));    
  }
  //Bin 3:
  TH1D* histbin3;TH1D* histMETAbin3;TH1D* histMETriggerbin3;
  if(Bin3){
    histbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histMETAbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETriggerbin3 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  }
  //Bin 4:
  TH1D* histbin4;TH1D* histMETAbin4;TH1D* histMETriggerbin4;
  if(Bin4){
    histbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histMETAbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETriggerbin4 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  }
  //Bin 5:
   TH1D* histbin5;TH1D* histMETAbin5;TH1D* histMETriggerbin5;
  if(Bin4){
    histbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histMETAbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETriggerbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  }
  //Bin 6:
  TH1D* histbin6;TH1D* histMETAbin6;TH1D* histMETriggerbin6;
  if(Bin6){
    histbin6		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6"		,histo->GetName())));
    histMETAbin6		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6META"	,histo->GetName())));
    histMETriggerbin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6METrigger"	,histo->GetName())));
  }
  //Bin 7:
  TH1D* histbin7;TH1D* histMETAbin7;TH1D* histMETriggerbin7;
  if(Bin7){
    histbin7		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7"		,histo->GetName())));
    histMETAbin7		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7META"	,histo->GetName())));
    histMETriggerbin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7METrigger"	,histo->GetName())));
  }  
  //Doubles to hold the values until they are put into the hists.

  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;
  Double_t BinContent = 0.0;Double_t BinError   = 0.0; Double_t BinContentMETA = 0.0;Double_t BinErrorMETA   = 0.0;
  Double_t BinContentMETrigger = 0.0;Double_t BinErrorMEtrigger   = 0.0; Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  Double_t BinContentbin6 = 0.0;Double_t BinErrorbin6   = 0.0;Double_t BinContentMETAbin6 = 0.0;Double_t BinErrorMETAbin6   = 0.0;Double_t BinContentMETriggerbin6 = 0.0;Double_t BinErrorMEtriggerbin6   = 0.0; 
  Double_t BinContentbin7 = 0.0;Double_t BinErrorbin7   = 0.0;Double_t BinContentMETAbin7 = 0.0;Double_t BinErrorMETAbin7   = 0.0;Double_t BinContentMETriggerbin7 = 0.0;Double_t BinErrorMEtriggerbin7   = 0.0; 

  //loop over the bins:
  for(int x=0;x<=histo->GetNbinsX()+1;x++){
    //loop over all Multiplicity-Vertex bins:
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
//       if((TString(directories->At(i)->GetName()).Contains("Z(-10.00)"))&&(TString(directories->At(i)->GetName()).Contains("Z(5.00)")))continue;

      //find the Multiplicity bin we are in
      for(int j = 1;j<multdirlist->GetEntries();j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
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
	else if(Mbin == 6){ 	BinContentbin6 += bincontl/(binerrorl*binerrorl); BinErrorbin6   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 7){ 	BinContentbin7 += bincontl/(binerrorl*binerrorl); BinErrorbin7   += 1.0/(binerrorl*binerrorl);}      }
      else{
	BinError += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorbin1   += 1.0;
	else if(Mbin == 2)	BinErrorbin2   += 1.0;
	else if(Mbin == 3)	BinErrorbin3   += 1.0;
	else if(Mbin == 4)	BinErrorbin4   += 1.0;	      
	else if(Mbin == 5)	BinErrorbin5   += 1.0;
	else if(Mbin == 6)	BinErrorbin6   += 1.0;
	else if(Mbin == 7)	BinErrorbin7   += 1.0;
      }
        if(bincontlMETA>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETA += bincontlMETA/(binerrorlMETA*binerrorlMETA);
	BinErrorMETA   += 1.0/(binerrorlMETA*binerrorlMETA);
	if(Mbin == 1){		BinContentMETAbin1 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin1   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 2){	BinContentMETAbin2 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin2   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 3){	BinContentMETAbin3 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin3   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 4){	BinContentMETAbin4 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin4   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 5){	BinContentMETAbin5 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin5   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 6){	BinContentMETAbin6 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin6   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 7){	BinContentMETAbin7 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin7   += 1.0/(binerrorlMETA*binerrorlMETA);}
	    }
      else{
	BinErrorMETA += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMETAbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMETAbin2   += 1.0;	      
	else if(Mbin == 3)	BinErrorMETAbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMETAbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMETAbin5   += 1.0;
	else if(Mbin == 6)	BinErrorMETAbin6   += 1.0;
	else if(Mbin == 7)	BinErrorMETAbin7   += 1.0;
      }
      if(bincontlMEtrigger>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETrigger += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger);
	BinErrorMEtrigger   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);
	if(Mbin == 1){		BinContentMETriggerbin1 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin1   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 2){	BinContentMETriggerbin2 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin2   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 3){	BinContentMETriggerbin3 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin3   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 4){	BinContentMETriggerbin4 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin4   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 5){	BinContentMETriggerbin5 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin5   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 6){	BinContentMETriggerbin6 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin6   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 7){	BinContentMETriggerbin7 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin7   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
      }
      else{
	BinErrorMEtrigger += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMEtriggerbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMEtriggerbin2   += 1.0;
	else if(Mbin == 3)	BinErrorMEtriggerbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMEtriggerbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMEtriggerbin5   += 1.0;
	else if(Mbin == 6)	BinErrorMEtriggerbin6   += 1.0;
	else if(Mbin == 7)	BinErrorMEtriggerbin7   += 1.0;
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
    if(BinErrorbin6>1.0e-10){		BinContentbin6 = BinContentbin6/BinErrorbin6;					BinErrorbin6 = 1.0/TMath::Sqrt(BinErrorbin6);			}
    else{				BinContentbin6=0.0;								BinErrorbin6=0.0;						}
    if(BinErrorbin7>1.0e-10){		BinContentbin7 = BinContentbin7/BinErrorbin7;					BinErrorbin7 = 1.0/TMath::Sqrt(BinErrorbin7);			}
    else{				BinContentbin7=0.0;								BinErrorbin7=0.0;						}
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
    if(BinErrorMETAbin6>1.0e-10){	BinContentMETAbin6 = BinContentMETAbin6/BinErrorMETAbin6;			BinErrorMETAbin6 = 1.0/TMath::Sqrt(BinErrorMETAbin6);		}
    else{				BinContentMETAbin6=0.0;								BinErrorMETAbin6=0.0;						}
    if(BinErrorMETAbin7>1.0e-10){	BinContentMETAbin7 = BinContentMETAbin7/BinErrorMETAbin7;			BinErrorMETAbin7 = 1.0/TMath::Sqrt(BinErrorMETAbin7);		}
    else{				BinContentMETAbin7=0.0;								BinErrorMETAbin7=0.0;						}
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
    if(BinErrorMEtriggerbin5>1.0e-10){	BinContentMETriggerbin5 = BinContentMETriggerbin5/BinErrorMEtriggerbin5;	BinErrorMEtriggerbin5 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin5);	}
    else{				BinContentMETriggerbin5=0.0;							BinErrorMEtriggerbin5=0.0;					}
    if(BinErrorMEtriggerbin6>1.0e-10){	BinContentMETriggerbin6 = BinContentMETriggerbin6/BinErrorMEtriggerbin6;	BinErrorMEtriggerbin6 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin6);	}
    else{				BinContentMETriggerbin6=0.0;							BinErrorMEtriggerbin6=0.0;					}
    if(BinErrorMEtriggerbin7>1.0e-10){	BinContentMETriggerbin7 = BinContentMETriggerbin7/BinErrorMEtriggerbin7;	BinErrorMEtriggerbin7 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin7);	}
    else{				BinContentMETriggerbin7=0.0;							BinErrorMEtriggerbin7=0.0;					}
    //Set the bin content and error in every histogram.

  if(BinContent>1.0e-10){
      hist->SetBinContent(x,BinContent);
      hist->SetBinError(x,BinError);
      histbin1->SetBinContent(x,BinContentbin1);
      histbin1->SetBinError(x,BinErrorbin1);
      if(Bin2){
	histbin2->SetBinContent(x,BinContentbin2);
	histbin2->SetBinError(x,BinErrorbin2);
      }
      if(Bin3){
	histbin3->SetBinContent(x,BinContentbin3);
	histbin3->SetBinError(x,BinErrorbin3);
      }
      if(Bin4){
	histbin4->SetBinContent(x,BinContentbin4);
	histbin4->SetBinError(x,BinErrorbin4);
      }
      if(Bin5){
	histbin5->SetBinContent(x,BinContentbin5);
	histbin5->SetBinError(x,BinErrorbin5);
      }
      if(Bin6){
	histbin6->SetBinContent(x,BinContentbin6);
	histbin6->SetBinError(x,BinErrorbin6);
      }
      if(Bin7){
	histbin7->SetBinContent(x,BinContentbin7);
	histbin7->SetBinError(x,BinErrorbin7);
      }
      histMETA->SetBinContent(x,BinContentMETA);
      histMETA->SetBinError(x,BinErrorMETA);
      histMETAbin1->SetBinContent(x,BinContentMETAbin1);
      histMETAbin1->SetBinError(x,BinErrorMETAbin1);
      if(Bin2){
	histMETAbin2->SetBinContent(x,BinContentMETAbin2);
	histMETAbin2->SetBinError(x,BinErrorMETAbin2);
      }
      if(Bin3){
	histMETAbin3->SetBinContent(x,BinContentMETAbin3);
	histMETAbin3->SetBinError(x,BinErrorMETAbin3);
      }
      if(Bin4){
	histMETAbin4->SetBinContent(x,BinContentMETAbin4);
	histMETAbin4->SetBinError(x,BinErrorMETAbin4);
      }
      if(Bin5){
	histMETAbin5->SetBinContent(x,BinContentMETAbin5);
	histMETAbin5->SetBinError(x,BinErrorMETAbin5);
      }
      if(Bin6){
	histMETAbin6->SetBinContent(x,BinContentMETAbin6);
	histMETAbin6->SetBinError(x,BinErrorMETAbin6);
      }
      if(Bin7){
	histMETAbin7->SetBinContent(x,BinContentMETAbin7);
	histMETAbin7->SetBinError(x,BinErrorMETAbin7);
      }
      histMETrigger->SetBinContent(x,BinContentMETrigger);
      histMETrigger->SetBinError(x,BinErrorMEtrigger);
      histMETriggerbin1->SetBinContent(x,BinContentMETriggerbin1);
      histMETriggerbin1->SetBinError(x,BinErrorMEtriggerbin1);
      if(Bin2){
      histMETriggerbin2->SetBinContent(x,BinContentMETriggerbin2);
      histMETriggerbin2->SetBinError(x,BinErrorMEtriggerbin2);
      }
      if(Bin3){
      histMETriggerbin3->SetBinContent(x,BinContentMETriggerbin3);
      histMETriggerbin3->SetBinError(x,BinErrorMEtriggerbin3);
      }
      if(Bin4){
      histMETriggerbin4->SetBinContent(x,BinContentMETriggerbin4);
      histMETriggerbin4->SetBinError(x,BinErrorMEtriggerbin4);
      }
      if(Bin5){
      histMETriggerbin5->SetBinContent(x,BinContentMETriggerbin5);
      histMETriggerbin5->SetBinError(x,BinErrorMEtriggerbin5);
      }
      if(Bin6){      
	histMETriggerbin6->SetBinContent(x,BinContentMETriggerbin6);
	histMETriggerbin6->SetBinError(x,BinErrorMEtriggerbin6);
      }
      if(Bin7){      
	histMETriggerbin7->SetBinContent(x,BinContentMETriggerbin7);
	histMETriggerbin7->SetBinError(x,BinErrorMEtriggerbin7);
      }      
    }
    //Reset all:
    BinContent = 0.0;BinError = 0.0;BinContentMETA = 0.0;BinErrorMETA = 0.0;BinContentMETrigger = 0.0;BinErrorMEtrigger = 0.0;
    BinContentbin1 = 0.0;BinErrorbin1 = 0.0;BinContentMETAbin1 = 0.0;BinErrorMETAbin1 = 0.0;BinContentMETriggerbin1 = 0.0;BinErrorMEtriggerbin1 = 0.0;
    BinContentbin2 = 0.0;BinErrorbin2 = 0.0;BinContentMETAbin2 = 0.0;BinErrorMETAbin2 = 0.0;BinContentMETriggerbin2 = 0.0;BinErrorMEtriggerbin2 = 0.0;
    BinContentbin3 = 0.0;BinErrorbin3 = 0.0;BinContentMETAbin3 = 0.0;BinErrorMETAbin3 = 0.0;BinContentMETriggerbin3 = 0.0;BinErrorMEtriggerbin3 = 0.0;
    BinContentbin4 = 0.0;BinErrorbin4 = 0.0;BinContentMETAbin4 = 0.0;BinErrorMETAbin4 = 0.0;BinContentMETriggerbin4 = 0.0;BinErrorMEtriggerbin4 = 0.0;
    BinContentbin5 = 0.0;BinErrorbin5 = 0.0;BinContentMETAbin5 = 0.0;BinErrorMETAbin5 = 0.0;BinContentMETriggerbin5 = 0.0;BinErrorMEtriggerbin5 = 0.0;	
    BinContentbin6 = 0.0;BinErrorbin6 = 0.0;BinContentMETAbin6 = 0.0;BinErrorMETAbin6 = 0.0;BinContentMETriggerbin6 = 0.0;BinErrorMEtriggerbin6 = 0.0;	
    BinContentbin7 = 0.0;BinErrorbin7 = 0.0;BinContentMETAbin7 = 0.0;BinErrorMETAbin7 = 0.0;BinContentMETriggerbin7 = 0.0;BinErrorMEtriggerbin7 = 0.0;	
  }//end binloop
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
  if(Bin2){
    Bin2->Samediv()->cd();
    histbin2->Write(histo->GetName());
    delete histbin2;
    Bin2->METAdiv()->cd();
    histMETAbin2->Write(histo->GetName());
    delete histMETAbin2;
    Bin2->METriggerdiv()->cd();
    histMETriggerbin2->Write(histo->GetName());
    delete histMETriggerbin2;  
  }
  if(Bin4){
    Bin3->Samediv()->cd();
    histbin3->Write(histo->GetName());
    delete histbin3;
    Bin3->METAdiv()->cd();
    histMETAbin3->Write(histo->GetName());
    delete histMETAbin3;
    Bin3->METriggerdiv()->cd();
    histMETriggerbin3->Write(histo->GetName());
    delete histMETriggerbin3;  
  }
  if(Bin4){
    Bin4->Samediv()->cd();
    histbin4->Write(histo->GetName());
    delete histbin4;
    Bin4->METAdiv()->cd();
    histMETAbin4->Write(histo->GetName());
    delete histMETAbin4;
    Bin4->METriggerdiv()->cd();
    histMETriggerbin4->Write(histo->GetName());
    delete histMETriggerbin4;
  }
  if(Bin5){
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
  if(Bin6){
    Bin6->Samediv()->cd();
    histbin6->Write(histo->GetName());
    delete histbin6;
    Bin6->METAdiv()->cd();
    histMETAbin6->Write(histo->GetName());
    delete histMETAbin6;
    Bin6->METriggerdiv()->cd();
    histMETriggerbin6->Write(histo->GetName());
    delete histMETriggerbin6; 
  }
  if(Bin7){
    Bin7->Samediv()->cd();
    histbin7->Write(histo->GetName());
    delete histbin7;
    Bin7->METAdiv()->cd();
    histMETAbin7->Write(histo->GetName());
    delete histMETAbin7;
    Bin7->METriggerdiv()->cd();
    histMETriggerbin7->Write(histo->GetName());
    delete histMETriggerbin7; 
  }
  
}
void CollectHist(TH2D* histo, TList * directories, TObjArray* multdirlist,bool pearsonserrors = false,bool collectdivfirst = false){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 =NULL;
  if(multdirlist->GetEntries()>2) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs* Bin3 = NULL;
  if(multdirlist->GetEntries()>3) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs* Bin4 = NULL;
  if(multdirlist->GetEntries()>4) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs* Bin5 = NULL;
  if(multdirlist->GetEntries()>5) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  BinDirs * Bin6 =NULL; 
  if(multdirlist->GetEntries()>6) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  BinDirs * Bin7=NULL;  
  if(multdirlist->GetEntries()>7) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(7));
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
    
  //all bins:
  TH2D* hist 			= dynamic_cast<TH2D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH2D* hista 			= dynamic_cast<TH2D*>( histo->Clone(Form("%sav"			,histo->GetName())));
  TH2D* histMETA 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH2D* histaMETA 		= dynamic_cast<TH2D*>( histo->Clone(Form("%savMETA"		,histo->GetName())));
  TH2D* histMETrigger 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  TH2D* histaMETrigger 		= dynamic_cast<TH2D*>( histo->Clone(Form("%savMETrigger"	,histo->GetName())));
  //Bin 1:
  TH2D* histbin1 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH2D* histMETAbin1 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH2D* histMETriggerbin1 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  //Bin 2:
  TH2D* histbin2;TH2D* histMETAbin2;TH2D* histMETriggerbin2;
  if(Bin2){
    histbin2		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histMETAbin2	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETriggerbin2 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));    
  }
  //Bin 3:
  TH2D* histbin3;TH2D* histMETAbin3;TH2D* histMETriggerbin3;
  if(Bin3){
    histbin3		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histMETAbin3	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETriggerbin3 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  }
  //Bin 4:
  TH2D* histbin4;TH2D* histMETAbin4;TH2D* histMETriggerbin4;
  if(Bin4){
    histbin4		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histMETAbin4	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETriggerbin4 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  }
  //Bin 5:
   TH2D* histbin5;TH2D* histMETAbin5;TH2D* histMETriggerbin5;
  if(Bin4){
    histbin5		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histMETAbin5	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETriggerbin5	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  }
  //Bin 6:
  TH2D* histbin6;TH2D* histMETAbin6;TH2D* histMETriggerbin6;
  if(Bin6){
    histbin6		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin6"		,histo->GetName())));
    histMETAbin6	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin6META"	,histo->GetName())));
    histMETriggerbin6	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin6METrigger"	,histo->GetName())));
  }
  //Bin 7:
  TH2D* histbin7;TH2D* histMETAbin7;TH2D* histMETriggerbin7;
  if(Bin7){
    histbin7		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin7"		,histo->GetName())));
    histMETAbin7	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin7META"	,histo->GetName())));
    histMETriggerbin7	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin7METrigger"	,histo->GetName())));
  }  
  //Doubles to hold the values until they are put into the hists.
  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;
  Double_t BinContent = 0.0;Double_t BinError   = 0.0; Double_t BinContentMETA = 0.0;Double_t BinErrorMETA   = 0.0;
  Double_t BinContentMETrigger = 0.0;Double_t BinErrorMEtrigger   = 0.0; Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  Double_t BinContentbin6 = 0.0;Double_t BinErrorbin6   = 0.0;Double_t BinContentMETAbin6 = 0.0;Double_t BinErrorMETAbin6   = 0.0;Double_t BinContentMETriggerbin6 = 0.0;Double_t BinErrorMEtriggerbin6   = 0.0; 
  Double_t BinContentbin7 = 0.0;Double_t BinErrorbin7   = 0.0;Double_t BinContentMETAbin7 = 0.0;Double_t BinErrorMETAbin7   = 0.0;Double_t BinContentMETriggerbin7 = 0.0;Double_t BinErrorMEtriggerbin7   = 0.0; 

  if(collectdivfirst){
    //in this case, first add all same_event/mixed_event, then scale by #triggers.
    Double_t scalingfactor=0.0;Double_t scalingfactorbin1=0.0;Double_t scalingfactorbin2=0.0;Double_t scalingfactorbin3=0.0;Double_t scalingfactorbin4=0.0;Double_t scalingfactorbin5=0.0;Double_t scalingfactorbin6=0.0;Double_t scalingfactorbin7=0.0;
    Double_t scalingfactorMETA=0.0;Double_t scalingfactorMETAbin1=0.0;Double_t scalingfactorMETAbin2=0.0;Double_t scalingfactorMETAbin3=0.0;Double_t scalingfactorMETAbin4=0.0;Double_t scalingfactorMETAbin5=0.0;Double_t scalingfactorMETAbin6=0.0;Double_t scalingfactorMETAbin7=0.0;
    Double_t scalingfactorMETrigger=0.0;Double_t scalingfactorMETriggerbin1=0.0;Double_t scalingfactorMETriggerbin2=0.0;Double_t scalingfactorMETriggerbin3=0.0;Double_t scalingfactorMETriggerbin4=0.0;Double_t scalingfactorMETriggerbin5=0.0;Double_t scalingfactorMETriggerbin6=0.0;Double_t scalingfactorMETriggerbin7=0.0;
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
//     if((TString(directories->At(i)->GetName()).Contains("Z(-10.00)"))&&(TString(directories->At(i)->GetName()).Contains("Z(5.00)")))continue;      
      //find the Multiplicity bin we are in
      for(int j = 1;j<multdirlist->GetEntries();j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/mixed_event",directories->At(i)->GetName()))->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.
      if(!dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/mixed_event",directories->At(i)->GetName()))->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.
      if(!dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/mixed_event",directories->At(i)->GetName()))->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.      
      TH2D* histnow = dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get(histo->GetName())->Clone("nowhist"));
      TH2D* histmixed = dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/mixed_event",directories->At(i)->GetName()))->Get(histo->GetName()));
      histnow->Divide(histmixed);
      hist->Add(histnow);
      if(Mbin == 1){histbin1->Add(histnow);scalingfactorbin1+= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 2){histbin2->Add(histnow);scalingfactorbin2+= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 3){histbin3->Add(histnow);scalingfactorbin3+= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 4){histbin4->Add(histnow);scalingfactorbin4+= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 5){histbin5->Add(histnow);scalingfactorbin5+= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 6){histbin6->Add(histnow);scalingfactorbin6+= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 7){histbin7->Add(histnow);scalingfactorbin7+= dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      scalingfactor += dynamic_cast<TH1D*>(All->Same()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();
      
      TH2D* histnowMETA = dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get(histo->GetName())->Clone("nowhistMETA"));
      histnowMETA->Divide(histmixed);
      histMETA->Add(histnowMETA);
      if(Mbin == 1){histMETAbin1->Add(histnowMETA);scalingfactorMETAbin1+= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 2){histMETAbin2->Add(histnowMETA);scalingfactorMETAbin2+= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 3){histMETAbin3->Add(histnowMETA);scalingfactorMETAbin3+= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 4){histMETAbin4->Add(histnowMETA);scalingfactorMETAbin4+= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 5){histMETAbin5->Add(histnowMETA);scalingfactorMETAbin5+= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 6){histMETAbin6->Add(histnowMETA);scalingfactorMETAbin6+= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 7){histMETAbin7->Add(histnowMETA);scalingfactorMETAbin7+= dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}      
      scalingfactorMETA += dynamic_cast<TH1D*>(All->META()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();
      
      TH2D* histnowMETrigger = dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get(histo->GetName())->Clone("nowhistMETrigger"));
      histnowMETrigger->Divide(histmixed);
      histMETrigger->Add(histnowMETrigger);
      if(Mbin == 1){histMETriggerbin1->Add(histnowMETrigger);scalingfactorMETriggerbin1+= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 2){histMETriggerbin2->Add(histnowMETrigger);scalingfactorMETriggerbin2+= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 3){histMETriggerbin3->Add(histnowMETrigger);scalingfactorMETriggerbin3+= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 4){histMETriggerbin4->Add(histnowMETrigger);scalingfactorMETriggerbin4+= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 5){histMETriggerbin5->Add(histnowMETrigger);scalingfactorMETriggerbin5+= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 6){histMETriggerbin6->Add(histnowMETrigger);scalingfactorMETriggerbin6+= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}
      else if(Mbin == 7){histMETriggerbin7->Add(histnowMETrigger);scalingfactorMETriggerbin7+= dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();}      
      scalingfactorMETrigger += dynamic_cast<TH1D*>(All->METrigger()->GetDirectory(Form("%s/same_event",directories->At(i)->GetName()))->Get("number_of_triggers"))->Integral();
    }


    if(!scalingfactor <1.0e-10)hist->Scale(1.0/scalingfactor);
    if(!scalingfactorbin1 <1.0e-10)histbin1->Scale(1.0/scalingfactorbin1);
    if(!scalingfactorbin2 <1.0e-10)histbin2->Scale(1.0/scalingfactorbin2);
    if(!scalingfactorbin3 <1.0e-10)histbin3->Scale(1.0/scalingfactorbin3);
    if(!scalingfactorbin4 <1.0e-10)histbin4->Scale(1.0/scalingfactorbin4);
    if(!scalingfactorbin5 <1.0e-10)histbin5->Scale(1.0/scalingfactorbin5);
    if(!scalingfactorbin6 <1.0e-10)histbin6->Scale(1.0/scalingfactorbin6);
    if(!scalingfactorbin7 <1.0e-10)histbin7->Scale(1.0/scalingfactorbin7);

    if(!scalingfactorMETA <1.0e-10)histMETA->Scale(1.0/scalingfactorMETA);
    if(!scalingfactorMETAbin1 <1.0e-10)histMETAbin1->Scale(1.0/scalingfactorMETAbin1);
    if(!scalingfactorMETAbin2 <1.0e-10)histMETAbin2->Scale(1.0/scalingfactorMETAbin2);
    if(!scalingfactorMETAbin3 <1.0e-10)histMETAbin3->Scale(1.0/scalingfactorMETAbin3);
    if(!scalingfactorMETAbin4 <1.0e-10)histMETAbin4->Scale(1.0/scalingfactorMETAbin4);
    if(!scalingfactorMETAbin5 <1.0e-10)histMETAbin5->Scale(1.0/scalingfactorMETAbin5);
    if(!scalingfactorMETAbin6 <1.0e-10)histMETAbin6->Scale(1.0/scalingfactorMETAbin6);
    if(!scalingfactorMETAbin7 <1.0e-10)histMETAbin7->Scale(1.0/scalingfactorMETAbin7);

    if(!scalingfactorMETrigger <1.0e-10)histMETrigger->Scale(1.0/scalingfactorMETrigger);
    if(!scalingfactorMETriggerbin1 <1.0e-10)histMETriggerbin1->Scale(1.0/scalingfactorMETriggerbin1);
    if(!scalingfactorMETriggerbin2 <1.0e-10)histMETriggerbin2->Scale(1.0/scalingfactorMETriggerbin2);
    if(!scalingfactorMETriggerbin3 <1.0e-10)histMETriggerbin3->Scale(1.0/scalingfactorMETriggerbin3);
    if(!scalingfactorMETriggerbin4 <1.0e-10)histMETriggerbin4->Scale(1.0/scalingfactorMETriggerbin4);
    if(!scalingfactorMETriggerbin5 <1.0e-10)histMETriggerbin5->Scale(1.0/scalingfactorMETriggerbin5);
    if(!scalingfactorMETriggerbin6 <1.0e-10)histMETriggerbin6->Scale(1.0/scalingfactorMETriggerbin6);
    if(!scalingfactorMETriggerbin7 <1.0e-10)histMETriggerbin7->Scale(1.0/scalingfactorMETriggerbin7);
  }
  
  else{
    //loop over the bins:
    for(int x=0;x<=histo->GetNbinsX()+1;x++){for(int y=0;y<=histo->GetNbinsY();y++){
      Double_t averageinbin =0.0;
      Double_t averageinbinMETA =0.0;
      Double_t averageinbinMETrigger =0.0;

      if(pearsonserrors){
	//find the average of the bin, and use it as an estimate of the true value.
	double norm = 0.0;
	for(int i=0;i<directories->GetEntries();i++){
	  if(!dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))){continue;}
	  averageinbin += dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
	  averageinbinMETA += dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
	  averageinbinMETrigger += dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
	  norm +=1.0;
	}
	averageinbin /=norm;
	averageinbinMETA /=norm;
	averageinbinMETrigger /=norm;
      }
      hista->SetBinContent(x,y,averageinbin);
      histaMETA->SetBinContent(x,y,averageinbinMETA);
      histaMETrigger->SetBinContent(x,y,averageinbinMETrigger);

      //loop over all Multiplicity-Vertex bins:
      for(int i=0;i<directories->GetEntries();i++){
	int Mbin = 0;
	//find the Multiplicity bin we are in
	for(int j = 1;j<multdirlist->GetEntries();j++)if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
	//test if the relevant histogram exists in this bin:
	if(!dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName())))
	{
  // 	cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
	  continue;}
	//extract bin content and error in the relevant bin:
	bincontl 		= dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
	binerrorl 		= dynamic_cast<TH2D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y);	
	bincontlMETA 		= dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
	binerrorlMETA 		= dynamic_cast<TH2D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y);	      
	bincontlMEtrigger 	= dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y);
	binerrorlMETrigger 	= dynamic_cast<TH2D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y); 

	if(pearsonserrors&& averageinbin >1.0e-10 && averageinbinMETA>1.0e-10&&averageinbinMETrigger>1.0e-10){
	  binerrorl = (bincontl-averageinbin)*(bincontl-averageinbin)/averageinbin;
	  binerrorlMETA = (bincontlMETA-averageinbinMETA)*(bincontlMETA-averageinbinMETA)/averageinbinMETA;
	  binerrorlMETrigger = (bincontlMEtrigger-averageinbinMETrigger)*(bincontlMEtrigger-averageinbinMETrigger)/averageinbinMETrigger;
	}    
	else if(pearsonserrors){
	  binerrorl = 0.0;
	  binerrorlMETA = 0.0;
	  binerrorlMETrigger = 0.0;
	}
	
	if(binerrorl>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	  BinContent += bincontl/(binerrorl*binerrorl);
	  BinError   += 2.0/(binerrorl*binerrorl);
	  if(Mbin == 1){		BinContentbin1 += bincontl/(binerrorl*binerrorl); BinErrorbin1   += 1.0/(binerrorl*binerrorl);}
	  else if(Mbin == 2){	BinContentbin2 += bincontl/(binerrorl*binerrorl); BinErrorbin2   += 1.0/(binerrorl*binerrorl);}	      
	  else if(Mbin == 3){	BinContentbin3 += bincontl/(binerrorl*binerrorl); BinErrorbin3   += 1.0/(binerrorl*binerrorl);}
	  else if(Mbin == 4){	BinContentbin4 += bincontl/(binerrorl*binerrorl); BinErrorbin4   += 1.0/(binerrorl*binerrorl);}	      
	  else if(Mbin == 5){ 	BinContentbin5 += bincontl/(binerrorl*binerrorl); BinErrorbin5   += 1.0/(binerrorl*binerrorl);}
	  else if(Mbin == 6){ 	BinContentbin6 += bincontl/(binerrorl*binerrorl); BinErrorbin6   += 1.0/(binerrorl*binerrorl);}
	  else if(Mbin == 7){ 	BinContentbin7 += bincontl/(binerrorl*binerrorl); BinErrorbin7   += 1.0/(binerrorl*binerrorl);}
	}
  //       else{
  // 	BinError += 1.0;//??? Still unshure about this.
  // 	if(Mbin == 1)		BinErrorbin1   += 1.0;
  // 	else if(Mbin == 2)	BinErrorbin2   += 1.0;
  // 	else if(Mbin == 3)	BinErrorbin3   += 1.0;
  // 	else if(Mbin == 4)	BinErrorbin4   += 1.0;	      
  // 	else if(Mbin == 5)	BinErrorbin5   += 1.0;
  // 	else if(Mbin == 6)	BinErrorbin6   += 1.0;
  // 	else if(Mbin == 7)	BinErrorbin7   += 1.0;
  //       }
	if(binerrorlMETA>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	  BinContentMETA += bincontlMETA/(binerrorlMETA*binerrorlMETA);
	  BinErrorMETA   += 1.0/(binerrorlMETA*binerrorlMETA);
	  if(Mbin == 1){		BinContentMETAbin1 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin1   += 1.0/(binerrorlMETA*binerrorlMETA);}
	  else if(Mbin == 2){	BinContentMETAbin2 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin2   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	  else if(Mbin == 3){	BinContentMETAbin3 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin3   += 1.0/(binerrorlMETA*binerrorlMETA);}
	  else if(Mbin == 4){	BinContentMETAbin4 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin4   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	  else if(Mbin == 5){	BinContentMETAbin5 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin5   += 1.0/(binerrorlMETA*binerrorlMETA);}
	  else if(Mbin == 6){	BinContentMETAbin6 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin6   += 1.0/(binerrorlMETA*binerrorlMETA);}
	  else if(Mbin == 7){	BinContentMETAbin7 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin7   += 1.0/(binerrorlMETA*binerrorlMETA);}
	      }
  //       else{
  // 	BinErrorMETA += 1.0;//??? Still unshure about this.
  // 	if(Mbin == 1)		BinErrorMETAbin1   += 1.0;
  // 	else if(Mbin == 2)	BinErrorMETAbin2   += 1.0;	      
  // 	else if(Mbin == 3)	BinErrorMETAbin3   += 1.0;
  // 	else if(Mbin == 4)	BinErrorMETAbin4   += 1.0;
  // 	else if(Mbin == 5)	BinErrorMETAbin5   += 1.0;
  // 	else if(Mbin == 6)	BinErrorMETAbin6   += 1.0;
  // 	else if(Mbin == 7)	BinErrorMETAbin7   += 1.0;
  //       }
	if(binerrorlMETrigger>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	  BinContentMETrigger += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger);
	  BinErrorMEtrigger   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);
	  if(Mbin == 1){		BinContentMETriggerbin1 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin1   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	  else if(Mbin == 2){	BinContentMETriggerbin2 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin2   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	  else if(Mbin == 3){	BinContentMETriggerbin3 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin3   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	  else if(Mbin == 4){	BinContentMETriggerbin4 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin4   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	  else if(Mbin == 5){	BinContentMETriggerbin5 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin5   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	  else if(Mbin == 6){	BinContentMETriggerbin6 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin6   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	  else if(Mbin == 7){	BinContentMETriggerbin7 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin7   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	}
  //       else{
  // 	BinErrorMEtrigger += 1.0;//??? Still unshure about this.
  // 	if(Mbin == 1)		BinErrorMEtriggerbin1   += 1.0;
  // 	else if(Mbin == 2)	BinErrorMEtriggerbin2   += 1.0;
  // 	else if(Mbin == 3)	BinErrorMEtriggerbin3   += 1.0;
  // 	else if(Mbin == 4)	BinErrorMEtriggerbin4   += 1.0;
  // 	else if(Mbin == 5)	BinErrorMEtriggerbin5   += 1.0;
  // 	else if(Mbin == 6)	BinErrorMEtriggerbin6   += 1.0;
  // 	else if(Mbin == 7)	BinErrorMEtriggerbin7   += 1.0;
  //       }
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
      if(BinErrorbin6>1.0e-10){		BinContentbin6 = BinContentbin6/BinErrorbin6;					BinErrorbin6 = 1.0/TMath::Sqrt(BinErrorbin6);			}
      else{				BinContentbin6=0.0;								BinErrorbin6=0.0;						}
      if(BinErrorbin7>1.0e-10){		BinContentbin7 = BinContentbin7/BinErrorbin7;					BinErrorbin7 = 1.0/TMath::Sqrt(BinErrorbin7);			}
      else{				BinContentbin7=0.0;								BinErrorbin7=0.0;						}
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
      if(BinErrorMETAbin6>1.0e-10){	BinContentMETAbin6 = BinContentMETAbin6/BinErrorMETAbin6;			BinErrorMETAbin6 = 1.0/TMath::Sqrt(BinErrorMETAbin6);		}
      else{				BinContentMETAbin6=0.0;								BinErrorMETAbin6=0.0;						}
      if(BinErrorMETAbin7>1.0e-10){	BinContentMETAbin7 = BinContentMETAbin7/BinErrorMETAbin7;			BinErrorMETAbin7 = 1.0/TMath::Sqrt(BinErrorMETAbin7);		}
      else{				BinContentMETAbin7=0.0;								BinErrorMETAbin7=0.0;						}
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
      if(!BinErrorMEtriggerbin6>1.0e-10){	BinContentMETriggerbin6 = BinContentMETriggerbin6/BinErrorMEtriggerbin6;	BinErrorMEtriggerbin6 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin6);	}
      else{				BinContentMETriggerbin6=0.0;							BinErrorMEtriggerbin6=0.0;					}
      if(!BinErrorMEtriggerbin7>1.0e-10){	BinContentMETriggerbin7 = BinContentMETriggerbin7/BinErrorMEtriggerbin7;	BinErrorMEtriggerbin7 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin7);	}
      else{				BinContentMETriggerbin7=0.0;							BinErrorMEtriggerbin7=0.0;					}
      //Set the bin content and error in every histogram.
      if(BinContent>1.0e-10){
	hist->SetBinContent(x,y,BinContent);
	hist->SetBinError(x,y,BinError);
	histbin1->SetBinContent(x,y,BinContentbin1);
	histbin1->SetBinError(x,y,BinErrorbin1);
	if(Bin2){
	  histbin2->SetBinContent(x,y,BinContentbin2);
	  histbin2->SetBinError(x,y,BinErrorbin2);
	}
	if(Bin3){
	  histbin3->SetBinContent(x,y,BinContentbin3);
	  histbin3->SetBinError(x,y,BinErrorbin3);
	}
	if(Bin4){
	  histbin4->SetBinContent(x,y,BinContentbin4);
	  histbin4->SetBinError(x,y,BinErrorbin4);
	}
	if(Bin5){
	  histbin5->SetBinContent(x,y,BinContentbin5);
	  histbin5->SetBinError(x,y,BinErrorbin5);
	}
	if(Bin6){
	  histbin6->SetBinContent(x,y,BinContentbin6);
	  histbin6->SetBinError(x,y,BinErrorbin6);
	}
	if(Bin7){
	  histbin7->SetBinContent(x,y,BinContentbin7);
	  histbin7->SetBinError(x,y,BinErrorbin7);
	}
	histMETA->SetBinContent(x,y,BinContentMETA);
	histMETA->SetBinError(x,y,BinErrorMETA);
	histMETAbin1->SetBinContent(x,y,BinContentMETAbin1);
	histMETAbin1->SetBinError(x,y,BinErrorMETAbin1);
	if(Bin2){
	  histMETAbin2->SetBinContent(x,y,BinContentMETAbin2);
	  histMETAbin2->SetBinError(x,y,BinErrorMETAbin2);
	}
	if(Bin3){
	  histMETAbin3->SetBinContent(x,y,BinContentMETAbin3);
	  histMETAbin3->SetBinError(x,y,BinErrorMETAbin3);
	}
	if(Bin4){
	  histMETAbin4->SetBinContent(x,y,BinContentMETAbin4);
	  histMETAbin4->SetBinError(x,y,BinErrorMETAbin4);
	}
	if(Bin5){
	  histMETAbin5->SetBinContent(x,y,BinContentMETAbin5);
	  histMETAbin5->SetBinError(x,y,BinErrorMETAbin5);
	}
	if(Bin6){
	  histMETAbin6->SetBinContent(x,y,BinContentMETAbin6);
	  histMETAbin6->SetBinError(x,y,BinErrorMETAbin6);
	}
	if(Bin7){
	  histMETAbin7->SetBinContent(x,y,BinContentMETAbin7);
	  histMETAbin7->SetBinError(x,y,BinErrorMETAbin7);
	}
	histMETrigger->SetBinContent(x,y,BinContentMETrigger);
	histMETrigger->SetBinError(x,y,BinErrorMEtrigger);
	histMETriggerbin1->SetBinContent(x,y,BinContentMETriggerbin1);
	histMETriggerbin1->SetBinError(x,y,BinErrorMEtriggerbin1);
	if(Bin2){
	  histMETriggerbin2->SetBinContent(x,y,BinContentMETriggerbin2);
	  histMETriggerbin2->SetBinError(x,y,BinErrorMEtriggerbin2);
	}
	if(Bin3){
	  histMETriggerbin3->SetBinContent(x,y,BinContentMETriggerbin3);
	  histMETriggerbin3->SetBinError(x,y,BinErrorMEtriggerbin3);
	}
	if(Bin4){
	  histMETriggerbin4->SetBinContent(x,y,BinContentMETriggerbin4);
	  histMETriggerbin4->SetBinError(x,y,BinErrorMEtriggerbin4);
	}
	if(Bin5){
	  histMETriggerbin5->SetBinContent(x,y,BinContentMETriggerbin5);
	  histMETriggerbin5->SetBinError(x,y,BinErrorMEtriggerbin5);
	}
	if(Bin6){      
	  histMETriggerbin6->SetBinContent(x,y,BinContentMETriggerbin6);
	  histMETriggerbin6->SetBinError(x,y,BinErrorMEtriggerbin6);
	}
	if(Bin7){      
	  histMETriggerbin7->SetBinContent(x,y,BinContentMETriggerbin7);
	  histMETriggerbin7->SetBinError(x,y,BinErrorMEtriggerbin7);
	}      
      }
      //Reset all:
      BinContent = 0.0;BinError = 0.0;BinContentMETA = 0.0;BinErrorMETA = 0.0;BinContentMETrigger = 0.0;BinErrorMEtrigger = 0.0;
      BinContentbin1 = 0.0;BinErrorbin1 = 0.0;BinContentMETAbin1 = 0.0;BinErrorMETAbin1 = 0.0;BinContentMETriggerbin1 = 0.0;BinErrorMEtriggerbin1 = 0.0;
      BinContentbin2 = 0.0;BinErrorbin2 = 0.0;BinContentMETAbin2 = 0.0;BinErrorMETAbin2 = 0.0;BinContentMETriggerbin2 = 0.0;BinErrorMEtriggerbin2 = 0.0;
      BinContentbin3 = 0.0;BinErrorbin3 = 0.0;BinContentMETAbin3 = 0.0;BinErrorMETAbin3 = 0.0;BinContentMETriggerbin3 = 0.0;BinErrorMEtriggerbin3 = 0.0;
      BinContentbin4 = 0.0;BinErrorbin4 = 0.0;BinContentMETAbin4 = 0.0;BinErrorMETAbin4 = 0.0;BinContentMETriggerbin4 = 0.0;BinErrorMEtriggerbin4 = 0.0;
      BinContentbin5 = 0.0;BinErrorbin5 = 0.0;BinContentMETAbin5 = 0.0;BinErrorMETAbin5 = 0.0;BinContentMETriggerbin5 = 0.0;BinErrorMEtriggerbin5 = 0.0;	
      BinContentbin6 = 0.0;BinErrorbin6 = 0.0;BinContentMETAbin6 = 0.0;BinErrorMETAbin6 = 0.0;BinContentMETriggerbin6 = 0.0;BinErrorMEtriggerbin6 = 0.0;	
      BinContentbin7 = 0.0;BinErrorbin7 = 0.0;BinContentMETAbin7 = 0.0;BinErrorMETAbin7 = 0.0;BinContentMETriggerbin7 = 0.0;BinErrorMEtriggerbin7 = 0.0;	
    }}//end binloop
  }
  //save the histograms in the relevant directories:
  All->Samediv()->cd();
  hist->Write(histo->GetName());
  hista->Write();
  delete hist;delete hista;
  All->METAdiv()->cd();
  histMETA->Write(histo->GetName());
  histaMETA->Write();
  delete histMETA; delete histaMETA;
  All->METriggerdiv()->cd();
  histMETrigger->Write(histo->GetName());
  histaMETrigger->Write();
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
  if(Bin2){
    Bin2->Samediv()->cd();
    histbin2->Write(histo->GetName());
    delete histbin2;
    Bin2->METAdiv()->cd();
    histMETAbin2->Write(histo->GetName());
    delete histMETAbin2;
    Bin2->METriggerdiv()->cd();
    histMETriggerbin2->Write(histo->GetName());
    delete histMETriggerbin2; 
  }
  if(Bin3){
    Bin3->Samediv()->cd();
    histbin3->Write(histo->GetName());
    delete histbin3;
    Bin3->METAdiv()->cd();
    histMETAbin3->Write(histo->GetName());
    delete histMETAbin3;
    Bin3->METriggerdiv()->cd();
    histMETriggerbin3->Write(histo->GetName());
    delete histMETriggerbin3;  
  }
  if(Bin4){
    Bin4->Samediv()->cd();
    histbin4->Write(histo->GetName());
    delete histbin4;
    Bin4->METAdiv()->cd();
    histMETAbin4->Write(histo->GetName());
    delete histMETAbin4;
    Bin4->METriggerdiv()->cd();
    histMETriggerbin4->Write(histo->GetName());
    delete histMETriggerbin4;  
  }
  if(Bin5){
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
  if(Bin6){
    Bin6->Samediv()->cd();
    histbin6->Write(histo->GetName());
    delete histbin6;
    Bin6->METAdiv()->cd();
    histMETAbin6->Write(histo->GetName());
    delete histMETAbin6;
    Bin6->METriggerdiv()->cd();
    histMETriggerbin6->Write(histo->GetName());
    delete histMETriggerbin6; 
  }
  if(Bin7){
    Bin7->Samediv()->cd();
    histbin7->Write(histo->GetName());
    delete histbin7;
    Bin7->METAdiv()->cd();
    histMETAbin7->Write(histo->GetName());
    delete histMETAbin7;
    Bin7->METriggerdiv()->cd();
    histMETriggerbin7->Write(histo->GetName());
    delete histMETriggerbin7; 
  }
}
void CollectHist(TH3D* histo, TList * directories, TObjArray* multdirlist){
  BinDirs * All  = dynamic_cast<BinDirs*>(multdirlist->At(0));
  BinDirs * Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(1));
  BinDirs * Bin2 =NULL;
  if(multdirlist->GetEntries()>2) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(2));
  BinDirs* Bin3 = NULL;
  if(multdirlist->GetEntries()>3) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(3));
  BinDirs* Bin4 = NULL;
  if(multdirlist->GetEntries()>4) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(4));
  BinDirs* Bin5 = NULL;
  if(multdirlist->GetEntries()>5) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(5));
  BinDirs * Bin6 =NULL; 
  if(multdirlist->GetEntries()>6) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  BinDirs * Bin7=NULL;  
  if(multdirlist->GetEntries()>7) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(7));
  //reset the histogram and create clones for all the types we want.
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
  TH3D* histbin2;TH3D* histMETAbin2;TH3D* histMETriggerbin2;
  if(Bin2){
    histbin2		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histMETAbin2	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETriggerbin2 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));
  }
  //Bin 3:
  TH3D* histbin3;TH3D* histMETAbin3;TH3D* histMETriggerbin3;
  if(Bin3){
    histbin3		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histMETAbin3	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETriggerbin3 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  }
  //Bin 4:
  TH3D* histbin4;TH3D* histMETAbin4;TH3D* histMETriggerbin4;
  if(Bin4){
    histbin4		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histMETAbin4	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETriggerbin4 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  }
  //Bin 5:
  TH3D* histbin5;TH3D* histMETAbin5;TH3D* histMETriggerbin5;
  if(Bin5){
    histbin5		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histMETAbin5		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETriggerbin5	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  }
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
//       if((TString(directories->At(i)->GetName()).Contains("Z(-10.00)"))&&(TString(directories->At(i)->GetName()).Contains("Z(5.00)")))continue;      
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
      if(Bin2){
	histbin2->SetBinContent(x,y,z,BinContentbin2);
	histbin2->SetBinError(x,y,z,BinErrorbin2);
	histbin3->SetBinContent(x,y,z,BinContentbin3);
	histbin3->SetBinError(x,y,z,BinErrorbin3);
	histbin4->SetBinContent(x,y,z,BinContentbin4);
	histbin4->SetBinError(x,y,z,BinErrorbin4);
	histbin5->SetBinContent(x,y,z,BinContentbin5);
	histbin5->SetBinError(x,y,z,BinErrorbin5);
      }
      histMETA->SetBinContent(x,y,z,BinContentMETA);
      histMETA->SetBinError(x,y,z,BinErrorMETA);
      histMETAbin1->SetBinContent(x,y,z,BinContentMETAbin1);
      histMETAbin1->SetBinError(x,y,z,BinErrorMETAbin1);
      if(Bin2){
	histMETAbin2->SetBinContent(x,y,z,BinContentMETAbin2);
	histMETAbin2->SetBinError(x,y,z,BinErrorMETAbin2);
	histMETAbin3->SetBinContent(x,y,z,BinContentMETAbin3);
	histMETAbin3->SetBinError(x,y,z,BinErrorMETAbin3);
	histMETAbin4->SetBinContent(x,y,z,BinContentMETAbin4);
	histMETAbin4->SetBinError(x,y,z,BinErrorMETAbin4);
	histMETAbin5->SetBinContent(x,y,z,BinContentMETAbin5);
	histMETAbin5->SetBinError(x,y,z,BinErrorMETAbin5);
      }
      histMETrigger->SetBinContent(x,y,z,BinContentMETrigger);
      histMETrigger->SetBinError(x,y,z,BinErrorMEtrigger);
      histMETriggerbin1->SetBinContent(x,y,z,BinContentMETriggerbin1);
      histMETriggerbin1->SetBinError(x,y,z,BinErrorMEtriggerbin1);
      if(Bin2){
	histMETriggerbin2->SetBinContent(x,y,z,BinContentMETriggerbin2);
	histMETriggerbin2->SetBinError(x,y,z,BinErrorMEtriggerbin2);
	histMETriggerbin3->SetBinContent(x,y,z,BinContentMETriggerbin3);
	histMETriggerbin3->SetBinError(x,y,z,BinErrorMEtriggerbin3);
	histMETriggerbin4->SetBinContent(x,y,z,BinContentMETriggerbin4);
	histMETriggerbin4->SetBinError(x,y,z,BinErrorMEtriggerbin4);
	histMETriggerbin5->SetBinContent(x,y,z,BinContentMETriggerbin5);
	histMETriggerbin5->SetBinError(x,y,z,BinErrorMEtriggerbin5);
      }
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
  if(Bin2){
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
}
void CollectHist(const char* histname,TList * directories, TObjArray* multdirlist,bool pearsonserrors = false,bool collecdivfirst = false){
//   if(TString(histname).CompareTo("DPhi_DEta")==0) return;
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(multdirlist->At(0))->Same()->GetDirectory(Form("%s/divided",directories->At(1)->GetName()))->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  TH2D* hist2d = dynamic_cast<TH2D*>(hist);
  TH3D* hist3d = dynamic_cast<TH3D*>(hist);
  if(hist1d){cout << "TH1D " <<histname<<endl;CollectHist(hist1d,directories,multdirlist);}
  if(hist2d){cout << "TH2D " <<histname<<endl;CollectHist(hist2d,directories,multdirlist, pearsonserrors,collecdivfirst);}
  if(hist3d){cout << "TH3D " <<histname<<endl;CollectHist(hist3d,directories,multdirlist);}
  canvasmaker(histname,multdirlist);
}
void CollectHistbinstats(const char* histname,TList * directories, TObjArray* multdirlist){
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(multdirlist->At(0))->Same()->GetDirectory(Form("%s/bin_stats",directories->At(1)->GetName()))->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  if(hist1d){cout << "TH1D " <<histname <<endl;CollectHistbs(hist1d,directories,multdirlist);}
  canvasmaker(histname,multdirlist);
}




void CollectHistbsfirst(TH1D* histo, TList * directories, BinDirs* All){
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH1D* hist 			= dynamic_cast<TH1D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH1D* histMETA 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH1D* histMETrigger 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  //Doubles to hold the values until they are put into the hists.
  TH1D* histtmp;TH1D* histMETAtmp;TH1D* histMETriggertmp;
  //loop over the bins:
  for(int i=0;i<directories->GetEntries();i++){
    BinDirs * thisbin = dynamic_cast<BinDirs*>(directories->At(i));    
    //test if the relevant histogram exists in this bin:
    if(!dynamic_cast<TH1D*>(thisbin->SameDir("bin_stats")->Get(histo->GetName())))continue;
    //test if All is a centrality bin:
    if(TString(All->Same()->GetName()).Contains("BinM")){
      //compare thisbin with the name:
      if(!TString(thisbin->Same()->GetName()).Contains(All->Same()->GetName()))continue;
    }
    //extract histogram in relevant bin:
    histtmp 		= dynamic_cast<TH1D*>(thisbin->SameDir("bin_stats")->Get(histo->GetName()));
    histMETAtmp 	= dynamic_cast<TH1D*>(thisbin->METADir("bin_stats")->Get(histo->GetName()));
    histMETriggertmp 	= dynamic_cast<TH1D*>(thisbin->METriggerDir("bin_stats")->Get(histo->GetName()));
    hist->Add(histtmp); histMETA->Add(histMETAtmp); histMETrigger->Add(histMETriggertmp);
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
}
void CollectHistfirst(TH1D* histo, TList * directories, BinDirs* All){
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH1D* hist 			= dynamic_cast<TH1D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH1D* histMETA 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH1D* histMETrigger 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  TH1D* histloc;TH1D*histMETAloc;TH1D*histMETriggerloc;
  //loop over all Multiplicity-Vertex bins:
  for(int i=0;i<directories->GetEntries();i++){
    BinDirs * thisbin = dynamic_cast<BinDirs*>(directories->At(i));    
    if(!dynamic_cast<TH1D*>(thisbin->SameDir("same_event")->Get(histo->GetName()))){continue;}
    //test if All is a centrality bin:
    if(TString(All->Same()->GetName()).Contains("BinM")){
      //compare thisbin with the name:
      if(!TString(thisbin->Same()->GetName()).Contains(All->Same()->GetName()))continue;
    }
    //get the hists:
    histloc = dynamic_cast<TH1D*>(thisbin->SameDir("same_event")->Get(histo->GetName()));
    histloc->SetBit(TH1::kIsAverage,kFALSE);
    hist->Add(histloc);
    histMETAloc = dynamic_cast<TH1D*>(thisbin->METADir("same_event")->Get(histo->GetName()));
    histMETAloc->SetBit(TH1::kIsAverage,kFALSE);    
    histMETA->Add(histMETAloc);
    histMETriggerloc = dynamic_cast<TH1D*>(thisbin->METriggerDir("same_event")->Get(histo->GetName()));
    histMETriggerloc->SetBit(TH1::kIsAverage,kFALSE);
    histMETrigger->Add(histMETriggerloc);
  }//end loop over M-V bins
  //save the histograms in the relevant directories:
  All->Same()->GetDirectory("same_event")->cd();
  hist->Write(histo->GetName());
  delete hist;
  All->META()->GetDirectory("same_event")->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  All->METrigger()->GetDirectory("same_event")->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
}
void CollectHistfirst(TH2D* histo, TList * directories, BinDirs* All){
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH2D* hist 			= dynamic_cast<TH2D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH2D* histMETA 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH2D* histMETrigger 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  TH2D* histloc;TH2D*histMETAloc;TH2D*histMETriggerloc;
  
  //loop over all Multiplicity-Vertex bins:
  for(int i=0;i<directories->GetEntries();i++){
    BinDirs * thisbin = dynamic_cast<BinDirs*>(directories->At(i));        
    //test if the relevant histogram exists in this bin:
    if(!dynamic_cast<TH2D*>(thisbin->SameDir("same_event")->Get(histo->GetName()))){continue;}
    //test if All is a centrality bin:
    if(TString(All->Same()->GetName()).Contains("BinM")){
      //compare thisbin with the name:
      if(!TString(thisbin->Same()->GetName()).Contains(All->Same()->GetName()))continue;
    }
    //get the hists:
    histloc = dynamic_cast<TH2D*>(thisbin->SameDir("same_event")->Get(histo->GetName()));
    histloc->SetBit(TH1::kIsAverage,kFALSE);
    hist->Add(histloc);
    histMETAloc = dynamic_cast<TH2D*>(thisbin->METADir("same_event")->Get(histo->GetName()));
    histMETAloc->SetBit(TH1::kIsAverage,kFALSE);    
    histMETA->Add(histMETAloc);
    histMETriggerloc = dynamic_cast<TH2D*>(thisbin->METriggerDir("same_event")->Get(histo->GetName()));
    histMETriggerloc->SetBit(TH1::kIsAverage,kFALSE);
    histMETrigger->Add(histMETriggerloc);
  }//end loop over M-V bins.
  //save the histograms in the relevant directories:
  All->Same()->GetDirectory("same_event")->cd();
  hist->Write(histo->GetName());
  delete hist;
  All->META()->GetDirectory("same_event")->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  All->METrigger()->GetDirectory("same_event")->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
}
void CollectHistfirst(TH3D* histo, TList * directories, BinDirs* All){
   //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH3D* hist 			= dynamic_cast<TH3D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH3D* histMETA 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH3D* histMETrigger 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  TH3D* histloc;TH3D*histMETAloc;TH3D*histMETriggerloc;
  
  //loop over all Multiplicity-Vertex bins:
  for(int i=0;i<directories->GetEntries();i++){
    BinDirs * thisbin = dynamic_cast<BinDirs*>(directories->At(i));        
    //test if the relevant histogram exists in this bin:
    if(!dynamic_cast<TH3D*>(thisbin->SameDir("same_event")->Get(histo->GetName()))){continue;}
    //test if All is a centrality bin:
    if(TString(All->Same()->GetName()).Contains("BinM")){
      //compare thisbin with the name:
      if(!TString(thisbin->Same()->GetName()).Contains(All->Same()->GetName()))continue;
    }
    //get the hists:
    histloc = dynamic_cast<TH3D*>(thisbin->SameDir("same_event")->Get(histo->GetName()));
    histloc->SetBit(TH1::kIsAverage,kFALSE);
    hist->Add(histloc);
    histMETAloc = dynamic_cast<TH3D*>(thisbin->METADir("same_event")->Get(histo->GetName()));
    histMETAloc->SetBit(TH1::kIsAverage,kFALSE);    
    histMETA->Add(histMETAloc);
    histMETriggerloc = dynamic_cast<TH3D*>(thisbin->METriggerDir("same_event")->Get(histo->GetName()));
    histMETriggerloc->SetBit(TH1::kIsAverage,kFALSE);
    histMETrigger->Add(histMETriggerloc);
  }//end loop over M-V bins.
  //save the histograms in the relevant directories:
  All->Same()->GetDirectory("same_event")->cd();
  hist->Write(histo->GetName());
  delete hist;
  All->META()->GetDirectory("same_event")->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  All->METrigger()->GetDirectory("same_event")->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
}
void CollectHistMEfirst(TH2D* histo, TList * directories, BinDirs* All){
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH2D* hist 			= dynamic_cast<TH2D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH2D* histMETA 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH2D* histMETrigger 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  TH2D* histloc;TH2D*histMETAloc;TH2D*histMETriggerloc;
  TH1D* histtrig;
  Double_t scalingf=0.0;
  
  //loop over all Multiplicity-Vertex bins:
  for(int i=0;i<directories->GetEntries();i++){
    BinDirs * thisbin = dynamic_cast<BinDirs*>(directories->At(i));
    //test if the relevant histogram exists in this bin:
    if(!dynamic_cast<TH2D*>(thisbin->SameDir("mixed_event")->Get(histo->GetName()))){continue;}
    //test if All is a centrality bin:
    if(TString(All->Same()->GetName()).Contains("BinM")){
      //compare thisbin with the name:
      if(!TString(thisbin->Same()->GetName()).Contains(All->Same()->GetName()))continue;
    }    
    //get the hists:
    histloc = dynamic_cast<TH2D*>(thisbin->SameDir("mixed_event")->Get(histo->GetName()));
    histloc->SetBit(TH1::kIsAverage,kFALSE);
    histtrig = dynamic_cast<TH1D*>(thisbin->SameDir("mixed_event")->Get("number_of_triggers"));
    hist->Add(histloc,histtrig->Integral());
    scalingf+=histtrig->Integral();
    histMETAloc = dynamic_cast<TH2D*>(thisbin->METADir("mixed_event")->Get(histo->GetName()));
    histMETAloc->SetBit(TH1::kIsAverage,kFALSE);    
    histMETA->Add(histMETAloc,histtrig->Integral());
    histMETriggerloc = dynamic_cast<TH2D*>(thisbin->METADir("mixed_event")->Get(histo->GetName()));
    histMETriggerloc->SetBit(TH1::kIsAverage,kFALSE);
    histMETrigger->Add(histMETriggerloc,histtrig->Integral());
  }//end loop over M-V bins.
  hist->Scale(1.0/scalingf);
  //save the histograms in the relevant directories:
  All->Same()->GetDirectory("mixed_event")->cd();
  hist->Write(histo->GetName());
  delete hist;
  histMETA->Scale(1.0/scalingf);
  All->META()->GetDirectory("mixed_event")->cd();
  histMETA->Write(histo->GetName());
  delete histMETA;
  histMETrigger->Scale(1.0/scalingf);
  All->METrigger()->GetDirectory("mixed_event")->cd();
  histMETrigger->Write(histo->GetName());
  delete histMETrigger;
}
void CollectHistfirst(const char* histname,TList * directories, BinDirs* bindirs){
//   if(TString(histname).CompareTo("DPhi_DEta")==0) return;
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("same_event")->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  TH2D* hist2d = dynamic_cast<TH2D*>(hist);
  TH3D* hist3d = dynamic_cast<TH3D*>(hist);
  if(hist1d){cout << "TH1D " <<histname<<endl;CollectHistfirst(hist1d,directories,bindirs);}
  if(hist2d){cout << "TH2D " <<histname<<endl;CollectHistfirst(hist2d,directories,bindirs);}
  if(hist3d){cout << "TH3D " <<histname<<endl;CollectHistfirst(hist3d,directories,bindirs);}
//   canvasmaker(histname,multdirlist);
}
void CollectHistMEfirst(const char* histname,TList * directories, BinDirs* bindirs){
//   if(TString(histname).CompareTo("DPhi_DEta")==0) return;
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("mixed_event")->Get(histname));
  TH2D* hist2d = dynamic_cast<TH2D*>(hist);
  if(hist2d){cout << "TH2D " <<histname<<endl;CollectHistMEfirst(hist2d,directories,bindirs);}
//   canvasmaker(histname,multdirlist);
}
void CollectHistDivfirst(const char* histname, BinDirs* bindirs){
//   if(TString(histname).CompareTo("DPhi_DEta")==0) return;
  
  TDirectory * divdir = bindirs->SameDir("divided");
  
  TH1D* norm = dynamic_cast<TH1D*>(bindirs->SameDir("same_event")->Get("number_of_triggers"));
  divdir->cd();
  if(!dynamic_cast<TH1D*>(divdir->Get("number_of_triggers")))norm->Write("number_of_triggers");
  TH2D* hist = dynamic_cast<TH2D*>(bindirs->SameDir("same_event")->Get(histname));
  TH2D* histME = dynamic_cast<TH2D*>(bindirs->SameDir("mixed_event")->Get(histname));
  if(norm&&hist&&histME){
//     cout << "1. "<< hist->GetName()<<endl;
    TH2D * histc = dynamic_cast<TH2D*>( hist->Clone(Form("%sdiv"		,hist->GetName())));
    histc->Divide(histME);
    histc->Scale(1.0/norm->Integral());
    divdir->cd();
//     cout << histname << " " << hist->GetName()<<endl;
    histc->Write(histname);
  }
  divdir = bindirs->METADir("divided");
  norm = dynamic_cast<TH1D*>(bindirs->METADir("same_event")->Get("number_of_triggers"));
  hist = dynamic_cast<TH2D*>(bindirs->METADir("same_event")->Get(histname));
  histME = dynamic_cast<TH2D*>(bindirs->METADir("mixed_event")->Get(histname));
  if(norm&&hist&&histME){
//     cout << "1. "<< hist->GetName()<<endl;
    TH2D * histc = dynamic_cast<TH2D*>( hist->Clone(Form("%sdiv"		,hist->GetName())));
    histc->Divide(histME);
    histc->Scale(1.0/norm->Integral());
    divdir->cd();
//     cout << histname << " " << hist->GetName()<<endl;
    histc->Write(histname);
  }
  divdir = bindirs->METriggerDir("divided");
  norm = dynamic_cast<TH1D*>(bindirs->METriggerDir("same_event")->Get("number_of_triggers"));
  hist = dynamic_cast<TH2D*>(bindirs->METriggerDir("same_event")->Get(histname));
  histME = dynamic_cast<TH2D*>(bindirs->METriggerDir("mixed_event")->Get(histname));
  if(norm&&hist&&histME){
//     cout << "1. "<< hist->GetName()<<endl;
    TH2D * histc = dynamic_cast<TH2D*>( hist->Clone(Form("%sdiv"		,hist->GetName())));
    histc->Divide(histME);
    histc->Scale(1.0/norm->Integral());
    divdir->cd();
//     cout << histname << " " << hist->GetName()<<endl;
    histc->Write(histname);
  }
  //   canvasmaker(histname,multdirlist);
}
void CollectHistbinstatsfirst(const char* histname,TList * directories, BinDirs* bindirs){
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("bin_stats")->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  if(hist1d){cout << "TH1D " <<histname <<endl;CollectHistbsfirst(hist1d,directories,bindirs);}
//   canvasmaker(histname,multdirlist);
}




///////////////////////////////////////


//functions to correct with ME types.
  
//function for minimizing the histogram edge:
Double_t Chi2(Double_t scaleMETA, Double_t scaleMETrigger){
  Double_t grad = 0.0;
  TH2D* testhist = dynamic_cast<TH2D*>(gSameEvent->Clone("testhist"));
  testhist->Add(gMETA,-1.0*scaleMETA);
  testhist->Add(gMETrigger,-1.0*scaleMETrigger); 
  int pihalfbin = testhist->GetXaxis()->FindBin(TMath::Pi()/2.0);
  for(int x=1;x<testhist->GetNbinsX();x++){
    Double_t gradloc =testhist->GetBinContent(x,1) - testhist->GetBinContent(x+1,1);
    grad += gradloc*gradloc;
    gradloc = testhist->GetBinContent(x,testhist->GetNbinsY()) - testhist->GetBinContent(x+1,testhist->GetNbinsY());
    grad += gradloc*gradloc;
    gradloc =testhist->GetBinContent(x,pihalfbin) - testhist->GetBinContent(x+1,pihalfbin);
    grad += gradloc*gradloc;
    gradloc = testhist->GetBinContent(x,pihalfbin) - testhist->GetBinContent(x+1,pihalfbin);
    grad += gradloc*gradloc;
//     gradloc = testhist->GetBinContent(x,testhist->GetXaxis()->FindBin(TMath::Pi()/2.0)) - testhist->GetBinContent(x+1,testhist->GetXaxis()->FindBin(TMath::Pi()/2.0));
//     grad += gradloc*gradloc;
//     gradloc = testhist->GetBinContent(x,testhist->GetXaxis()->FindBin(TMath::Pi()/2.0)-1) - testhist->GetBinContent(x+1,testhist->GetXaxis()->FindBin(TMath::Pi()/2.0)-1);
//     grad += gradloc*gradloc;
//     gradloc = testhist->GetBinContent(x,testhist->GetXaxis()->FindBin(TMath::Pi()/2.0)+1) - testhist->GetBinContent(x+1,testhist->GetXaxis()->FindBin(TMath::Pi()/2.0)+1);
//     grad += gradloc*gradloc;
//     gradloc = testhist->GetBinContent(x,2) - testhist->GetBinContent(x+1,2);
//     grad += gradloc*gradloc;
//     gradloc = testhist->GetBinContent(x,testhist->GetNbinsY()-1) - testhist->GetBinContent(x+1,testhist->GetNbinsY()-1);
//     grad += gradloc*gradloc;
//     gradloc = testhist->GetBinContent(x,3) - testhist->GetBinContent(x+1,3);
//     grad += gradloc*gradloc;
//     gradloc = testhist->GetBinContent(x,testhist->GetNbinsY()-2) - testhist->GetBinContent(x+1,testhist->GetNbinsY()-2);
//     grad += gradloc*gradloc;    
  }
//   Double_t gradEta = 0.0;
//   TH2D* testhist2 = dynamic_cast<TH2D*>(gSameEventETA->Clone("testhisteta"));  
//   testhist2->Add(gMETAETA,-1.0*scaleMETA);
//   testhist2->Add(gMETriggerETA,-1.0*scaleMETrigger);    
//   TH1D* deta12ss = testhist2->ProjectionX("testhistetaproj");deta12ss->Reset();
//   fillwithbinnr(testhist2,deta12ss,1);
//   gEtasigRange = 1.0;
//   TF1* fitbg =  new TF1("fgs",CBKGPol0,-gkEtaFitRange,gkEtaFitRange,2) ;
//   fitbg->SetParNames("B", "peakpos" ) ;
//   fitbg->SetParameters(0.0,0.0);
//   fitbg->SetParLimits(1,-0.1,0.1);
//   TF1* rembg = new TF1("fgs",CBKGPol0p,-gkEtaFitRange,gkEtaFitRange,1) ;
//   rembg->SetParNames("B", "peakpos" ) ;
//   rembg->SetParameters(0.0,0.0);    
// 
//   TFitResultPtr bkgresult = deta12ss->Fit(fitbg,"SQ","",-gkEtaFitRange,gkEtaFitRange);
//     if(int(bkgresult)!=4000)//if "error", try again with more:
//       bkgresult = deta12ss->Fit(fitbg,"MSQ","",-gkEtaFitRange,gkEtaFitRange);    
//   rembg->SetParameter(0,fitbg->GetParameter(0));
//   deta12ss->Add(rembg,-1.0);
//   Double_t integral = deta12ss->Integral(deta12ss->GetXaxis()->FindBin(-gEtasigRange),deta12ss->GetXaxis()->FindBin(gEtasigRange));
//   gradEta = integral*integral;
  delete testhist;//delete testhist2;delete fitbg;delete rembg;
  return 100000*grad;
}
//____________________________________________________________________________________
void FcnFitScale(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t) {
  //
  // Chi2 function
  // par[0] - scale factor for the META 	background histogram
  // par[0] - scale factor for the METrigger 	background histogram
  //
  Double_t chi2 = Chi2(par[0],par[1]);
//   cout << chi2<< " " << Chi2(par[0],par[1]+0.1)<<endl;
  f = chi2;
} 

  
void GetScalingFactors(BinDirs * BinDir,const char* in, const char* out,Double_t * METAscale, Double_t * METriggerscale,bool ispp, bool isPbPb, int iter)
{
  //Function to find the relevant scaling factors.
  //Clean and create the output directory
  resultsdirectory(BinDir->Same(),out);
  if(!BinDir->SameDir(in)) return;
  bool issmallptbin =  TString(BinDir->Same()->GetPath()).Contains("0_1");
  gislowpTbin = issmallptbin;
  
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * SameHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METAHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_META")));
  tmp = BinDir->METriggerdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METRiggerHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_METrigger")));
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DEta_12_SameSide");if(!tmp){return;}
  TH2D * SameHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DEta_12_SameSide");if(!tmp){return;}  
  TH2D * METAHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_META")));
   tmp = BinDir->METriggerdiv()->Get("DPhi_1_DEta_12_SameSide");if(!tmp){return;}
  TH2D * METRiggerHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_METrigger")));  
  delete tmp;

  gSameEvent 	= SameHistPhiPhi;
  gMETA 	= METAHistPhiPhi;
  gMETrigger 	= METRiggerHistPhiPhi;
  gSameEventETA	= SameHistPhiEta;
  gMETAETA 	= METAHistPhiEta;
  gMETriggerETA	= METRiggerHistPhiEta;  
  

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
  step[0] = 0.0001; 
  vstart[1] = *METriggerscale;
//   if(issmallptbin){
//     Double_t scalingMETrigger = 0.5*SameHistPhiEta->GetBinContent(SameHistPhiEta->GetXaxis()->FindBin(0.0)+1,1);//bin at det12 =0 in the first bin in dphi
//     scalingMETrigger += 0.5*SameHistPhiEta->GetBinContent(SameHistPhiEta->GetXaxis()->FindBin(0.0)-1,1);//bin at det12 =0 in the first bin in dphi
//     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(1,1);//bin at the edge of det12 in the first bin in dphi
//     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(2,1);//bin at the edge of det12 in the first bin in dphi
//     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX(),1);//bin at the edge of det12 in the first bin in dphi
//     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX()-1,1);//bin at the edge of det12 in the first bin in dphi
// //     scalingMETrigger += SameHistPhiEta->GetBinContent(SameHistPhiEta->GetXaxis()->FindBin(0.0),2);//bin at det12 =0 in the first bin in dphi
// //     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(1,2);//bin at the edge of det12 in the first bin in dphi
// //     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(2,2);//bin at the edge of det12 in the first bin in dphi
// //     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX(),2);//bin at the edge of det12 in the first bin in dphi
// //     scalingMETrigger -= 0.25*SameHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX()-1,2);//bin at the edge of det12 in the first bin in dphi
//     Double_t Divby = 0.5*METRiggerHistPhiEta->GetBinContent(SameHistPhiEta->GetXaxis()->FindBin(0.0)+1,1);//bin at det12 =0 in the first bin in dphi
//     Divby += 0.5*METRiggerHistPhiEta->GetBinContent(SameHistPhiEta->GetXaxis()->FindBin(0.0)-1,1);//bin at det12 =0 in the first bin in dphi
//     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(1,1);//bin at the edge of det12 in the first bin in dphi
//     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(2,1);//bin at the edge of det12 in the first bin in dphi
//     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX(),1);//bin at the edge of det12 in the first bin in dphi
//     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX()-1,1);//bin at the edge of det12 in the first bin in dphi
// //     Divby = METRiggerHistPhiEta->GetBinContent(SameHistPhiEta->GetXaxis()->FindBin(0.0),2);//bin at det12 =0 in the first bin in dphi
// //     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(1,2);//bin at the edge of det12 in the first bin in dphi
// //     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(2,2);//bin at the edge of det12 in the first bin in dphi
// //     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX(),2);//bin at the edge of det12 in the first bin in dphi
// //     Divby -= 0.25*METRiggerHistPhiEta->GetBinContent(SameHistPhiEta->GetNbinsX()-1,2);//bin at the edge of det12 in the first bin in dphi
//     vstart[1] = scalingMETrigger/Divby;
//   }
  step[1] = 0.0001;  
  
  minuit->mnparm(0, "METAscale", 	vstart[0], step[0], 0.0, 1000000.0, ierflg);
  minuit->mnparm(1, "METriggerscale", 	vstart[1], step[1], 0.0, 1000000.0, ierflg);
//   if(issmallptbin) minuit->FixParameter(1);
  
  arglist[0] = 500;
  arglist[1] = 1.;

  minuit->SetMaxIterations(10000);
  minuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
  minuit->mnexcm("MIGRAD" , arglist ,2,ierflg); 
//   if(issmallptbin){
//     Double_t t1 = 0.0;
//     Double_t te = 0.0;
//     minuit->GetParameter(0,t1,te);
//     cout << "parameters are: "<< t1<<" and ";
//     minuit->GetParameter(1,t1,te);
//     cout <<t1 <<". After fixing p1: ";
//     minuit->FixParameter(1);
//     minuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
//     minuit->GetParameter(0,t1,te);
//     cout << t1<<" and ";
//     minuit->GetParameter(1,t1,te);
//     cout <<t1<<". After fixing p2: ";
//     minuit->Release(1);
//     minuit->FixParameter(0);
//     minuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
//     minuit->GetParameter(0,t1,te);
//     cout << t1<<" and ";
//     minuit->GetParameter(1,t1,te);
//     cout <<t1<<endl;
//   }
  
// //   if(ispp){
// //     if(iter ==1){
// //       Double_t ScalingMETA =SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
// //       ScalingMETA +=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
// //       ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       Double_t Divideby = METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
// //       Divideby +=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
// //       Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       if(Divideby >0.0) *METAscale = ScalingMETA/Divideby;    
// //       else *METAscale = 0.0;
// //       if(*METAscale<0.0)*METAscale = 0.0;
// //       
// //       Double_t ScalingMETrigger = SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(0.0),SameHistPhiEta->GetXaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
// //       ScalingMETrigger -=0.5* SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 + 0.2));
// //       ScalingMETrigger -=0.5* SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 + 0.2));
// //       Double_t DividebyMETrigger = METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(0.0),METRiggerHistPhiEta->GetXaxis()->FindBin(0.0),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
// //       DividebyMETrigger -=0.5* METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(-1.0),METRiggerHistPhiEta->GetXaxis()->FindBin(-1.0),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
// //       DividebyMETrigger -=0.5* METRiggerHistPhiEta->Integral(METRiggerHistPhiEta->GetXaxis()->FindBin(1.0),METRiggerHistPhiEta->GetXaxis()->FindBin(1.0),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0 - 0.2),METRiggerHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0+0.2));
// //       if(DividebyMETrigger>0.0) *METriggerscale = ScalingMETrigger/DividebyMETrigger;
// //       else *METriggerscale = 0.0;
// //       if(*METriggerscale<0.0)*METriggerscale = 0.0;
// //     }
// //   }
// //   if(isPbPb){
// //     if(iter ==1||iter==2){
// //       Double_t ScalingMETA =SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
// //       ScalingMETA +=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(0.0),SameHistPhiEta->GetYaxis()->FindBin(0.0));
// //       ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(-1.5),SameHistPhiEta->GetXaxis()->FindBin(-1.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       ScalingMETA -=SameHistPhiEta->Integral(SameHistPhiEta->GetXaxis()->FindBin(1.0),SameHistPhiEta->GetXaxis()->FindBin(1.5),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),SameHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       Double_t Divideby = METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
// //       Divideby +=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(0.0),METAHistPhiEta->GetYaxis()->FindBin(0.0));
// //       Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(-1.5),METAHistPhiEta->GetXaxis()->FindBin(-1.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       Divideby -=METAHistPhiEta->Integral(METAHistPhiEta->GetXaxis()->FindBin(1.0),METAHistPhiEta->GetXaxis()->FindBin(1.5),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0),METAHistPhiEta->GetYaxis()->FindBin(TMath::Pi()/2.0));
// //       if(Divideby >0.0) *METAscale = ScalingMETA/Divideby;    
// //       else *METAscale = 0.0;
// //       if(*METAscale<0.0)*METAscale = 0.0;
// // 
// //       Double_t ScalingMETrigger =SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetNbinsX(),1);
// //       ScalingMETrigger +=SameHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetNbinsX());
// //       ScalingMETrigger -=SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-1.0),SameHistPhiPhi->GetXaxis()->FindBin(4.0));
// //       ScalingMETrigger -=SameHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(4.0),SameHistPhiPhi->GetXaxis()->FindBin(-1.0));
// //       Double_t DividebyMETrigger = METRiggerHistPhiPhi->GetBinContent(1,SameHistPhiPhi->GetNbinsX()); 
// //       DividebyMETrigger += METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetNbinsX(),1);
// //       DividebyMETrigger -=METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(-1.0),SameHistPhiPhi->GetXaxis()->FindBin(4.0));
// //       DividebyMETrigger -=METRiggerHistPhiPhi->GetBinContent(SameHistPhiPhi->GetXaxis()->FindBin(4.0),SameHistPhiPhi->GetXaxis()->FindBin(-1.0));
// //       if(DividebyMETrigger>0.0) *METriggerscale = ScalingMETrigger/DividebyMETrigger;
// //       else *METriggerscale = 0.0;
// //       if(*METriggerscale<0.0)*METriggerscale = 0.0;
// //     }    
// //   }
  Double_t temp =0.0;
  Double_t tempe = 0.0;
  minuit->GetParameter(0, temp, tempe); 
  *METAscale = temp;
  minuit->GetParameter(1, temp, tempe); 
  *METriggerscale = temp;
  cout << *METAscale <<" " <<*METriggerscale<<endl;
  delete SameHistPhiPhi;delete METAHistPhiPhi;delete METRiggerHistPhiPhi;delete SameHistPhiEta;delete METAHistPhiEta;delete METRiggerHistPhiEta;

  
}


void Correct(const char* histname, BinDirs * BinDir,const char* in,const char* out, Double_t  METAscale, Double_t  METriggerscale)
{
//   cout << histname << " "<<BinDir->Same()->GetName()<<"  "<< METAscale<<" "<< BinDir->META()->GetName()<< "  "<< METAscale<<" "<< BinDir->METrigger()->GetName()<<endl;
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get(histname);if(!tmp){return;}
  TH2D * SameHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname)));
  tmp = BinDir->METAdiv()->Get(histname);if(!tmp){return;}
  TH2D * METAHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname)));
  tmp = BinDir->METriggerdiv()->Get(histname);if(!tmp){return;}
  TH2D * METriggerHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname)));  
  delete tmp;
//   cout << "Signal 0,edge ="<<SameHist->GetBinContent(5,SameHist->GetYaxis()->FindBin(0.0))<< " META 0,edge:"<< METAHist->GetBinContent(5,SameHist->GetYaxis()->FindBin(0.0))*METAscale<<" METrigger 0,edge:"<< METriggerHist->GetBinContent(5,SameHist->GetYaxis()->FindBin(0.0))*METriggerscale<<endl;
  //Scale the Mixed hists with the correct scales:
  METAHist->Scale(METAscale);
  METriggerHist->Scale(METriggerscale);
 //perform the removal:  
  AddSigBins(SameHist,METAHist,METriggerHist);
//   SameHist->Add(METAHist,-1.0*METAscale);SameHist->Add(METriggerHist,-1.0*METriggerscale);
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
  //perform the removal:  
  METAHist1->Scale(METAscale);
  METriggerHist1->Scale(METriggerscale);
  AddSigBins(SameHist1,METAHist1,METriggerHist1);
//   SameHist1->Add(METAHist1,-1.0*METAscale);SameHist1->Add(METriggerHist1,-1.0*METriggerscale);
  METAHist2->Scale(METAscale);
  METriggerHist2->Scale(METriggerscale);
  AddSigBins(SameHist2,METAHist2,METriggerHist2);
//   SameHist2->Add(METAHist2,-1.0*METAscale);SameHist2->Add(METriggerHist2,-1.0*METriggerscale);
  METAHist3->Scale(METAscale);
  METriggerHist3->Scale(METriggerscale);
  AddSigBins(SameHist3,METAHist3,METriggerHist3);
//   SameHist3->Add(METAHist3,-1.0*METAscale);SameHist3->Add(METriggerHist3,-1.0*METriggerscale);
  METAHist4->Scale(METAscale);
  METriggerHist4->Scale(METriggerscale);
  AddSigBins(SameHist4,METAHist4,METriggerHist4);
//   SameHist4->Add(METAHist4,-1.0*METAscale);SameHist4->Add(METriggerHist4,-1.0*METriggerscale);
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
//   resultsdirectory(BinDir->Same(),"iteration1p1");
//   resultsdirectory(BinDir->Same(),"iteration1m1");
  cout << BinDir->path()<<" META:"<<METAScale<<" METriggerScale:"<<METriggerScale<<endl;
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,METriggerScale,"iteration1");
//   savedircontent(BinDir,METAScale+0.01,METriggerScale+0.01,"iteration1p1");
//   savedircontent(BinDir,METAScale-0.01,METriggerScale-0.01,"iteration1m1");

  Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","iteration1",METAScale,METriggerScale);

//   Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","iteration1p1",METAScale+0.01,METriggerScale+0.01);
//   Correct("DPhi_1_DEta_12",BinDir,"divided","iteration1p1",METAScale+0.01,METriggerScale+0.01);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","iteration1p1",METAScale+0.01,METriggerScale+0.01);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","iteration1p1",METAScale+0.01,METriggerScale+0.01);
//   Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","iteration1p1",METAScale+0.01,METriggerScale+0.01);
// 
//   Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","iteration1m1",METAScale,METriggerScale-0.25);
//   Correct("DPhi_1_DEta_12",BinDir,"divided","iteration1m1",METAScale,METriggerScale-0.25);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","iteration1m1",METAScale,METriggerScale-0.25);
//   Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","iteration1m1",METAScale,METriggerScale-0.25);
//   Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","iteration1m1",METAScale,METriggerScale-0.25);  
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

}   
void CorrectScan(BinDirs* BinDir)
{
  //If the bin is empty, skip the entire bin:
  if(!dynamic_cast<TH1D*>(BinDir->Samediv()->Get("number_of_triggers"))){cout << "skipping empty bin "<<BinDir->path()<<endl;return;}
  //create the scaling doubles.
  Double_t METAScale 		= 0.0;
  Double_t METriggerScale 	= 0.0;
  //First iteration: 
  GetScalingFactors(BinDir,"divided","scan",&METAScale,&METriggerScale,false,true,1);
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,METriggerScale,"scan");
  TDirectory * Basedir = BinDir->Same();
  for(int i =-10;i<=10;i++){
    for(int j =-10;j<=10;j++){
      cout << i<<" " << j << endl;
      Double_t METAScaleloc   = METAScale+METAScale/100.0*i;
      Double_t METriggerScaleloc = METriggerScale+METriggerScale/100.0*j;
      TString direct;
      if(i<0&&j<0) direct = TString(Form("META_min_%i_METrigger_min_%i",i,j));
      if(i<0&&j>=0) direct = TString(Form("META_min_%i_METrigger_pl_%i",i,j));
      if(i>=0&&j<0) direct = TString(Form("META_pl_%i_METrigger_min_%i",i,j));
      if(i>=0&&j>=0) direct = TString(Form("META_pl_%i_METrigger_pl_%i",i,j));
//       TDirectory * locdir =
      resultsdirectory(Basedir,direct.Data());
      Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided",direct.Data(),METAScaleloc,METriggerScaleloc);
      Correct("DPhi_1_DEta_12",BinDir,"divided",direct.Data(),METAScaleloc,METriggerScaleloc);
      Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided",direct.Data(),METAScaleloc,METriggerScaleloc);
      Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided",direct.Data(),METAScaleloc,METriggerScaleloc);
      Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided",direct.Data(),METAScaleloc,METriggerScaleloc);  
    }
  }
  
  


}  
////////////////////////////////////



void CollectMVbins(bool CDivfirst = false){
  //take the histograms from the bins and average over them.
  TFile * outfile =  TFile::Open("results.root","UPDATE");
  TList * folderlist = outfile->GetListOfKeys();
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    if(!TString(folder).Contains("ThreePartTracks"))continue;
    TObjArray * multdirlist = new TObjArray(6);
    BinDirs * divsame = new BinDirs(outfile->GetDirectory(Form("%s/",folder)),outfile->GetDirectory(Form("%s/META",folder)),outfile->GetDirectory(Form("%s/METrigger",folder)),true);

    multdirlist->Add(divsame);
    //List of directories for multiplicity bins:
    TList * directories = GetMZDirectories(divsame);
    //Go through the list and find the multiplicity/centrality binning in order to sum over them and add them to the TObjArray:
    TString  s = TString("");
    for(int i=0; i<directories->GetEntries();i++){
      if(TString(directories->At(i)->GetName()).Contains("Z")){
	TString * bin = new TString(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName());
	if(bin->CompareTo(s.Data())!=0&&bin->Contains("BinM")){
	  BinDirs * dir = new BinDirs(outfile->GetDirectory(Form("%s/",folder)),bin->Data(),true);
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
    while(histtokens.NextToken()){CollectHist(histtokens.Data(),directories,multdirlist,false,CDivfirst);}

  }
  outfile->Close();
  delete outfile;
}

void CollectMVbinsFirst(const char* options){
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
  
  //take the histograms from the bins and average over them.
  TFile * outfile = new TFile("results.root","UPDATE");
  TList * folderlist = outfile->GetListOfKeys();
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    cout << "current folder: "<< folder<<endl;
    if(TString(folder).Contains("ThreePartTracks")&&(TString(folder).Contains("8_16")||TString(folder).Contains("4_8"))){
      if(ispp){
	TDirectory* collectdir = resultsdirectory(outfile->GetDirectory(Form("%s/",folder)),"Collected");
	TDirectory* collectdirMETA = resultsdirectory(outfile->GetDirectory(Form("%s/META",folder)),"Collected");
	TDirectory* collectdirMETrigger = resultsdirectory(outfile->GetDirectory(Form("%s/METrigger",folder)),"Collected");
	BinDirs * divsame = new BinDirs(collectdir,collectdirMETA,collectdirMETrigger,true);
	//make the correct directories for saving:
	resultsdirectory(divsame->Same(),"same_event");resultsdirectory(divsame->META(),"same_event");resultsdirectory(divsame->METrigger(),"same_event");
	resultsdirectory(divsame->Same(),"mixed_event");resultsdirectory(divsame->META(),"mixed_event");resultsdirectory(divsame->METrigger(),"mixed_event");
	//List of directories for multiplicity bins:
	TList * directories = GetMZDirectories(outfile->GetDirectory(Form("%s/",folder)));
	//Get Tokens for all histograms:
	TStringToken histtokensbinstats = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("bin_stats"));
	while(histtokensbinstats.NextToken()){CollectHistbinstatsfirst(histtokensbinstats.Data(),directories,divsame);}
	TStringToken histtokens = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("same_event"));
	while(histtokens.NextToken()){CollectHistfirst(histtokens.Data(),directories,divsame);}
	TStringToken histtokens2 = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("mixed_event"));
	while(histtokens2.NextToken()){CollectHistMEfirst(histtokens2.Data(),directories,divsame);}      
	TStringToken histtokens3 = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("same_event"));
	while(histtokens3.NextToken()){CollectHistDivfirst(histtokens3.Data(),divsame);}            
      }
      if(isPbPb){
	TObjArray * multdirlist = new TObjArray(6);
	BinDirs * divsame = new BinDirs(outfile->GetDirectory(Form("%s/",folder)),outfile->GetDirectory(Form("%s/META",folder)),outfile->GetDirectory(Form("%s/METrigger",folder)),true);
	TDirectory* collectdir = resultsdirectory(divsame->Same(),"Collected");
	TDirectory* collectdirMETA = resultsdirectory(divsame->META(),"Collected");
	TDirectory* collectdirMETrigger = resultsdirectory(divsame->METrigger(),"Collected");
	//List of directories for multiplicity bins:
	TList * directories = GetMZDirectories(outfile->GetDirectory(Form("%s/",folder)));
	//Go through the list and find the multiplicity/centrality binning in order to sum over them and add them to the TObjArray:
	TString  s = TString("");
	for(int i=0; i<directories->GetEntries();i++){
	  if(TString(dynamic_cast<BinDirs*>(directories->At(i))->Same()->GetName()).Contains("Z(")){
	    TString * bin = new TString(TString(dynamic_cast<BinDirs*>(directories->At(i))->Same()->GetName()).Tokenize("Z")->At(0)->GetName());
	    if(bin->CompareTo(s.Data())!=0&&bin->Contains("BinM")){
	      BinDirs * dir = new BinDirs(resultsdirectory(collectdir,bin->Data()),resultsdirectory(collectdirMETA,bin->Data()),resultsdirectory(collectdirMETrigger,bin->Data()));
	      resultsdirectory(dir->Same(),"same_event");resultsdirectory(dir->META(),"same_event");resultsdirectory(dir->METrigger(),"same_event");
	      resultsdirectory(dir->Same(),"mixed_event");resultsdirectory(dir->META(),"mixed_event");resultsdirectory(dir->METrigger(),"mixed_event");
	      multdirlist->Add(dir);
	      TStringToken histtokensbinstats = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("bin_stats"));
	      while(histtokensbinstats.NextToken()){CollectHistbinstatsfirst(histtokensbinstats.Data(),directories,dir);}	      
	      TStringToken histtokens = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("same_event"));
	      while(histtokens.NextToken()){CollectHistfirst(histtokens.Data(),directories,dir);}
	      TStringToken histtokens2 = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("mixed_event"));
	      while(histtokens2.NextToken()){CollectHistMEfirst(histtokens2.Data(),directories,dir);}      
	      TStringToken histtokens3 = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("same_event"));
	      while(histtokens3.NextToken()){CollectHistDivfirst(histtokens3.Data(),dir);}      
	      s.Clear();
	      s.Append(bin->Data());
	      
	    }
	  }
	}


      }
    }
  }
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
  TList * folderlist = rfile->GetListOfKeys();
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();  
    TDirectory * folderdir = rfile->GetDirectory(folder);
    //Find a list over all directories beginning with BinM and create a BinDirs object for each:
    TList * dirs = folderdir->GetListOfKeys();
    TObjArray * BinDirAr = new TObjArray();
    for(int i = 0;i<dirs->GetEntries();i++)
    {
//       if(TString(dirs->At(i)->GetName()).BeginsWith("BinM")&&!TString(dirs->At(i)->GetName()).Contains("Z")&&isPbPb){
      if(TString(dirs->At(i)->GetName()).BeginsWith("BinM")&&isPbPb){
	BinDirs * dir = new BinDirs(folderdir,dirs->At(i)->GetName(),false);
	BinDirAr->Add(dir);
      }
       if(TString(dirs->At(i)->GetName()).CompareTo("divided")==0&&ispp){
	BinDirs * dir = new BinDirs(folderdir->GetDirectory(""),folderdir->GetDirectory("META"),folderdir->GetDirectory("METrigger"),false);
	BinDirAr->Add(dir);
      }
//       else if(TString(dirs->At(i)->GetName()).CompareTo("Collected")==0&&ispp){
// 	TDirectory * dividir = folderdir->GetDirectory("Collected");
// 	TDirectory * METAdir = folderdir->GetDirectory("META")->GetDirectory("Collected");
// 	TDirectory * METriggerdir = folderdir->GetDirectory("METrigger")->GetDirectory("Collected");
// 	BinDirs * dir = new BinDirs(dividir,METAdir,METriggerdir,false);
// 	BinDirAr->Add(dir);
//       }
//       else if(TString(dirs->At(i)->GetName()).CompareTo("Collected")==0&&isPbPb){
// 	TDirectory * dividir = folderdir->GetDirectory("Collected");
// 	TDirectory * METAdir = folderdir->GetDirectory("META")->GetDirectory("Collected");
// 	TDirectory * METriggerdir = folderdir->GetDirectory("METrigger")->GetDirectory("Collected");
// 	TList * underdirs = dividir->GetListOfKeys();
// 	for(int k=0;k<underdirs->GetEntries();k++){
// 	  TDirectory * same = dividir->GetDirectory(underdirs->At(k)->GetName());
// 	  TDirectory * META = METAdir->GetDirectory(underdirs->At(k)->GetName());
// 	  TDirectory * METrigger = METriggerdir->GetDirectory(underdirs->At(k)->GetName());
// 	  BinDirs* dir = new BinDirs(same,META,METrigger);
// 	  BinDirAr->Add(dir);	  
// 	}
//       }
    }
    for(int i=0;i<BinDirAr->GetEntriesFast();i++){
      if(ispp)	Correctpp(	dynamic_cast<BinDirs*>(BinDirAr->At(i)));
      if(isPbPb)CorrectPbPb(	dynamic_cast<BinDirs*>(BinDirAr->At(i)));
    }
  }
  rfile->Close();
}

void ScanCorrections(const char* options)
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
  TList * folderlist = rfile->GetListOfKeys();
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();  
    if(!TString(folder).Contains("4_8_0_1"))continue;
    TDirectory * folderdir = rfile->GetDirectory(folder);
    //Find a list over all directories beginning with BinM and create a BinDirs object for each:
    TList * dirs = folderdir->GetListOfKeys();
    TObjArray * BinDirAr = new TObjArray();
    for(int i = 0;i<dirs->GetEntries();i++)
    {
      if(TString(dirs->At(i)->GetName()).BeginsWith("BinM(0.00")&&!TString(dirs->At(i)->GetName()).Contains("Z")){
	BinDirs * dir = new BinDirs(folderdir,dirs->At(i)->GetName(),false);
	BinDirAr->Add(dir);
      }

    }
    for(int i=0;i<BinDirAr->GetEntriesFast();i++){
      if(isPbPb)CorrectScan(	dynamic_cast<BinDirs*>(BinDirAr->At(i)));
      if(ispp) continue;
    }
  }
  rfile->Close();
}

void yield(const char* options){
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
  if(!(ispp||isPbPb))return;
  TFile * rfile = TFile::Open("results.root","UPDATE");
  TObjArray * dirarray = new TObjArray();
  TObjArray * yieldarray = new TObjArray();
  TList * folderlist = rfile->GetListOfKeys();
  for(int l = 0;l<folderlist->GetEntries();l++){
    const char* folder = folderlist->At(l)->GetName();  
    if(!TString(folder).Contains("ThreePart"))continue;
    TDirectory * folderdir = rfile->GetDirectory(folder);
    TList * dirs = folderdir->GetListOfKeys();
    for(int i = 0;i<dirs->GetEntries();i++){
      if(TString(dirs->At(i)->GetName()).BeginsWith("BinM")&&!TString(dirs->At(i)->GetName()).Contains("Z")){
	unsigned int j=1;
	cout << dirs->At(i)->GetName()<<endl;
	while(j<2){
	  TDirectory * tmp1 = dynamic_cast<TDirectory*>(folderdir->Get(Form("%s/iteration%u",dirs->At(i)->GetName(),j)));
	  if(!tmp1){j = 11;continue;}
	  TDirectory * tmp2 = dynamic_cast<TDirectory*>(folderdir->Get(Form("%s/iteration%u",dirs->At(i)->GetName(),j+1)));
	  if(tmp1&&!tmp2){
	    dirarray->Add(tmp1);
// 	    tmp1->pwd();
	    yieldarray->Add(resultsdirectory(folderdir->GetDirectory(Form("%s",dirs->At(i)->GetName())),"yield"));
	    j = 11;
	  }
	  if(tmp1&&tmp2){j +=1;}
	}
      }
    }
    unsigned int k = 1;
    while(k<2){
      TDirectory * tmp1 = dynamic_cast<TDirectory*>(folderdir->Get(Form("iteration%u",k)));
      if(!tmp1){k = 11;continue;}
      TDirectory * tmp2 = dynamic_cast<TDirectory*>(folderdir->Get(Form("iteration%u",k+1)));
      if(tmp1&&!tmp2){
	dirarray->Add(folderdir->GetDirectory(Form("iteration%u",k)));
	yieldarray->Add(resultsdirectory(folderdir->GetDirectory(""),"yield"));
	
	k = 11;
      }
      if(tmp1&&tmp2){k +=1;}
    }
//     TDirectory * tmp = dynamic_cast<TDirectory*>(folderdir->GetDirectory("Collected")->GetDirectory("iteration1"));
//     if(tmp){ dirarray->Add(tmp);
// 	     yieldarray->Add(resultsdirectory(folderdir->GetDirectory("Collected"),"yield"));
//     }
  }
  Double_t etalimit = 0.0;
  if(ispp) etalimit = 1.4;
  if(isPbPb) etalimit = 1.0;
  TDirectory * dir;
  TDirectory * yielddir;
  for(int i=0;i<dirarray->GetEntriesFast();i++){
    dir = dynamic_cast<TDirectory*>(dirarray->At(i));
    dir->pwd();
    yielddir = dynamic_cast<TDirectory*>(yieldarray->At(i));
    extractbinyield(dir,yielddir,etalimit);
  }
  rfile->Close();
}
void draw(const char* option = ""){
  cout << option << endl;
  TFile * rfile = TFile::Open("results.root","UPDATE");
  TList * folderlist = rfile->GetListOfKeys();
  TObjArray * Bins = new TObjArray();
  if(TString("pp").CompareTo(option)==0){
    for(int l = folderlist->GetEntries()-1;l>=0;l--){
      TString folder = TString(folderlist->At(l)->GetName());  
      if(folder.Contains("ThreePartTracks")){
	Bins->Add(new ResHists(folder,rfile,false));
      }
    }
    TDirectory * results = resultsdirectory(rfile->GetDirectory("/"),"Results");
    results->cd();
    TCanvas * canvas = new TCanvas("yieldcanvas");
    TCanvas * canvasfit = new TCanvas("yieldfitcanvas");
    TCanvas * canvasdiv = new TCanvas("dividedyieldcanvas");
    canvasdiv->Divide(2,4);
    
    TLegend * legend = new TLegend(0.56,0.7,0.9,0.9);
    TLegend * legendfit = new TLegend(0.56,0.7,0.9,0.9);

    Double_t bins[4] = {0.5,1.0,2.0,4.0};
    TH1D * near48 = new TH1D("nearsidey48","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c.",3,bins);
    near48->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    TH1D * near48fit = new TH1D("nearsidey48fit","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c.",3,bins);
    near48fit->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    TH1D * away48 = new TH1D("awaysidey48","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c.",3,bins);
    away48->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    TH1D * away48fit = new TH1D("awaysidey48fit","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c.",3,bins);
    away48fit->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    Double_t bins816[5] = {0.5,1.0,2.0,4.0,8.0};
    TH1D * near816 = new TH1D("nearsidey816","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c.",4,bins816);
    near816->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    TH1D * near816fit = new TH1D("nearsidey816fit","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c.",4,bins816);
    near816fit->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    TH1D * away816 = new TH1D("awaysidey816","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c.",4,bins816);
    away816->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    TH1D * away816fit = new TH1D("awaysidey816fit","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c.",4,bins816);
    away816fit->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12}} (rad)^{-1}");
    for(int i = 0;i<Bins->GetEntries();i++){
      ResHists * now = dynamic_cast<ResHists*>(Bins->At(i));
      const char* lable = Form("Trigger pT from %1.1f-%1.1f GeV/c and Associated pT from %1.1f-%1.1f GeV/c",now->GetMinTpT(),now->GetMaxTpT(),now->GetMinApT(),now->GetMaxApT());
      TH1D* hist = now->GetYield();
      TH1D* histfit = now->GetYieldfit();
      hist->SetTitle(lable);
      int padnumber = 0;

      if(now->GetMinApT()==0.5){hist->SetLineColor(2);histfit->SetLineColor(2);padnumber=1;}
      if(now->GetMinApT()==1.0){hist->SetLineColor(3);histfit->SetLineColor(3);padnumber=3;}
      if(now->GetMinApT()==2.0){hist->SetLineColor(4);histfit->SetLineColor(4);padnumber=5;}
      if(now->GetMinApT()==4.0){hist->SetLineColor(4);histfit->SetLineColor(4);padnumber=7;}

      hist->SetMarkerColor(1); histfit->SetMarkerColor(1);
      hist->SetMarkerSize(0.5); histfit->SetMarkerSize(1);
      if(now->GetMinTpT()==4.0){hist->SetMarkerStyle(20);histfit->SetMarkerStyle(20);}
      if(now->GetMinTpT()==8.0){hist->SetMarkerStyle(21);histfit->SetMarkerStyle(21);padnumber +=1;}
      //hist->GetYaxis()->SetRange(0.0,0.045);
      canvas->cd();
      hist->SetStats(kFALSE);histfit->SetStats(kFALSE);
      hist->Draw("Esame");
      legend->AddEntry(hist,lable);
      canvasfit->cd();
      histfit->Draw("Esame");
      legendfit->AddEntry(hist,lable);
      if(padnumber !=0){
	canvasdiv->cd(padnumber);
	hist->Draw("E");
      }
      Double_t errn=0.0;
      Double_t neary = hist->IntegralAndError(1,hist->FindBin(TMath::Pi()*0.5),errn,"width");
      Double_t erra=0.0;
      Double_t awayy = hist->IntegralAndError(hist->FindBin(TMath::Pi()*0.5),hist->GetNbinsX(),erra,"width");
      Double_t errnfit=0.0;
      Double_t nearyfit = histfit->IntegralAndError(1,histfit->FindBin(TMath::Pi()*0.5),errnfit,"width");
      Double_t errafit=0.0;
      Double_t awayyfit = histfit->IntegralAndError(histfit->FindBin(TMath::Pi()*0.5),histfit->GetNbinsX(),errafit,"width");
      if(now->GetMinTpT()==4.0){
	if(now->GetMinApT()==0.5){near48->SetBinContent(1,neary);near48->SetBinError(1,errn);away48->SetBinContent(1,awayy);away48->SetBinError(1,erra);near48fit->SetBinContent(1,nearyfit);near48fit->SetBinError(1,errnfit);away48fit->SetBinContent(1,awayyfit);away48fit->SetBinError(1,errafit);}
	if(now->GetMinApT()==1.0){near48->SetBinContent(2,neary);near48->SetBinError(2,errn);away48->SetBinContent(2,awayy);away48->SetBinError(2,erra);near48fit->SetBinContent(2,nearyfit);near48fit->SetBinError(2,errnfit);away48fit->SetBinContent(2,awayyfit);away48fit->SetBinError(2,errafit);}
	if(now->GetMinApT()==2.0){near48->SetBinContent(3,neary);near48->SetBinError(3,errn);away48->SetBinContent(3,awayy);away48->SetBinError(3,erra);near48fit->SetBinContent(3,nearyfit);near48fit->SetBinError(3,errnfit);away48fit->SetBinContent(3,awayyfit);away48fit->SetBinError(3,errafit);}
      }
      if(now->GetMinTpT()==8.0){
	if(now->GetMinApT()==0.5){near816->SetBinContent(1,neary);near816->SetBinError(1,errn);away816->SetBinContent(1,awayy);away816->SetBinError(1,erra);near816fit->SetBinContent(1,nearyfit);near816fit->SetBinError(1,errnfit);away816fit->SetBinContent(1,awayyfit);away816fit->SetBinError(1,errafit);}
	if(now->GetMinApT()==1.0){near816->SetBinContent(2,neary);near816->SetBinError(2,errn);away816->SetBinContent(2,awayy);away816->SetBinError(2,erra);near816fit->SetBinContent(2,nearyfit);near816fit->SetBinError(2,errnfit);away816fit->SetBinContent(2,awayyfit);away816fit->SetBinError(2,errafit);}
	if(now->GetMinApT()==2.0){near816->SetBinContent(3,neary);near816->SetBinError(3,errn);away816->SetBinContent(3,awayy);away816->SetBinError(3,erra);near816fit->SetBinContent(3,nearyfit);near816fit->SetBinError(3,errnfit);away816fit->SetBinContent(3,awayyfit);away816fit->SetBinError(3,errafit);}
	if(now->GetMinApT()==4.0){near816->SetBinContent(4,neary);near816->SetBinError(4,errn);away816->SetBinContent(4,awayy);away816->SetBinError(4,erra);near816fit->SetBinContent(4,nearyfit);near816fit->SetBinError(4,errnfit);away816fit->SetBinContent(4,awayyfit);away816fit->SetBinError(4,errafit);}
      }    
    }
    canvasdiv->Write();
    canvas->cd();
    legend->Draw("same");
    canvas->Write();
    canvasfit->cd();
    legendfit->Draw("same");
    canvasfit->Write();
    near48->Write();
    near48fit->Write();
    away48->Write();
    away48fit->Write();
    near816->Write();
    near816fit->Write();
    away816->Write();
    away816fit->Write();
    away48->Divide(near48);
    away48->SetTitle("Away side yield divided by near side yield for Trigger pT= 4.0-8.0 GeV/c");
    away48->GetYaxis()->SetTitle("");
    away48->Write("Ratio48");
    away48fit->Divide(near48fit);
    away48fit->SetTitle("Away side yield divided by near side yield for Trigger pT= 4.0-8.0 GeV/c");
    away48fit->GetYaxis()->SetTitle("");
    away48fit->Write("Ratio48fit");
    away816->Divide(near816);
    away816->SetTitle("Away side yield divided by near side yield for Trigger pT= 4.0-8.0 GeV/c");
    away816->GetYaxis()->SetTitle("");
    away816->Write("Ratio816");     
    away816fit->Divide(near816fit);
    away816fit->SetTitle("Away side yield divided by near side yield for Trigger pT= 4.0-8.0 GeV/c");
    away816fit->GetYaxis()->SetTitle("");
    away816fit->Write("Ratio816fit");     
    
  }
}
void Checkmixed(const char* option){
  cout << option<<endl;

  TFile * rfile = TFile::Open("results.root","UPDATE");
  TList * folderlist = rfile->GetListOfKeys();
  if(TString("pp").CompareTo(option)==0){
    for(int l = folderlist->GetEntries()-1;l>=0;l--){
      TString folder = TString(folderlist->At(l)->GetName());  
      if(folder.Contains("ThreePartTracks")){
	TObjArray * HistList = new TObjArray();
	TObjArray * CanvasList = new TObjArray();
	TObjArray * CanvasListVertex = new TObjArray();

	TH2D * collectedhist = NULL;
	TH2D * collectedhistV1 = NULL;
	TH2D * collectedhistV2 = NULL;
	TH2D * collectedhistV3 = NULL;
	TH2D * collectedhistV4 = NULL;
	TH2D * collectedhistV5 = NULL;

	TDirectory * foldernow = rfile->GetDirectory(folder.Data());
	TDirectory * savefolder = resultsdirectory(foldernow, "Mixedcheck");
	TList* subfolderlist = foldernow->GetListOfKeys();
	for(int k = 0; k <subfolderlist->GetEntries();k++){
	  if(TString(subfolderlist->At(k)->GetName()).Contains("Z")){
	    if(TString(subfolderlist->At(k)->GetName()).Contains("Z(-10.00)")){
	      TObjString * name = dynamic_cast<TObjString*>(TString(subfolderlist->At(k)->GetName()).Tokenize("Z")->At(0));
	      TCanvas * canvas = new TCanvas(name->GetString().Data());
	      canvas->Divide(2,3);
	      CanvasList->Add(canvas);
	    }
	    if(TString(subfolderlist->At(k)->GetName()).Contains("M(0.")){
	      TObjString * name = dynamic_cast<TObjString*>(TString(subfolderlist->At(k)->GetName()).Tokenize("Z")->At(1));
	      TCanvas * canvas = new TCanvas(name->GetString().Data());
	      canvas->Divide(2,3);
	      CanvasListVertex->Add(canvas);
	    }	    
	    TDirectory * MEdir = foldernow->GetDirectory(subfolderlist->At(k)->GetName())->GetDirectory("mixed_event");
	    TH2D * MEthisbin = dynamic_cast<TH2D*>(MEdir->Get("DPhi_1_DEta_12_SameSide"));
	    if(!collectedhist) collectedhist = dynamic_cast<TH2D*>(MEthisbin->Clone("MEsumoverall"));
	    else collectedhist->Add(MEthisbin);
	    if(TString(subfolderlist->At(k)->GetName()).Contains("Z(-10")&&!collectedhistV1) collectedhistV1 = dynamic_cast<TH2D*>(MEthisbin->Clone("MEsumV1"));
	    else if (TString(subfolderlist->At(k)->GetName()).Contains("Z(-10")&&collectedhistV1) collectedhistV1->Add(MEthisbin);
	    if(TString(subfolderlist->At(k)->GetName()).Contains("Z(-5")&&!collectedhistV2) collectedhistV2 = dynamic_cast<TH2D*>(MEthisbin->Clone("MEsumV2"));
	    else if (TString(subfolderlist->At(k)->GetName()).Contains("Z(-5")&&collectedhistV2) collectedhistV2->Add(MEthisbin);
	    if(TString(subfolderlist->At(k)->GetName()).Contains("Z(-2")&&!collectedhistV3) collectedhistV3 = dynamic_cast<TH2D*>(MEthisbin->Clone("MEsumV3"));
	    else if (TString(subfolderlist->At(k)->GetName()).Contains("Z(-2")&&collectedhistV3) collectedhistV3->Add(MEthisbin);	    
	    if(TString(subfolderlist->At(k)->GetName()).Contains("Z(2")&&!collectedhistV4) collectedhistV4 = dynamic_cast<TH2D*>(MEthisbin->Clone("MEsumV4"));
	    else if (TString(subfolderlist->At(k)->GetName()).Contains("Z(2")&&collectedhistV3) collectedhistV4->Add(MEthisbin);
	    if(TString(subfolderlist->At(k)->GetName()).Contains("Z(5")&&!collectedhistV5) collectedhistV5 = dynamic_cast<TH2D*>(MEthisbin->Clone("MEsumV5"));
	    else if (TString(subfolderlist->At(k)->GetName()).Contains("Z(5")&&collectedhistV3) collectedhistV5->Add(MEthisbin);	    
	    HistList->Add(dynamic_cast<TH2D*>(MEthisbin->Clone(Form("ME_%s",subfolderlist->At(k)->GetName()))));
	  }
	}
	collectedhist->Scale(1.0/collectedhist->GetBinContent(collectedhist->GetXaxis()->FindBin(0.0),collectedhist->GetYaxis()->FindBin(0.0)));
	collectedhistV1->Scale(1.0/collectedhistV1->GetBinContent(collectedhist->GetXaxis()->FindBin(0.0),collectedhist->GetYaxis()->FindBin(0.0)));
	collectedhistV2->Scale(1.0/collectedhistV2->GetBinContent(collectedhist->GetXaxis()->FindBin(0.0),collectedhist->GetYaxis()->FindBin(0.0)));
	collectedhistV3->Scale(1.0/collectedhistV3->GetBinContent(collectedhist->GetXaxis()->FindBin(0.0),collectedhist->GetYaxis()->FindBin(0.0)));
	collectedhistV4->Scale(1.0/collectedhistV4->GetBinContent(collectedhist->GetXaxis()->FindBin(0.0),collectedhist->GetYaxis()->FindBin(0.0)));
	collectedhistV5->Scale(1.0/collectedhistV5->GetBinContent(collectedhist->GetXaxis()->FindBin(0.0),collectedhist->GetYaxis()->FindBin(0.0)));

	savefolder->cd();

	  
	TCanvas * verteces = new TCanvas("Verteces");
	verteces->Divide(2,3);
	verteces->cd(1);
	collectedhistV1->Draw("colz");
	verteces->cd(2);
	collectedhistV2->Draw("colz");
	verteces->cd(3);
	collectedhistV3->Draw("colz");
	verteces->cd(4);
	collectedhistV4->Draw("colz");
	verteces->cd(5);
	collectedhistV5->Draw("colz");
	verteces->Write();delete verteces;
	
	TCanvas * vertecesdiv = new TCanvas("Verteces");
	vertecesdiv->Divide(2,3);
	vertecesdiv->cd(1);
	TH2D* collectedhistV1div =dynamic_cast<TH2D*>(collectedhistV1->Clone("v1"));
	collectedhistV1div->Divide(collectedhistV1);
	collectedhistV1div->Draw("colz");
	vertecesdiv->cd(2);
	TH2D* collectedhistV2div =dynamic_cast<TH2D*>(collectedhistV2->Clone("v2"));
	collectedhistV2div->Divide(collectedhistV1);
	collectedhistV2div->Draw("colz");
	vertecesdiv->cd(3);
	TH2D* collectedhistV3div =dynamic_cast<TH2D*>(collectedhistV3->Clone("v1"));
	collectedhistV3div->Divide(collectedhistV1);
	collectedhistV3div->Draw("colz");
	vertecesdiv->cd(4);
	TH2D* collectedhistV4div =dynamic_cast<TH2D*>(collectedhistV4->Clone("v1"));
	collectedhistV4div->Divide(collectedhistV1);
	collectedhistV4div->Draw("colz");
	vertecesdiv->cd(5);
	TH2D* collectedhistV5div =dynamic_cast<TH2D*>(collectedhistV5->Clone("v1"));
	collectedhistV5div->Divide(collectedhistV1);
	collectedhistV5div->Draw("colz");
	vertecesdiv->Write();delete vertecesdiv;
	
	for(int i=0;i<HistList->GetEntries();i++){
	  TH2D* MEhist = dynamic_cast<TH2D*>(HistList->At(i));
	  TH2D* MEhist2 = dynamic_cast<TH2D*>(MEhist->Clone(Form("%s_2",MEhist->GetName())));

	  MEhist->Divide(collectedhist);
	  MEhist->SetTitle(dynamic_cast<TObjString*>(TString(MEhist->GetName()).Tokenize("_")->At(1))->GetString().Data());
	  MEhist->Write();
	  MEhist->SetStats(false);
	  TCanvas * canvas;TCanvas * canvasvert;
	  if(TString(MEhist->GetName()).Contains("BinM(0.00)")){canvas = dynamic_cast<TCanvas*>(CanvasList->At(0));}
	  if(TString(MEhist->GetName()).Contains("BinM(44.44)"))canvas = dynamic_cast<TCanvas*>(CanvasList->At(1));
	  if(TString(MEhist->GetName()).Contains("BinM(88"))canvas = dynamic_cast<TCanvas*>(CanvasList->At(2));
	  if(TString(MEhist->GetName()).Contains("BinM(133"))canvas = dynamic_cast<TCanvas*>(CanvasList->At(3));
	  if(TString(MEhist->GetName()).Contains("BinM(177"))canvas = dynamic_cast<TCanvas*>(CanvasList->At(4));
	  
	  if(TString(MEhist->GetName()).Contains("Z(-1")){canvas->cd(1);canvasvert = dynamic_cast<TCanvas*>(CanvasListVertex->At(0));MEhist2->Divide(collectedhistV1);}
	  if(TString(MEhist->GetName()).Contains("Z(-5")){canvas->cd(2);canvasvert = dynamic_cast<TCanvas*>(CanvasListVertex->At(1));MEhist2->Divide(collectedhistV2);}
	  if(TString(MEhist->GetName()).Contains("Z(-2")){canvas->cd(3);canvasvert = dynamic_cast<TCanvas*>(CanvasListVertex->At(2));MEhist2->Divide(collectedhistV3);}
	  if(TString(MEhist->GetName()).Contains("Z(2.")){canvas->cd(4);canvasvert = dynamic_cast<TCanvas*>(CanvasListVertex->At(3));MEhist2->Divide(collectedhistV4);}
	  if(TString(MEhist->GetName()).Contains("Z(5.")){canvas->cd(5);canvasvert = dynamic_cast<TCanvas*>(CanvasListVertex->At(4));MEhist2->Divide(collectedhistV5);}
	  MEhist->Draw("colz");
	  if(TString(MEhist->GetName()).Contains("BinM(0.00)")){canvasvert->cd(1);}
	  if(TString(MEhist->GetName()).Contains("BinM(44.44)")){canvasvert->cd(2);}
	  if(TString(MEhist->GetName()).Contains("BinM(88")){canvasvert->cd(3);}
	  if(TString(MEhist->GetName()).Contains("BinM(133")){canvasvert->cd(4);}
	  if(TString(MEhist->GetName()).Contains("BinM(177")){canvasvert->cd(5);}	  
	  MEhist2->Draw("colz");

	  
	}
	collectedhist->Write();
	for(int i = 0;i<CanvasList->GetEntries();i++){dynamic_cast<TCanvas*>(CanvasList->At(i))->Write();delete dynamic_cast<TCanvas*>(CanvasList->At(i));}
	for(int i = 0;i<CanvasListVertex->GetEntries();i++){dynamic_cast<TCanvas*>(CanvasListVertex->At(i))->Write();delete dynamic_cast<TCanvas*>(CanvasListVertex->At(i));}

	
      }
    }
  }
}

TCanvas * Periodcanvas(TH1 * hist, int ncol){
  TCanvas * Canvas = new TCanvas(Form("%sCanvas",hist->GetName()));
//   TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,hist->GetTitle());
//   title->Draw();
//   TPad* graphPad = new TPad("Graphs","Graphs",0.01,0.05,0.95,0.95);
//   graphPad->Draw();
//   graphPad->cd();
  Canvas->Divide(ncol,ncol);
  return Canvas;
}

TObjArray * PadArray(TCanvas * canvas, int ntrig, int nass,TObjArray* trignames, TObjArray* assnames, const char* title){
  canvas->cd();
  if(ntrig!=trignames->GetEntries()||nass!=assnames->GetEntries()){ cout << "wrong dimensions on the names"<<endl; return 0x0;}
  TObjArray * arrayofpads = new TObjArray();
  
  TPaveLabel* overtitle = new TPaveLabel(0.15,0.95,0.9,0.99,title,"nb");
  overtitle->Draw();
  TPaveLabel* asstitle = new TPaveLabel(0.0,0.8,0.15,0.99,"p_T^{ass}:","nb");
  asstitle->Draw();
  
  Double_t xdiv = 0.75/ntrig;
  Double_t ydiv = 0.85/nass;
  
  
  //make trigger lables:
  for(int nx = 1;nx<=ntrig;nx++){
    canvas->cd();
    Double_t xmin = 0.15+(nx-1)*xdiv;
    Double_t xmax = xmin + xdiv;
    Double_t ymin = 0.85;
    Double_t ymax = 0.95;
    TPaveLabel* Triglable = new TPaveLabel(xmin,ymin,xmax,ymax,trignames->At(nx-1)->GetName(),"nb");
    Triglable->Draw();
  }
  //so the lables for associated:
  for(int ny = 1;ny<=nass;ny++){
    canvas->cd();
    Double_t xmin = 0.0;
    Double_t xmax = 0.15;
    Double_t ymax = 0.85-(ny-1)*ydiv;
    Double_t ymin = ymax-ydiv;
    TPaveText* Asslable = new TPaveText(xmin,ymin,xmax,ymax);
    Asslable->SetOption("nb");
    TString string = dynamic_cast<TObjString*>( assnames->At(ny-1))->GetString();
    Asslable->AddText(string.Data());
//     Asslable->AddText(dynamic_cast<TObjString*>(string.Tokenize("NL")->At(1))->GetString().Data());
//     Asslable->AddText(dynamic_cast<TObjString*>(string.Tokenize("NL")->At(2))->GetString().Data());
    Asslable->Draw();
  }
  
  for(int nx = 1;nx<=ntrig;nx++){for(int ny=1;ny<=nass;ny++){
    canvas->cd();
    Double_t xmin = 0.15+(nx-1)*xdiv;
    Double_t xmax = xmin + xdiv;
    Double_t ymax = 0.85-(ny-1)*ydiv;
    Double_t ymin = ymax-ydiv;
    TPad * thispad = new TPad(Form("padno_%i_%i",nx,ny),"",xmin,ymin,xmax,ymax);
    arrayofpads->Add(thispad);
    thispad->Draw();
  }}
  
  return arrayofpads;
  
  
}

TVirtualPad * cdppcanvas(TCanvas * canvas, int i){
  TPad * pad = dynamic_cast<TPad*>(canvas->cd(i));
  return pad;//->cd(i);
}

void Periods(){
  //go over periods and collect the relevant histograms in a new results file.
  TFile * rfile = TFile::Open("results.root","UPDATE");
  TString * path = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  TString * path10b = new TString("LHC10b/Train");
  TString * path10c = new TString("LHC10c/Train");
  TString * path10d = new TString("LHC10d/Traineff");
  TString * path10e = new TString("LHC10e/Train");
  TString * path11a = new TString("LHC11a/Train");
  TString * path10h = new TString("LHC10h/TrainEff");
  TString * path10h2 = new TString("LHC10h/TrainEff/highpt");
  
  TDirectory * basedir = resultsdirectory(rfile->GetDirectory("/"),"Periods");
  TObjArray * filearray = new TObjArray();
//   filearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10b->Data()),"READ"));
//   filearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10c->Data()),"READ"));
  filearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10d->Data()),"READ"));
//   filearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10e->Data()),"READ"));
  filearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path11a->Data()),"READ"));

  TObjArray * PbPbfilearray = new TObjArray();
  PbPbfilearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10h->Data()),"READ"));
  PbPbfilearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10h2->Data()),"READ"));

  
  TCanvas *  testcanvas = new TCanvas("YieldsinPeriods");
//   TCanvas * corryields10b = new TCanvas("10bCorrYields");
//   TCanvas * corryields10c = new TCanvas("10cCorrYields");
  TCanvas * corryields10d = new TCanvas("10dCorrYields");
//   TCanvas * corryields10e = new TCanvas("10eCorrYields");
  TCanvas * corryields11a = new TCanvas("11aCorrYields");
//   TCanvas * corryields10h = new TCanvas("10hCorrYields");

  //For PbPb: in centrality bins:
  TCanvas * corryields10h05 = new TCanvas("10h05CorrYields");
  TCanvas * corryields10h510 = new TCanvas("10h510CorrYields");
  TCanvas * corryields10h1020 = new TCanvas("10h1020CorrYields");
  TCanvas * corryields10h2040 = new TCanvas("10h2040CorrYields");
  TCanvas * corryields10h4060 = new TCanvas("10h4060CorrYields");
  TCanvas * corryields10h6080 = new TCanvas("10h6080CorrYields");

//   TCanvas * corryieldsphi10b = new TCanvas("10bCorrYieldsphi");
//   TCanvas * corryieldsphi10c = new TCanvas("10cCorrYieldsphi");
  TCanvas * corryieldsphi10d = new TCanvas("10dCorrYieldsphi");
//   TCanvas * corryieldsphi10e = new TCanvas("10eCorrYieldsphi");
  TCanvas * corryieldsphi11a = new TCanvas("11aCorrYieldsphi");
//   TCanvas * corryieldsphi10h = new TCanvas("10hCorrYieldsphi");

  //For PbPb: in centrality bins:
  TCanvas * corryieldsphi10h05 = new TCanvas("10h05CorrYieldsphi");
  TCanvas * corryieldsphi10h510 = new TCanvas("10h510CorrYieldsphi");
  TCanvas * corryieldsphi10h1020 = new TCanvas("10h1020CorrYieldsphi");
  TCanvas * corryieldsphi10h2040 = new TCanvas("10h2040CorrYieldsphi");
  TCanvas * corryieldsphi10h4060 = new TCanvas("10h4060CorrYieldsphi");
  TCanvas * corryieldsphi10h6080 = new TCanvas("10h6080CorrYieldsphi");
  
//   TCanvas * yieldsphi10b = new TCanvas("10bCorrYieldsphi");
//   TCanvas * yieldsphi10c = new TCanvas("10cCorrYieldsphi");
  TCanvas * yieldsphi10d = new TCanvas("10dYieldsphi");
//   TCanvas * yieldsphi10e = new TCanvas("10eCorrYieldsphi");
  TCanvas * yieldsphi11a = new TCanvas("11aYieldsphi");
//   TCanvas * corryieldsphi10h = new TCanvas("10hCorrYieldsphi");

  //For PbPb: in centrality bins:
  TCanvas * yieldsphi10h05 = new TCanvas("10h05Yieldsphi");
  TCanvas * yieldsphi10h510 = new TCanvas("10h510Yieldsphi");
  TCanvas * yieldsphi10h1020 = new TCanvas("10h1020Yieldsphi");
  TCanvas * yieldsphi10h2040 = new TCanvas("10h2040Yieldsphi");
  TCanvas * yieldsphi10h4060 = new TCanvas("10h4060Yieldsphi");
  TCanvas * yieldsphi10h6080 = new TCanvas("10h6080Yieldsphi");
  
  //For PbPb: Yields for a given pT bin:
  TCanvas* yields816_816_10h = new TCanvas("10hYields816816");
  TLegend* leg816_816 = new TLegend(0.1,0.7,0.48,0.9);
  TCanvas* yields816_68_10h = new TCanvas("10hYields81668");
  TLegend* leg816_68 = new TLegend(0.1,0.7,0.48,0.9);
  TCanvas* yields48_68_10h = new TCanvas("10hYields4868");
  TLegend* leg48_68 = new TLegend(0.1,0.7,0.48,0.9);
  
  TObjArray * trigarray = new TObjArray();
  trigarray->Add(new TObjString("4 GeV/c <= p_{T}^{trigger}<=8 GeV/c"));
  trigarray->Add(new TObjString("8 GeV/c <= p_{T}^{trigger}<=16 GeV/c"));
  TObjArray * assarray  = new TObjArray();
  assarray->Add(new TObjString("0.5-1.0 GeV/c"));
  assarray->Add(new TObjString("1.0-2.0 GeV/c"));
  assarray->Add(new TObjString("2.0-3.0 GeV/c"));
  assarray->Add(new TObjString("3.0-4.0 GeV/c"));
  assarray->Add(new TObjString("4.0-6.0 GeV/c"));
  assarray->Add(new TObjString("6.0-8.0 GeV/c"));
  assarray->Add(new TObjString("8.0-16.0 GeV/c"));
  
  TObjArray* arrayofpads = PadArray(testcanvas,2,7,trigarray,assarray,"Yields from bin counting in different pT ranges");
//   TObjArray* arrayofpadscor10b = PadArray(corryields10b,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10b in different pT ranges");
//   TObjArray* arrayofpadscor10c = PadArray(corryields10c,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10c in different pT ranges");
  TObjArray* arrayofpadscor10d = PadArray(corryields10d,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10d in different pT ranges");
//   TObjArray* arrayofpadscor10e = PadArray(corryields10e,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10e in different pT ranges");  
  TObjArray* arrayofpadscor11a = PadArray(corryields11a,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 11a in different pT ranges");  
//   TObjArray* arrayofpadscor10h = PadArray(corryields10h,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h in different pT ranges");  

  TObjArray* arrayofpadscor10h05 = PadArray(corryields10h05,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 0%-5% events in different pT ranges");
  TObjArray* arrayofpadscor10h510 = PadArray(corryields10h510,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 5%-10% events in different pT ranges");
  TObjArray* arrayofpadscor10h1020 = PadArray(corryields10h1020,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 10%-20% events in different pT ranges");
  TObjArray* arrayofpadscor10h2040 = PadArray(corryields10h2040,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 20%-40% events in different pT ranges");
  TObjArray* arrayofpadscor10h4060 = PadArray(corryields10h4060,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 40%-60% events in different pT ranges");
  TObjArray* arrayofpadscor10h6080 = PadArray(corryields10h6080,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 60%-80% events in different pT ranges");

  
//   TObjArray* arrayofpadscorphi10b = PadArray(corryieldsphi10b,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10b in different pT ranges");
//   TObjArray* arrayofpadscorphi10c = PadArray(corryieldsphi10c,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10c in different pT ranges");
  TObjArray* arrayofpadscorphi10d = PadArray(corryieldsphi10d,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10d in different pT ranges");
//   TObjArray* arrayofpadscorphi10e = PadArray(corryieldsphi10e,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10e in different pT ranges"); 
  TObjArray* arrayofpadscorphi11a = PadArray(corryieldsphi11a,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 11a in different pT ranges"); 
//   TObjArray* arrayofpadscorphi10h = PadArray(corryieldsphi10h,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10h in different pT ranges"); 

  TObjArray* arrayofpadscorphi10h05 = PadArray(corryieldsphi10h05,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10h 0%-5% events in different pT ranges");
  TObjArray* arrayofpadscorphi10h510 = PadArray(corryieldsphi10h510,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10h 5%-10% events in different pT ranges");
  TObjArray* arrayofpadscorphi10h1020 = PadArray(corryieldsphi10h1020,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10h 10%-20% events in different pT ranges");
  TObjArray* arrayofpadscorphi10h2040 = PadArray(corryieldsphi10h2040,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10h 20%-40% events in different pT ranges");
  TObjArray* arrayofpadscorphi10h4060 = PadArray(corryieldsphi10h4060,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10h 40%-60% events in different pT ranges");
  TObjArray* arrayofpadscorphi10h6080 = PadArray(corryieldsphi10h6080,2,7,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#Phi_{2} in 10h 60%-80% events in different pT ranges");


  
//   TObjArray* arrayofpadsyield10b = PadArray(corryieldsphi10b,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10b in different pT ranges");
//   TObjArray* arrayofpadsyield10c = PadArray(corryieldsphi10c,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10c in different pT ranges");
  TObjArray* arrayofpadsyield10d = PadArray(yieldsphi10d,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10d in different pT ranges");
//   TObjArray* arrayofpadsyield10e = PadArray(corryieldsphi10e,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10e in different pT ranges"); 
  TObjArray* arrayofpadsyield11a = PadArray(yieldsphi11a,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1}  in 11a in different pT ranges"); 
//   TObjArray* arrayofpadsyield10h = PadArray(corryieldsphi10h,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10h in different pT ranges"); 

  TObjArray* arrayofpadsyield10h05 = PadArray(yieldsphi10h05,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10h 0%-5% events in different pT ranges");
  TObjArray* arrayofpadsyield10h510 = PadArray(yieldsphi10h510,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10h 5%-10% events in different pT ranges");
  TObjArray* arrayofpadsyield10h1020 = PadArray(yieldsphi10h1020,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10h 10%-20% events in different pT ranges");
  TObjArray* arrayofpadsyield10h2040 = PadArray(yieldsphi10h2040,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10h 20%-40% events in different pT ranges");
  TObjArray* arrayofpadsyield10h4060 = PadArray(yieldsphi10h4060,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10h 40%-60% events in different pT ranges");
  TObjArray* arrayofpadsyield10h6080 = PadArray(yieldsphi10h6080,2,7,trigarray,assarray,"Yield vs #Delta#Phi_{1} in 10h 60%-80% events in different pT ranges");

  
  
  //get the relevant bins:
  TObjArray * binarray= new TObjArray();
  TObjArray * Resultsbinarray = new TObjArray();
  TList * binlist = dynamic_cast<TDirectory*>(filearray->At(0))->GetListOfKeys();
  for(int i = 0;i<binlist->GetEntries();i++){
    TObjString * bin = new TObjString(binlist->At(i)->GetName());
    if(bin->GetString().Contains("ThreePartTracks")){
      binarray->Add(bin);
      Resultsbinarray->Add(resultsdirectory(basedir,bin->GetString().Data()));
    }
    else delete bin;
  }

  //get the relevant bins:
  TObjArray * PbPbbinarray= new TObjArray();
  TList * PbPbbinlist = dynamic_cast<TDirectory*>(PbPbfilearray->At(0))->GetListOfKeys();
  for(int i = 0;i<PbPbbinlist->GetEntries();i++){
    TObjString * bin = new TObjString(PbPbbinlist->At(i)->GetName());
    if(bin->GetString().Contains("ThreePartTracks")){
      PbPbbinarray->Add(bin);
    }
    else delete bin;
  }
  TList * PbPbbinlist2 = dynamic_cast<TDirectory*>(PbPbfilearray->At(1))->GetListOfKeys();
  for(int i = 0;i<PbPbbinlist2->GetEntries();i++){
    TObjString * bin = new TObjString(PbPbbinlist2->At(i)->GetName());
    if(bin->GetString().Contains("ThreePartTracks")){
      PbPbbinarray->Add(bin);
    }
    else delete bin;
  }
  
  //get the relevant histograms:
  TObjArray * histlist= new TObjArray();
  TObjArray * histposlist= new TObjArray();
  {
    TObjString * iteration = new TObjString("iteration1");
    TObjString * yield = new TObjString("yield/original_binning/GP0");
    TObjString * dphi1dphi2 = new TObjString("DPhi_1_DPHI_2");
    histlist->Add(dphi1dphi2);
    histposlist->Add(iteration);
    TObjString * dphi1detass = new TObjString("DPhi_1_DEta_12_SameSide");
    histlist->Add(dphi1detass);
    histposlist->Add(iteration);
    TObjString * yieldh = new TObjString("dphiyieldbc");
    histlist->Add(yieldh);
    histposlist->Add(yield);
  }
  
  //get to the relevant centrality bin:
  TObjArray * cbinarray= new TObjArray();
  TList * cbinlist = dynamic_cast<TDirectory*>(PbPbfilearray->At(0))->GetDirectory(PbPbbinarray->At(0)->GetName())->GetListOfKeys();
  for(int i = 0;i<cbinlist->GetEntries();i++){
    TObjString * bin = new TObjString(cbinlist->At(i)->GetName());
    if(bin->GetString().Contains("BinM")&&!(bin->GetString().Contains("Z"))){
      cbinarray->Add(bin);
    }
    else delete bin;
  }
  
  
  TH1D * yieldbcperiodcollected = dynamic_cast<TH1D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(dynamic_cast<TObjString*>(binarray->At(0))->GetString().Data())->GetDirectory("yield/original_binning/GP0")->Get("dphiyieldbc"));
  yieldbcperiodcollected->Reset();
  
  Double_t bins[7] = {0.5,1.0,2.0,3.0,4.0,6.0,8.0};
  TH1D * near48 = new TH1D("nearsidey48","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c.",6,bins);
  near48->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away48 = new TH1D("awaysidey48","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c.",6,bins);
  away48->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  Double_t bins816[8] = {0.5,1.0,2.0,3.0,4.0,6.0,8.0,16.0};
  TH1D * near816 = new TH1D("nearsidey816","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c.",7,bins816);
  near816->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away816 = new TH1D("awaysidey816","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c.",7,bins816);
  away816->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  
  TH1D * near48_10h_05 = new TH1D("nearsidey48_10h_05","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c in 0% - 5% PbPb events.",6,bins);
  near48_10h_05->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away48_10h_05 = new TH1D("awaysidey48_10h_05","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c in 0% - 5% PbPb events.",6,bins);
  away48_10h_05->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near816_10h_05 = new TH1D("nearsidey816_10h_05","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c in 0% - 5% PbPb events.",7,bins816);
  near816_10h_05->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away816_10h_05 = new TH1D("awaysidey816_10h_05","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c in 0% - 5% PbPb events.",7,bins816);
  away816_10h_05->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near48_10h_510 = new TH1D("nearsidey48_10h_510","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c in 5% - 10% PbPb events.",6,bins);
  near48_10h_510->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away48_10h_510 = new TH1D("awaysidey48_10h_510","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c in 5% - 10% PbPb events.",6,bins);
  away48_10h_510->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near816_10h_510 = new TH1D("nearsidey816_10h_510","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c in 5% - 10% PbPb events.",7,bins816);
  near816_10h_510->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away816_10h_510 = new TH1D("awaysidey816_10h_510","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c in 5% - 10% PbPb events.",7,bins816);
  away816_10h_510->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
   
  TH1D * near48_10h_1020 = new TH1D("nearsidey48_10h_1020","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c in 10% - 20% PbPb events.",6,bins);
  near48_10h_1020->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away48_10h_1020 = new TH1D("awaysidey48_10h_1020","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c in 10% - 20% PbPb events.",6,bins);
  away48_10h_1020->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near816_10h_1020 = new TH1D("nearsidey816_10h_1020","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c in 10% - 20% PbPb events.",7,bins816);
  near816_10h_1020->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away816_10h_1020 = new TH1D("awaysidey816_10h_1020","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c in 10% - 20% PbPb events.",7,bins816);
  away816_10h_1020->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near48_10h_2040 = new TH1D("nearsidey48_10h_2040","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c in 20% - 40% PbPb events.",6,bins);
  near48_10h_2040->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away48_10h_2040 = new TH1D("awaysidey48_10h_2040","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c in 20% - 40% PbPb events.",6,bins);
  away48_10h_2040->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near816_10h_2040 = new TH1D("nearsidey816_10h_2040","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c in 20% - 40% PbPb events.",7,bins816);
  near816_10h_2040->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away816_10h_2040 = new TH1D("awaysidey816_10h_2040","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c in 20% - 40% PbPb events.",7,bins816);
  away816_10h_2040->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near48_10h_4060 = new TH1D("nearsidey48_10h_4060","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c in 40% - 60% PbPb events.",6,bins);
  near48_10h_4060->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away48_10h_4060 = new TH1D("awaysidey48_10h_4060","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c in 40% - 60% PbPb events.",6,bins);
  away48_10h_4060->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * near816_10h_4060 = new TH1D("nearsidey816_10h_4060","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c in 40% - 60% PbPb events.",7,bins816);
  near816_10h_4060->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  TH1D * away816_10h_4060 = new TH1D("awaysidey816_10h_4060","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c in 40% - 60% PbPb events.",7,bins816);
  away816_10h_4060->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
  
//   TCanvas * yieldp = new TCanvas("YieldsinPeriods");
//   yieldp->Divide(2,4);
  
  TObjArray * CanvasList = new TObjArray();
  for(int i = 0;i<histlist->GetEntries();i++){
    TDirectory * dir = dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(dynamic_cast<TObjString*>(binarray->At(0))->GetString().Data())->GetDirectory(dynamic_cast<TObjString*>(histposlist->At(i))->GetString().Data());
    TH1 * hist = dynamic_cast<TH1*>(dir->Get(dynamic_cast<TObjString*>(histlist->At(i))->GetString().Data()));
    TCanvas * can = Periodcanvas(hist,2);
    CanvasList->Add(can);
  }
  for(int i = 0;i<binarray->GetEntries();i++){
    TDirectory * savedir = dynamic_cast<TDirectory*>(Resultsbinarray->At(i));
    for(int x = 0;x<=yieldbcperiodcollected->GetNbinsX();x++){
      Double_t S = 0.0;
      Double_t error = 0.0;
      for(int k = 0;k<filearray->GetEntries();k++){
	TH1D * hist = dynamic_cast<TH1D*>(dynamic_cast<TFile*>(filearray->At(k))->GetDirectory(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Data())->GetDirectory("yield/original_binning/GP0")->Get("dphiyieldbc"));

	  
	Double_t loccont = hist->GetBinContent(x);
	Double_t locerror = hist->GetBinError(x);
	if(TMath::Abs(locerror)>1.0E-12){
	  S += loccont/(locerror*locerror);
	  error += 1.0/(locerror*locerror);
	}
      }
      if(error>1.0E-12){
	yieldbcperiodcollected->SetBinContent(x,S/error);
	yieldbcperiodcollected->SetBinError(x,TMath::Sqrt(1.0/error));
      }
    }
    savedir->cd();
    yieldbcperiodcollected->SetTitle(Form("Yield for %s",binarray->At(i)->GetName()));
    yieldbcperiodcollected->Write();
    Double_t near = 0.0; Double_t neare=0.0;Double_t away = 0.0; Double_t awaye=0.0;
    near = yieldbcperiodcollected->IntegralAndError(1,yieldbcperiodcollected->FindBin(TMath::Pi()/2.0),neare,"width");
    away = yieldbcperiodcollected->IntegralAndError(yieldbcperiodcollected->FindBin(TMath::Pi()/2.0),yieldbcperiodcollected->GetNbinsX(),awaye,"width");

    if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_0_1")){near48->SetBinContent(1,near/0.5);near48->SetBinError(1,neare/0.5);away48->SetBinContent(1,away/0.5);away48->SetBinError(1,awaye/0.5);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_1_2")){near48->SetBinContent(2,near);near48->SetBinError(2,neare);away48->SetBinContent(2,away);away48->SetBinError(2,awaye);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_2_3")){near48->SetBinContent(3,near/1.0);near48->SetBinError(3,neare/1.0);away48->SetBinContent(3,away/1.0);away48->SetBinError(3,awaye/1.0);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_3_4")){near48->SetBinContent(4,near/1.0);near48->SetBinError(4,neare/1.0);away48->SetBinContent(4,away/1.0);away48->SetBinError(4,awaye/1.0);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_4_6")){near48->SetBinContent(5,near/2.0);near48->SetBinError(5,neare/2.0);away48->SetBinContent(5,away/2.0);away48->SetBinError(5,awaye/2.0);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_6_8")){near48->SetBinContent(6,near/2.0);near48->SetBinError(6,neare/2.0);away48->SetBinContent(6,away/2.0);away48->SetBinError(6,awaye/2.0);}

    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_0_1")){near816->SetBinContent(1,near/0.5);near816->SetBinError(1,neare/0.5);away816->SetBinContent(1,away/0.5);away816->SetBinError(1,awaye/0.5);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_1_2")){near816->SetBinContent(2,near);near816->SetBinError(2,neare);away816->SetBinContent(2,away);away816->SetBinError(2,awaye);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_2_3")){near816->SetBinContent(3,near/1.0);near816->SetBinError(3,neare/1.0);away816->SetBinContent(3,away/1.0);away816->SetBinError(3,awaye/1.0);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_3_4")){near816->SetBinContent(4,near/1.0);near816->SetBinError(4,neare/1.0);away816->SetBinContent(4,away/1.0);away816->SetBinError(4,awaye/1.0);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_4_6")){near816->SetBinContent(5,near/2.0);near816->SetBinError(5,neare/2.0);away816->SetBinContent(5,away/2.0);away816->SetBinError(5,awaye/2.0);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_6_8")){near816->SetBinContent(6,near/2.0);near816->SetBinError(6,neare/2.0);away816->SetBinContent(6,away/2.0);away816->SetBinError(6,awaye/2.0);}
    else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_8_16")){near816->SetBinContent(7,near/8.0);near816->SetBinError(7,neare/8.0);away816->SetBinContent(7,away/8.0);away816->SetBinError(7,awaye/8.0);}


    for(int k = 0;k<filearray->GetEntries();k++){
      TObjArray* arr;
      TH1D * hist = dynamic_cast<TH1D*>(dynamic_cast<TFile*>(filearray->At(k))->GetDirectory(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Data())->GetDirectory("yield/original_binning/GP0")->Get("dphiyieldbc"));
      if(k==1) arr = arrayofpadsyield11a;
// 	if(k==1) arr = arrayofpadscor10c;
      if(k==0) arr = arrayofpadsyield10d;
// 	if(k==3) arr = arrayofpadscor10e;
      if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_0_1")){dynamic_cast<TPad*>(arr->At(0))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_1_2")){dynamic_cast<TPad*>(arr->At(1))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_2_3")){dynamic_cast<TPad*>(arr->At(2))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_3_4")){dynamic_cast<TPad*>(arr->At(3))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_4_6")){dynamic_cast<TPad*>(arr->At(4))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_6_8")){dynamic_cast<TPad*>(arr->At(5))->cd();}

      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_0_1")){dynamic_cast<TPad*>(arr->At(7))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_1_2")){dynamic_cast<TPad*>(arr->At(8))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_2_3")){dynamic_cast<TPad*>(arr->At(9))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_3_4")){dynamic_cast<TPad*>(arr->At(10))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_4_6")){dynamic_cast<TPad*>(arr->At(11))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_6_8")){dynamic_cast<TPad*>(arr->At(12))->cd();}
      else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_8_16")){dynamic_cast<TPad*>(arr->At(13))->cd();}
      hist->SetStats(false);
      hist->GetXaxis()->SetLabelSize(0.08);
      hist->GetXaxis()->SetTitleSize(0.08);
      hist->GetXaxis()->SetTitleOffset(0.5);
      hist->GetYaxis()->SetLabelSize(0.08);
      hist->GetYaxis()->SetTitleSize(0.065);
      hist->GetYaxis()->SetTitleOffset(0.5);
      hist->SetTitle("");
      hist->Draw("E");
      
    }
    
    
    for(int j = 0;j<histlist->GetEntries();j++){
      TCanvas * canvnow = dynamic_cast<TCanvas*>(CanvasList->At(j));
      for(int k = 0;k<filearray->GetEntries();k++){
	TH1 * hist = dynamic_cast<TH1*>(dynamic_cast<TFile*>(filearray->At(k))->GetDirectory(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Data())->GetDirectory(dynamic_cast<TObjString*>(histposlist->At(j))->GetString().Data())->Get(histlist->At(j)->GetName()));
	
	
	

	TH1D * hist1d = dynamic_cast<TH1D*>(hist);
	TH2D * hist2d = dynamic_cast<TH2D*>(hist);
	
	if(TString(hist->GetName()).Contains("DPhi_1_DEta_12_SameSide")){
	  TObjArray* arr;
	  if(k==1) arr = arrayofpadscor11a;
// 	  if(k==1) arr = arrayofpadscor10c;
	  if(k==0) arr = arrayofpadscor10d;
// 	  if(k==3) arr = arrayofpadscor10e;
	  if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_0_1")){dynamic_cast<TPad*>(arr->At(0))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_1_2")){dynamic_cast<TPad*>(arr->At(1))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_2_3")){dynamic_cast<TPad*>(arr->At(2))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_3_4")){dynamic_cast<TPad*>(arr->At(3))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_4_6")){dynamic_cast<TPad*>(arr->At(4))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_6_8")){dynamic_cast<TPad*>(arr->At(5))->cd();}

	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_0_1")){dynamic_cast<TPad*>(arr->At(7))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_1_2")){dynamic_cast<TPad*>(arr->At(8))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_2_3")){dynamic_cast<TPad*>(arr->At(9))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_3_4")){dynamic_cast<TPad*>(arr->At(10))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_4_6")){dynamic_cast<TPad*>(arr->At(11))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_6_8")){dynamic_cast<TPad*>(arr->At(12))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_8_16")){dynamic_cast<TPad*>(arr->At(13))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_8_16")){dynamic_cast<TPad*>(arr->At(13))->cd();}
	  gStyle->SetPalette(1);
	  hist2d->SetStats(false);
	  hist2d->GetXaxis()->SetLabelSize(0.08);
	  hist2d->GetXaxis()->SetTitleSize(0.09);
	  hist2d->GetXaxis()->SetTitleOffset(0.3);
	  hist2d->GetYaxis()->SetLabelSize(0.08);
	  hist2d->GetYaxis()->SetTitleSize(0.09);
	  hist2d->GetYaxis()->SetTitleOffset(0.2);
	  hist2d->SetTitle("");
	  hist2d->Draw("COLZ");
	  gPad->Update();
	  TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
	  if(palette) palette->SetLabelSize(0.08);
	}
	if(TString(hist->GetName()).Contains("DPhi_1_DPHI_2")){
	  TObjArray* arr;
	  if(k==1) arr = arrayofpadscorphi11a;
// 	  if(k==1) arr = arrayofpadscorphi10c;
	  if(k==0) arr = arrayofpadscorphi10d;
// 	  if(k==3) arr = arrayofpadscorphi10e;
	  if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_0_1")){dynamic_cast<TPad*>(arr->At(0))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_1_2")){dynamic_cast<TPad*>(arr->At(1))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_2_3")){dynamic_cast<TPad*>(arr->At(2))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_3_4")){dynamic_cast<TPad*>(arr->At(3))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_4_6")){dynamic_cast<TPad*>(arr->At(4))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("4_8_6_8")){dynamic_cast<TPad*>(arr->At(5))->cd();}

	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_0_1")){dynamic_cast<TPad*>(arr->At(7))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_1_2")){dynamic_cast<TPad*>(arr->At(8))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_2_3")){dynamic_cast<TPad*>(arr->At(9))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_3_4")){dynamic_cast<TPad*>(arr->At(10))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_4_6")){dynamic_cast<TPad*>(arr->At(11))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_6_8")){dynamic_cast<TPad*>(arr->At(12))->cd();}
	  else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_8_16")){dynamic_cast<TPad*>(arr->At(13))->cd();}
	  hist2d->SetStats(false);
	  hist2d->GetXaxis()->SetLabelSize(0.08);
	  hist2d->GetXaxis()->SetTitleSize(0.09);
	  hist2d->GetXaxis()->SetTitleOffset(0.3);
	  hist2d->GetYaxis()->SetLabelSize(0.08);
	  hist2d->GetYaxis()->SetTitleSize(0.09);
	  hist2d->GetYaxis()->SetTitleOffset(0.2);
	  hist2d->SetTitle("");
	  hist2d->Draw("colz");
	  gPad->Update();
	  TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
	  if(palette) palette->SetLabelSize(0.08);
	}
	
	cdppcanvas(canvnow,k+1)->cd();
	
// 	if(k==1) hist->SetTitle("LHC11a");
// 	if(k==1) hist->SetTitle("LHC10c");
// 	if(k==0) hist->SetTitle("LHC10d");
// 	if(k==3) hist->SetTitle("LHC10e");
	
	if(hist1d) hist1d->Draw("E");
	if(hist2d) hist2d->Draw("colz");

	
	
      }
      savedir->cd();
      canvnow->Write();
    }
  }
  
  for(int i=0;i<binarray->GetEntries();i++){
    TDirectory * savedir = dynamic_cast<TDirectory*>(Resultsbinarray->At(i));
    TH1D * hist = dynamic_cast<TH1D*>(savedir->Get("dphiyieldbc")->Clone(Form("dphiyieldbcclonenr%i",i)));
    if(TString(binarray->At(i)->GetName()).Contains("4_8_0_1"))dynamic_cast<TPad*>(arrayofpads->At(0))->cd();//yieldp->cd(1);
    if(TString(binarray->At(i)->GetName()).Contains("4_8_1_2"))dynamic_cast<TPad*>(arrayofpads->At(1))->cd();//yieldp->cd(3);
    if(TString(binarray->At(i)->GetName()).Contains("4_8_2_3"))dynamic_cast<TPad*>(arrayofpads->At(2))->cd();//yieldp->cd(5);
    if(TString(binarray->At(i)->GetName()).Contains("4_8_3_4"))dynamic_cast<TPad*>(arrayofpads->At(3))->cd();//yieldp->cd(5);
    if(TString(binarray->At(i)->GetName()).Contains("4_8_4_6"))dynamic_cast<TPad*>(arrayofpads->At(4))->cd();//yieldp->cd(5);
    if(TString(binarray->At(i)->GetName()).Contains("4_8_6_8"))dynamic_cast<TPad*>(arrayofpads->At(5))->cd();//yieldp->cd(5);   
    
    if(TString(binarray->At(i)->GetName()).Contains("8_16_0_1"))dynamic_cast<TPad*>(arrayofpads->At(7))->cd();//yieldp->cd(2);
    if(TString(binarray->At(i)->GetName()).Contains("8_16_1_2"))dynamic_cast<TPad*>(arrayofpads->At(8))->cd();//yieldp->cd(4);
    if(TString(binarray->At(i)->GetName()).Contains("8_16_2_3"))dynamic_cast<TPad*>(arrayofpads->At(9))->cd();//yieldp->cd(6);
    if(TString(binarray->At(i)->GetName()).Contains("8_16_3_4"))dynamic_cast<TPad*>(arrayofpads->At(10))->cd();//yieldp->cd(8);
    if(TString(binarray->At(i)->GetName()).Contains("8_16_4_6"))dynamic_cast<TPad*>(arrayofpads->At(11))->cd();//yieldp->cd(8);
    if(TString(binarray->At(i)->GetName()).Contains("8_16_6_8"))dynamic_cast<TPad*>(arrayofpads->At(12))->cd();//yieldp->cd(8);
    if(TString(binarray->At(i)->GetName()).Contains("8_16_8_16"))dynamic_cast<TPad*>(arrayofpads->At(13))->cd();//yieldp->cd(8);

    hist->SetStats(false);
    hist->GetXaxis()->SetLabelSize(0.08);
    hist->GetXaxis()->SetTitleSize(0.08);
    hist->GetXaxis()->SetTitleOffset(0.5);
    hist->GetYaxis()->SetLabelSize(0.08);
    hist->GetYaxis()->SetTitleSize(0.065);
    hist->GetYaxis()->SetTitleOffset(0.5);
    
    hist->SetTitle("");
    hist->Draw("E");
  }

  
  for(int i = 0;i<PbPbbinarray->GetEntries();i++){
    for(int l = 0;l<cbinarray->GetEntries();l++){
      int centindex = 0;
      if(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Contains("(0.00)->(5.00)")) centindex = 1;
      else if(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Contains("(5.00)->(10.00)")) centindex = 2;
      else if(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Contains("(10.00)->(20.00)")) centindex = 3;
      else if(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Contains("(20.00)->(40.00)")) centindex = 4;
      else if(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Contains("(40.00)->(60.00)")) centindex = 5;
      else if(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Contains("(60.00)->(80.00)")) centindex = 6;
      for(int k = 0;k<PbPbfilearray->GetEntries();k++){
	if(!dynamic_cast<TFile*>(PbPbfilearray->At(k))->GetDirectory(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Data())) continue;
	TH1D * hist = dynamic_cast<TH1D*>(dynamic_cast<TFile*>(PbPbfilearray->At(k))->GetDirectory(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Data())->GetDirectory(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Data())->GetDirectory("yield/original_binning/GP0")->Get("dphiyieldbc"));
	TObjArray* arr;
	TH1D * near48tmp = NULL;
	TH1D * away48tmp = NULL;
	TH1D * near816tmp = NULL;
	TH1D * away816tmp = NULL;
	if(k==0||k==1){ 

	  if(centindex ==1){ arr = arrayofpadsyield10h05;near48tmp = near48_10h_05;away48tmp=away48_10h_05;near816tmp=near816_10h_05;away816tmp=away816_10h_05;}
	  else if(centindex ==2){ arr = arrayofpadsyield10h510;near48tmp = near48_10h_510;away48tmp=away48_10h_510;near816tmp=near816_10h_510;away816tmp=away816_10h_510;}
	  else if(centindex ==3){ arr = arrayofpadsyield10h1020;near48tmp = near48_10h_1020;away48tmp=away48_10h_1020;near816tmp=near816_10h_1020;away816tmp=away816_10h_1020;}
	  else if(centindex ==4){ arr = arrayofpadsyield10h2040;near48tmp = near48_10h_2040;away48tmp=away48_10h_2040;near816tmp=near816_10h_2040;away816tmp=away816_10h_2040;}
	  else if(centindex ==5){ arr = arrayofpadsyield10h4060;near48tmp = near48_10h_4060;away48tmp=away48_10h_4060;near816tmp=near816_10h_4060;away816tmp=away816_10h_4060;}
	  else if(centindex ==6) arr = arrayofpadsyield10h6080;	    
	}
	Double_t near = 0.0; Double_t neare=0.0;Double_t away = 0.0; Double_t awaye=0.0;
	near = hist->IntegralAndError(1,hist->FindBin(TMath::Pi()/2.0),neare,"width");
	away = hist->IntegralAndError(hist->FindBin(TMath::Pi()/2.0),hist->GetNbinsX(),awaye,"width");
	
	if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_0_1")){dynamic_cast<TPad*>(arr->At(0))->cd();if(near48tmp&&away48tmp){near48tmp->SetBinContent(1,near/0.5);away48tmp->SetBinContent(1,away/0.5);near48tmp->SetBinError(1,neare/0.5);away48tmp->SetBinError(1,awaye/0.5);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_1_2")){dynamic_cast<TPad*>(arr->At(1))->cd();if(near48tmp&&away48tmp){near48tmp->SetBinContent(2,near/1.0);away48tmp->SetBinContent(2,away/1.0);near48tmp->SetBinError(2,neare/1.0);away48tmp->SetBinError(2,awaye/1.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_2_3")){dynamic_cast<TPad*>(arr->At(2))->cd();if(near48tmp&&away48tmp){near48tmp->SetBinContent(3,near/1.0);away48tmp->SetBinContent(3,away/1.0);near48tmp->SetBinError(3,neare/1.0);away48tmp->SetBinError(3,awaye/1.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_3_4")){dynamic_cast<TPad*>(arr->At(3))->cd();if(near48tmp&&away48tmp){near48tmp->SetBinContent(4,near/1.0);away48tmp->SetBinContent(4,away/1.0);near48tmp->SetBinError(4,neare/1.0);away48tmp->SetBinError(4,awaye/1.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_4_6")){dynamic_cast<TPad*>(arr->At(4))->cd();if(near48tmp&&away48tmp){near48tmp->SetBinContent(5,near/2.0);away48tmp->SetBinContent(5,away/2.0);near48tmp->SetBinError(5,neare/2.0);away48tmp->SetBinError(5,awaye/2.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_6_8")){dynamic_cast<TPad*>(arr->At(5))->cd();if(near48tmp&&away48tmp){near48tmp->SetBinContent(6,near/2.0);away48tmp->SetBinContent(6,away/2.0);near48tmp->SetBinError(6,neare/2.0);away48tmp->SetBinError(6,awaye/2.0);}}

	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_0_1")){dynamic_cast<TPad*>(arr->At(7))->cd();if(near816tmp&&away816tmp){near816tmp->SetBinContent(1,near/0.5);away816tmp->SetBinContent(1,away/0.5);near816tmp->SetBinError(1,neare/0.5);away816tmp->SetBinError(1,awaye/0.5);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_1_2")){dynamic_cast<TPad*>(arr->At(8))->cd();if(near816tmp&&away816tmp){near816tmp->SetBinContent(2,near/1.0);away816tmp->SetBinContent(2,away/1.0);near816tmp->SetBinError(2,neare/1.0);away816tmp->SetBinError(2,awaye/1.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_2_3")){dynamic_cast<TPad*>(arr->At(9))->cd();if(near816tmp&&away816tmp){near816tmp->SetBinContent(3,near/1.0);away816tmp->SetBinContent(3,away/1.0);near816tmp->SetBinError(3,neare/1.0);away816tmp->SetBinError(3,awaye/1.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_3_4")){dynamic_cast<TPad*>(arr->At(10))->cd();if(near816tmp&&away816tmp){near816tmp->SetBinContent(4,near/1.0);away816tmp->SetBinContent(4,away/1.0);near816tmp->SetBinError(4,neare/1.0);away816tmp->SetBinError(4,awaye/1.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_4_6")){dynamic_cast<TPad*>(arr->At(11))->cd();if(near816tmp&&away816tmp){near816tmp->SetBinContent(5,near/2.0);away816tmp->SetBinContent(5,away/2.0);near816tmp->SetBinError(5,neare/2.0);away816tmp->SetBinError(5,awaye/2.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_6_8")){dynamic_cast<TPad*>(arr->At(12))->cd();if(near816tmp&&away816tmp){near816tmp->SetBinContent(6,near/2.0);away816tmp->SetBinContent(6,away/2.0);near816tmp->SetBinError(6,neare/2.0);away816tmp->SetBinError(6,awaye/2.0);}}
	else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_8_16")){dynamic_cast<TPad*>(arr->At(13))->cd();if(near816tmp&&away816tmp){near816tmp->SetBinContent(7,near/8.0);away816tmp->SetBinContent(7,away/8.0);near816tmp->SetBinError(7,neare/8.0);away816tmp->SetBinError(7,awaye/8.0);}}
	hist->SetStats(false);
	hist->GetXaxis()->SetLabelSize(0.08);
	hist->GetXaxis()->SetTitleSize(0.08);
	hist->GetXaxis()->SetTitleOffset(0.5);
	hist->GetYaxis()->SetLabelSize(0.08);
	hist->GetYaxis()->SetTitleSize(0.065);
	hist->GetYaxis()->SetTitleOffset(0.5);
	hist->SetTitle("");
	hist->Draw("E");
	
	TCanvas* can=NULL;
	TLegend* leg=NULL;
	if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_6_8")){
	  can = yields48_68_10h;
	  leg = leg48_68;
	}
	if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_6_8")){
	  can = yields816_68_10h;
	  leg = leg816_68;
	}
	if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_8_16")){
	  can = yields816_816_10h;
	  leg = leg816_816;
	}
	if(can&&leg){
	  can->cd();
// 	  if(centindex==1){hist->SetLineColor(kGreen);leg->AddEntry(hist, "C = (0.00) % ->(5.00) %");}
// 	  if(centindex==2){hist->SetLineColor(kRed);leg->AddEntry(hist,   "C = (5.00) % ->(10.00) %");}
// 	  if(centindex==3){hist->SetLineColor(kBlue);leg->AddEntry(hist,  "C = (10.00) % ->(20.00) %");}
// 	  if(centindex==4){hist->SetLineColor(kBlack);leg->AddEntry(hist, "C = (20.00) % ->(40.00) %");}
// 	  if(centindex==5){hist->SetLineColor(kViolet);leg->AddEntry(hist,"C = (40.00) % ->(60.00) %");}
	  hist->Draw("SAMEE");
	}
	
	for(int j = 0;j<histlist->GetEntries();j++){
	  TH1 * hist = dynamic_cast<TH1*>(dynamic_cast<TFile*>(PbPbfilearray->At(k))->GetDirectory(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Data())->GetDirectory(dynamic_cast<TObjString*>(cbinarray->At(l))->GetString().Data())->GetDirectory(dynamic_cast<TObjString*>(histposlist->At(j))->GetString().Data())->Get(histlist->At(j)->GetName()));


	  TH1D * hist1d = dynamic_cast<TH1D*>(hist);
	  TH2D * hist2d = dynamic_cast<TH2D*>(hist);
	
	  if(TString(hist->GetName()).Contains("DPhi_1_DEta_12_SameSide")&&centindex!=0){
	    TObjArray* arr;
	    if(k==0||k==1){ 
	      if(centindex ==1) arr = arrayofpadscor10h05;
	      else if(centindex ==2) arr = arrayofpadscor10h510;
	      else if(centindex ==3) arr = arrayofpadscor10h1020;
	      else if(centindex ==4) arr = arrayofpadscor10h2040;
	      else if(centindex ==5) arr = arrayofpadscor10h4060;
	      else if(centindex ==6) arr = arrayofpadscor10h6080;	    
	    }
	    
	    if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_0_1")){dynamic_cast<TPad*>(arr->At(0))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_1_2")){dynamic_cast<TPad*>(arr->At(1))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_2_3")){dynamic_cast<TPad*>(arr->At(2))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_3_4")){dynamic_cast<TPad*>(arr->At(3))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_4_6")){dynamic_cast<TPad*>(arr->At(4))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_6_8")){dynamic_cast<TPad*>(arr->At(5))->cd();}

	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_0_1")){dynamic_cast<TPad*>(arr->At(7))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_1_2")){dynamic_cast<TPad*>(arr->At(8))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_2_3")){dynamic_cast<TPad*>(arr->At(9))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_3_4")){dynamic_cast<TPad*>(arr->At(10))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_4_6")){dynamic_cast<TPad*>(arr->At(11))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_6_8")){dynamic_cast<TPad*>(arr->At(12))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_8_16")){dynamic_cast<TPad*>(arr->At(13))->cd();}
	    hist2d->SetStats(false);
	    hist2d->GetXaxis()->SetLabelSize(0.08);
	    hist2d->GetXaxis()->SetTitleSize(0.09);
	    hist2d->GetXaxis()->SetTitleOffset(0.3);
	    hist2d->GetYaxis()->SetLabelSize(0.08);
	    hist2d->GetYaxis()->SetTitleSize(0.09);
	    hist2d->GetYaxis()->SetTitleOffset(0.2);
	    hist2d->SetTitle("");
	    hist2d->Draw("colz");
	    gPad->Update();
	    TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
	    if(palette) palette->SetLabelSize(0.08);
	  }
	  if(TString(hist->GetName()).Contains("DPhi_1_DPHI_2")){
	    TObjArray* arr;
	    if(k==0||k==1){ 
	      if(centindex ==1) arr = arrayofpadscorphi10h05;
	      else if(centindex ==2) arr = arrayofpadscorphi10h510;
	      else if(centindex ==3) arr = arrayofpadscorphi10h1020;
	      else if(centindex ==4) arr = arrayofpadscorphi10h2040;
	      else if(centindex ==5) arr = arrayofpadscorphi10h4060;
	      else if(centindex ==6) arr = arrayofpadscorphi10h6080;	    
	    }
	    if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_0_1")){dynamic_cast<TPad*>(arr->At(0))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_1_2")){dynamic_cast<TPad*>(arr->At(1))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_2_3")){dynamic_cast<TPad*>(arr->At(2))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_3_4")){dynamic_cast<TPad*>(arr->At(3))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_4_6")){dynamic_cast<TPad*>(arr->At(4))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("4_8_6_8")){dynamic_cast<TPad*>(arr->At(5))->cd();}

	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_0_1")){dynamic_cast<TPad*>(arr->At(7))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_1_2")){dynamic_cast<TPad*>(arr->At(8))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_2_3")){dynamic_cast<TPad*>(arr->At(9))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_3_4")){dynamic_cast<TPad*>(arr->At(10))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_4_6")){dynamic_cast<TPad*>(arr->At(11))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_6_8")){dynamic_cast<TPad*>(arr->At(12))->cd();}
	    else if(dynamic_cast<TObjString*>(PbPbbinarray->At(i))->GetString().Contains("8_16_8_16")){dynamic_cast<TPad*>(arr->At(13))->cd();}
	    hist2d->SetStats(false);
	    hist2d->GetXaxis()->SetLabelSize(0.08);
	    hist2d->GetXaxis()->SetTitleSize(0.09);
	    hist2d->GetXaxis()->SetTitleOffset(0.3);
	    hist2d->GetYaxis()->SetLabelSize(0.08);
	    hist2d->GetYaxis()->SetTitleSize(0.09);
	    hist2d->GetYaxis()->SetTitleOffset(0.2);
	    hist2d->SetTitle("");
	    hist2d->Draw("colz");
	    gPad->Update();
	    TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
	    if(palette) palette->SetLabelSize(0.08);
	  }
	}
      }
    }
  }
  
  
  
  
  basedir->cd();
  testcanvas->Update();
  testcanvas->Write();
// 
//   corryields10b->Update();
//   corryields10b->Write();
//   corryields10c->Update();
//   corryields10c->Write();
  corryields10d->Update();
  corryields10d->Write();
//   corryields10e->Update();
//   corryields10e->Write();
  corryields11a->Update();
  corryields11a->Write();  
  corryields10h05->Update();
  corryields10h05->Write();
  corryields10h510->Update();
  corryields10h510->Write();
  corryields10h1020->Update();
  corryields10h1020->Write();
  corryields10h2040->Update();
  corryields10h2040->Write();
  corryields10h4060->Update();
  corryields10h4060->Write();
  corryields10h6080->Update();
  corryields10h6080->Write();

//   corryieldsphi10b->Update();
//   corryieldsphi10b->Write();
//   corryieldsphi10c->Update();
//   corryieldsphi10c->Write();
  corryieldsphi10d->Update();
  corryieldsphi10d->Write();
//   corryieldsphi10e->Update();
//   corryieldsphi10e->Write();
  corryieldsphi11a->Update();
  corryieldsphi11a->Write();
  corryieldsphi10h05->Update();
  corryieldsphi10h05->Write();
  corryieldsphi10h510->Update();
  corryieldsphi10h510->Write();
  corryieldsphi10h1020->Update();
  corryieldsphi10h1020->Write();
  corryieldsphi10h2040->Update();
  corryieldsphi10h2040->Write();
  corryieldsphi10h4060->Update();
  corryieldsphi10h4060->Write();
  corryieldsphi10h6080->Update();
  corryieldsphi10h6080->Write();
  
  yieldsphi10d->Update();
  yieldsphi10d->Write();
  yieldsphi11a->Update();
  yieldsphi11a->Write();
  yieldsphi10h05->Update();
  yieldsphi10h05->Write();
  yieldsphi10h510->Update();
  yieldsphi10h510->Write();
  yieldsphi10h1020->Update();
  yieldsphi10h1020->Write();
  yieldsphi10h2040->Update();
  yieldsphi10h2040->Write();
  yieldsphi10h4060->Update();
  yieldsphi10h4060->Write();
  yieldsphi10h6080->Update();
  yieldsphi10h6080->Write();
  
  delete testcanvas;delete corryields10d;delete corryieldsphi10d;delete corryields11a;delete corryieldsphi11a;delete corryields10h05;delete corryieldsphi10h05;delete corryields10h1020;delete corryieldsphi10h1020;delete corryields10h2040;delete corryieldsphi10h2040;delete corryields10h4060;delete corryieldsphi10h4060;delete corryields10h510;delete corryieldsphi10h510;delete corryields10h6080;delete corryieldsphi10h6080;
  delete yieldsphi10d;  delete yieldsphi11a;  delete yieldsphi10h05; delete yieldsphi10h510;delete yieldsphi10h1020;delete yieldsphi10h2040;delete yieldsphi10h4060;delete yieldsphi10h6080;
  yields48_68_10h->cd();
  leg48_68->Draw("SAME");
  yields48_68_10h->Update();
  yields48_68_10h->Write();
  yields816_68_10h->cd();
  leg816_68->Draw("SAME");
  yields816_68_10h->Update();
  yields816_68_10h->Write();
  
  near48->Write();
  away48->Write();
  near816->Write();
  away816->Write();
  away48->Divide(near48);
  away48->SetTitle("Away side yield divided by near side yield for 4.0GeV/c <pT_{trig}< 8.0GeV/c.");
  away48->GetYaxis()->SetTitle("#frac{dN_{pairs-away}}{dN_{pairs - near}}");
  away48->Write("awaydivnear48");
  away816->Divide(near816);
  away816->SetTitle("Away side yield divided by near side yield for 8.0GeV/c <pT_{trig}<16.0GeV/c.");
  away816->GetYaxis()->SetTitle("#frac{dN_{pairs-away}}{dN_{pairs - near}}");
  away816->Write("awaydivnear816");
  
  near48_10h_05->Write();
  away48_10h_05->Write();
//   away48_10h_05->Divide(near48_10h_05);
//   away48_10h_05->SetTitle("Away side yield divided by near side yield for 4.0GeV/c <pT_{trig}< 8.0GeV/c in 0% - 5% PbPb events.");
//   away48_10h_05->Write("awaydivnear48_05");
  near816_10h_05->Write();
  away816_10h_05->Write();
//   away816_10h_05->Divide(near816_10h_05);
//   away816_10h_05->SetTitle("Away side yield divided by near side yield for 8.0GeV/c <pT_{trig}<16.0GeV/c in 0% - 5% PbPb events.");
//   away816_10h_05->Write("awaydivnear816_05");
  near48_10h_510->Write();
  away48_10h_510->Write();
//   away48_10h_510->Divide(near48_10h_510);
//   away48_10h_510->SetTitle("Away side yield divided by near side yield for 4.0GeV/c <pT_{trig}< 8.0GeV/c in 5% - 10% PbPb events.");
//   away48_10h_510->Write("awaydivnear48_510");
  near816_10h_510->Write();
  away816_10h_510->Write();
//   away816_10h_510->Divide(near816_10h_510);
//   away816_10h_510->SetTitle("Away side yield divided by near side yield for 8.0GeV/c <pT_{trig}<16.0GeV/c in 5% - 10% PbPb events.");
//   away816_10h_510->Write("awaydivnear816_510");
  near48_10h_1020->Write();
  away48_10h_1020->Write();
//   away48_10h_1020->Divide(near48_10h_1020);
//   away48_10h_1020->SetTitle("Away side yield divided by near side yield for 4.0GeV/c <pT_{trig}< 8.0GeV/c in 5% - 10% PbPb events.");
//   away48_10h_1020->Write("awaydivnear48_1020");
  near816_10h_1020->Write();
  away816_10h_1020->Write();
//   away816_10h_1020->Divide(near816_10h_1020);
//   away816_10h_1020->SetTitle("Away side yield divided by near side yield for 8.0GeV/c <pT_{trig}<16.0GeV/c in 5% - 10% PbPb events.");
//   away816_10h_1020->Write("awaydivnear816_1020");
  near48_10h_2040->Write();
  away48_10h_2040->Write();
//   away48_10h_2040->Divide(near48_10h_2040);
//   away48_10h_2040->SetTitle("Away side yield divided by near side yield for 4.0GeV/c <pT_{trig}< 8.0GeV/c in 5% - 10% PbPb events.");
//   away48_10h_2040->Write("awaydivnear48_2040");
  near816_10h_2040->Write();
  away816_10h_2040->Write();
//   away816_10h_2040->Divide(near816_10h_2040);
//   away816_10h_2040->SetTitle("Away side yield divided by near side yield for 8.0GeV/c <pT_{trig}<16.0GeV/c in 5% - 10% PbPb events.");
//   away816_10h_2040->Write("awaydivnear816_2040");
  near48_10h_4060->Write();
  away48_10h_4060->Write();
//   away48_10h_4060->Divide(near48_10h_4060);
//   away48_10h_4060->SetTitle("Away side yield divided by near side yield for 4.0GeV/c <pT_{trig}< 8.0GeV/c in 5% - 10% PbPb events.");
//   away48_10h_4060->Write("awaydivnear48_4060");
  near816_10h_4060->Write();
  away816_10h_4060->Write();
//   away816_10h_4060->Divide(near816_10h_4060);
//   away816_10h_4060->SetTitle("Away side yield divided by near side yield for 8.0GeV/c <pT_{trig}<16.0GeV/c in 5% - 10% PbPb events.");
//   away816_10h_4060->Write("awaydivnear816_4060");
  near48_10h_05->Divide(near48_10h_4060);
  near48_10h_05->SetTitle("Near side yield in 0% - 5% divided by near side yield in 40% - 60 % in PbPb collisions with 4.0GeV/c <pT_{trig}< 8.0GeV/c ");
  near48_10h_05->GetYaxis()->SetTitle("I_CP");
  near48_10h_05->Write("ICP_48_near");
  away48_10h_05->Divide(away48_10h_4060);
  away48_10h_05->SetTitle("Away side yield in 0% - 5% divided by away side yield in 40% - 60 % in PbPb collisions with 4.0GeV/c <pT_{trig}< 8.0GeV/c ");  
  away48_10h_05->GetYaxis()->SetTitle("I_CP");
  away48_10h_05->Write("ICP_48_away");

  near816_10h_05->Divide(near816_10h_4060);
  near816_10h_05->SetTitle("Near side yield in 0% - 5% divided by near side yield in 40% - 60 % in PbPb collisions with 8.0GeV/c <pT_{trig}< 16.0GeV/c ");
  near816_10h_05->GetYaxis()->SetTitle("I_CP");
  near816_10h_05->Write("ICP_816_near");
  away816_10h_05->Divide(away816_10h_4060);
  away816_10h_05->SetTitle("Away side yield in 0% - 5% divided by away side yield in 40% - 60 % in PbPb collisions with 8.0GeV/c <pT_{trig}< 16.0GeV/c ");  
  away816_10h_05->GetYaxis()->SetTitle("I_CP");
  away816_10h_05->Write("ICP_816_away");

  for(int i=0;i<filearray->GetEntries();i++){
    dynamic_cast<TFile*>(filearray->At(i))->Close();
  }
  delete filearray;
  for(int i=0;i<PbPbfilearray->GetEntries();i++){
    dynamic_cast<TFile*>(PbPbfilearray->At(i))->Close();
  }
  delete PbPbfilearray;
  rfile->Close();
  CanvasList->SetOwner(true); CanvasList->Clear();delete CanvasList;
  binarray->SetOwner(true);binarray->Clear();delete binarray;
  delete Resultsbinarray;
  histlist->SetOwner(true);histlist->Clear();delete histlist;
  delete histposlist;
//   delete near48;delete near816;delete away48;delete away816;
}

void Centralities(){
  //go over periods and collect the relevant histograms in a new results file.
  TFile * rfile = TFile::Open("results.root","UPDATE");
  TString * path = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  TString * path10h48 = new TString("LHC10h/Train48");
  TString * path10h816 = new TString("LHC10h/Train816");

  
  TDirectory * basedir = resultsdirectory(rfile->GetDirectory("/"),"Cents");
  TObjArray * filearray = new TObjArray();
  filearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10h48->Data()),"READ"));
  filearray->Add(TFile::Open(Form("%s%s/results.root",path->Data(),path10h816->Data()),"READ"));

  
  TCanvas * c05p = new TCanvas("05perc");
  TCanvas * c510p = new TCanvas("510perc");
  TCanvas * c1020p = new TCanvas("1020perc");
  TCanvas * c2040p = new TCanvas("2040perc");
  TCanvas * c4060p = new TCanvas("4060perc");
  TCanvas * c6080p = new TCanvas("6080perc");

  
  TObjArray * trigarray = new TObjArray();
  trigarray->Add(new TObjString("4 GeV/c <= p_{T}^{trigger}<=8 GeV//c"));
  trigarray->Add(new TObjString("8 GeV/c <= p_{T}^{trigger}<=16 GeV//c"));
  TObjArray * assarray  = new TObjArray();
  assarray->Add(new TObjString("0.5-1.0 GeV//c"));
  assarray->Add(new TObjString("1.0-2.0 GeV//c"));
  assarray->Add(new TObjString("2.0-4.0 GeV//c"));
  assarray->Add(new TObjString("4.0-8.0 GeV//c"));
  
  TObjArray* array05 = PadArray(c05p,2,4,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 0%%-5%% centrality in different pT ranges");
  TObjArray* array510 = PadArray(c510p,2,4,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 5%%-10%%centrality in different pT ranges");
  TObjArray* array1020 = PadArray(c1020p,2,4,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 10%%-20%%centrality in different pT ranges");
  TObjArray* array2040 = PadArray(c2040p,2,4,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 20%%-40%%centrality in different pT ranges");  
  TObjArray* array4060p = PadArray(c4060p,2,4,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 40%%-60%%centrality in different pT ranges");  
  TObjArray* array6080p = PadArray(c6080p,2,4,trigarray,assarray,"Corrected correlation function #Delta#Phi_{1} vs #Delta#eta_{12} in 10h 60%%-80%%centrality in different pT ranges");  

  TList * list48 = dynamic_cast<TFile*>(filearray->At(0))->GetListOfKeys();
  for(int i = 0;i<list48->GetEntries();i++){
    int npad = -1;
    if(TString(list48->At(i)->GetName()).Contains("4_8_0_1")){npad = 0;}
    else if(TString(list48->At(i)->GetName()).Contains("4_8_1_2")){npad = 1;}
    else if(TString(list48->At(i)->GetName()).Contains("4_8_2_4")){npad = 2;}
    if(npad !=-1){
      dynamic_cast<TPad*>(array05->At(npad))->cd();
      TH2D * hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(list48->At(i)->GetName())->GetDirectory("BinM(0.00)->(5.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array510->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(list48->At(i)->GetName())->GetDirectory("BinM(5.00)->(10.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array1020->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(list48->At(i)->GetName())->GetDirectory("BinM(10.00)->(20.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array2040->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(list48->At(i)->GetName())->GetDirectory("BinM(20.00)->(40.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");        
      dynamic_cast<TPad*>(array4060p->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(list48->At(i)->GetName())->GetDirectory("BinM(40.00)->(60.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array6080p->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(list48->At(i)->GetName())->GetDirectory("BinM(60.00)->(80.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");        
    }
  }
  
  TList * list816 = dynamic_cast<TFile*>(filearray->At(1))->GetListOfKeys();
  for(int i = 0;i<list816->GetEntries();i++){
    int npad = -1;
    if(TString(list816->At(i)->GetName()).Contains("8_16_0_1")){npad = 4;}
    else if(TString(list816->At(i)->GetName()).Contains("8_16_1_2")){npad = 5;}
    else if(TString(list816->At(i)->GetName()).Contains("8_16_2_4")){npad = 6;}
    else if(TString(list816->At(i)->GetName()).Contains("8_16_4_8")){npad = 7;}
    if(npad !=-1){
      dynamic_cast<TPad*>(array05->At(npad))->cd();
      TH2D * hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(1))->GetDirectory(list816->At(i)->GetName())->GetDirectory("BinM(0.00)->(5.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array510->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(1))->GetDirectory(list816->At(i)->GetName())->GetDirectory("BinM(5.00)->(10.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array1020->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(1))->GetDirectory(list816->At(i)->GetName())->GetDirectory("BinM(10.00)->(20.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array2040->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(1))->GetDirectory(list816->At(i)->GetName())->GetDirectory("BinM(20.00)->(40.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");        
      dynamic_cast<TPad*>(array4060p->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(1))->GetDirectory(list816->At(i)->GetName())->GetDirectory("BinM(40.00)->(60.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
      dynamic_cast<TPad*>(array6080p->At(npad))->cd();
      hist = dynamic_cast<TH2D*>(dynamic_cast<TFile*>(filearray->At(1))->GetDirectory(list816->At(i)->GetName())->GetDirectory("BinM(60.00)->(80.00)/iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      if(hist) hist->Draw("colz");    
    }
  }
  basedir->cd();
  c05p->Update();
  c05p->Write();
  c510p->Update();
  c510p->Write();
  c1020p->Update();
  c1020p->Write();
  c2040p->Update();  
  c2040p->Write();  
  c4060p->Update();  
  c4060p->Write(); 
  c6080p->Update();  
  c6080p->Write(); 
  rfile->Close();
//       else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_0_1")){dynamic_cast<TPad*>(arr->At(4))->cd();}
//     else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_1_2")){dynamic_cast<TPad*>(arr->At(5))->cd();}
//     else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_2_4")){dynamic_cast<TPad*>(arr->At(6))->cd();}
//     else if(dynamic_cast<TObjString*>(binarray->At(i))->GetString().Contains("8_16_4_8")){dynamic_cast<TPad*>(arr->At(7))->cd();}
    
  
//   //get the relevant bins:
//   TObjArray * binarray= new TObjArray();
//   TObjArray * Resultsbinarray = new TObjArray();
//   TList * binlist = dynamic_cast<TDirectory*>(filearray->At(0))->GetListOfKeys();
//   for(int i = 0;i<binlist->GetEntries();i++){
//     TObjString * bin = new TObjString(binlist->At(i)->GetName());
//     if(bin->GetString().Contains("ThreePartTracks")){
//       binarray->Add(bin);
//       Resultsbinarray->Add(resultsdirectory(basedir,bin->GetString().Data()));
//     }
//     else delete bin;
//   }
/*
  //get the relevant histograms:
  TObjArray * histlist= new TObjArray();
  TObjArray * histposlist= new TObjArray();
  {
    TObjString * iteration = new TObjString("iteration1");
    TObjString * yield = new TObjString("yield/original_binning/GP0");
    TObjString * dphi1dphi2 = new TObjString("DPhi_1_DPHI_2");
    histlist->Add(dphi1dphi2);
    histposlist->Add(iteration);
    TObjString * dphi1detass = new TObjString("DPhi_1_DEta_12_SameSide");
    histlist->Add(dphi1detass);
    histposlist->Add(iteration);
    TObjString * yieldh = new TObjString("dphiyieldbc");
    histlist->Add(yieldh);
    histposlist->Add(yield);
  }*/
  
//   TH1D * yieldbcperiodcollected = dynamic_cast<TH1D*>(dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(dynamic_cast<TObjString*>(binarray->At(0))->GetString().Data())->GetDirectory("yield/original_binning/GP0")->Get("dphiyieldbc"));
//   yieldbcperiodcollected->Reset();
  
//   Double_t bins[4] = {0.5,1.0,2.0,4.0};
//   TH1D * near48 = new TH1D("nearsidey48","Yield on the near side with Trigger pT= 4.0-8.0 GeV/c.",3,bins);
//   near48->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
//   TH1D * away48 = new TH1D("awaysidey48","Yield on the away side with Trigger pT= 4.0-8.0 GeV/c.",3,bins);
//   away48->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
//   Double_t bins816[5] = {0.5,1.0,2.0,4.0,8.0};
//   TH1D * near816 = new TH1D("nearsidey816","Yield on the near side with Trigger pT= 8.0-16.0 GeV/c.",4,bins816);
//   near816->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
//   TH1D * away816 = new TH1D("awaysidey816","Yield on the away side with Trigger pT= 8.0-16.0 GeV/c.",4,bins816);
//   away816->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}d#eta_{12} dpT_{A}} (rad)^{-1}");
//   
//   
//   TCanvas * yieldp = new TCanvas("YieldsinPeriods");
//   yieldp->Divide(2,4);
  
//   TObjArray * CanvasList = new TObjArray();
// //   for(int i = 0;i<histlist->GetEntries();i++){
//     TDirectory * dir = dynamic_cast<TFile*>(filearray->At(0))->GetDirectory(dynamic_cast<TObjString*>(binarray->At(0))->GetString().Data())->GetDirectory(dynamic_cast<TObjString*>(histposlist->At(i))->GetString().Data());
//     TH1 * hist = dynamic_cast<TH1*>(dir->Get(dynamic_cast<TObjString*>(histlist->At(i))->GetString().Data()));
//     TCanvas * can = Periodcanvas(hist,2);
//     CanvasList->Add(can);
// //   }
//   delete near48;delete near816;delete away48;delete away816;
}

void CountRawNumbers(){
   //take the same event histograms and compare to different mixed.
  TFile * infile = TFile::Open("results.root","READ");
  TFile * outfile =  TFile::Open("CompareMixed.root","RECREATE");
  TList * folderlist = infile->GetListOfKeys();
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    if(!TString(folder).Contains("ThreePartTracks"))continue;
   
    BinDirs * divsame = new BinDirs(infile->GetDirectory(Form("%s/",folder)),infile->GetDirectory(Form("%s/META",folder)),infile->GetDirectory(Form("%s/METrigger",folder)),true);
    TDirectory* outdir =  outfile->mkdir(folder);
    
    //List of directories for multiplicity bins:
    TList * directories = GetMZDirectories(divsame);
    //Go through the list and find the multiplicity/centrality binning in order to sum over them and add them to the TObjArray:
    TString  s = TString("");
    for(int i=0; i<directories->GetEntries();i++){
      if(TString(directories->At(i)->GetName()).Contains("Z")){
	TDirectory* outdir2= outdir->mkdir(directories->At(i)->GetName());
	TDirectory* SAMEDir = divsame->SameDir(directories->At(i)->GetName());
	TDirectory* METADir = divsame->METADir(directories->At(i)->GetName());
	TDirectory* METriggerDir = divsame->METriggerDir(directories->At(i)->GetName());
	TH3F * Same3d = dynamic_cast<TH3F*>(SAMEDir->GetDirectory("same_event")->Get("DPhi_1_DPhi_2_DEta_12")->Clone("CorSame3d"));
	TH1D * SameNtriggers = dynamic_cast<TH1D*>(SAMEDir->GetDirectory("same_event")->Get("number_of_triggers")->Clone("NTriggersSame"));
// 	if(SameNtriggers->Integral()>0.1)cout << Same3d->Integral()<< " "<<SameNtriggers->Integral()<<" "<< Same3d->Integral()/SameNtriggers->Integral() <<endl;
	if(SameNtriggers->Integral()>0.1)Same3d->Scale(1.0/SameNtriggers->Integral());
	else Same3d->Scale(0.0);
// 	cout << "s"<< Same3d->Integral()<<endl;
	Same3d->SetTitle(Form("Same Event, integral = %f",Same3d->Integral()));
	TH3F * META3d = dynamic_cast<TH3F*>(METADir->GetDirectory("same_event")->Get("DPhi_1_DPhi_2_DEta_12")->Clone("CorMETA3d"));
	TH1D * METANtriggers = dynamic_cast<TH1D*>(METADir->GetDirectory("same_event")->Get("number_of_triggers")->Clone("NTriggersMETA"));
// 	if(METANtriggers->Integral()>0.1)cout << META3d->Integral()<< " "<<METANtriggers->Integral()<<" "<< META3d->Integral()/METANtriggers->Integral() <<endl;
	if(METANtriggers->Integral()>0.1)META3d->Scale(1.0/METANtriggers->Integral());
	else META3d->Scale(0.0);
// 	cout<<"o" << META3d->Integral()<<endl;
	
	META3d->SetTitle(Form("Mixed Event Trigger - Associated, integral = %f",META3d->Integral()));
	TH3F * METrigger3d = dynamic_cast<TH3F*>(METriggerDir->GetDirectory("same_event")->Get("DPhi_1_DPhi_2_DEta_12")->Clone("CorMETrigger3d"));
	TH1D * METriggerNtriggers = dynamic_cast<TH1D*>(METriggerDir->GetDirectory("same_event")->Get("number_of_triggers")->Clone("NTriggersMETrigger"));
// 	if(METriggerNtriggers->Integral()>0.1)cout << METrigger3d->Integral()<< " "<<METriggerNtriggers->Integral()<<" "<< METrigger3d->Integral()/METriggerNtriggers->Integral() <<endl;
	if(METriggerNtriggers->Integral()>0.1)METrigger3d->Scale(1.0/METriggerNtriggers->Integral());
	else METrigger3d->Scale(0.0);
// 	cout<<"a" << METrigger3d->Integral()<<endl;
	METrigger3d->SetTitle(Form("Mixed Event Associated - Associated, integral = %f",METrigger3d->Integral()));
	TCanvas * can = new TCanvas("3dintegrals");
	can->Divide(2,2);
	can->cd(1);
	Same3d->Draw();
	can->cd(2);
	META3d->Draw();
	can->cd(3);
	METrigger3d->Draw();
	can->cd(4);
	TH3F* Substr = dynamic_cast<TH3F*>(Same3d->Clone("Substracted"));
	Substr->Add(META3d,-1.0/3.0);
	Substr->Add(METrigger3d,-1.0/3.0);
	Substr->SetTitle(Form("Substracted, integral = %f",Substr->Integral()));
	Substr->Draw();
	
	outdir2->cd();
	can->Write();
	
	TH2D * SameDphiDphi = dynamic_cast<TH2D*>(SAMEDir->GetDirectory("same_event")->Get("DPhi_1_DPHI_2")->Clone("CorSamePhiPhi"));
	if(SameNtriggers->Integral()>0.1)SameDphiDphi->Scale(1.0/SameNtriggers->Integral());
	else SameDphiDphi->Scale(0.0);
	TH2D * METADphiDphi = dynamic_cast<TH2D*>(METADir->GetDirectory("same_event")->Get("DPhi_1_DPHI_2")->Clone("CorMETAPhiPhi"));
	if(METANtriggers->Integral()>0.1)METADphiDphi->Scale(1.0/METANtriggers->Integral());
	else METADphiDphi->Scale(0.0);
	TH2D * METriggerDphiDphi = dynamic_cast<TH2D*>(METriggerDir->GetDirectory("same_event")->Get("DPhi_1_DPHI_2")->Clone("CorMETriggerPhiPhi"));
	if(METriggerNtriggers->Integral()>0.1)METriggerDphiDphi->Scale(1.0/METriggerNtriggers->Integral());
	else METriggerDphiDphi->Scale(0.0);
		
	
	can->cd(1);
	SameDphiDphi->Draw("colz");
	can->cd(2);
	METADphiDphi->Draw("colz");
	can->cd(3);
	METriggerDphiDphi->Draw("colz");
	can->cd(4);
	TH2D* SubstrPhiPhi = dynamic_cast<TH2D*>(SameDphiDphi->Clone("SubstractedPhiPhi"));
	SubstrPhiPhi->Add(METADphiDphi,-1.0);
	SubstrPhiPhi->Add(METriggerDphiDphi,-1.0);
	SubstrPhiPhi->SetTitle(Form("Substracted, integral = %f",SubstrPhiPhi->Integral()));
	SubstrPhiPhi->Draw("colz");
	can->Write("2dintegrals");
	
	
	TH2D * SameDphiDphid = dynamic_cast<TH2D*>(SAMEDir->GetDirectory("divided")->Get("DPhi_1_DPHI_2")->Clone("CorSamePhiPhid"));
// 	if(SameNtriggers->Integral()>0.1)SameDphiDphi->Scale(1.0/SameNtriggers->Integral());
// 	else SameDphiDphi->Scale(0.0);
	TH2D * METADphiDphid = dynamic_cast<TH2D*>(METADir->GetDirectory("divided")->Get("DPhi_1_DPHI_2")->Clone("CorMETAPhiPhid"));
// 	if(METANtriggers->Integral()>0.1)METADphiDphi->Scale(1.0/METANtriggers->Integral());
// 	else METADphiDphi->Scale(0.0);
	TH2D * METriggerDphiDphid = dynamic_cast<TH2D*>(METriggerDir->GetDirectory("divided")->Get("DPhi_1_DPHI_2")->Clone("CorMETriggerPhiPhid"));
// 	if(METriggerNtriggers->Integral()>0.1)METriggerDphiDphi->Scale(1.0/METriggerNtriggers->Integral());
// 	else METriggerDphiDphi->Scale(0.0);
		
	
	can->cd(1);
	SameDphiDphid->Draw("colz");
	can->cd(2);
	METADphiDphid->Draw("colz");
	can->cd(3);
	METriggerDphiDphid->Draw("colz");
	can->cd(4);
	TH2D* SubstrPhiPhid = dynamic_cast<TH2D*>(SameDphiDphi->Clone("SubstractedPhiPhid"));
	SubstrPhiPhid->Add(METADphiDphi,-1.0);
	SubstrPhiPhid->Add(METriggerDphiDphi,-1.0);
	SubstrPhiPhid->SetTitle(Form("Substracted, integral = %f",SubstrPhiPhid->Integral()));
	SubstrPhiPhid->Draw("colz");
	can->Write("2dintegralsd");
	
	
	TH2D * SameDphiDetad = dynamic_cast<TH2D*>(SAMEDir->GetDirectory("divided")->Get("DPhi_1_DEta_12_SameSide")->Clone("CorSamePhiEtad"));
	TH2D * MEDphiDetad   = dynamic_cast<TH2D*>(SAMEDir->GetDirectory("mixed_event")->Get("DPhi_1_DEta_12_SameSide")->Clone("CorMEPhiEtad"));
	TParameter<double> * par = dynamic_cast<TParameter<double>*>(SAMEDir->GetDirectory("mixed_event")->Get("DPhi_1_DEta_12_SameSide_scale")->Clone("CorMEPhiEtad_scale"));
	TH1D * MENtriggers = dynamic_cast<TH1D*>(SAMEDir->GetDirectory("mixed_event")->Get("number_of_triggers")->Clone("NTriggersME"));
	MEDphiDetad->Divide(MEDphiDetad);
	MEDphiDetad->Scale(par->GetVal());
	if(MENtriggers->Integral()>0.1) MEDphiDetad->Scale(1.0/MENtriggers->Integral());
	else MEDphiDetad->Scale(0.0);
	// 	if(SameNtriggers->Integral()>0.1)SameDphiDphi->Scale(1.0/SameNtriggers->Integral());
// 	else SameDphiDphi->Scale(0.0);
	TH2D * METADphiDetad = dynamic_cast<TH2D*>(METADir->GetDirectory("divided")->Get("DPhi_1_DEta_12_SameSide")->Clone("CorMETAPhiEtad"));
// 	if(METANtriggers->Integral()>0.1)METADphiDphi->Scale(1.0/METANtriggers->Integral());
// 	else METADphiDphi->Scale(0.0);
	TH2D * METriggerDphiDetad = dynamic_cast<TH2D*>(METriggerDir->GetDirectory("divided")->Get("DPhi_1_DEta_12_SameSide")->Clone("CorMETriggerPhiEtad"));
// 	if(METriggerNtriggers->Integral()>0.1)METriggerDphiDphi->Scale(1.0/METriggerNtriggers->Integral());
// 	else METriggerDphiDphi->Scale(0.0);
		
	
	can->cd(1);
	SameDphiDetad->Draw("colz");
	can->cd(2);
	METADphiDetad->Draw("colz");
	can->cd(3);
	METriggerDphiDetad->Draw("colz");
	can->cd(4);
	TH2D* SubstrPhiEtad = dynamic_cast<TH2D*>(SameDphiDetad->Clone("SubstractedPhiEtad"));
	SubstrPhiEtad->Add(METADphiDetad,-1.0);
	SubstrPhiEtad->Add(METriggerDphiDetad,-1.0);
	SubstrPhiEtad->Add(MEDphiDetad,1.0);
	SubstrPhiEtad->SetTitle(Form("Substracted, integral = %f",SubstrPhiEtad->Integral()));
	SubstrPhiEtad->Draw("colz");
	can->Write("2dintegralsd2");
	
	
      }
    }
    
    //Get Tokens for all histograms:
//     TStringToken histtokensbinstats = GetHistTokens(divsame->Same()->GetDirectory(Form("%s/bin_stats",directories->At(1)->GetName())));
//     while(histtokensbinstats.NextToken()){CollectHistbinstats(histtokensbinstats.Data(),directories,multdirlist);}
//     TStringToken histtokens = GetHistTokens(divsame->Same()->GetDirectory(Form("%s/divided",directories->At(1)->GetName())));
//     while(histtokens.NextToken()){CollectHist(histtokens.Data(),directories,multdirlist,false,CDivfirst);}

  }
  infile->Close();
  delete infile;
  outfile->Close();
  delete outfile;
}