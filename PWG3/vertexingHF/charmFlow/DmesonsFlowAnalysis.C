#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TASImage.h>
#include <TPad.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TLatex.h>

#include <AliHFMassFitter.h>
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"

#endif

//methods for the analysis of AliAnalysisTaskSEHFv2 output
//Authors: Chiara Bianchin cbianchi@pd.infn.it
//         Giacomo Ortona  ortona@to.infn.it
//         Francesco Prino prino@to.infn.it

//global variables to be set
// const Int_t nptbinsnew=6;
// Float_t ptbinsnew[nptbinsnew+1]={2,3,4,5,8,12,999.};
const Int_t nptbinsnew=3;
Float_t ptbinsnew[nptbinsnew+1]={2,5,8,16.};
//const Int_t nptbinsnew=4;
//Float_t ptbinsnew[nptbinsnew+1]={2,4,5,8,999.};
//Float_t ptbinsnew[nptbinsnew+1]={2,3,5,8,999.};
// const Int_t nptbinsnew=5;
// Float_t ptbinsnew[nptbinsnew+1]={2,3,4,5,8,999.};
// const Int_t nptbinsnew=2;
// Float_t ptbinsnew[nptbinsnew+1]={2,8,999.};
Int_t fittype=0;
//Int_t rebin[nptbinsnew]={4,4,4,5};
Int_t rebin[nptbinsnew]={4,4,4};
Double_t nsigma=3;

//methods
Bool_t ReadFile(TList* &list,TH1F* &hstat,AliRDHFCuts* &cutobj,TString listname,TString partname,TString path="./",TString filename="AnalysisResults.root");
void WriteCanvas(TCanvas* cv,TString text);
Int_t FindPtBin(Int_t nbins, Float_t* array,Float_t value);
void InOutPic(TVirtualPad *c,Int_t inout=0,TString where="tr");//inout: 0=IN, 1=OUT 
void PhiBinPic(TVirtualPad *c,Int_t angle,TString where);
Double_t FindChi(Double_t res, Int_t k);
Double_t Pol(Double_t x, Int_t k);
Double_t ResolK1(Double_t x);
Double_t ComputeResol(Double_t resSub, Int_t k);


//methods implementation

Bool_t ReadFile(TList* &list,TH1F* &hstat,AliRDHFCuts* &cutsobj,TString listname,TString partname,TString path,TString filename){

  TString hstatname="hEventsInfo",dirname="PWG3_D2H_HFv2";
  filename.Prepend(path);
  listname+=partname;
  hstatname+=partname;
  // TString tmpsuff="NoCos"; //"Nod0d0"
  // listname+=tmpsuff;
  // hstatname+=tmpsuff;

  TFile* f=new TFile(filename.Data());
  if(!f){
    cout<<filename.Data()<<" not found"<<endl;
    return kFALSE;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname);
  if(!f){
    cout<<dirname.Data()<<" not found  in "<<filename.Data()<<endl;
    return kFALSE;
  }

  list=(TList*)dir->Get(listname);
  if(!list){
    cout<<"List "<<listname.Data()<<" not found"<<endl;
    dir->ls();
    return kFALSE;
  }

  hstat=(TH1F*)dir->Get(hstatname);
  if(!hstat){
    cout<<hstatname.Data()<<" not found"<<endl;
    return kFALSE;
  }
  
  if(partname.Contains("D0")) cutsobj=((AliRDHFCutsD0toKpi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
  if(partname.Contains("Dplus")) cutsobj=((AliRDHFCutsDplustoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
  if(partname.Contains("Dstar")) cutsobj=((AliRDHFCutsDStartoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
		
  if(!cutsobj){
    cout<<"Cut object not found! check position of the key!"<<endl;
    return kFALSE;
  }
  
  return kTRUE;
}

Int_t FindPtBin(Int_t nbins, Float_t* array,Float_t value){
  for (Int_t i=0;i<nbins;i++){
    //cout<<value<<" from "<<array[i]<<" to "<<array[i+1]<<"?"<<endl;
    if(value>=array[i] && value<array[i+1]){
      return i;
    }
  }
  cout<<value<< " out of range "<<array[0]<<", "<<array[nbins]<<endl;
  return -1;

}

void WriteCanvas(TCanvas* cv,TString text){
  cv->SaveAs(Form("%s%s.png",cv->GetName(),text.Data()));
}

void InOutPic(TVirtualPad *c,Int_t inout,TString where){
  TPad *myPadLogo = 0x0;
  if(where=="tr"){
    myPadLogo = new TPad("myPadPic", "Pad for In/Out pic",0.66,0.68,0.86,0.88);
  }

  myPadLogo->SetFillColor(0); 
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage(Form("/home/cbianchi/macros/sketch%s.png",inout==0 ? "IN" : "OUT"));
  myAliceLogo->Draw();
  c->cd();
  TPaveText* t1=0x0;//=new TPaveText(0.62,0.63,0.89,0.71,"NDC");
  if(where=="tr"){
    t1=new TPaveText(0.65,0.58,0.85,0.67,"NDC");
  }
  t1->SetFillStyle(0);
  t1->AddText(Form("%s-PLANE",inout==0 ? "IN" : "OUT-OF"));
  if(inout==0)t1->SetTextColor(kBlue);
  else t1->SetTextColor(kRed);
  t1->Draw();
}

void PhiBinPic(TVirtualPad *c,Int_t angle,TString where){

  TString picname="angles";
  if(angle==0) picname.Append("0pi4");
  if(angle==1) picname.Append("pi4pi2");
  if(angle==2) picname.Append("pi232pi");
  if(angle==3) picname.Append("32pipi");
  picname.Append(".png");

  TPad *myPadLogo = 0x0;
  if(where=="tr"){
    myPadLogo = new TPad("myPadPic", "Pad for bin of phi pic",0.5,0.75,0.8,0.8);
  }

  myPadLogo->SetFillColor(0); 
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage(Form("/home/cbianchi/macros/%s",picname.Data()));
  c->cd();
  myAliceLogo->Draw();
  

  /*
  TPaveText* t1=0x0;//=new TPaveText(0.62,0.63,0.89,0.71,"NDC");
  if(where=="tr"){
    t1=new TPaveText(0.65,0.58,0.85,0.67,"NDC");
  }
  t1->SetFillStyle(0);
  t1->AddText(Form("%s-PLANE",inout==0 ? "IN" : "OUT-OF"));
  if(inout==0)t1->SetTextColor(kBlue);
  else t1->SetTextColor(kRed);
  t1->Draw();
  */
}


void DmesonsSignalExtraction(TString partname="D0",Int_t mincentr=30,Int_t maxcentr=70,TString textleg="",TString path="./",Bool_t official=kFALSE){

  //read the output of HFv2task and extract the signal in pt bins in plane and out of plane and in bins of phi pt integrated

  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  Int_t info=1;
  if (official) info=0;

  if(textleg==""){
    textleg=Form("centr%d-%d",mincentr,maxcentr);
  }
  //  Double_t pi=TMath::Pi();
  TString listname="coutputv2";

  TList* list;
  TH1F * hstat;
  AliRDHFCuts* cutobj;

  Bool_t isRead=ReadFile(list,hstat,cutobj,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat || !cutobj){
    cout<<":-( null pointers..."<<endl;
    return;
  }
  //uncomment when problems with the cut obj
  /* 
  const Int_t nptbins=12;
  Float_t ptbins[nptbins+1]={0,1,2,3,4,5,6,8,12,16,999};
  */

  Float_t* ptbins=cutobj->GetPtBinLimits();
  const Int_t nptbins=cutobj->GetNPtBins();

  TH1F* hnphibins=(TH1F*)list->FindObject("hPhiBins");
  Int_t nphibins=hnphibins->GetNbinsX();
  Double_t phibins[nphibins],dx[nphibins];
  cout<<"(";
  for(Int_t i=0;i<nphibins;i++){
    Double_t width=hnphibins->GetBinWidth(i+1);
    phibins[i]=hnphibins->GetBinLowEdge(i+1)+width/2.;
    dx[i]=width/2.;
    cout<<phibins[i]<<",";
  }
  cout<<")"<<endl;
  cout<<"Number of phi bins = "<<nphibins<<endl;

  cout<<"Number of pt bins = "<<cutobj->GetNPtBins()<<endl;
  Double_t deltapt[nptbinsnew],pt[nptbinsnew],ptledges[nptbinsnew+1];
  for(Int_t j=0;j<nptbinsnew;j++){
    deltapt[j]=(ptbinsnew[j+1]-ptbinsnew[j])/2.;
    pt[j]=ptbinsnew[j]+deltapt[j];
    ptledges[j]=ptbinsnew[j];

    if(j==nptbinsnew-1 && deltapt[j]>10) {
      pt[j]=ptbinsnew[j]+10;
      deltapt[j]=10;
    }
  }
  ptledges[nptbinsnew]=ptledges[nptbinsnew-1]+deltapt[nptbinsnew-1];

  //invariant mass pt integrated, in bins of phi
  TH1F** hmassptint=new TH1F*[nphibins];
  //invariant mass in bins of pt and in/out plane
  TH1F*** hmassnonintphi=new TH1F**[2];

  //init
  for(Int_t j=0;j<2;j++){
    hmassnonintphi[j]=new TH1F*[nptbinsnew];
    for(Int_t i=0;i<nptbinsnew;i++){
      hmassnonintphi[j][i]=0x0;
      cout<<"Init "<<hmassnonintphi[j][i]<<endl;
    }
  }
  for(Int_t j=0;j<nphibins;j++){
    hmassptint[j]=0x0;
  }  


  //output file
  TFile* fout=new TFile(Form("HistoInputv2Calc%s.root",textleg.Data()),"recreate");
  TClonesArray ptblim("TParameter<float>",nptbinsnew+2); //the first is the number of bins
  new(ptblim[0])TParameter<float>("nptbins",nptbinsnew);
  for(Int_t ival=0;ival<nptbinsnew+1;ival++){
    TString name=Form("ptbin%d",ival);
    new(ptblim[ival+1])TParameter<float>(name.Data(),ptbinsnew[ival]);
  }
  fout->cd();
  ptblim.Write();
  cutobj->Write();


  //reading the hostograms and filling hmassnonintphi and hmassptint

  Bool_t isInPlane=kFALSE, isOutOfPlane=kFALSE;
  Int_t indx=-1;

  for(Int_t i=0;i<nphibins;i++){ //loop on phi
    cout<<" -- phi bin "<<i<<"/"<<nphibins-1<<endl;

    for(Int_t j=0;j<nptbins;j++){ //loop on pt
      cout<<"   ** pt bin "<<j<<"/"<<nptbins-1<<endl;
  
      isInPlane=kFALSE, isOutOfPlane=kFALSE;

      TH1F* h=(TH1F*)list->FindObject(Form("hMass_pt%dphi%dcentr%d_%d",j,i,mincentr,mincentr+5));
      if(!h){
  	cout<<"qui Histogram pt "<<j<<", phi "<<i<<", centr "<<mincentr<<"-"<<mincentr+5<<" not found"<<endl;
	//list->ls();
  	continue;
	//return;
      }
      for(Int_t icentr=mincentr+5;icentr<maxcentr;icentr=icentr+5){
	TH1F* hcentr05=(TH1F*)list->FindObject(Form("hMass_pt%dphi%dcentr%d_%d",j,i,icentr,icentr+5));
	if(!hcentr05){
	  cout<<"Histogram pt "<<j<<", phi "<<i<<", centr "<<icentr<<"-"<<icentr+5<<" not found"<<endl;
	  return;
	}
	h->Add(hcentr05);
      }
      Int_t kptnew=FindPtBin(nptbinsnew,ptbinsnew,ptbins[j]);
      cout<<"----*** pt"<<kptnew<<endl;
      if(kptnew==-1) continue;
      //Filling in plane and out of plane
      if((phibins[i]>0 && phibins[i]<= 0.79) || (phibins[i]> 2.36 && phibins[i]<3.14)) {isInPlane=kTRUE; indx=0;}
      if(phibins[i]> 0.79 && phibins[i]<= 2.36) {isOutOfPlane=kTRUE; indx=1;cout<<" IS OUT of plane"<<endl;}
      if(!hmassnonintphi[indx][kptnew]) {

	cout<<"\tpt"<<j<<", phi"<<i<<" in phi "<<i<<" ptnew "<<kptnew<<endl;
	hmassnonintphi[indx][kptnew]=(TH1F*)h->Clone(Form("hMass_ptnw%dphi%d",kptnew,isInPlane? 0 : 1));
	hmassnonintphi[indx][kptnew]->SetTitle(Form("Mass ptnw%d phi%s",kptnew,isInPlane ? "In plane" : "Out of plane"));

      }
      else {
	cout<<"\tadding pt"<<j<<", phi"<<i<<" in phi "<<i<<" ptnew "<<kptnew<<endl;
	hmassnonintphi[indx][kptnew]->Add(h);
      }


      //pt integrated mass plots
      if(!hmassptint[i]) {
	cout<<"Filling histo ptint["<<i<<"] with hMass_pt"<<j<<"phi"<<i<<endl;
	hmassptint[i]=(TH1F*)h->Clone(Form("hphi_%d",i));
	hmassptint[i]->SetTitle(Form("Invariant Mass (Phi %d, Pt integr)",i));
	//cout<<"Entries "<<hmassptint[i]->GetEntries()<<endl;
      }
      else {
	cout<<"Adding histo hMass_pt"<<j<<"phi"<<i<<" to ptint["<<i<<"]"<<endl;
	hmassptint[i]->Add(h);
      }
    } //end pt loop
  } //end phi loop

  //In-plane Out-of-plane anisotropy

  //canvas
  //mass in bins of pt
  TCanvas* cvnewptbins=new TCanvas("cvnewptbins","Some bins in pt, in-plane/out-of-plane",1600,800);
  cvnewptbins->Divide(nptbinsnew,2);
  //chi square per ptbin
  TCanvas* cvchinewptbins=new TCanvas("cvchinewptbins", "Chi vs pt");

  //legend IN-PLANE OUT-OF-PLANE
  TLegend* legsig=new TLegend(0.7,0.6,0.9,0.8,"");
  legsig->SetBorderSize(0);
  legsig->SetFillStyle(0);

  TPaveText* txtoff=0x0;
  if(official){
    gStyle->SetOptTitle(0);
    gStyle->SetFrameBorderMode(0);
    txtoff=new TPaveText(0.1,0.53,0.55,0.83,"NDC");
    txtoff->SetBorderSize(0);
    txtoff->SetFillStyle(0);
    txtoff->AddText(Form("D^{0}#rightarrow K^{-}#pi^{+} %.0f#times10^{6} events",hstat->Integral(0,1)/1e6));
    txtoff->AddText("Pb-Pb @ #sqrt{s_{NN}} = 2.76TeV");
    txtoff->AddText(Form("Centrality %d-%d%s",mincentr,maxcentr,"%"));
    txtoff->SetTextColor(kCyan+3);
  }

  Double_t signal[2][nptbinsnew],background[2][nptbinsnew],significance[2][nptbinsnew];
  Double_t errsgn[2][nptbinsnew],errbkg[2][nptbinsnew],errsignificance[2][nptbinsnew];
  Int_t k=0;
  //TString textinoutpl[2]={"0<#Delta#phi<#pi/4 && 3/4#pi<#Delta#phi<#pi","#pi/4<#Delta#phi<3/4#pi"};


  //histograms
  //chi square
  TH1F** hchivsptnew=new TH1F*[2];


  for(Int_t i=0;i<2;i++){ //in plane/out of plane
    printf("Drawing %s\n",i==0 ? "in plane" : "out of plane");
    hchivsptnew[i]=new TH1F(Form("hchivsptnew%d",i), "Chi square from fit;p_{t} (GeV/c);#chi^{2}",nptbinsnew,ptledges);
    hchivsptnew[i]->SetMarkerColor(i+1);
    hchivsptnew[i]->SetMarkerStyle(20);

     for(Int_t j=0;j<nptbinsnew;j++){
      TPaveText* pvbinphipt=new TPaveText(0.1,0.8,0.7,0.9,"NDC");
      pvbinphipt->SetBorderSize(0);
      pvbinphipt->SetFillStyle(0);
      if(j!=nptbinsnew-1) pvbinphipt->AddText(Form("%.0f<p_{t}<%.0f GeV/c",ptbinsnew[j],ptbinsnew[j+1]));
      else pvbinphipt->AddText(Form("p_{t}>%.0f GeV/c",ptbinsnew[j]));
      k++;

      if(!hmassnonintphi[i][j]) {
	printf("hist %s pt %d not found\n",i==0 ? "In-plane" : "out-of-plane",j);
	continue;
      }

      fout->cd();
      hmassnonintphi[i][j]->Write();

      cvnewptbins->cd(k);

      //Fit mass histograms in bins of pt in-plane or out-of-plane
      AliHFMassFitter fitter(hmassnonintphi[i][j],hmassnonintphi[i][j]->GetBinLowEdge(1),hmassnonintphi[i][j]->GetBinLowEdge(hmassnonintphi[i][j]->GetNbinsX()+1),rebin[j],fittype);

      Bool_t ok=fitter.MassFitter(kFALSE);
      if(ok) {
	fitter.DrawHere(cvnewptbins->cd(k),nsigma,info);

	fitter.Signal(nsigma,signal[i][j],errsgn[i][j]);
	fitter.Background(nsigma,background[i][j],errbkg[i][j]);
	fitter.Significance(nsigma,significance[i][j],errsignificance[i][j]);

	//fill chi square
	hchivsptnew[i]->SetBinContent(j,fitter.GetChiSquare());

      }
      else { //fit failed
	hmassnonintphi[i][j]->Draw();
	hchivsptnew[i]->SetBinContent(j,-1);
      }

      //draw mass
      pvbinphipt->Draw();
     }

     //draw chi square
     cvchinewptbins->cd();
     if(i==0){
       hchivsptnew[i]->Draw("ptext");
       legsig->AddEntry(hchivsptnew[i],"In-plane","p");
     }
     else {
       hchivsptnew[i]->Draw("psamestext");
       legsig->AddEntry(hchivsptnew[i],"Out-on-plane","p");
     }
  }


  //write png of mass as a function of pt in plane and out of plane and chi square
  cvchinewptbins->cd();
  legsig->Draw();
  WriteCanvas(cvchinewptbins,textleg);
  WriteCanvas(cvnewptbins,Form("inout%dptbinscentr%d_%d",nptbinsnew,mincentr,maxcentr));

  //write canvas chi square
  fout->cd();
  cvchinewptbins->Write();

  //graphs of yields
  TGraphErrors** gr=new TGraphErrors*[2];
  TGraphErrors** grb=new TGraphErrors*[2];
  TGraphErrors** grsgf=new TGraphErrors*[2];
  TGraph** grerr=new TGraph*[2];
  TGraph** grberr=new TGraph*[2];
  //canvas for these graphs
  TCanvas* cvsig=new TCanvas("cvsig","Signal vs pt");
  TCanvas* cvsgnf=new TCanvas("cvsgnf","Significance vs pt");
  TCanvas* cvbkg=new TCanvas("cvbkg","Background vs pt");
  TCanvas* cvsiger=new TCanvas("cvsiger","Error on Signal vs pt");
  TCanvas* cvbkger=new TCanvas("cvbkger","Error on Background vs pt");

  for(Int_t i=0;i<2;i++) { //in/out
    gr[i]=new TGraphErrors(0);
    grb[i]=new TGraphErrors(0);
    grsgf[i]=new TGraphErrors(0);
    grerr[i]=new TGraph(0);
    grberr[i]=new TGraph(0);

    for(Int_t j=0;j<nptbinsnew;j++){ //loop on ptbins

      //signal
      gr[i]->SetPoint(j,pt[j],signal[i][j]);
      gr[i]->SetPointError(j,deltapt[j],errsgn[i][j]);
      cout<<j<<" pt "<<pt[j]<<" pm "<<deltapt[j]<<", sgn "<<signal[i][j]<<" pm "<<errsgn[i][j]<<endl;

      //background
      grb[i]->SetPoint(j,pt[j],background[i][j]);
      grb[i]->SetPointError(j,deltapt[j],errbkg[i][j]);
      grsgf[i]->SetPoint(j,pt[j],significance[i][j]);
      grsgf[i]->SetPointError(j,deltapt[j],errsignificance[i][j]);

      //error on signal
      grerr[i]->SetPoint(j,pt[j],errsgn[i][j]/signal[i][j]);
      grberr[i]->SetPoint(j,pt[j],errbkg[i][j]/background[i][j]);
    }

    //titles etc.
    gr[i]->SetTitle("Signal vs p_{t}");
    gr[i]->SetName(Form("sig%s", i==0 ? "inplane" : "outofplane"));
    gr[i]->GetXaxis()->SetTitle("p_{t}");
    gr[i]->GetYaxis()->SetTitle("Signal");
    gr[i]->SetMinimum(0);
    gr[i]->SetMarkerStyle(20);
    gr[i]->SetMarkerColor(i+1);
    
    grb[i]->SetTitle("Background vs p_{t}");
    grb[i]->SetName(Form("bkg%s", i==0 ? "inplane" : "outofplane"));
    grb[i]->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    grb[i]->GetYaxis()->SetTitle("Background");
    grb[i]->SetMinimum(0);
    grb[i]->SetMarkerStyle(20);
    grb[i]->SetMarkerColor(i+1);

    grsgf[i]->SetTitle("Significance vs p_{t};p_{t} (GeV/c);Significance");
    grsgf[i]->SetName(Form("sgf%s", i==0 ? "inplane" : "outofplane"));
    grsgf[i]->SetMinimum(0);
    grsgf[i]->SetMarkerStyle(20);
    grsgf[i]->SetMarkerColor(i+1);

    grerr[i]->SetTitle("#epsilon_{r} Signal");
    grerr[i]->SetName(Form("sigerr%s", i==0 ? "inplane" : "outofplane"));
    grerr[i]->GetYaxis()->SetTitle("#epsilon_{r} Signal");
    grerr[i]->GetXaxis()->SetTitle("p_{t}");
    grerr[i]->SetMinimum(0);
    grerr[i]->SetMarkerStyle(20);
    grerr[i]->SetMarkerColor(i+1);

    grberr[i]->SetTitle("#epsilon_{r} Background");
    grberr[i]->SetName(Form("bkgerr%s", i==0 ? "inplane" : "outofplane"));
    grberr[i]->GetYaxis()->SetTitle("#epsilon_{r} Background");
    grberr[i]->GetXaxis()->SetTitle("p_{t}");
    grberr[i]->SetMinimum(0);
    grberr[i]->SetMarkerStyle(20);
    grberr[i]->SetMarkerColor(i+1);

    //draw signal, bkg, significance, errors
   if(i==0) {
      cvsig->cd();
      gr[i]->Draw("AP");

      cvbkg->cd();
      grb[i]->Draw("AP");

      cvsgnf->cd();
      grsgf[i]->Draw("AP");

      cvsiger->cd();
      grerr[i]->Draw("AP");

      cvbkger->cd();
      grberr[i]->Draw("AP");
    }
    else {
      cvsig->cd();
      gr[i]->Draw("P");

      cvbkg->cd();
      grb[i]->Draw("P");

      cvsgnf->cd();
      grsgf[i]->Draw("P");

      cvsiger->cd();
      grerr[i]->Draw("P");

      cvbkger->cd();
      grberr[i]->Draw("P");
   }

   //write in output file
   fout->cd();
   gr[i]->Write();
   grb[i]->Write();
   grerr[i]->Write();
   grsgf[i]->Write();
  }

  cvsig->cd();
  legsig->Draw();
  cvbkg->cd();
  legsig->Draw();
  cvsgnf->cd();
  legsig->Draw();
  cvsiger->cd();
  legsig->Draw();
  cvbkger->cd();
  legsig->Draw();

  //write png
  WriteCanvas(cvsig,textleg);
  WriteCanvas(cvbkg,textleg);
  WriteCanvas(cvsiger,textleg);
  WriteCanvas(cvbkger,textleg);
  WriteCanvas(cvsgnf,textleg);

  //fit again with sigma obtained from total inv mass (in-plane+out-of-plane)
  Double_t sigmafullmass[nptbinsnew];

  //canvas
  TCanvas* cvmasstotnewptbins=new TCanvas("cvmasstotnewptbins", "Mass in some pt bins in+out",1600,400);
  cvmasstotnewptbins->Divide(nptbinsnew,1);
  TCanvas* cvnewptbins2=new TCanvas("cvnewptbins2", "Mass in some pt bins in/out of plane fixed sigma",1600,800);
  cvnewptbins2->Divide(nptbinsnew,2);

  for(Int_t i=0;i<nptbinsnew;i++){//pt bins

    //pave text
    TPaveText* pvbinpt=new TPaveText(0.1,0.8,0.7,0.9,"NDC");
    pvbinpt->SetBorderSize(0);
    pvbinpt->SetFillStyle(0);
    if(i!= nptbinsnew-1) pvbinpt->AddText(Form("%.0f<p_{t}<%.0f GeV/c",ptbinsnew[i],ptbinsnew[i+1]));
    else if(ptbinsnew[i+1]-ptbinsnew[i] > 10) pvbinpt->AddText(Form("p_{t}>%.0f GeV/c",ptbinsnew[i]));

    //mass histogram integrated in phi
    TH1F* hmassfull=(TH1F*)hmassnonintphi[0][i]->Clone();
    hmassfull->Add(hmassnonintphi[1][i]);

    AliHFMassFitter fitter(hmassfull,hmassfull->GetBinLowEdge(1),hmassfull->GetBinLowEdge(hmassfull->GetNbinsX()+1),rebin[i],fittype);
    Bool_t ok=fitter.MassFitter(kFALSE);
    if(ok){
      fitter.DrawHere(cvmasstotnewptbins->cd(i+1),nsigma,info);
      sigmafullmass[i]=fitter.GetSigma();
    }else{
      cvmasstotnewptbins->cd(i+1);
      hmassfull->Draw();
      sigmafullmass[i]=-2;
    }
    cvmasstotnewptbins->cd(i+1);
    pvbinpt->Draw();


    //in-plane
    AliHFMassFitter fitter0(hmassnonintphi[0][i],hmassnonintphi[0][i]->GetBinLowEdge(1),hmassnonintphi[0][i]->GetBinLowEdge(hmassnonintphi[0][i]->GetNbinsX()+1),rebin[i],fittype);
    if(sigmafullmass[i]!=-2) fitter0.SetFixGaussianSigma(sigmafullmass[i]);
    ok=fitter0.MassFitter(kFALSE);
    if(ok){
      fitter0.DrawHere(cvnewptbins2->cd(i+1),nsigma,info);
      fitter0.Signal(nsigma,signal[0][i],errsgn[0][i]);
      fitter0.Background(nsigma,background[0][i],errbkg[0][i]);
      fitter0.Significance(nsigma,significance[0][i],errsignificance[0][i]);

    }else{
      cvnewptbins2->cd(i+1);
      hmassnonintphi[0][i]->Draw();
    }
    cvnewptbins2->cd(i+1);
    pvbinpt->Draw();


    if(official) {
      cvnewptbins2->cd(i+1);
      TPaveText* txtsignif=new TPaveText(0.1,0.1,0.6,0.2,"NDC");
      txtsignif->SetBorderSize(0);
      txtsignif->SetFillStyle(0);
      txtsignif->SetTextColor(kGreen+3);
      txtsignif->AddText(Form("Sgnf(%.0f#sigma) = %.1f#pm%.1f",nsigma,significance[0][i],errsignificance[0][i]));
 

      txtsignif->Draw();
      TH1F* htmp=(TH1F*)cvnewptbins2->cd(i+1)->FindObject("fhistoInvMass");
      htmp->SetYTitle(Form("Entries/%.3f GeV/c^{2}",htmp->GetBinWidth(3)));
      htmp->SetXTitle("D^{0} invariant mass (GeV/c^{2})");
    }

    if(official) InOutPic(cvnewptbins2->cd(i+1),i);

    //out-of-plane
    AliHFMassFitter fitter1(hmassnonintphi[1][i],hmassnonintphi[1][i]->GetBinLowEdge(1),hmassnonintphi[1][i]->GetBinLowEdge(hmassnonintphi[1][i]->GetNbinsX()+1),rebin[i],fittype);
    if(sigmafullmass[i]!=-2) fitter1.SetFixGaussianSigma(sigmafullmass[i]);
    ok=fitter1.MassFitter(kFALSE);
    if(ok){
      fitter1.DrawHere(cvnewptbins2->cd(nptbinsnew+i+1),nsigma,info);
      fitter1.Signal(nsigma,signal[1][i],errsgn[1][i]);
      fitter1.Background(nsigma,background[1][i],errbkg[1][i]);
      fitter1.Significance(nsigma,significance[1][i],errsignificance[1][i]);
    }else{
      cvnewptbins2->cd(nptbinsnew+i+1);
      hmassnonintphi[1][i]->Draw();
    }
    cvnewptbins2->cd(nptbinsnew+i+1);
    pvbinpt->Draw();

    if(official) {
      cvnewptbins2->cd(nptbinsnew+i+1);
      TPaveText* txtsignif=new TPaveText(0.1,0.1,0.6,0.2,"NDC");
      txtsignif->SetBorderSize(0);
      txtsignif->SetFillStyle(0);
      txtsignif->SetTextColor(kGreen+3);
      txtsignif->AddText(Form("Sgnf(%.0f#sigma) = %.1f#pm%.1f",nsigma,significance[1][i],errsignificance[1][i]));
 

      txtsignif->Draw();
      TH1F* htmp=(TH1F*)cvnewptbins2->cd(nptbinsnew+i+1)->FindObject("fhistoInvMass");
      htmp->SetYTitle(Form("Entries/%.3f GeV/c^{2}",htmp->GetBinWidth(3)));
      htmp->SetXTitle("D^{0} invariant mass (GeV/c^{2})");
    }
    if(official) InOutPic(cvnewptbins2->cd(nptbinsnew+i+1),i);

  }

  //other drawing details
  if(official){
    cvnewptbins2->cd(2);
    txtoff->Draw();
    gROOT->LoadMacro("/home/cbianchi/macros/ALICEPerformance.C");
    gROOT->ProcessLine(Form("ALICEPerformance((TVirtualPad*)%p,\"today\",\"br\")",cvnewptbins2->cd(1)));
  }


  //png mass histograms
  WriteCanvas(cvnewptbins2,Form("inout%dptbinscentr%d_%d",nptbinsnew,mincentr,maxcentr));
  WriteCanvas(cvmasstotnewptbins,Form("totmass%dptbinscentr%d_%d",nptbinsnew,mincentr,maxcentr));
  //graphs of yields
  TGraphErrors** gr2=new TGraphErrors*[2];
  TGraphErrors** grb2=new TGraphErrors*[2];
  TGraphErrors** grsgf2=new TGraphErrors*[2];
  TGraph** grerr2=new TGraph*[2];
  TGraph** grberr2=new TGraph*[2];

  //canvas for these graphs
  TCanvas* cvsig2=new TCanvas("cvsig2","Signal vs pt (fixed sigma)");
  TCanvas* cvbkg2=new TCanvas("cvbkg2","Background vs pt (fixed sigma)");
  TCanvas* cvsgnf2=new TCanvas("cvsgnf2","Significance vs pt (fixed sigma)");
  TCanvas* cvsiger2=new TCanvas("cvsiger2","Error on Signal vs pt (fixed sigma)");
  TCanvas* cvbkger2=new TCanvas("cvbkger2","Error on Background vs pt (fixed sigma)");

  for(Int_t i=0;i<2;i++) { //in/out

    gr2[i]=new TGraphErrors(0);
    grb2[i]=new TGraphErrors(0);
    grsgf2[i]=new TGraphErrors(0);
    grerr2[i]=new TGraph(0);
    grberr2[i]=new TGraph(0);

    for(Int_t j=0;j<nptbinsnew;j++){ //bins of pt

      if(significance[i][j]>3.){
	gr2[i]->SetPoint(j,pt[j],signal[i][j]);
	gr2[i]->SetPointError(j,deltapt[j],errsgn[i][j]);
	grsgf2[i]->SetPoint(j,pt[j],significance[i][j]);
	grsgf2[i]->SetPointError(j,deltapt[j],errsignificance[i][j]);
	grerr2[i]->SetPoint(j,pt[j],errsgn[i][j]/signal[i][j]);
      }
      cout<<j<<" pt "<<pt[j]<<" pm "<<deltapt[j]<<", sgn "<<signal[i][j]<<" pm "<<errsgn[i][j]<<endl;
      grb2[i]->SetPoint(j,pt[j],background[i][j]);
      grb2[i]->SetPointError(j,deltapt[j],errbkg[i][j]);
      grberr2[i]->SetPoint(j,pt[j],errbkg[i][j]/background[i][j]);
    }

    gr2[i]->SetTitle("Signal vs p_{t} (#sigma fixed)");
    gr2[i]->SetName(Form("sig%sfs", i==0 ? "inplane" : "outofplane"));
    gr2[i]->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    gr2[i]->GetYaxis()->SetTitle("Signal");
    gr2[i]->SetMinimum(0);
    gr2[i]->SetMarkerStyle(20);
    gr2[i]->SetMarkerColor(i+1);
    
    grb2[i]->SetTitle("Background vs p_{t} (#sigma fixed)");
    grb2[i]->SetName(Form("bkg%sfs", i==0 ? "inplane" : "outofplane"));
    grb2[i]->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    grb2[i]->GetYaxis()->SetTitle("Background");
    grb2[i]->SetMinimum(0);
    grb2[i]->SetMarkerStyle(20);
    grb2[i]->SetMarkerColor(i+1);

    grsgf2[i]->SetTitle("Significance vs p_{t} (#sigma fixed);p_{t} (GeV/c);Significance");
    grsgf2[i]->SetName(Form("sgf%sfs", i==0 ? "inplane" : "outofplane"));
    grsgf2[i]->SetMinimum(0);
    grsgf2[i]->SetMarkerStyle(20);
    grsgf2[i]->SetMarkerColor(i+1);

    grerr2[i]->SetTitle("#epsilon_{r} Signal (#sigma fixed)");
    grerr2[i]->SetName(Form("sgnerr%sfs", i==0 ? "inplane" : "outofplane"));
    grerr2[i]->GetYaxis()->SetTitle("#epsilon_{r} Signal");
    grerr2[i]->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    grerr2[i]->SetMinimum(0);
    grerr2[i]->SetMarkerStyle(20);
    grerr2[i]->SetMarkerColor(i+1);

    grberr2[i]->SetTitle("#epsilon_{r} Background (#sigma fixed)");
    grberr2[i]->SetName(Form("bkgerr%sfs", i==0 ? "inplane" : "outofplane"));
    grberr2[i]->GetYaxis()->SetTitle("#epsilon_{r} Background");
    grberr2[i]->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    grberr2[i]->SetMinimum(0);
    grberr2[i]->SetMarkerStyle(20);
    grberr2[i]->SetMarkerColor(i+1);

   
    if(i==0) {
      cvsig2->cd();
      gr2[i]->Draw("AP");

      cvbkg2->cd();
      grb2[i]->Draw("AP");

      cvsgnf2->cd();
      grsgf2[i]->Draw("AP");

      cvsiger2->cd();
      grerr2[i]->Draw("AP");

      cvbkger2->cd();
      grberr2[i]->Draw("AP");
    }
    else {
      cvsig2->cd();
      gr2[i]->Draw("P");

      cvbkg2->cd();
      grb2[i]->Draw("P");

      cvsgnf2->cd();
      grsgf2[i]->Draw("P");

      cvsiger2->cd();
      grerr2[i]->Draw("P");

      cvbkger2->cd();
      grberr2[i]->Draw("P");
    }

    //write in output file
    fout->cd();
    gr2[i]->Write();
    grb2[i]->Write();
    grerr2[i]->Write();
    grsgf2[i]->Write();
  }

  cvsig2->cd();
  legsig->Draw();
  cvbkg2->cd();
  legsig->Draw();
  cvsgnf2->cd();
  legsig->Draw();
  cvsiger2->cd();
  legsig->Draw();
  cvbkger2->cd();
  legsig->Draw();

  //write png
  WriteCanvas(cvsig2,textleg);
  WriteCanvas(cvbkg2,textleg);
  WriteCanvas(cvsgnf2,textleg);
  WriteCanvas(cvsiger2,textleg);
  WriteCanvas(cvbkger2,textleg);


  //invariant mass in bins of phi, pt integrated

  Double_t sigvsphiptint[nphibins],errsigvsphiptint[nphibins],sgnfvsphiptint[nphibins],errsgnfvsphiptint[nphibins];

  //canvas
  TCanvas* cvptint=new TCanvas("cvptint","Pt integrated mass plots",1600,600);
  cvptint->Divide(nphibins,1);
  TCanvas* cvphierr=new TCanvas("cvphierr","Error on Signal phi bins");
  TCanvas* cvchiphi=new TCanvas("cvchiphi","Chi as a function of phi");
  TCanvas* cvsgnvsphiptint=new TCanvas("cvsgnvsphiptint","Sigmal vs #Delta #phi");

  //chi square and yields
  TH1F* hchivsphi=(TH1F*)hnphibins->Clone("hchivsphi");
  hchivsphi->SetTitle("#chi^{2} vs #phi;#phi (rad);#chi^{2}");
  hchivsphi->SetMarkerStyle(20);
  TGraphErrors* grphi=new TGraphErrors(0);
  TGraph* grphierr=new TGraph(0);

  for (Int_t i=0;i<nphibins;i++){ //bins of phi
    cvptint->cd(i+1);

    if(hmassptint[i]) {
      AliHFMassFitter fitter(hmassptint[i],hmassptint[i]->GetBinLowEdge(1),hmassptint[i]->GetBinLowEdge(hmassptint[i]->GetNbinsX()+1),rebin[0],fittype);
      Bool_t ok=fitter.MassFitter(kFALSE);

      TPaveText* txtsignif=new TPaveText(0.1,0.1,0.6,0.2,"NDC");
      txtsignif->SetBorderSize(0);
      txtsignif->SetFillStyle(0);
      txtsignif->SetTextColor(kGreen+3);

      if(ok) {
	fitter.DrawHere(cvptint->cd(i+1),nsigma,info);
	fitter.Signal(nsigma,sigvsphiptint[i],errsigvsphiptint[i]);

	fitter.Significance(nsigma,sgnfvsphiptint[i],errsgnfvsphiptint[i]);
	txtsignif->AddText(Form("Sgnf(%.0f#sigma) = %.1f#pm%.1f",nsigma,sgnfvsphiptint[i],errsgnfvsphiptint[i]));

	//fill chi square
	hchivsphi->SetBinContent(i,fitter.GetChiSquare());
	//fill yields vs phi
	grphi->SetPoint(i,phibins[i],sigvsphiptint[i]);
	grphi->SetPointError(i,dx[i],errsigvsphiptint[i]);
	grphierr->SetPoint(i, phibins[i],errsigvsphiptint[i]/sigvsphiptint[i]);
      }
      else {
	//retry fit if it fails
	fitter.Reset();
	fitter.SetHisto(hmassptint[i]);
	fitter.SetRangeFit(hmassptint[i]->GetBinLowEdge(3),hmassptint[i]->GetBinLowEdge(hmassptint[i]->GetNbinsX()-2));
	fitter.RebinMass(rebin[0]*2);
	fitter.SetType(fittype,0);
	ok=fitter.MassFitter(kFALSE);
	if(ok){
	  fitter.DrawHere(cvptint->cd(i+1),nsigma,info);
	  TPaveText* pvnote=new TPaveText(0.1,0.7,0.8,0.9,"NDC");
	  pvnote->SetBorderSize(0);
	  pvnote->SetFillStyle(0);
	  pvnote->AddText("Refitted w/ different range and rebin");
	  cvptint->cd(i+1);
	  pvnote->Draw();

	  fitter.Signal(nsigma,sigvsphiptint[i],errsigvsphiptint[i]);

	  //fill chi square
	  hchivsphi->SetBinContent(i,fitter.GetChiSquare());
	  //fill yields vs phi
	  grphi->SetPoint(i,phibins[i],sigvsphiptint[i]);
	  grphi->SetPointError(i,dx[i],errsigvsphiptint[i]);
	  grphierr->SetPoint(i, phibins[i],errsigvsphiptint[i]/sigvsphiptint[i]);
	} else hmassptint[i]->Draw(); //fit didn't work at all
      }
      cvptint->cd(i+1);

      if(official){
	txtsignif->Draw();
	//PhiBinPic(cvptint->cd(i+1),i,"tr");
      }

    }

  }

  cvchiphi->cd();
  hchivsphi->Draw("ptext");

  grphi->SetTitle("Signal vs #phi");
  grphi->SetName("gsigvsphi");
  grphi->GetXaxis()->SetTitle("#phi (rad)");
  grphi->GetYaxis()->SetTitle("Signal");
  grphi->SetMinimum(0);
  grphi->SetMarkerStyle(20);
  cvsgnvsphiptint->cd();
  grphi->Draw("AP");

  grphierr->SetTitle("#epsilon_{r} Signal vs #phi");
  grphierr->SetName("gsigerrvsphi");
  grphierr->GetXaxis()->SetTitle("#phi (rad)");
  grphierr->GetYaxis()->SetTitle("#epsilon_{r} Signal");
  grphierr->SetMinimum(0);
  grphierr->SetMarkerStyle(20);
  cvphierr->cd();
  grphierr->Draw("AP");

  //fit dN/dDeltaphi
  TF1 *flowFunc = new TF1("flow","[0]*(1.+2.*[1]*TMath::Cos(2.*x))");
  grphi->Fit(flowFunc);

  TPaveText* pvv2=new TPaveText(0.1,0.1,0.6,0.2,"NDC");
  pvv2->SetBorderSize(0);
  pvv2->SetFillStyle(0);
  pvv2->AddText(Form("v'_{2} = %.3f#pm%.3f",flowFunc->GetParameter(1),flowFunc->GetParError(1)));

  cvsgnvsphiptint->cd();
  pvv2->Draw();

  //write chi square, yields and errors in output file
  fout->cd();
  cvchiphi->Write();
  grphierr->Write();
  grphi->Write();

  //other drawing details
  if(official){
    cvptint->cd(1);
    txtoff->Draw();
    gROOT->LoadMacro("/home/cbianchi/macros/ALICEPerformance.C");
    gROOT->ProcessLine(Form("ALICEPerformance((TVirtualPad*)%p,\"today\",\"br\")",cvptint->cd(1)));
  }

  //write png mass histos yields and error
  WriteCanvas(cvptint,textleg);
  WriteCanvas(cvphierr,textleg);
  WriteCanvas(cvsgnvsphiptint,textleg);

  //statistics
  TCanvas* cst=new TCanvas("cst","Stat");
  cst->SetGridy();
  cst->cd();
  hstat->Draw("htext0");
  WriteCanvas(cst,textleg);


  //save event plane info in the output file
  //fix to the latest version of the task!!!!!
  TH1F* hresos=(TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",mincentr,mincentr+5));
  /*
  TH2F* hMphis=(TH2F*)list->FindObject(Form("hMphicentr%d_%d",mincentr,mincentr+5));
  TH2F* hMc2phis=(TH2F*)list->FindObject(Form("hMc2phicentr%d_%d",mincentr,mincentr+5));
  */
  TString buildname=Form("centr%d_%d",mincentr,maxcentr);
  hresos->SetName(Form("hEvPlaneReso%s",buildname.Data()));
  /*
  hMphis->SetName(Form("hMphi%s",buildname.Data()));
  hMc2phis->SetName(Form("hMc2phi%s",buildname.Data()));
  */
  for(Int_t icentr=mincentr+5;icentr<=maxcentr;icentr=icentr+5){
    TH1F* hreso=(TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",icentr,icentr+5));
    hresos->Add(hreso);
    // for(Int_t i=0;i<nptbins;i++){
    // TH2F* hMphi=(TH2F*)list->FindObject(Form("hMphicentr_pt%d%d_%d",icentr,icentr+5));
    // TH2F* hMc2phi=(TH2F*)list->FindObject(Form("hMc2phi_pt%dcentr%d_%d",icentr,icentr+5));
    // hMphis->Add(hMphi);
    // hMc2phis->Add(hMc2phi);
    // }
  }

  fout->cd();
  hresos->Write();
  /*
  hMphis->Write();
  hMc2phis->Write();
  */
}

void Dmesonsv2InOutAnisotropy(Int_t mincentr,Int_t maxcentr,Bool_t fixedsigma=kTRUE){

  //Compute v2 with the method of the anisotropy between in-plane and out-of plane. Needs the output file from the previous method
  //Author: Francesco Prino prino@to.infn.it

  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptTitle(0);

  TString name=Form("centr%d-%d",mincentr,maxcentr),
    nameh=Form("centr%d_%d",mincentr,maxcentr);

  TFile* fil=new TFile(Form("HistoInputv2Calc%s.root",name.Data()));
  if(!fil->IsOpen()){
    cout<<"Input file not found"<<endl;
    return;
  }
  
 
  //histogram resolution with subevents
  TH1F* hResolSubAB=(TH1F*)fil->Get(Form("hEvPlaneReso%s",nameh.Data()));
  hResolSubAB->GetXaxis()->SetTitle("cos[2(#psi_{A}-#psi_{B})]");

  //graphs yield
  TGraphErrors* gsigin=(TGraphErrors*)fil->Get(Form("siginplane%s",fixedsigma ? "fs" : ""));
  TGraphErrors* gsigout=(TGraphErrors*)fil->Get(Form("sigoutofplane%s",fixedsigma ? "fs" : ""));

  TCanvas* c1=new TCanvas("c1","EvPlaneResol");
  hResolSubAB->Draw();
  Double_t resolSub=TMath::Sqrt(hResolSubAB->GetMean());
  printf("Resolution on sub event  = %.4f\n",resolSub);
  Double_t chisub=FindChi(resolSub,1);
  printf("Chi Subevent             = %.4f\n",chisub);
  Double_t chifull=chisub*TMath::Sqrt(2);
  printf("Chi Full Event           = %.4f\n",chifull);
  Double_t resolFull=ComputeResol(resolSub,1);
  printf("Resolution on full event = %.4f\n",resolFull);
  TLatex* tres=new TLatex(0.15,0.7,Form("Resolution on full event = %.4f\n",resolFull));
  tres->SetNDC();
  tres->Draw();
  c1->SaveAs(Form("EvPlaneResol%s.png",nameh.Data()));

  TCanvas* c2=new TCanvas("c2",Form("Yields %s",fixedsigma ? "fixed sigma" : ""));
  c2->cd();
  gsigin->Draw("AP");
  gsigout->Draw("Psame");

  Int_t nPtBins=gsigin->GetN();

  //graphs results
  TGraphErrors *gAnis=new TGraphErrors(0);
  gAnis->SetName(Form("gAnisotropy%s",nameh.Data()));
  gAnis->SetTitle("(Nin-Nout)/(Nin+Nout)");

  TGraphErrors *gv2obs=new TGraphErrors(0);
  gv2obs->SetName(Form("gv2obs%s",nameh.Data()));
  gv2obs->SetTitle("v2obs");

  TGraphErrors *gv2=new TGraphErrors(0);
  gv2->SetName(Form("gv2%s",nameh.Data()));
  gv2->SetTitle("v2");

  TGraph *gv2err=new TGraph(0);
  gv2err->SetName(Form("gv2err%s",nameh.Data()));

  for(Int_t iBin=0; iBin<nPtBins; iBin++){
    Double_t pt,nIn,ept,enIn;
    Double_t nOut,enOut;
    gsigin->GetPoint(iBin,pt,nIn);
    ept=gsigin->GetErrorX(iBin);
    enIn=gsigin->GetErrorY(iBin);
    gsigout->GetPoint(iBin,pt,nOut);
    enOut=gsigout->GetErrorY(iBin);
    if(nIn<1e-6 || nOut<1e-6) continue;

    Double_t anis=(nIn-nOut)/(nIn+nOut);
    Double_t eAnis=TMath::Sqrt(enIn*enIn+enOut*enOut)/(nIn+nOut);
    gAnis->SetPoint(iBin,pt,anis);
    gAnis->SetPointError(iBin,ept,eAnis);

    Double_t v2obs=anis*TMath::Pi()/4.;
    Double_t ev2obs=eAnis*TMath::Pi()/4.;
    gv2obs->SetPoint(iBin,pt,v2obs);
    gv2obs->SetPointError(iBin,ept,ev2obs);

    Double_t v2=v2obs/resolFull;
    Double_t ev2=ev2obs/resolFull;
    gv2->SetPoint(iBin,pt,v2);
    gv2->SetPointError(iBin,ept,ev2);

    gv2err->SetPoint(iBin,pt,ev2/TMath::Abs(v2));
  }

  TCanvas* c3=new TCanvas("c3","Anisotropy");
  gAnis->SetMarkerStyle(20);
  gAnis->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  gAnis->GetYaxis()->SetTitle("(N_{IN-PLANE}-N_{OUT-OF-PLANE})/(N_{IN-PLANE}+N_{OUT-OF-PLANE})");
  gAnis->GetYaxis()->SetTitleOffset(1.1);
  gAnis->Draw("AP");
  gAnis->Draw("Psame");
  c3->SaveAs(Form("Anisotropy%s%s.png",name.Data(),fixedsigma ? "fs" : ""));

  TCanvas* c4=new TCanvas("c4","v2obs");
  gv2obs->SetMarkerStyle(20);
  gv2obs->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  gv2obs->GetYaxis()->SetTitle("v_{2}^{obs}");
  gv2obs->GetYaxis()->SetTitleOffset(1.1);
  gv2obs->Draw("AP");
  gv2obs->Draw("Psame");
  c4->SaveAs(Form("v2obsD%s%s.png",name.Data(),fixedsigma ? "fs" : ""));

  TCanvas* c5=new TCanvas("c5","v2");
  gv2->SetMarkerStyle(20);
  gv2->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  gv2->GetYaxis()->SetTitle("v_{2}");
  gv2->GetYaxis()->SetTitleOffset(1.1);
  gv2->Draw("AP");
  gv2->Draw("Psame");
  c5->SaveAs(Form("v2D%s%s.png",name.Data(),fixedsigma ? "fs" : ""));

  TCanvas* c6=new TCanvas("c6","error v2");
  gv2err->SetTitle("Relative error v2;p_{t} (GeV/c);#epsilon_{r} v2");
  gv2err->SetMarkerStyle(20);
  gv2err->GetYaxis()->SetTitleOffset(1.1);
  c6->cd();
  gv2err->Draw("AP");
  c6->SaveAs(Form("Errv2%s%s.png",name.Data(),fixedsigma ? "fs" : ""));

  TFile* fout=new TFile(Form("Outputv2%s.root",name.Data()),"recreate");

  fout->cd();
  gAnis->Write();
  gv2obs->Write();
  gv2->Write();
  gv2err->Write();
  c1->Write(); //canvas resolution
}

//methods needed inside Dmesonsv2InOutAnisotropy
//Author: Francesco Prino prino@to.infn.it

Double_t ComputeResol(Double_t resSub, Int_t k){
  Double_t chisub=FindChi(resSub,k);
  Double_t chifull=chisub*TMath::Sqrt(2);
  if(k==1) return ResolK1(chifull);
  else if(k==2) return Pol(chifull,2);
  else{
    printf("k should be <=2\n");
    return 1.;
  }
}

Double_t FindChi(Double_t res, Int_t k){
  Double_t x1=0;
  Double_t x2=15;
  Double_t y1,y2;
  if(k==1){
    y1=ResolK1(x1)-res;
    y2=ResolK1(x2)-res;
  }
  else if(k==2){
    y1=Pol(x1,2)-res;
    y2=Pol(x2,2)-res;
  }
  else return -1;

  if(y1*y2>0.) return -1;
  if(y1==0.) return y1;
  if(y2==0.) return y2;
  Double_t xmed,ymed;
  Int_t jiter=0;
  while((x2-x1)>0.0001){
    xmed=0.5*(x1+x2);
    if(k==1){
      y1=ResolK1(x1)-res;
      y2=ResolK1(x2)-res;
      ymed=ResolK1(xmed)-res;
    }
    else if(k==2){
      y1=Pol(x1,2)-res;
      y2=Pol(x2,2)-res;
      ymed=Pol(xmed,2)-res;
    }
    else return -1;
    if((y1*ymed)<0) x2=xmed;
    if((y2*ymed)<0)x1=xmed;
    if(ymed==0) return xmed;
    jiter++;
  }
  return 0.5*(x1+x2);
}

Double_t Pol(Double_t x, Int_t k){
  Double_t c[5];
  if(k==1){ 
    c[0]=0.626657; c[1]=0.; c[2]=-0.09694; c[3]=0.02754; c[4]=-0.002283;
  }
  else if(k==2){
    c[0]=0.; c[1]=0.25; c[2]=-0.011414; c[3]=-0.034726; c[4]=0.006815;
  }
  else return -1;
  return c[0]*x+c[1]*x*x+c[2]*x*x*x+c[3]*x*x*x*x+c[4]*x*x*x*x*x;
}

Double_t ResolK1(Double_t x){
  return TMath::Sqrt(TMath::Pi()/8)*x*TMath::Exp(-x*x/4)*(TMath::BesselI0(x*x/4)+TMath::BesselI1(x*x/4));
}


void DrawEventPlane(Int_t mincentr=0,Int_t maxcentr=0,TString partname="D0",TString path="./"){

 //draw the histograms correlated to the event plane

  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  //gStyle->SetPalette(1);
  gStyle->SetOptFit(1);

  TString listname="coutputv2";
  TList* list;
  TH1F * hstat;
  AliRDHFCuts* cutobj;

  Bool_t isRead=ReadFile(list,hstat,cutobj,listname,partname,path);
  if(!isRead) return;
  if(!list || !hstat || !cutobj){
    cout<<":-( null pointers..."<<endl;
    return;
  }
 
  if(!(mincentr==0 && maxcentr==0)){ //draw the total in mincentr-maxcentr
    TString suffixcentr=Form("centr%d_%d",mincentr,maxcentr);
    TH1F* hevpls=(TH1F*)list->FindObject(Form("hEvPlanecentr%d_%d",mincentr,mincentr+5));
    hevpls->SetName(Form("hEvPlane%s",suffixcentr.Data()));
    hevpls->SetTitle(Form("Event Plane angle %s",suffixcentr.Data()));
    TH1F* hevplresos=(TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",mincentr,mincentr+5));
    hevplresos->SetName(Form("hEvPlaneReso%s",suffixcentr.Data()));
    hevplresos->SetTitle(Form("Event Plane Resolution %s",suffixcentr.Data()));

    for(Int_t icentr=mincentr+5;icentr<=maxcentr;icentr=icentr+5){
      TH1F* h=(TH1F*)list->FindObject(Form("hEvPlanecentr%d_%d",icentr,icentr+5));
      if(h)hevpls->Add(h);
      else cout<<"skipping ev plane "<<icentr<<"_"<<icentr+5<<endl;
      h=(TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",icentr,icentr+5));
      if(h)hevplresos->Add(h);
      else cout<<"skipping ev pl reso "<<icentr<<"_"<<icentr+5<<endl;
    }

    TCanvas* cvtotevpl=new TCanvas("cvtotevpl",Form("Ev plane %s",suffixcentr.Data()));
    cvtotevpl->cd();
    hevpls->Draw();
    hevpls->Fit("pol0");

    TCanvas* cvtotevplreso=new TCanvas("cvtotevplreso",Form("Ev plane Resolution %s",suffixcentr.Data()));
    cvtotevplreso->cd();
    hevplresos->Draw();
    Double_t resolSub=TMath::Sqrt(hevplresos->GetMean());
    Double_t resolFull=ComputeResol(resolSub,1);

    TPaveText* pvreso=new TPaveText(0.1,0.1,0.6,0.2,"NDC");
    pvreso->SetBorderSize(0);
    pvreso->SetFillStyle(0);
    pvreso->AddText(Form("Resolution on full event = %.4f\n",resolFull));
    pvreso->Draw();

    TFile* fout=new TFile(Form("EvPlanecentr%d-%d.root",mincentr,maxcentr),"recreate");
    fout->cd();
    hevpls->Write();
    hevplresos->Write();
  }
  else{ //draw all in 5% centrality bins
    for(Int_t i=0;i<list->GetEntries();i++){
      TH1F* h=(TH1F*)list->At(i);
      if(!h){
	cout<<"Histogram "<<i<<" not found"<<endl;
	continue;
      }
      TString hname=h->GetName();
      if(hname.Contains("EvPlane")){
	TCanvas* cv=new TCanvas(Form("cv%s",hname.Data()),hname.Data());
	cv->cd();
	h->Draw();

	if(!hname.Contains("Reso")){
	  h->Fit("pol0");
	  
	}else{
	  Double_t resolSub=TMath::Sqrt(h->GetMean());
	  Double_t resolFull=ComputeResol(resolSub,1);

	  TPaveText* pvreso=new TPaveText(0.1,0.1,0.6,0.2,"NDC");
	  pvreso->SetBorderSize(0);
	  pvreso->SetFillStyle(0);
	  pvreso->AddText(Form("Resolution on full event = %.4f\n",resolFull));
	  pvreso->Draw();

	}
      }
    }
  }

}


void DmesonsFlowAnalysisvsCentrality(Int_t mincentr, Int_t maxcentr, Int_t step=10){
  //needs outputs from the other methods HistoInputv2Calccentr%d-%d.root,EvPlanecentr%d-%d.root, Outputv2centr%d-%d.root

  //min(max)centr= min(max) value of centrality among those files - note that if there are "holes" they are skipped
  //step= step in centrality between mincentr and maxcentr - default 10%

  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //number of files expected
  Int_t n=(maxcentr-mincentr)/step;
  cout<<"Trying to read "<<n<<" files"<<endl;
  
  TString filename="",gname="";
  Int_t icentr=0;

  //v2
  TCanvas* cvv2=new TCanvas("cvv2","v_{2}");
  TCanvas* cvv2obs=new TCanvas("cvv2obs","v_{2}'");

  //Event plane distribution and resolution
  TCanvas* cvevplres=new TCanvas("cvevplres","Resolution");
  cvevplres->SetLogy();
  TH1F* hresolFull=new TH1F("hresolFull","Resolution vs centrality;centrality;resolution",n,mincentr,maxcentr);
  TCanvas* cvres=new TCanvas("cvres", "Resolution vs centrality");
  TCanvas* cvevplane=new TCanvas("cvevplane", "Event plane vs centrality");

  //v2 in phi bins 
  TCanvas* cvfitvsphi=new TCanvas("cvfitvsphi","Fit flow vs phi");
  TH1F* hv2fcorr=new TH1F("hv2fcorr","v2 from fit corrected;centrality;v_{2}",n,mincentr,maxcentr);
  hv2fcorr->SetMarkerStyle(30);
  hv2fcorr->SetMarkerSize(2);
  TH1F* hv2corr=new TH1F("hv2corr","v2 from asym corrected;centrality;v_{2}",n,mincentr,maxcentr);
  hv2corr->SetMarkerStyle(29);
  hv2corr->SetMarkerSize(1.7);
  hv2corr->SetMarkerColor(2);
  hv2corr->SetLineColor(2);
 
  TCanvas* cvv2cor=new TCanvas("v2corr", "v2 from asymmetry corrected vs centrality");

  TLegend* leg=new TLegend(0.7,0.6,0.9,0.8,"Centrality %");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  TLegend* leg2=new TLegend(0.5,0.6,0.9,0.9,"");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  TLegend* legresoval=new TLegend(0.2,0.6,0.5,0.9,"Full resolution");
  legresoval->SetBorderSize(0);
  legresoval->SetFillStyle(0);

  Int_t nread=0;
  for(Int_t i=0;i<n;i++){
    icentr=mincentr+i*step;
    filename=Form("Outputv2centr%d-%d.root",icentr,icentr+step);
    TFile* f=new TFile(filename);
    if(!f->IsOpen()){
      cout<<filename.Data()<<" not found"<<endl;
      continue;
    }
    nread++;
    gname=Form("gv2centr%d_%d",icentr,icentr+step);
    TGraphErrors* grv2=(TGraphErrors*)f->Get(gname);
    if(!grv2){
      cout<<gname<<" not found in "<<filename<<endl;
      continue;
    }
    grv2->GetYaxis()->SetRangeUser(-0.2,0.4);
    grv2->SetMarkerColor(i+1);
    grv2->SetLineColor(i+1);
    cvv2->cd();
    if(nread==1)grv2->Draw("APL");
    else grv2->Draw("PL");

    leg->AddEntry(grv2,Form("%d-%d",icentr,icentr+step),"p");

    gname=Form("gv2obscentr%d_%d",icentr,icentr+step);
    TGraphErrors* grv2obs=(TGraphErrors*)f->Get(gname);
    if(!grv2obs){
      cout<<gname<<" not found in "<<filename<<endl;
      continue;
    }
    grv2obs->GetYaxis()->SetRangeUser(-0.2,0.4);
    grv2obs->SetMarkerColor(i+1);
    grv2obs->SetLineColor(i+1);
    cvv2obs->cd();
    if(nread==1)grv2obs->Draw("APL");
    else grv2obs->Draw("PL");

    //draw ev plane reso
    filename=Form("EvPlanecentr%d-%d.root",icentr,icentr+step);
    TFile* f2=new TFile(filename);
    if(!f2->IsOpen()){
      cout<<filename.Data()<<" not found"<<endl;
      continue;
    }
    TString hname=Form("hEvPlaneResocentr%d_%d",icentr,icentr+step);
    TH1F* hevplres=(TH1F*)f2->Get(hname);
    if(!hevplres){
      cout<<hname<<" not found in "<<filename<<endl;
      continue;
    }
    cvevplres->cd();
    hevplres->SetLineColor(i+1);
    if(nread==1) hevplres->Draw();
    else hevplres->Draw("sames");

    Double_t resolSub=TMath::Sqrt(hevplres->GetMean());
    Double_t resolFull=ComputeResol(resolSub,1);
    hresolFull->SetBinContent(i+1,resolFull);
    TLegendEntry* lentry=legresoval->AddEntry(hevplres,Form("%.2f",resolFull),"");
    lentry->SetTextColor(i+1);
    cvevplres->cd();
    legresoval->Draw();

    hname=Form("hEvPlanecentr%d_%d",icentr,icentr+step);
    TH1F* hevpl=(TH1F*)f2->Get(hname);
    if(!hevpl){
      cout<<hname<<" not found in "<<filename<<endl;
      continue;
    }
    hevpl->Scale(1./hevpl->Integral());
    cvevplane->cd();
    hevpl->SetLineColor(i+1);
    if(nread==1) hevpl->Draw();
    else hevpl->Draw("sames");

 
    filename=Form("HistoInputv2Calccentr%d-%d.root",icentr,icentr+step);
    TFile* f3=new TFile(filename);
    if(!f3->IsOpen()){
      cout<<filename.Data()<<" not found"<<endl;
      continue;
    }
    gname="gsigvsphi";
    TGraphErrors* grdNdphi=(TGraphErrors*)f3->Get(gname);
    if(!grdNdphi){
      cout<<gname<<" not found in "<<filename<<endl;
      continue;
    }
    TF1* ffit=(TF1*)grdNdphi->FindObject("flow");
    ffit->SetName(Form("flow%d_%d",icentr,icentr+step));
    ffit->SetLineColor(i+1);
    // Double_t xmin,xmax;
    // ffit->GetRange(xmin,xmax);
    //ffit->Scale(1./ffit->Integral(xmin,xmax));
    cvfitvsphi->cd();
    if(nread==1){
      ffit->SetMinimum(0);
      ffit->Draw();
    }
    else ffit->Draw("sames");

    Double_t v2obsf=ffit->GetParameter(1);
    Double_t v2obsferr=ffit->GetParError(1);
    Double_t v2fcorr=v2obsf/resolFull;
    Double_t v2fcorrerr=v2obsferr/resolFull;
    hv2fcorr->SetBinContent(i+1,v2fcorr);
    hv2fcorr->SetBinError(i+1,v2fcorrerr);

    gname="siginplane";
    TGraphErrors* gsigin=(TGraphErrors*)f3->Get(gname);
    gname="sigoutofplane";
    TGraphErrors* gsigout=(TGraphErrors*)f3->Get(gname);
    if(!gsigin || !gsigout){
      cout<<"graphs in/out not found"<<endl;
      continue;
    }
    Int_t nbins=((TParameter<float>*)f3->Get("nptbins"))->GetVal();
    Double_t countin=0,countout=0,e2in=0,e2out=0;
    for(Int_t k=0;k<nbins;k++){
      Double_t x,y,ey;
      gsigin->GetPoint(k,x,y);
      ey=gsigin->GetErrorY(k);
      e2in+=ey*ey;
      countin+=y;
      gsigout->GetPoint(k,x,y);
      ey=gsigin->GetErrorY(k);
      e2out+=ey*ey;
      countout+=y;
    }
    if(countin<1e-6 || countout<1e-6) continue;
    Double_t ein=TMath::Sqrt(e2in),eout=TMath::Sqrt(e2out);
    Double_t anis=(countin-countout)/(countin+countout);
    Double_t eAnis=TMath::Sqrt(ein*ein+eout*eout)/(countin+countout);
    Double_t v2obs=anis*TMath::Pi()/4.;
    Double_t ev2obs=eAnis*TMath::Pi()/4.;
    Double_t v2corr=v2obs/resolFull;
    Double_t v2correrr=ev2obs/resolFull;
    hv2corr->SetBinContent(i+1,v2corr);
    hv2corr->SetBinError(i+1,v2correrr);

  }

  cvv2->cd();
  leg->Draw();
  cvv2obs->cd();
  leg->Draw();
  cvevplres->cd();
  leg->Draw();
  cvevplane->cd();
  leg->Draw();

  cvres->cd();
  hresolFull->Draw();

  cvfitvsphi->cd();
  leg->Draw();

  cvv2cor->cd();
  hv2corr->SetLineWidth(2);
  hv2fcorr->SetLineWidth(2);
  hv2corr->GetYaxis()->SetRangeUser(-0.2,0.4);
  hv2corr->Draw();
  hv2fcorr->Draw("sames");
  leg2->AddEntry(hv2fcorr,"v_{2} from dN/d#phi","l");
  leg2->AddEntry(hv2corr,"v_{2} from asymmetry","l");
  leg2->Draw();

  cvv2->SaveAs(Form("v2Dcentr%d-%d.png",mincentr, maxcentr));
  cvv2obs->SaveAs(Form("v2obsDcentr%d-%d.png",mincentr, maxcentr));
  cvevplres->SaveAs(Form("evplaneresocentr%d-%d.png",mincentr, maxcentr));
  cvevplane->SaveAs(Form("evplanedistrcentr%d-%d.png",mincentr, maxcentr));
  cvres->SaveAs(Form("evplresovscentr%d-%d.png",mincentr, maxcentr));
  cvfitvsphi->SaveAs(Form("flowvsphicentr%d-%d.png",mincentr, maxcentr));
  cvv2cor->SaveAs(Form("cmpv2centr%d-%d.png",mincentr, maxcentr));


}

void ReflectdNdphiPoints(Int_t mincentr,Int_t maxcentr){

  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetOptTitle(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0);

  TString textleg=Form("centr%d-%d",mincentr,maxcentr);
  TFile* fin=new TFile(Form("HistoInputv2Calc%s.root",textleg.Data()));
  if(!fin->IsOpen()){
    printf("HistoInputv2Calc%s.root not found",textleg.Data());
    return;
  }
  TGraphErrors* grphi=(TGraphErrors*)fin->Get("gsigvsphi");
  if(!grphi){
    cout<<"gsigvsphi not found in "<<fin->GetName()<<endl;
    return;
  }

  TGraphErrors* grphidouble=new TGraphErrors(0);
  Int_t n=grphi->GetN();
  //Int_t doublen=n*2;
  Double_t xn=0,yn=0;
  grphi->GetPoint(n-1,xn,yn);
  xn+=grphi->GetErrorX(n-1);
  cout<<"xn = "<<xn<<endl;
  for(Int_t i=0;i<n;i++){
    Double_t x,y;
    grphi->GetPoint(i,x,y);
    //x+=xn;
    x=(-1)*x;
    grphidouble->SetPoint(i,x,y);
    grphidouble->SetPointError(i,grphi->GetErrorX(i),grphi->GetErrorY(i));
  }
  grphidouble->SetMarkerStyle(24);
  grphidouble->SetName(Form("%sReflected",grphi->GetName()));
  grphidouble->SetTitle(grphi->GetTitle());
  grphidouble->GetXaxis()->SetRangeUser(0,6.5);

  TMultiGraph* grtot=new TMultiGraph();
  grtot->Add(grphi,"P");
  grtot->Add(grphidouble,"P");
  grtot->SetNameTitle("grphi2pi","dN/d#Delta#phi;#Delta#phi (rad);dN/d#Delta#phi");
  TPaveText* pv=new TPaveText(0.3,0.7,0.7,0.9,"NDC");
  pv->SetBorderSize(0);
  pv->SetFillStyle(0);

  TPaveText* pvcentr=new TPaveText(0.1,0.7,0.5,0.8,"NDC");
  pvcentr->SetBorderSize(0);
  pvcentr->SetFillStyle(0);
  pvcentr->AddText(Form("Centrality %d-%d",mincentr,maxcentr));

  TPaveText* pvwp=new TPaveText(0.5,0.1,0.9,0.3,"NDC");
  pvwp->SetBorderSize(0);
  pvwp->SetFillStyle(0);
  pvwp->SetTextColor(kRed);
  pvwp->SetTextFont(65);
  pvwp->AddText("ALICE Work in progress");

  pv->AddText("Open points are reflected");
  TCanvas *cvphidouble=new TCanvas("cvphidouble","dN/dphi in 2*pi");
  cvphidouble->cd();
  grtot->Draw("AP");
  pv->Draw();
    
  TCanvas *cvphi=new TCanvas("cvphi","dN/dphi");
  cvphi->cd();
  grphi->SetNameTitle(grphi->GetName(),"dN/d#Delta#phi;#Delta#phi (rad);dN/d#Delta#phi");
  grphi->Draw("AP");
  pvcentr->Draw();
  pvwp->Draw();
}
