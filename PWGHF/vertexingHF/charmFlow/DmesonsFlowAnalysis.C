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
#include <TGraphAsymmErrors.h>

#include <AliHFMassFitter.h>
#include "AliRDHFCutsD0toKpi.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include <AliVertexingHFUtils.h>

#endif

//methods for the analysis of AliAnalysisTaskSEHFv2 output
//Authors: Chiara Bianchin cbianchi@pd.infn.it
//         Giacomo Ortona  ortona@to.infn.it
//         Francesco Prino prino@to.infn.it

//global variables to be set
Bool_t gnopng=kTRUE; //don't save in png format (only root and eps)
const Int_t nptbinsnew=3;
Float_t ptbinsnew[nptbinsnew+1]={3,5,8,12.};
Int_t fittype=0;
Int_t rebin[nptbinsnew]={4,4,4};
Double_t nsigma=3;
const Int_t nphibins=4;
Float_t phibinslim[nphibins+1]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};
Int_t minPtBin[nptbinsnew]={-1,-1,-1};
Int_t maxPtBin[nptbinsnew]={-1,-1,-1};
Double_t mass;
//methods
//Bool_t ReadFile(TList* &list,TH1F* &hstat,AliRDHFCuts* &cutobj,TString listname,TString partname,TString path="./",TString filename="AnalysisResults.root");
Int_t FindPtBin(Int_t nbins, Float_t* array,Float_t value);
// void InOutPic(TVirtualPad *c,Int_t inout=0,TString where="tr");//inout: 0=IN, 1=OUT 
// void PhiBinPic(TVirtualPad *c,Int_t angle,TString where);
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs,Bool_t inoutanis);
void DrawEventPlane(Int_t mincentr=0,Int_t maxcentr=0,TString filename="AnalysisResults.root",TString dirname="PWG3_D2H_HFv2",TString listname="coutputv2");
void DrawEventPlane(TList *list,Int_t mincentr=0,Int_t maxcentr=0);
//Aggiungere <pt> method

//methods implementation
//________________________________________________________________________________
TList *LoadMassHistos(TList *inputlist,Int_t minCent,Int_t maxCent,Bool_t inoutanis){
  // printf("Start load histos\n");
  //  const Int_t nptbins=cutobj->GetNPtBins();
  TList *outlist = new TList();
  outlist->SetName("azimuthalhistoslist");

  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;

  //Create 2D histogram in final pt bins
  for(Int_t iFinalPtBin=0; iFinalPtBin<nptbinsnew; iFinalPtBin++){
    for(Int_t iphi=0;iphi<nphi;iphi++){
      TH1F *hMass=0x0;//=new TH1F();
      for(Int_t iPtBin=minPtBin[iFinalPtBin]; iPtBin<=maxPtBin[iFinalPtBin];iPtBin++){
	for(Int_t iHisC=minCent; iHisC<=maxCent-5; iHisC+=5){    
	  TString hisname=Form("hMphi_pt%dcentr%d_%d",iPtBin,iHisC,iHisC+5);
	  TH2F* htmp=(TH2F*)inputlist->FindObject(hisname.Data());
	  Int_t startX=htmp->FindBin(phibinslim[iphi]);
	  Int_t endX=htmp->FindBin(phibinslim[iphi+1]);
	  TH1F *h1tmp;
	  if(inoutanis){
	    if(iphi==0){
	      h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi0",iPtBin),htmp->FindBin(0),htmp->FindBin(TMath::Pi()/4.));
	      h1tmp->Add((TH1F*)htmp->ProjectionY(Form("hMass%d",iPtBin),htmp->FindBin(3.*TMath::Pi()/4.),htmp->FindBin(TMath::Pi())));
	    }else{
	      h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi1",iPtBin),htmp->FindBin(TMath::Pi()/4.),htmp->FindBin(3.*TMath::Pi()/4.));
	    }
	  }else h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi%d",iPtBin,iphi),startX,endX);
	  if(hMass==0)hMass=(TH1F*)h1tmp->Clone();
	  else hMass->Add((TH1F*)h1tmp->Clone());
	}
      }
      hMass->SetTitle(Form("hMass_pt%d_phi%d",iFinalPtBin,iphi));
      hMass->SetName(Form("hMass_pt%d_phi%d",iFinalPtBin,iphi));
      outlist->Add(hMass->Clone());
      //      hMass->DrawClone();
      delete hMass;
      hMass=0x0;
    }
  }
  return outlist;
}
//______________________________________________________________
Bool_t DefinePtBins(AliRDHFCuts *cutobj){
  Int_t nPtBinsCuts=cutobj->GetNPtBins();
  Float_t *ptlimsCuts=cutobj->GetPtBinLimits();
  //  for(Int_t iPt=0; iPt<nPtBinsCuts; iPt++) printf(" %d %f-%f\n",iPt,ptlimsCuts[iPt],ptlimsCuts[iPt+1]);
  for(Int_t iPtCuts=0; iPtCuts<nPtBinsCuts; iPtCuts++){
    for(Int_t iFinalPtBin=0; iFinalPtBin<nptbinsnew; iFinalPtBin++){  
      if(TMath::Abs(ptlimsCuts[iPtCuts]-ptbinsnew[iFinalPtBin])<0.0001){ 
	minPtBin[iFinalPtBin]=iPtCuts;
	if(iFinalPtBin>0) maxPtBin[iFinalPtBin-1]=iPtCuts-1;
      }
    }
    if(TMath::Abs(ptlimsCuts[iPtCuts]-ptbinsnew[nptbinsnew])<0.0001) maxPtBin[nptbinsnew-1]=iPtCuts-1;
  }
  if(TMath::Abs(ptbinsnew[nptbinsnew]-ptlimsCuts[nPtBinsCuts])<0.0001) maxPtBin[nptbinsnew-1]=nPtBinsCuts-1;
  for(Int_t iFinalPtBin=0; iFinalPtBin<nptbinsnew; iFinalPtBin++){
    printf("Pt bins to be merged: %d %d\n",minPtBin[iFinalPtBin],maxPtBin[iFinalPtBin]);
    if(minPtBin[iFinalPtBin]<0 || maxPtBin[iFinalPtBin]<0) return kFALSE;
  }

  return kTRUE;
}
//______________________________________________________________
Int_t GetPadNumber(Int_t ix,Int_t iy){
  return (iy)*nptbinsnew+ix+1;
}
//________________________________________________________________________________
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs,Bool_t inoutanis){

  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;
  
  //Canvases for drawing histograms
  TCanvas *cDeltaPhi = new TCanvas("cinvmassdeltaphi","Invariant mass distributions",1200,700);
  TCanvas *cDeltaPhifs = new TCanvas("cinvmassdeltaphifs","Invariant mass distributions - fit with fixed sigma",1200,700);
  TCanvas *cPhiInteg = new TCanvas("cinvmass","Invariant mass distributions - #phi integrated",1200,350);
  cDeltaPhi->Divide(nptbinsnew,nphi);
  cDeltaPhifs->Divide(nptbinsnew,nphi);
  cPhiInteg->Divide(nptbinsnew,1);

  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){    
    TH1F *histtofitfullsigma=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi0",ipt))->Clone();
    for(Int_t iphi=0;iphi<nphi;iphi++){
      Int_t ipad=GetPadNumber(ipt,iphi);
      Double_t signal=0,esignal=0;
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      if(iphi>0)histtofitfullsigma->Add((TH1F*)histtofit->Clone());
      if(!histtofit){
	gSignal[ipt]->SetPoint(iphi,iphi,signal);
	gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
	return;
      }
      histtofit->SetTitle(Form("%.1f<p_{t}<%.1f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      AliHFMassFitter fitter(histtofit,histtofit->GetBinLowEdge(2),histtofit->GetBinLowEdge(histtofit->GetNbinsX()-2),rebin[ipt]);
      fitter.SetInitialGaussianMean(mass);
      //      fitter.SetInitialGaussianSigma(0.012);
      Bool_t ok=fitter.MassFitter(kFALSE);
      if(ok){
	fitter.DrawHere(cDeltaPhi->cd(ipad),3,1);
	fitter.Signal(3,signal,esignal);
      }
      gSignal[ipt]->SetPoint(iphi,iphi,signal);
      gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
    }
    //fit for fixed sigma
    histtofitfullsigma->SetTitle(Form("%.1f<p_{t}<%.1f",ptbinsnew[ipt],ptbinsnew[ipt+1]));
    AliHFMassFitter fitter(histtofitfullsigma,histtofitfullsigma->GetBinLowEdge(2),histtofitfullsigma->GetBinLowEdge(histtofitfullsigma->GetNbinsX()-2),rebin[ipt]);
    fitter.SetInitialGaussianMean(mass);
    //      fitter.SetInitialGaussianSigma(0.012);
    Bool_t ok=fitter.MassFitter(kFALSE);
    if(ok){
      fitter.DrawHere(cPhiInteg->cd(ipt+1),3,1);
    }
    Double_t sigma=fitter.GetSigma();
    for(Int_t iphi=0;iphi<nphi;iphi++){
      Int_t ipad=GetPadNumber(ipt,iphi);
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      histtofit->SetTitle(Form("%.1f<p_{t}<%.1f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      AliHFMassFitter fitter2(histtofit,histtofit->GetBinLowEdge(2),histtofit->GetBinLowEdge(histtofit->GetNbinsX()-2),rebin[ipt]);
      fitter2.SetInitialGaussianMean(mass);
      fitter2.SetFixGaussianSigma(sigma);
      Bool_t ok2=fitter2.MassFitter(kFALSE);
      Double_t signal=0,esignal=0;
      if(ok2){
	fitter2.DrawHere(cDeltaPhifs->cd(ipad),3,1);
	fitter2.Signal(3,signal,esignal);
      }
      gSignalfs[ipt]->SetPoint(iphi,iphi,signal);
      gSignalfs[ipt]->SetPointError(iphi,0,0,esignal,esignal);

    }
  }//end loop on pt bin

  cDeltaPhi->SaveAs("InvMassDeltaPhi.eps");
  cDeltaPhifs->SaveAs("InvMassDeltaPhi_fs.eps");
  cPhiInteg->SaveAs("InvMassfullphi.eps");

  if(!gnopng){
    cDeltaPhi->SaveAs("InvMassDeltaPhi.png");
    cDeltaPhifs->SaveAs("InvMassDeltaPhi_fs.png");
    cPhiInteg->SaveAs("InvMassfullphi.png");
  }
  //cDeltaPhifs->DrawClone();
  cDeltaPhifs->Close();
 
}
//______________________________________________________________
void DmesonsFlowAnalysis(Bool_t inoutanis=kTRUE,Int_t minC=30,Int_t maxC=50,TString partname="Dplus"){
  TString filename="AnalysisResults.root";
  TString dirname="PWG3_D2H_HFv2";
  TString listname="coutputv2";

  AliRDHFCuts *cutsobj=0x0;
  //Load input data from AliAnalysisTaskSEHFv2
  TFile *f = TFile::Open(filename.Data());
  if(!f){
    printf("file %s not found, please check file name\n",filename.Data());return;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname.Data());
  if(!dir){
    printf("Directory %s not found, please check dir name\n",dirname.Data());return;
  }
  if(partname.Contains("D0")) {
    cutsobj=((AliRDHFCutsD0toKpi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    mass=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  if(partname.Contains("Dplus")){
    cutsobj=((AliRDHFCutsDplustoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    mass=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
  if(partname.Contains("Dstar")) {
    cutsobj=((AliRDHFCutsDStartoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    mass=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass()); 
  }

  TList *list =(TList*)dir->Get(listname.Data());
  if(!list){
    printf("list %s not found in file, please check list name\n",listname.Data());return;
  }
  if(!cutsobj){
    printf("cut object not found in file, please check keylist number\n");return;
  }
  //Define new pt bins
  if(!DefinePtBins(cutsobj)){
    printf("cut not define pt bins\n");return;
  }
  //Load mass histograms corresponding to the required centrality, pt range and phi binning
  printf("Load mass histos \n");
  TList *histlist=LoadMassHistos(list,minC,maxC,inoutanis);
  TString anisstr="";
  if(inoutanis)anisstr+="anis";
  histlist->SaveAs(Form("v2Histograms_%d_%d_%s.root",minC,maxC,anisstr.Data()),"RECREATE");

  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;

  //EP resolution
  TH1F* hevplresos=(TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",minC,minC+5));
  for(Int_t icent=minC+5;icent<maxC;icent++)hevplresos->Add((TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",icent,icent+5)));
  //Double_t resolSub=TMath::Sqrt(hevplresos->GetMean());
  Double_t resol=AliVertexingHFUtils::GetFullEvResol(hevplresos);//ComputeResol(resolSub,1);

  printf("average pt for pt bin \n");
  //average pt for pt bin
  AliVertexingHFUtils *utils=new AliVertexingHFUtils();
  TH2F* hmasspt=(TH2F*)list->FindObject(Form("hMPtCandcentr%d_%d",minC,minC+5));
  for(Int_t icent=minC+5;icent<maxC;icent=icent+5)hmasspt->Add((TH2F*)list->FindObject(Form("hMPtCandCentr%d_%d",icent,icent+5)));
  Float_t averagePt[nptbinsnew];
  Float_t errorPt[nptbinsnew];
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    TH1F *histtofit = (TH1F*)hmasspt->ProjectionY("_py",hmasspt->FindBin(ptbinsnew[ipt]),hmasspt->FindBin(ptbinsnew[ipt+1]));
    AliHFMassFitter fitter(histtofit,histtofit->GetBinLowEdge(2),histtofit->GetBinLowEdge(histtofit->GetNbinsX()-2),1);
    fitter.MassFitter(kFALSE);
    Float_t massFromFit=fitter.GetMean();
    Float_t sigmaFromFit=fitter.GetSigma();
    TF1* funcB2=fitter.GetBackgroundRecalcFunc();
    utils->AveragePt(averagePt[ipt],errorPt[ipt],ptbinsnew[ipt],ptbinsnew[ipt+1],hmasspt,massFromFit,sigmaFromFit,funcB2);
  }
  printf("Fill TGraphs for signal \n");
  //Fill TGraphs for signal
  TGraphAsymmErrors *gSignal[nptbinsnew];
  TGraphAsymmErrors *gSignalfs[nptbinsnew];
  for(Int_t i=0;i<nptbinsnew;i++){
    gSignal[i]=new TGraphAsymmErrors(nphi);
    gSignal[i]->SetName(Form("gasigpt%d",i));
    gSignal[i]->SetTitle(Form("Signal %.1f<p_{t}<%.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignal[i]->SetMarkerStyle(25);
    gSignalfs[i]=new TGraphAsymmErrors(nphi);
    gSignalfs[i]->SetName(Form("gasigfspt%d",i));
    gSignalfs[i]->SetTitle(Form("Signal (fixed sigma) %.1f<p_{t}<%.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignalfs[i]->SetMarkerStyle(21);
  }
  FillSignalGraph(histlist,gSignal,gSignalfs,inoutanis);

  printf("Compute v2 \n");
  //compute v2
  TCanvas *cv2 =new TCanvas("v2","v2");
  TCanvas *cv2fs =new TCanvas("v2_fs","v2 Fit with fixed sigma");
  TGraphAsymmErrors *gv2=new TGraphAsymmErrors(nptbinsnew);
  TGraphAsymmErrors *gv2fs=new TGraphAsymmErrors(nptbinsnew);
  gv2->SetName("gav2");
  gv2->SetTitle("v_{2}(p_{t})");
  gv2fs->SetName("gav2fs");
  gv2fs->SetTitle("v_{2}(p_{t}) (fit with fixed sigma)");

  //Prepare output file
  TFile *fout=new TFile(Form("v2Output_%d_%d_%s.root",minC,maxC,anisstr.Data()),"RECREATE");

  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    fout->cd();
    gSignal[ipt]->Write();
    gSignalfs[ipt]->Write();

    if(inoutanis){
      //v2 from in-out anisotropy
      Double_t *y,*yfs;
      y=gSignal[ipt]->GetY();
      yfs=gSignalfs[ipt]->GetY();
      Double_t nIn=y[0];
      Double_t nOut=y[1];
      Double_t enIn=gSignal[ipt]->GetErrorY(0);
      Double_t enOut=gSignal[ipt]->GetErrorY(1);
      Double_t anis=(nIn-nOut)/(nIn+nOut);
      Double_t eAnis=TMath::Sqrt(enIn*enIn+enOut*enOut)/(nIn+nOut);
      Double_t v2=anis*TMath::Pi()/4./resol;
      Double_t ev2=eAnis*TMath::Pi()/4./resol;
      gv2->SetPoint(ipt,averagePt[ipt],v2);
      gv2->SetPointError(ipt,averagePt[ipt]-ptbinsnew[ipt],ptbinsnew[ipt+1]-averagePt[ipt],ev2,ev2);

      nIn=yfs[0];
      nOut=yfs[1];
      enIn=gSignal[ipt]->GetErrorY(0);
      enOut=gSignal[ipt]->GetErrorY(1);
      anis=(nIn-nOut)/(nIn+nOut);
      eAnis=TMath::Sqrt(enIn*enIn+enOut*enOut)/(nIn+nOut);
      v2=anis*TMath::Pi()/4./resol;
      ev2=eAnis*TMath::Pi()/4./resol;
      gv2fs->SetPoint(ipt,averagePt[ipt],v2);
      gv2fs->SetPointError(ipt,averagePt[ipt]-ptbinsnew[ipt],ptbinsnew[ipt+1]-averagePt[ipt],ev2,ev2);
    }else{
      TF1 *flowFunc = new TF1("flow","[0]*(1.+2.*[1]*TMath::Cos(2.*x))");
      //v2 from fit to Deltaphi distribution
      gSignal[ipt]->Fit(flowFunc);
      Double_t v2 = flowFunc->GetParameter(1)/resol;
      Double_t ev2=flowFunc->GetParError(1)/resol;
      gv2->SetPoint(ipt,averagePt[ipt],v2);
      gv2->SetPointError(ipt,averagePt[ipt]-ptbinsnew[ipt],ptbinsnew[ipt+1]-averagePt[ipt],ev2,ev2);

      gSignalfs[ipt]->Fit(flowFunc);
      v2 = flowFunc->GetParameter(1)/resol;
      ev2=flowFunc->GetParError(1)/resol;
      gv2fs->SetPoint(ipt,averagePt[ipt],v2);
      gv2fs->SetPointError(ipt,averagePt[ipt]-ptbinsnew[ipt],ptbinsnew[ipt+1]-averagePt[ipt],ev2,ev2);
    }
  }
  gv2->SetMarkerStyle(20);
  gv2fs->SetMarkerStyle(20);
  cv2->cd();
  gv2->Draw("AP");
  gv2->GetYaxis()->SetTitle("v_{2}");
  gv2->GetXaxis()->SetTitle("p_{t} [GeV/c]");
  cv2fs->cd();
  gv2fs->Draw("AP");
  gv2fs->GetYaxis()->SetTitle("v_{2}");
  gv2fs->GetXaxis()->SetTitle("p_{t} [GeV/c]");

  fout->cd();
  gv2->Write();
  gv2fs->Write();
  printf("Event plane resolution %f\n",resol);
  printf("Average pt\n");
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++) printf("%f +- %f\n",averagePt[ipt],errorPt[ipt]);
}
//_______________________________________________________________________________________________________________________________
Int_t FindPtBin(Int_t nbins, Float_t* array,Float_t value){
  for (Int_t i=0;i<nbins;i++){
    if(value>=array[i] && value<array[i+1]){
      return i;
    }
  }
  cout<<value<< " out of range "<<array[0]<<", "<<array[nbins]<<endl;
  return -1;
}

//_______________________________________________________________________________________________________________________________
void DrawEventPlane(Int_t mincentr,Int_t maxcentr,TString filename,TString dirname,TString listname){
  TFile *f = TFile::Open(filename.Data());
  if(!f){
    printf("file %s not found, please check file name\n",filename.Data());return;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname.Data());
  if(!dir){
    printf("Directory %s not found, please check dir name\n",dirname.Data());return;
  }
  TList *list =(TList*)dir->Get(listname.Data());
  if(!list){
    printf("list %s not found in file, please check list name\n",listname.Data());return;
  }
  DrawEventPlane(list,mincentr,maxcentr);
}
//_______________________________________________________________________________________________________________________________
void DrawEventPlane(TList *list,Int_t mincentr,Int_t maxcentr){

 //draw the histograms correlated to the event plane

  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  //gStyle->SetPalette(1);
  gStyle->SetOptFit(1);
 
  if(!(mincentr==0 && maxcentr==0)){ //draw the total in mincentr-maxcentr
    TString suffixcentr=Form("centr%d_%d",mincentr,maxcentr);
    TH2F* hevpls=(TH2F*)list->FindObject(Form("hEvPlanecentr%d_%d",mincentr,mincentr+5));
    hevpls->SetName(Form("hEvPlane%s",suffixcentr.Data()));
    hevpls->SetTitle(Form("Event Plane angle %s",suffixcentr.Data()));
    TH1F* hevplresos=(TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",mincentr,mincentr+5));
    hevplresos->SetName(Form("hEvPlaneReso%s",suffixcentr.Data()));
    hevplresos->SetTitle(Form("Event Plane Resolution %s",suffixcentr.Data()));

    for(Int_t icentr=mincentr+5;icentr<maxcentr;icentr=icentr+5){
      TH2F* h=(TH2F*)list->FindObject(Form("hEvPlanecentr%d_%d",icentr,icentr+5));
      if(h)hevpls->Add(h);
      else cout<<"skipping ev plane "<<icentr<<"_"<<icentr+5<<endl;
      TH1F *hr=(TH1F*)list->FindObject(Form("hEvPlaneResocentr%d_%d",icentr,icentr+5));
      if(hr)hevplresos->Add(hr);
      else cout<<"skipping ev pl reso "<<icentr<<"_"<<icentr+5<<endl;
    }

    TCanvas* cvtotevpl=new TCanvas("cvtotevpl",Form("Ev plane %s",suffixcentr.Data()),1200,400);
    cvtotevpl->Divide(3,1);
    cvtotevpl->cd(1);
    hevpls->Draw("COLZ");
    TH1F* htpc = (TH1F*)hevpls->ProjectionX();
    cvtotevpl->cd(2);
    htpc->Draw();
    htpc->Fit("pol0");
    TH1F* hv0 = (TH1F*)hevpls->ProjectionY();
    cvtotevpl->cd(3);
    hv0->Draw();
    hv0->Fit("pol0");

    TCanvas* cvtotevplreso=new TCanvas("cvtotevplreso",Form("Ev plane Resolution %s",suffixcentr.Data()));
    cvtotevplreso->cd();
    hevplresos->Draw();
    //AliVertexingHFUtils::SetSubEventHisto(hevplresos);
    //Double_t resolSub=TMath::Sqrt(hevplresos->GetMean());
    Double_t resolFull=AliVertexingHFUtils::GetFullEvResol(hevplresos);//ComputeResol(resolSub,1);

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

    TGraph* gresovscentr=new TGraph(0);
    gresovscentr->SetName("gresovscentr");
    gresovscentr->SetTitle(Form("Resolution vs Centrality;centrality (%s);Resolution","%"));
    Int_t ic=0;
    for(Int_t i=0;i<list->GetEntries();i++){
       TClass* objtype=list->At(i)->IsA();
       TString tpname=objtype->GetName();
       if(tpname=="TH2F"){
	 TH2F* h=(TH2F*)list->At(i);
	 // if(!h){
	 //   cout<<"Histogram "<<i<<" not found"<<endl;
	 //   continue;
	 // }
	 TString hname=h->GetName();
	 if(hname.Contains("EvPlane")){
	   TString fixedname="hEvPlanecentr";
	   TString scentr=hname.Copy();
	   //printf("%s\n",scentr.Data());
	   scentr.Remove(0,fixedname.Length());
	   printf("%s\n",scentr.Data());
	   TCanvas* cv=new TCanvas(Form("cv%s",hname.Data()),hname.Data());
	   cv->cd();
	   TH1F* htpc = (TH1F*)h->ProjectionX();
	   htpc->SetTitle(Form("TPC event plane (%s)",scentr.Data()));
	   htpc->Draw();
	   htpc->Fit("pol0");
	 }
       }else{ //it is a TH1F
	 TH1F* h=(TH1F*)list->At(i);
	 TString hname=h->GetName();
	 if(hname.Contains("Reso")){
	 //Double_t resolSub=TMath::Sqrt(h->GetMean());
	   Double_t resolFull=AliVertexingHFUtils::GetFullEvResol(h);//ComputeResol(resolSub,1);
	   TString smaxc=hname.Copy();
	   smaxc.Remove(0,hname.Length()-2);
	   Int_t maxc=smaxc.Atoi();
	   gresovscentr->SetPoint(ic,maxc-2.5,resolFull);
	   ic++;
	   TPaveText* pvreso=new TPaveText(0.1,0.1,0.6,0.2,"NDC");
	   pvreso->SetBorderSize(0);
	   pvreso->SetFillStyle(0);
	   pvreso->AddText(Form("Resolution on full event = %.4f\n",resolFull));
	   pvreso->Draw();
	}
      }
    }//loop on centrality
    TCanvas* cresovscentr=new TCanvas("cresovscentr", "Resolution vs Centrality");
    cresovscentr->cd();
    gresovscentr->SetMarkerStyle(20);
    gresovscentr->Draw("AP");
    if(!gnopng) gresovscentr->SaveAs(Form("%s.png",gresovscentr->GetName()));
    gresovscentr->SaveAs(Form("%s.eps",gresovscentr->GetName()));
    TFile* fout=new TFile(Form("%s.root",gresovscentr->GetName()),"recreate");
    fout->cd();
    gresovscentr->Write();
  }

}
//D^{0}#rightarrow K^{-}#pi^{+}
void DrawPaveStat(TVirtualPad* c,Int_t nev,Int_t mincentr,Int_t maxcentr,TString string,TString string2){

  TPaveText* txtoff=new TPaveText(0.2,0.3,0.6,0.5,"NDC");
  txtoff->SetBorderSize(0);
  txtoff->SetFillStyle(0);
  //txtoff->AddText(Form("%.0f#times10^{6} events",nev/1e6));
  txtoff->AddText("Pb-Pb, #sqrt{s_{NN}} = 2.76TeV");
  printf("Centrality range %d%s, nev=%d, corrected=%f \n",maxcentr-mincentr,"%",nev,(Float_t)nev/(Float_t)(maxcentr-mincentr)/100.);
  txtoff->AddText(Form("%.0f#times10^{6} evts in %d-%d%s centr cl",nev/1e6*(Float_t)(maxcentr-mincentr)/100.,mincentr,maxcentr,"%"));
  txtoff->AddText(string.Data());
  if(string2!="") {
    txtoff->AddText(string2.Data());
    txtoff->SetY1(0.3);
  }
  txtoff->SetTextColor(kAzure-5);
  c->cd();
  txtoff->Draw();
}

void ALICEPerformance2(TVirtualPad *c,TString today,Double_t* where,TString type){
 
 //date must be in the form: 04/05/2010
  if(today=="today"){
    TDatime startt;                                                        
    int date=startt.GetDate();
    int y=date/10000;
    int m=(date%10000)/100;
    int d=date%100;


    today="";
    today+=d;
    if(m<10)
      today.Append("/0");
    else today.Append("/");
    today+=m;
    today.Append("/");
    today+=y;  
    
  }
  cout<<"here"<<endl;
  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",where[0],where[1],where[2],where[3]);
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
  TASImage *myAliceLogo = new TASImage("/home/cbianchi/macros/alice_logo.png");
  myAliceLogo->Draw();
  printf("Draw logo\n");
  c->cd();
  TPaveText* t1=new TPaveText(where[0]-0.03,where[1]-0.11,where[2]+0.01,where[1]-0.01,"NDC");

  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  //t1->AddText(0.,0.,Form("%s%s",type.Data(),(type=="Performace") ? today.Data() : ""));
t1->AddText(0.,0.,Form("%s",type.Data()));
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  
  if(type=="Preliminary") {
    cout<<"Preliminary"<<endl;
    //t1->Print();
    t1->Draw();
  }
  else {
    printf("performance\n");
    t1->AddText(0.,0.,today.Data());
    t1->SetY1NDC(where[1]-0.17);
    t1->SetY2NDC(where[1]-0.05);
    t1->Draw();
    printf("Drawn date\n");
  }
  
}
