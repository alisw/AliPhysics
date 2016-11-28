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

#include "AliHFMassFitter.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliEventPlaneResolutionHandler.h"
#include "AliVertexingHFUtils.h"

#endif

/* $Id$ */

//methods for the analysis of AliAnalysisTaskSEHFv2 output
//Authors: Chiara Bianchin cbianchi@pd.infn.it
//         Giacomo Ortona  ortona@to.infn.it
//         Francesco Prino prino@to.infn.it

//_________________________________________________________________
//GLOBAL VARIABLES TO BE SET
//input file
const TString filename="AnalysisResults.root";
//const TString filename="$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/HR_bunch1-3/AnalysisResultsv2_EP_pass0_diffdet_topod0cut_3050.root";
const TString suffix="_Topod0Cut_pass0_QoverM_VZERO_EP";
const TString partname="Dplus";
const Int_t minCent=30;
const Int_t maxCent=50;

const TString outputdir = ".";

//EP resolution
//kTPCFullEta, kTPCPosEta,kVZERO,kVZEROA,kVZEROC
const Int_t evPlane=AliEventPlaneResolutionHandler::kVZERO;
//resolution flag fromAliEventPlaneResolutionHandler:
//kTwoRandSub,kTwoChargeSub,kTwoEtaSub,kThreeSub,kThreeSubTPCGap
const Bool_t useAliHandlerForRes=kFALSE;
Int_t evPlaneRes=AliEventPlaneResolutionHandler::kThreeSub;
const Bool_t useNcollWeight=kFALSE;

// pt and phi binning
const Int_t nptbinsnew=10;
const Double_t ptbinsnew[nptbinsnew+1]={2,3,4,5,6.,7.,8,10,12,16,24};

const Int_t nphibins=4;
const Double_t phibinslim[nphibins+1]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};

// mass fit configuration
const Int_t rebin[]={2,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4};
const Int_t typeb=AliHFMassFitter::kExpo;  // Background: 0=expo, 1=linear, 2=pol2
const Bool_t fixAlsoMass=kFALSE;
const Double_t minMassForFit=1.69;
const Double_t maxMassForFit=2.02;
const Double_t nSigmaForCounting=3.5;

//not to be set
Int_t minPtBin[nptbinsnew]={-1};
Int_t maxPtBin[nptbinsnew]={-1};
const Double_t effInOverEffOut=1.03;
Double_t massD;

const Int_t colors[] = {kRed+1,kBlack,kBlue+1,kGreen+2,kOrange+7,kBlue-7};
const Int_t markers[] = {kFullSquare,kFullCircle,kFullTriangleUp,kFullDiamond,kOpenSquare,kOpenCircle,kOpenTriangleUp,kOpenDiamond};

//_________________________________________________________________
//METHODS PROTOTYPES
void DmesonsFlowAnalysis(Bool_t inoutanis=kTRUE);
TList *LoadMassHistos(TList *inputlist,Bool_t inoutanis);
TList* LoadResolutionHistos(TList *inputlist);
Int_t FindPtBin(Int_t nbins, Double_t* array,Double_t value);
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2,TH1F **hSigmaFree,TH1F **hSigmaFixed, TH1F **hMean,TH1F **hMeanfs, TH1F **hChiSquare, TH1F **hChiSquarefs, Bool_t inoutanis, Int_t bkgfunc, Int_t minfit, Int_t maxfit, const Int_t *rebin);
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Bool_t inoutanis, TGraphAsymmErrors *gRelSystEff);
void DrawEventPlane();
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3);
void SetStyle(Int_t optfit=0);

//_________________________________________________________________
//METHODS IMPLEMENTATION
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3){
  Double_t resolFull=1.;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kTwoRandSub ||
     evPlaneRes==AliEventPlaneResolutionHandler::kTwoChargeSub){
    resolFull=AliVertexingHFUtils::GetFullEvResol(hsubev1);
    error = TMath::Abs(resolFull-AliVertexingHFUtils::GetFullEvResolLowLim(hsubev1));
  }else if(evPlaneRes==AliEventPlaneResolutionHandler::kTwoEtaSub){
    if(evPlane==AliEventPlaneResolutionHandler::kTPCFullEta){
      resolFull=AliVertexingHFUtils::GetFullEvResol(hsubev1);
      error = TMath::Abs(resolFull-AliVertexingHFUtils::GetFullEvResolLowLim(hsubev1));
    }else if(evPlane==AliEventPlaneResolutionHandler::kTPCPosEta){
      resolFull=AliVertexingHFUtils::GetSubEvResol(hsubev1);
      error = TMath::Abs(resolFull-AliVertexingHFUtils::GetSubEvResolLowLim(hsubev1));      
    }
  }else{
    Double_t resolSub[3];
    Double_t errors[3];
    TH1F* hevplresos[3];
    hevplresos[0]=hsubev1;
    hevplresos[1]=hsubev2;
    hevplresos[2]=hsubev3;
    TString namereso1 = hevplresos[0]->GetName();
    for(Int_t ires=0;ires<3;ires++){
      resolSub[ires]=hevplresos[ires]->GetMean();
      errors[ires]=hevplresos[ires]->GetMeanError();
    }
    Double_t lowlim[3];
    for(Int_t ie=0;ie<3;ie++) {lowlim[ie]=TMath::Abs(resolSub[ie]-errors[ie]);}
    
    if(!namereso1.Contains("hEvPlaneReso1")) {
      if(evPlane==AliEventPlaneResolutionHandler::kVZEROC ||
         evPlane==AliEventPlaneResolutionHandler::kVZERO){
        resolFull=TMath::Sqrt(resolSub[1]*resolSub[2]/resolSub[0]);
        error=resolFull-TMath::Sqrt(lowlim[2]*lowlim[1]/lowlim[0]);
      }
      else if(evPlane==AliEventPlaneResolutionHandler::kVZEROA){
        resolFull=TMath::Sqrt(resolSub[0]*resolSub[2]/resolSub[1]);
        error=resolFull-TMath::Sqrt(lowlim[2]*lowlim[0]/lowlim[1]);
      }
      else if(evPlane==AliEventPlaneResolutionHandler::kTPCFullEta ||
              evPlane==AliEventPlaneResolutionHandler::kTPCPosEta){
        resolFull=TMath::Sqrt(resolSub[0]*resolSub[1]/resolSub[2]);
        error=resolFull-TMath::Sqrt(lowlim[0]*lowlim[1]/lowlim[2]);
      }
    }
    else {
      resolFull=TMath::Sqrt(resolSub[1]*resolSub[2]/resolSub[0]);
      error=resolFull-TMath::Sqrt(lowlim[2]*lowlim[1]/lowlim[0]);
    }
  }
  return resolFull;
}

//____________________________________________________________________
TList* LoadResolutionHistos(TList *inputlist){

  TList *outlist = new TList();
  outlist->SetName("eventplanehistlist");

  const Int_t nBins=20;
  Double_t ncoll[nBins]={1790.77,1578.44,1394.82,1236.17
    ,1095.08,969.86,859.571,759.959,669.648,589.588,516.039
    ,451.409,392.853,340.493,294.426,252.385,215.484,183.284
    ,155.101,130.963};

  Int_t minCentTimesTen=minCent*10;
  Int_t maxCentTimesTen=maxCent*10;
  TGraphErrors* gResolVsCent=new TGraphErrors(0);
  Int_t iPt=0;
  Int_t nSubRes=1;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
     evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
  TString namereso[3]={"Reso","Reso2","Reso3"};
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
  TH1F* htestversion=(TH1F*)inputlist->FindObject(Form("hEvPlaneResocentr%d_%d",minCentTimesTen,minCentTimesTen+25));
  if(htestversion){
    printf("Old version of the task\n");
  }else{
    printf("New version of the task\n");
    namereso[0]="Reso1";
  }
  TH2F* hevpls=(TH2F*)inputlist->FindObject(Form("hEvPlanecentr%d_%d",minCentTimesTen,minCentTimesTen+25));
  hevpls->SetName(Form("hEvPlane%s",suffixcentr.Data()));
  hevpls->SetTitle(Form("Event Plane angle %s",suffixcentr.Data()));
  TH1F* hevplresos[3];
  Int_t ncBin=minCentTimesTen/25;
  
  for(Int_t ires=0;ires<nSubRes;ires++){
    hevplresos[ires]=(TH1F*)inputlist->FindObject(Form("hEvPlane%scentr%d_%d",namereso[ires].Data(),minCentTimesTen,minCentTimesTen+25));
    if(hevplresos[ires]){
      hevplresos[ires]->SetName(Form("hEvPlane%s%s",namereso[ires].Data(),suffixcentr.Data()));
      hevplresos[ires]->SetTitle(Form("Event Plane Resolution %s%s",namereso[ires].Data(),suffixcentr.Data()));
      if(useNcollWeight){
	printf("Centr %d Bin %d  Ncoll %f\n",minCentTimesTen,ncBin,ncoll[ncBin]);
	hevplresos[ires]->Scale(ncoll[ncBin]);
      }
    }
  }
  Double_t error;
  Double_t lowestRes=1;
  Double_t highestRes=0;
  Double_t resolBin=GetEventPlaneResolution(error,hevplresos[0],hevplresos[1],hevplresos[2]);
  if(resolBin<lowestRes) lowestRes=resolBin;
  if(resolBin>highestRes) highestRes=resolBin;

  Double_t binHalfWid=25./20.;
  Double_t binCentr=(Double_t)minCentTimesTen/10.+binHalfWid;
  gResolVsCent->SetPoint(iPt,binCentr,resolBin);
  gResolVsCent->SetPointError(iPt,binHalfWid,error);
  ++iPt;
  
  for(Int_t icentr=minCentTimesTen+25;icentr<maxCentTimesTen;icentr=icentr+25){
    TH2F* h=(TH2F*)inputlist->FindObject(Form("hEvPlanecentr%d_%d",icentr,icentr+25));
    if(h)hevpls->Add(h);
    else cout<<"skipping ev plane "<<icentr<<"_"<<icentr+5<<endl;
    TH1F* htmpresos[3];
    for(Int_t ires=0;ires<nSubRes;ires++){
      htmpresos[ires]=(TH1F*)inputlist->FindObject(Form("hEvPlane%scentr%d_%d",namereso[ires].Data(),icentr,icentr+25));
      if(!htmpresos[ires])cout<<"skipping ev pl reso "<<icentr<<"_"<<icentr+25<<endl;
    }
    resolBin=GetEventPlaneResolution(error,htmpresos[0],htmpresos[1],htmpresos[2]);
    if(resolBin<lowestRes) lowestRes=resolBin;
    if(resolBin>highestRes) highestRes=resolBin;
    binCentr=(Double_t)icentr/10.+binHalfWid;
    gResolVsCent->SetPoint(iPt,binCentr,resolBin);
    gResolVsCent->SetPointError(iPt,binHalfWid,error);
    ++iPt;
    ncBin=icentr/25;
    for(Int_t ires=0;ires<nSubRes;ires++){
      if(htmpresos[ires]){
	if(useNcollWeight){
	  printf("Centr %d Bin %d  Ncoll %f\n",icentr,ncBin,ncoll[ncBin]);
	  htmpresos[ires]->Scale(ncoll[ncBin]);
	}
	hevplresos[ires]->Add(htmpresos[ires]);
      }
    }
  }
  outlist->Add(hevpls->Clone());
  for(Int_t ires=0;ires<nSubRes;ires++){
    if(hevplresos[ires]) outlist->Add(hevplresos[ires]->Clone());
  }
  gResolVsCent->SetName("gResolVsCent");
  gResolVsCent->SetTitle("Resolution vs. Centrality");
  outlist->Add(gResolVsCent->Clone());
  return outlist;
}

//__________________________________________________________
TList *LoadMassHistos(TList *inputlist,Bool_t inoutanis){
  // printf("Start load histos\n");
  //  const Int_t nptbins=cutobj->GetNPtBins();
  TList *outlist = new TList();
  outlist->SetName("azimuthalhistoslist");
  
  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;
  Int_t minCentTimesTen=minCent*10;
  Int_t maxCentTimesTen=maxCent*10;
  
  //Create 2D histogram in final pt bins
  for(Int_t iFinalPtBin=0; iFinalPtBin<nptbinsnew; iFinalPtBin++){
    for(Int_t iphi=0;iphi<nphi;iphi++){
      TH1F *hMass=0x0;//=new TH1F();
      for(Int_t iPtBin=minPtBin[iFinalPtBin]; iPtBin<=maxPtBin[iFinalPtBin];iPtBin++){
	for(Int_t iHisC=minCentTimesTen; iHisC<=maxCentTimesTen-25; iHisC+=25){    
	  TString hisname=Form("hMdeltaphi_pt%dcentr%d_%d",iPtBin,iHisC,iHisC+25);
	  TH2F* htmp=(TH2F*)inputlist->FindObject(hisname.Data());
	  printf("---> Histogram: %s\n",htmp->GetName());
	  Int_t startX=htmp->FindBin(phibinslim[iphi]);
	  Int_t endX=htmp->FindBin(phibinslim[iphi+1]-0.0001); // -0.0001 to be sure that the upper limit of the bin is properly set
	  TH1F *h1tmp;
	  if(inoutanis){
	    if(iphi==0){
	      Int_t firstBin=htmp->FindBin(0);
	      Int_t lastBin=htmp->FindBin(TMath::Pi()/4.-0.0001); // -0.0001 to be sure that the upper limit of the bin is pi/4
	      h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi0",iPtBin),firstBin,lastBin);
	      printf("In-plane, Range: bins %d-%d -> phi %f - %f\n",firstBin,lastBin,htmp->GetXaxis()->GetBinLowEdge(firstBin),htmp->GetXaxis()->GetBinUpEdge(lastBin));
	      firstBin=htmp->FindBin(3.*TMath::Pi()/4.);
	      lastBin=htmp->FindBin(TMath::Pi()-0.0001); // -0.0001 to be sure that the upper limit of the bin is pi
	      h1tmp->Add((TH1F*)htmp->ProjectionY(Form("hMass%d",iPtBin),firstBin,lastBin));
	      printf("In-plane, Range: bins %d-%d -> phi %f - %f\n",firstBin,lastBin,htmp->GetXaxis()->GetBinLowEdge(firstBin),htmp->GetXaxis()->GetBinUpEdge(lastBin));
	    }else{
	      Int_t firstBin=htmp->FindBin(TMath::Pi()/4.);
	      Int_t lastBin=htmp->FindBin(3.*TMath::Pi()/4.-0.0001); // -0.0001 to be sure that the upper limit of the bin is pi/4
	      h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi1",iPtBin),firstBin,lastBin);
	      printf("Out-of-plane, Range: bins %d-%d -> phi %f - %f\n",firstBin,lastBin,htmp->GetXaxis()->GetBinLowEdge(firstBin),htmp->GetXaxis()->GetBinUpEdge(lastBin));
	    }
	  }else{
	    h1tmp=(TH1F*)htmp->ProjectionY(Form("hMass%d_phi%d",iPtBin,iphi),startX,endX);
	  }
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
  Float_t *ptlimsCuts=(Float_t*)cutobj->GetPtBinLimits();
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
//______________________________________________________________
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2,TH1F **hSigmaFree,TH1F **hSigmaFixed, TH1F **hMean,TH1F **hMeanfs, TH1F **hChiSquare, TH1F **hChiSquarefs, Bool_t inoutanis, Int_t bkgfunc, Int_t minfit, Int_t maxfit, const Int_t *rebin) {
  
  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;
  
  //Canvases for drawing histograms
  TCanvas *cDeltaPhi = new TCanvas("cinvmassdeltaphi","Invariant mass distributions",1920,1080);
  TCanvas *cDeltaPhifs = new TCanvas("cinvmassdeltaphifs","Invariant mass distributions - fit with fixed sigma",1920,1080);
  TCanvas *cPhiInteg = new TCanvas("cinvmass","Invariant mass distributions - #phi integrated",1920,1080);
  cDeltaPhi->Divide(nptbinsnew,nphi);
  cDeltaPhifs->Divide(nptbinsnew,nphi);
  if(nptbinsnew%2==0) {cPhiInteg->Divide(nptbinsnew/2,2);}
  else {cPhiInteg->Divide(nptbinsnew/2+1,2);}
  Int_t nMassBins;
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    TH1F *histtofitfullsigma=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi0",ipt))->Clone();
    for(Int_t iphi=0;iphi<nphi;iphi++){
      Int_t ipad=GetPadNumber(ipt,iphi);
      Double_t signal=0,esignal=0;
      Double_t sigma=0, esigma=0;
      Double_t mean=0, emean=0;
      Double_t chisquare=0;
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      if(iphi>0)histtofitfullsigma->Add((TH1F*)histtofit->Clone());
      if(!histtofit){
        gSignal[ipt]->SetPoint(iphi,iphi,signal);
        gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
        return;
      }
      histtofit->SetTitle(Form("%.0f < #it{p}_{T} < %.0f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      nMassBins=histtofit->GetNbinsX();
      histtofit->Rebin(rebin[ipt]);
      AliHFMassFitter fitter(histtofit,minfit,maxfit,1,bkgfunc);
      fitter.SetInitialGaussianMean(massD);
      fitter.SetInitialGaussianSigma(0.012);
      Bool_t ok=fitter.MassFitter(kFALSE);
      Double_t sigmaforcounting=0;
      Double_t meanforcounting=0;
      if(ok){
        fitter.DrawHere(cDeltaPhi->cd(ipad),3,1);
        signal=fitter.GetRawYield();
        esignal=fitter.GetRawYieldError();
        sigma=fitter.GetSigma();
        esigma=fitter.GetSigmaUncertainty();
        sigmaforcounting=sigma;
        mean=fitter.GetMean();
        emean=fitter.GetMeanUncertainty();
        meanforcounting=mean;
        chisquare=fitter.GetReducedChiSquare();
      }
      gSignal[ipt]->SetPoint(iphi,iphi,signal);
      gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
      hSigmaFree[iphi]->SetBinContent(ipt+1,sigma);
      hSigmaFree[iphi]->SetBinError(ipt+1,esigma);
      hMean[iphi]->SetBinContent(ipt+1,mean);
      hMean[iphi]->SetBinError(ipt+1,emean);
      hChiSquare[iphi]->SetBinContent(ipt+1,chisquare);
      hChiSquare[iphi]->SetBinError(ipt+1,1.e-15);
      TF1* fB1=fitter.GetBackgroundFullRangeFunc();
      TF1* fB2=fitter.GetBackgroundRecalcFunc();
      Double_t minBinSum=histtofit->FindBin(meanforcounting-nSigmaForCounting*sigmaforcounting);
      Double_t maxBinSum=histtofit->FindBin(meanforcounting+nSigmaForCounting*sigmaforcounting);
      Double_t cntSig1=0.;
      Double_t cntSig2=0.;
      Double_t cntErr=0.;
      for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
        Double_t bkg1=fB1 ? fB1->Eval(histtofit->GetBinCenter(iMB)) : 0;
        Double_t bkg2=fB2 ? fB2->Eval(histtofit->GetBinCenter(iMB)) : 0;
        cntSig1+=(histtofit->GetBinContent(iMB)-bkg1);
        cntSig2+=(histtofit->GetBinContent(iMB)-bkg2);
        cntErr+=(histtofit->GetBinContent(iMB));
      }
      cntErr=TMath::Sqrt(cntErr);
      gSignalBC1[ipt]->SetPoint(iphi,iphi,cntSig1);
      gSignalBC1[ipt]->SetPointError(iphi,0,0,cntErr,cntErr);
      gSignalBC2[ipt]->SetPoint(iphi,iphi,cntSig2);
      gSignalBC2[ipt]->SetPointError(iphi,0,0,cntErr,cntErr);
    }
    //fit for fixed sigma
    histtofitfullsigma->SetTitle(Form("%.0f < #it{p}_{T} < %.0f GeV/c",ptbinsnew[ipt],ptbinsnew[ipt+1]));
    histtofitfullsigma->GetXaxis()->SetTitle("M_{K#pi#pi} (GeV/c^{2})");
    histtofitfullsigma->GetXaxis()->SetTitleSize(0.05);
    nMassBins=histtofitfullsigma->GetNbinsX();
    histtofitfullsigma->Rebin(rebin[ipt]);
    AliHFMassFitter fitter(histtofitfullsigma,minfit,maxfit,1,bkgfunc);
    fitter.SetInitialGaussianMean(massD);
    Bool_t ok=fitter.MassFitter(kFALSE);
    Double_t sigmatot=0;
    Double_t massFromFit=0;
    if(ok){
      fitter.DrawHere(cPhiInteg->cd(ipt+1),3,1);
      sigmatot=fitter.GetSigma();
      massFromFit=fitter.GetMean();
    }
    for(Int_t iphi=0;iphi<nphi;iphi++){
      Int_t ipad=GetPadNumber(ipt,iphi);
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      histtofit->SetTitle(Form("%.1f<#it{p}_{T}<%.1f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      nMassBins=histtofit->GetNbinsX();
      histtofit->Rebin(rebin[ipt]);
      AliHFMassFitter fitter2(histtofit,minfit,maxfit,1,bkgfunc);
      fitter2.SetInitialGaussianMean(massD);
      fitter2.SetFixGaussianSigma(sigmatot);
      if(fixAlsoMass) fitter2.SetFixGaussianMean(massFromFit);
      Bool_t ok2=fitter2.MassFitter(kFALSE);
      Double_t signal=0,esignal=0;
      Double_t sigma=0, esigma=0;
      Double_t mean=0, emean=0;
      Double_t chisquare=0;
      if(ok2){
        fitter2.DrawHere(cDeltaPhifs->cd(ipad),3,1);
        signal=fitter2.GetRawYield();
        esignal=fitter2.GetRawYieldError();
        sigma=fitter.GetSigma();
        esigma=fitter.GetSigmaUncertainty();
        mean=fitter2.GetMean();
        emean=fitter2.GetMeanUncertainty();
        chisquare=fitter2.GetReducedChiSquare();
      }
      gSignalfs[ipt]->SetPoint(iphi,iphi,signal);
      gSignalfs[ipt]->SetPointError(iphi,0,0,esignal,esignal);
      hSigmaFixed[iphi]->SetBinContent(ipt+1,sigma);
      hSigmaFixed[iphi]->SetBinError(ipt+1,esigma);
      hMeanfs[iphi]->SetBinContent(ipt+1,mean);
      hMeanfs[iphi]->SetBinError(ipt+1,emean);
      hChiSquarefs[iphi]->SetBinContent(ipt+1,chisquare);
      hChiSquarefs[iphi]->SetBinError(ipt+1,1.e-15);
    }
  }//end loop on pt bin
  
  cDeltaPhi->SaveAs(Form("%s/InvMassDeltaPhi%s.pdf",outputdir.Data(),suffix.Data()));
  cDeltaPhi->SaveAs(Form("%s/InvMassDeltaPhi%s.root",outputdir.Data(),suffix.Data()));
  cDeltaPhifs->SaveAs(Form("%s/InvMassDeltaPhi_fs%s.pdf",outputdir.Data(),suffix.Data()));
  cDeltaPhifs->SaveAs(Form("%s/InvMassDeltaPhi_fs%s.root",outputdir.Data(),suffix.Data()));
  cPhiInteg->SaveAs(Form("%s/InvMassfullphi%s.pdf",outputdir.Data(),suffix.Data()));
}

//______________________________________________________________
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Bool_t inoutanis, TGraphAsymmErrors *gRelSystEff) {
  
  TGraphAsymmErrors* gv2 = new TGraphAsymmErrors(nptbinsnew);
  
  if(inoutanis) {
    for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
      Double_t *y=gSignal[iPt]->GetY();
      Double_t nIn=y[0];
      Double_t nOut=y[1];
      Double_t enIn=gSignal[iPt]->GetErrorY(0);
      Double_t enOut=gSignal[iPt]->GetErrorY(1);
      Double_t anis=(nIn-nOut)/(nIn+nOut);
      Double_t eAnis=2./((nIn+nOut)*(nIn+nOut))*TMath::Sqrt(nIn*nIn*enOut*enOut+nOut*nOut*enIn*enIn);
      Double_t v2=anis*TMath::Pi()/4./resol;
      Double_t ev2=eAnis*TMath::Pi()/4./resol;
      gv2->SetPoint(iPt,averagePt[iPt],v2);
      gv2->SetPointError(iPt,averagePt[iPt]-ptbinsnew[iPt],ptbinsnew[iPt+1]-averagePt[iPt],ev2,ev2);
      if(gRelSystEff) {
        //systematic uncertainty for in-out efficiency
        Double_t anis1=(nIn-nOut*effInOverEffOut)/(nIn+nOut*effInOverEffOut);
        Double_t anis2=(nIn*effInOverEffOut-nOut)/(nIn*effInOverEffOut+nOut);
        Double_t systEffUp=0.,systEffDown=0.;
        if(anis1>anis && anis1>anis2) systEffUp=anis1/anis;
        if(anis2>anis && anis2>anis1) systEffUp=anis2/anis;
        if(anis1<anis && anis1<anis2) systEffDown=anis1/anis;
        if(anis2<anis && anis2<anis1) systEffDown=anis2/anis;
        cout << Form(" Bin %d <pt>=%.3f  v2=%f+-%f systEff=%f %f\n",iPt,averagePt[iPt],v2,ev2,systEffUp*v2,systEffDown*v2)<<endl;
        gRelSystEff->SetPoint(iPt,averagePt[iPt],v2);
        gRelSystEff->SetPointError(iPt,0.4,0.4,v2*(1-systEffDown),v2*(systEffUp-1));
      }
    }
    return gv2;
  }
  else {
    TF1 *flowFunc = new TF1("flow","[0]*(1.+2.*[1]*TMath::Cos(2.*x))");
    for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
      //v2 from fit to Deltaphi distribution
      gSignal[iPt]->Fit(flowFunc);
      Double_t v2 = flowFunc->GetParameter(1)/resol;
      Double_t ev2=flowFunc->GetParError(1)/resol;
      gv2->SetPoint(iPt,averagePt[iPt],v2);
      gv2->SetPointError(iPt,averagePt[iPt]-ptbinsnew[iPt],ptbinsnew[iPt+1]-averagePt[iPt],ev2,ev2);
    }
    return gv2;
  }
}

//______________________________________________________________
void DmesonsFlowAnalysis(Bool_t inoutanis){
  
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s",partname.Data(),suffix.Data());

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
  if(partname.Contains("Dzero")) {
    cutsobj=((AliRDHFCutsD0toKpi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  if(partname.Contains("Dplus")){
    cutsobj=((AliRDHFCutsDplustoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
  if(partname.Contains("Dstar")) {
    cutsobj=((AliRDHFCutsDStartoKpipi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass()); 
  }
  if(partname.Contains("Ds")) {
    cutsobj=((AliRDHFCutsDstoKKpi*)dir->Get(dir->GetListOfKeys()->At(2)->GetName()));
    massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
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
  TList *histlist=LoadMassHistos(list,inoutanis);
  TString aniss="";
  if(inoutanis)aniss+="anis";
  TString fine="";
  if(nptbinsnew>=15) fine="_fineptbin";
  histlist->SaveAs(Form("%s/v2Histograms_%d_%d_%s%s%s.root",outputdir.Data(),minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()),"RECREATE");

  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;

  printf("average pt for pt bin \n");
  //average pt for pt bin
  AliVertexingHFUtils *utils=new AliVertexingHFUtils();
  Int_t minCentTimesTen=minCent*10;
  Int_t maxCentTimesTen=maxCent*10;
  TH2F* hmasspt=(TH2F*)list->FindObject(Form("hMPtCandcentr%d_%d",minCentTimesTen,minCentTimesTen+25));
  for(Int_t icent=minCentTimesTen+25;icent<maxCentTimesTen;icent=icent+25)hmasspt->Add((TH2F*)list->FindObject(Form("hMPtCandcentr%d_%d",icent,icent+25)));
  Float_t averagePt[nptbinsnew];
  Float_t errorPt[nptbinsnew];
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    Int_t binMin=hmasspt->FindBin(ptbinsnew[ipt]);
    Int_t binMax=hmasspt->FindBin(ptbinsnew[ipt+1]-0.001);
    if(TMath::Abs(hmasspt->GetXaxis()->GetBinLowEdge(binMin)-ptbinsnew[ipt])>0.001 || 
       TMath::Abs(hmasspt->GetXaxis()->GetBinUpEdge(binMax)-ptbinsnew[ipt+1])>0.001){
      printf("Error in pt bin limits for projection!\n");
      return;
    }
    TH1F *histtofit = (TH1F*)hmasspt->ProjectionY("_py",binMin,binMax);
    Int_t nMassBins=histtofit->GetNbinsX();
    Double_t hmin=histtofit->GetBinLowEdge(2); // need wide range for <pt>
    Double_t hmax=histtofit->GetBinLowEdge(nMassBins-2); // need wide range for <pt>
    AliHFMassFitter fitter(histtofit,hmin,hmax,1);
    fitter.MassFitter(kFALSE);
    Double_t massFromFit=fitter.GetMean();
    Double_t sigmaFromFit=fitter.GetSigma();
    TF1* funcB2=fitter.GetBackgroundRecalcFunc();
    utils->AveragePt(averagePt[ipt],errorPt[ipt],ptbinsnew[ipt],ptbinsnew[ipt+1],hmasspt,massFromFit,sigmaFromFit,funcB2,2.5,4.5,0.,3.,1);
    if(averagePt[ipt]==0) {averagePt[ipt]=(ptbinsnew[ipt]+ptbinsnew[ipt+1])/2;}
  }
  printf("Average pt\n");
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++) printf("%f +- %f\n",averagePt[ipt],errorPt[ipt]);

  printf("Fill TGraphs for signal \n");
  //Fill TGraphs for signal
  TGraphAsymmErrors *gSignal[nptbinsnew];
  TGraphAsymmErrors *gSignalfs[nptbinsnew];
  TGraphAsymmErrors *gSignalBC1[nptbinsnew];
  TGraphAsymmErrors *gSignalBC2[nptbinsnew];
  TH1F *hSigmaFree[nphi];
  TH1F *hSigmaFixed[nphi];
  TH1F *hMean[nphi];
  TH1F *hMeanfs[nphi];
  TH1F *hChiSquare[nphi];
  TH1F *hChiSquarefs[nphi];
  for(Int_t i=0;i<nptbinsnew;i++){
    gSignal[i]=new TGraphAsymmErrors(nphi);
    gSignal[i]->SetName(Form("gasigpt%d",i));
    gSignal[i]->SetTitle(Form("Signal %.0f < #it{p}_{T} < %.0f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignal[i]->SetMarkerStyle(25);
    gSignalfs[i]=new TGraphAsymmErrors(nphi);
    gSignalfs[i]->SetName(Form("gasigfspt%d",i));
    gSignalfs[i]->SetTitle(Form("Signal (fixed sigma) %.0f < #it{p}_{T} < %.0f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignalfs[i]->SetMarkerStyle(21);
    gSignalBC1[i]=new TGraphAsymmErrors(nphi);
    gSignalBC1[i]->SetName(Form("gasigBC1pt%d",i));
    gSignalBC1[i]->SetTitle(Form("Signal (BC1) %.0f < #it{p}_{T} < %.0f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignalBC2[i]=new TGraphAsymmErrors(nphi);
    gSignalBC2[i]->SetName(Form("gasigBC2pt%d",i));
    gSignalBC2[i]->SetTitle(Form("Signal (BC2) %.0f < #it{p}_{T} < %.0f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
  }
  for(Int_t iphi=0; iphi<nphi; iphi++) {
    hSigmaFree[iphi] = new TH1F(Form("hSigmaFree_phi%d",iphi),Form("Sigma free - phi%d",iphi),nptbinsnew,ptbinsnew);
    hSigmaFree[iphi]->SetMarkerStyle(markers[0]);
    hSigmaFree[iphi]->SetMarkerColor(colors[0]);
    hSigmaFree[iphi]->SetLineColor(colors[0]);
    hSigmaFixed[iphi] = new TH1F(Form("hSigmaFixed_phi%d",iphi),Form("Sigma fixed - phi%d",iphi),nptbinsnew,ptbinsnew);
    hSigmaFixed[iphi]->SetMarkerStyle(markers[1]);
    hSigmaFixed[iphi]->SetMarkerColor(colors[1]);
    hSigmaFixed[iphi]->SetLineColor(colors[1]);

    hMean[iphi] = new TH1F(Form("hMean_phi%d",iphi),Form("Mean (sigma free) - phi%d",iphi),nptbinsnew,ptbinsnew);
    hMean[iphi]->SetMarkerStyle(markers[0]);
    hMean[iphi]->SetMarkerColor(colors[0]);
    hMean[iphi]->SetLineColor(colors[0]);
    hMeanfs[iphi] = new TH1F(Form("hMeanfs_phi%d",iphi),Form("Mean (sigma fixed) - phi%d",iphi),nptbinsnew,ptbinsnew);
    hMeanfs[iphi]->SetMarkerStyle(markers[1]);
    hMeanfs[iphi]->SetMarkerColor(colors[1]);
    hMeanfs[iphi]->SetLineColor(colors[1]);

    hChiSquare[iphi] = new TH1F(Form("hChiSquare_phi%d",iphi),Form("#chi^{2} (sigma free) - phi%d",iphi),nptbinsnew,ptbinsnew);
    hChiSquare[iphi]->SetMarkerStyle(markers[0]);
    hChiSquare[iphi]->SetMarkerColor(colors[0]);
    hChiSquare[iphi]->SetLineColor(colors[0]);
    hChiSquarefs[iphi] = new TH1F(Form("hChiSquarefs_phi%d",iphi),Form("#chi^{2} (sigma fixed) - phi%d",iphi),nptbinsnew,ptbinsnew);
    hChiSquarefs[iphi]->SetMarkerStyle(markers[1]);
    hChiSquarefs[iphi]->SetMarkerColor(colors[1]);
    hChiSquarefs[iphi]->SetLineColor(colors[1]);
  }
  
  FillSignalGraph(histlist,gSignal,gSignalfs,gSignalBC1,gSignalBC2,hSigmaFree,hSigmaFixed,hMean,hMeanfs,hChiSquare,hChiSquarefs,inoutanis,typeb,minMassForFit,maxMassForFit,rebin);

  //EP resolution
  Double_t resol=-1.;
  Double_t errorres=-1.;

  if(useAliHandlerForRes) {
    AliEventPlaneResolutionHandler* epres=new AliEventPlaneResolutionHandler();
    epres->SetEventPlane(evPlane);
    epres->SetResolutionOption(evPlaneRes);
    if(useNcollWeight)
      epres->SetUseNcollWeights();
    resol=epres->GetEventPlaneResolution(minCent,maxCent);
    delete epres;
  }
  else {
    TList* resolhist=LoadResolutionHistos(list);
    TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
    TH1F* hevplresos[3];
    TString namereso[3]={"Reso","Reso2","Reso3"};
    TH1F* htestversion=(TH1F*)resolhist->FindObject(Form("hEvPlaneReso%s",suffixcentr.Data()));
    cout<<htestversion<<endl;
    if(htestversion){
      printf("Old version of the task\n");
    }else{
      printf("New version of the task\n");
      namereso[0]="Reso1";
    }
    Int_t nSubRes=1;
    if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
       evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
    for(Int_t ires=0;ires<nSubRes;ires++){
      hevplresos[ires]=(TH1F*)resolhist->FindObject(Form("hEvPlane%s%s",namereso[ires].Data(),suffixcentr.Data()));
      }
    resol=GetEventPlaneResolution(errorres,hevplresos[0],hevplresos[1],hevplresos[2]);
  }

  printf("Event plane resolution %f\n",resol);
  printf("Compute v2 \n");
  //compute v2
  SetStyle();
  
  TCanvas *cv2fs =new TCanvas("v2_fs","v2 Fit with fixed sigma",1920,1080);
  TCanvas *cv2 =new TCanvas("v2","v2 - systematic on yield extraction",1920,1080);
  
  TGraphAsymmErrors *grelSystEff=new TGraphAsymmErrors(nptbinsnew);
  TGraphAsymmErrors *gv2=Computev2(gSignal,resol,averagePt,inoutanis,0x0);
  TGraphAsymmErrors *gv2fs=Computev2(gSignalfs,resol,averagePt,inoutanis,grelSystEff);
  TGraphAsymmErrors *gv2BC1=Computev2(gSignalBC1,resol,averagePt,inoutanis,0x0);
  TGraphAsymmErrors *gv2BC2=Computev2(gSignalBC2,resol,averagePt,inoutanis,0x0);
 
  gv2->SetName("gav2");
  gv2->GetYaxis()->SetTitle("v_{2} {EP}");
  gv2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gv2->SetTitle("");
  gv2->SetLineColor(colors[0]);
  gv2->SetMarkerColor(colors[0]);
  gv2->SetLineWidth(2);
  gv2->SetMarkerStyle(markers[0]);
  gv2->SetMarkerSize(1.5);
  gv2fs->SetName("gav2fs");
  gv2fs->GetYaxis()->SetTitle("v_{2} {EP}");
  gv2fs->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gv2fs->SetTitle("");
  gv2fs->SetLineColor(colors[1]);
  gv2fs->SetMarkerColor(colors[1]);
  gv2fs->SetLineWidth(2);
  gv2fs->SetMarkerStyle(markers[1]);
  gv2fs->SetMarkerSize(1.5);
  gv2BC1->SetName("gav2BC1");
  gv2BC1->GetYaxis()->SetTitle("v_{2} {EP}");
  gv2BC1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gv2BC1->SetTitle("");
  gv2BC1->SetLineColor(colors[2]);
  gv2BC1->SetMarkerColor(colors[2]);
  gv2BC1->SetLineWidth(2);
  gv2BC1->SetMarkerStyle(markers[2]);
  gv2BC1->SetMarkerSize(1.5);
  gv2BC2->SetName("gav2BC2");
  gv2BC2->GetYaxis()->SetTitle("v_{2} {EP}");
  gv2BC2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  gv2BC2->SetTitle("");
  gv2BC2->SetLineColor(colors[3]);
  gv2BC2->SetMarkerColor(colors[3]);
  gv2BC2->SetLineWidth(2);
  gv2BC2->SetMarkerStyle(markers[3]);
  gv2BC2->SetMarkerSize(1.5);
  grelSystEff->SetName("grelSystEff");
  grelSystEff->SetTitle("");
  grelSystEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  grelSystEff->GetYaxis()->SetTitle("Efficiency Syst.");
  grelSystEff->SetLineColor(colors[0]);
  grelSystEff->SetMarkerColor(colors[0]);
  grelSystEff->SetLineWidth(2);
  grelSystEff->SetMarkerStyle(markers[0]);
  grelSystEff->SetMarkerSize(1.5);
  
  //Prepare output file
  TFile *fout=new TFile(Form("%s/v2Output_%d_%d_%s%s%s.root",outputdir.Data(),minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()),"RECREATE");

  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    fout->cd();
    gSignal[ipt]->Write();
    gSignalfs[ipt]->Write();
    gSignalBC1[ipt]->Write();
    gSignalBC2[ipt]->Write();
  }
  for(Int_t iphi=0; iphi<nphi; iphi++) {
    hSigmaFree[iphi]->Write();
    hSigmaFixed[iphi]->Write();
    hMean[iphi]->Write();
    hMeanfs[iphi]->Write();
    hChiSquare[iphi]->Write();
    hChiSquarefs[iphi]->Write();
  }
  
  cv2->cd();
  gv2fs->Draw("AP");
  gv2fs->SetMinimum(-0.2);
  gv2fs->SetMaximum(0.6);
  gv2->Draw("PSAME");
  gv2BC1->Draw("PSAME");
  gv2BC2->Draw("PSAME");
  TLegend* leg=new TLegend(0.55,0.7,0.89,0.89);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->AddEntry(gv2fs,"Fixed sigma","lpe");
  leg->AddEntry(gv2,"Free sigma","lpe");
  leg->AddEntry(gv2BC1,"Bin counting 1","lpe");
  leg->AddEntry(gv2BC2,"Bin counting 2","lpe");
  leg->Draw();
  cv2fs->cd();
  gv2fs->Draw("AP");

  fout->cd();
  gv2->Write();
  gv2fs->Write();
  gv2BC1->Write();
  gv2BC2->Write();
  grelSystEff->Write();

  cv2->SaveAs(Form("%s/v2Output_%d_%d_%s%s%s.pdf",outputdir.Data(),minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()),"RECREATE");
}
//___________________________________________________________
Int_t FindPtBin(Int_t nbins, Double_t* array,Double_t value){
  for (Int_t i=0;i<nbins;i++){
    if(value>=array[i] && value<array[i+1]){
      return i;
    }
  }
  cout<<value<< " out of range "<<array[0]<<", "<<array[nbins]<<endl;
  return -1;
}

//__________________________________________________________
void DrawEventPlane(){
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s",partname.Data(),suffix.Data());
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
  if(evPlane==AliEventPlaneResolutionHandler::kTPCPosEta &&
     evPlaneRes==AliEventPlaneResolutionHandler::kTwoEtaSub){
    printf("Wrong setting of event plane resolution forced it to kTwoSub\n");
    evPlaneRes=AliEventPlaneResolutionHandler::kTwoRandSub;
  }
  TList* resolhist=LoadResolutionHistos(list);
  TGraphErrors* gResolVsCent=(TGraphErrors*)resolhist->FindObject("gResolVsCent");
  Double_t lowestRes=1;
  Double_t highestRes=0;
  for(Int_t ipt=0; ipt<gResolVsCent->GetN(); ipt++){    
    Double_t x,resolBin;
    gResolVsCent->GetPoint(ipt,x,resolBin);
    if(resolBin<lowestRes) lowestRes=resolBin;
    if(resolBin>highestRes) highestRes=resolBin;
  }
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
  TH2F* hevpls=(TH2F*)resolhist->FindObject(Form("hEvPlane%s",suffixcentr.Data()));
  TH1F* hevplresos[3];
  TString namereso[3]={"Reso","Reso2","Reso3"};
  TH1F* htestversion=(TH1F*)resolhist->FindObject(Form("hEvPlaneReso%s",suffixcentr.Data()));
  if(htestversion){
    printf("Old version of the task\n");
  }else{
    printf("New version of the task\n");
    namereso[0]="Reso1";
  }

  Int_t nSubRes=1;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
     evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
  for(Int_t ires=0;ires<nSubRes;ires++){
    hevplresos[ires]=(TH1F*)resolhist->FindObject(Form("hEvPlane%s%s",namereso[ires].Data(),suffixcentr.Data()));
  }

  SetStyle();
  
  TH1F* htpc = (TH1F*)hevpls->ProjectionX();
  TH1F* hv0 = (TH1F*)hevpls->ProjectionY();
  TH1F* hUsedPlane;
  if(evPlane==AliEventPlaneResolutionHandler::kVZERO ||
     evPlane==AliEventPlaneResolutionHandler::kVZEROA ||
     evPlane==AliEventPlaneResolutionHandler::kVZEROC) hUsedPlane=(TH1F*)hv0->Clone("hUsedPlane");
  else hUsedPlane=(TH1F*)htpc->Clone("hUsedPlane");

  TCanvas* cvtotevpl=new TCanvas(Form("cvtotevpl%s",suffixcentr.Data()),Form("Ev plane %s",suffixcentr.Data()),1920,1080);
  cvtotevpl->Divide(3,1);
  cvtotevpl->cd(1);
  hevpls->SetTitleFont(42);
  hevpls->Draw("COLZ");
  cvtotevpl->cd(2);
  htpc->Draw();
  htpc->Fit("pol0");
  cvtotevpl->cd(3);
  hv0->Draw();
  hv0->Fit("pol0");

  TCanvas* cevpl=new TCanvas(Form("cevpl%s",suffixcentr.Data()),"",1920,1080);
  htpc->Draw();
  htpc->Fit("pol0");
  cevpl->SaveAs(Form("%s/EvPlaneDistTPC.pdf",outputdir.Data()));
  
  TCanvas* cvfitevpl=new TCanvas(Form("cvditevpl%s",suffixcentr.Data()),Form("Fit to Ev plane %s",suffixcentr.Data()));
  cvfitevpl->cd();
  hUsedPlane->Draw();
  TF1* four2=new TF1("four2","[0]*(1+2*[1]*cos(2*x))",hUsedPlane->GetXaxis()->GetXmin(),hUsedPlane->GetXaxis()->GetXmax());
  TF1* four2s=new TF1("four2s","[0]*(1+2*[1]*sin(2*x))",hUsedPlane->GetXaxis()->GetXmin(),hUsedPlane->GetXaxis()->GetXmax());
  hUsedPlane->Fit(four2s);
  hUsedPlane->Fit(four2);
  four2->SetLineColor(1);
  four2s->SetLineColor(4);
  four2->Draw("same");
  four2s->Draw("same");
  hUsedPlane->SetMaximum(hUsedPlane->GetMaximum()*1.07);
  TLatex* tsin=new TLatex(0.15,0.84,Form("1+2*(%.4f)*sin(2*#Phi_{EP})",four2s->GetParameter(1)));
  tsin->SetNDC();
  tsin->SetTextColor(4);
  tsin->Draw();
  TLatex* tcos=new TLatex(0.15,0.77,Form("1+2*(%.4f)*cos(2*#Phi_{EP})",four2->GetParameter(1)));
  tcos->SetNDC();
  tcos->SetTextColor(1);
  tcos->Draw();
  
  // Compute the second Fourier component for sine and cosine
  Double_t aveCos2Phi=0.;
  Double_t aveSin2Phi=0.;
  Double_t counts=0.;
  for(Int_t i=1; i<=hUsedPlane->GetNbinsX(); i++){
    Double_t co=hUsedPlane->GetBinContent(i);
    Double_t phi=hUsedPlane->GetBinCenter(i);
    aveCos2Phi+=co*TMath::Cos(2*phi);
    aveSin2Phi+=co*TMath::Sin(2*phi);
    counts+=co;
  }
  if(counts>0){ 
    aveCos2Phi/=counts;
    aveSin2Phi/=counts;
  }
  printf("\n------ Fourier 2nd components of EP distribution ------\n");
  printf("<cos(2*Psi_EP)>=%f\n",aveCos2Phi);
  printf("<sin(2*Psi_EP)>=%f\n",aveSin2Phi);
  printf("\n");


  TCanvas* cvtotevplreso=new TCanvas(Form("cvtotevplreso%s",suffixcentr.Data()),Form("Ev plane Resolution %s",suffixcentr.Data()));
  cvtotevplreso->cd();
  hevplresos[0]->SetTitle("TPC cos2(#Psi_{A}-#Psi_{B})");
  hevplresos[0]->Draw();
  gStyle->SetOptStat(0);
  if(nSubRes>1){
    hevplresos[1]->SetLineColor(2);
    hevplresos[2]->SetLineColor(3);
    hevplresos[0]->SetTitle("cos2(#Psi_{TPC}-#Psi_{V0A})");
    hevplresos[1]->SetTitle("cos2(#Psi_{TPC}-#Psi_{V0C})");
    hevplresos[2]->SetTitle("cos2(#Psi_{V0A}-#Psi_{V0C})");
    hevplresos[1]->Draw("SAME");
    hevplresos[2]->Draw("SAME");
  }
  TLegend *leg = cvtotevplreso->BuildLegend();
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  Double_t error;
  Double_t resolFull=GetEventPlaneResolution(error,hevplresos[0],hevplresos[1],hevplresos[2]);
  
  TPaveText* pvreso=new TPaveText(0.1,0.35,0.6,0.48,"NDC");
  pvreso->SetBorderSize(0);
  pvreso->SetFillStyle(0);
  pvreso->AddText(Form("Number of events = %.0f\n",hevplresos[0]->GetEntries()));
  pvreso->AddText(Form("Resolution on full event = %.4f#pm%.4f\n",resolFull,error));
  pvreso->Draw();

  gResolVsCent->SetMarkerStyle(20);
  gResolVsCent->GetXaxis()->SetTitle("Centrality");
  gResolVsCent->GetYaxis()->SetTitle("EP Resolution");
 
  Int_t colorave=colors[5];
  if(useNcollWeight) colorave=colors[4];

  TCanvas* cresvscent=new TCanvas("cresvscent","ResolVsCent");
  cresvscent->cd();
  gResolVsCent->Draw("AP");
  TLine* lin=new TLine(gResolVsCent->GetXaxis()->GetXmin(),resolFull,gResolVsCent->GetXaxis()->GetXmax(),resolFull);
  lin->SetLineColor(colorave);
  lin->SetLineWidth(2);
  lin->Draw();
  TLatex* taveres=new TLatex(gResolVsCent->GetXaxis()->GetXmin()+0.5,resolFull*1.01,Form("<R_{2}>=%.4f",resolFull));
  taveres->SetTextColor(colorave);
  taveres->Draw();
  printf("\n\n---- Syst from centrality dependence of resolution ----\n");
  printf("MinRes=%f  MaxRes=%f AveRes=%f ---> Syst: -%f%% +%f%%\n",lowestRes,highestRes,resolFull,(resolFull-lowestRes)/resolFull,(highestRes-resolFull)/resolFull);

  TFile* fout=new TFile(Form("%s/EvPlanecentr%d-%d.root",outputdir.Data(),minCent,maxCent),"recreate");
  fout->cd();
  hevpls->Write();
  for(Int_t ires=0;ires<nSubRes;ires++)hevplresos[ires]->Write();
  gResolVsCent->Write();

  cresvscent->SaveAs(Form("%s/EvPlaneResolutionVsCentrality.pdf",outputdir.Data()));
}

void SetStyle(Int_t optfit) {
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(optfit);
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetStatFont(42);
  gStyle->SetStatY(0.89);
  gStyle->SetStatX(0.89);
  gStyle->SetTitleFont(42,"xyzg");
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
}
