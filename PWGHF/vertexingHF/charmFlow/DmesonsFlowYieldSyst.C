#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TASImage.h>
#include <TPad.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>

#include "AliHFMassFitter.h"
#include "AliHFMassFitterVAR.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliEventPlaneResolutionHandler.h"
#include "AliVertexingHFUtils.h"

#endif

//methods for the extraction yield systematics for the v2 computed with the Event Plane method
//Author: Fabrizio Grosa, INFN Turin grosa@to.infn.it

//*************************************************//
//                                                 //
//      Main Function: DmesonsFlowYieldSyst()      //
//                                                 //
//*************************************************//

//_________________________________________________________________
//GLOBAL VARIABLES TO BE SET
//input file
const TString filename="$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/AnalysisResults_v2_3050_step2_EP_VZERO.root";
const TString suffix="_Topod0Cut_VZERO_EP";
const TString partname="Dplus";
const Int_t minCent=30;
const Int_t maxCent=50;

const TString outputdir = "Cent3050/v2";

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
const Double_t ptbinsnew[nptbinsnew+1]={2.,3.,4.,5.,6.,7.,8.,10.,12.,16.,24.};
//const Int_t nptbinsnew=16;
//const Double_t ptbinsnew[nptbinsnew+1]={2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,10.,12.,16.,24.};
const Int_t nphibins=4;
const Double_t phibinslim[nphibins+1]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};

// mass fit configuration
const Int_t nReb=5;
const Int_t rebin[nptbinsnew]={2,3,4,5,6};
enum {kGaus=0, kDoubleGaus, kReflTempl};
const Int_t types=kGaus;
const Int_t nBkgFcn=3;
const Int_t typeb[nBkgFcn]={AliHFMassFitter::kExpo,AliHFMassFitter::kPol2,AliHFMassFitter::kLin};
const Int_t nMins=5;
const Double_t minMassForFit[nMins]={1.66,1.68,1.70,1.72,1.76};
const Int_t nMaxs=5;
const Double_t maxMassForFit[nMaxs]={2.08,2.06,2.04,2.02,2.00};
const Double_t nSigmaForCounting=3.5;
const Bool_t fixAlsoMass=kFALSE;
const Double_t maxchi=2.0;
Bool_t useTemplD0Refl=kFALSE;
TString rflFitType="DoubleGaus";
TString fileNameMCD0refl="../reflections/reflections_fitted_DoubleGaus.root";

//not to be set
Int_t minPtBin[nptbinsnew]={-1,-1,-1,-1};
Int_t maxPtBin[nptbinsnew]={-1,-1,-1,-1};
const Double_t effInOverEffOut=1.03;
Double_t massD;

const Int_t colors[] = {kRed+1,kBlack,kBlue+1,kGreen+2,kOrange+7,kBlue-7};
const Int_t markers[] = {kFullSquare,kFullCircle,kFullTriangleUp,kFullDiamond,kOpenSquare,kOpenCircle,kOpenTriangleUp,kOpenDiamond};

//_________________________________________________________________
//METHODS PROTOTYPES
Int_t DmesonsFlowYieldSyst(Bool_t inoutanis=kTRUE);
TList *LoadMassHistos(TList *inputlist,Bool_t inoutanis);
TList* LoadResolutionHistos(TList *inputlist);
Int_t FindPtBin(Int_t nbins, Double_t* array,Double_t value);
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2,TGraphAsymmErrors **gSigmaFree,TGraphAsymmErrors **gSigmaFixed, TGraphAsymmErrors **gChiSquareFree, TGraphAsymmErrors **gChiSquareFixed, Bool_t inoutanis, Int_t bkgfunc, Double_t minfit, Double_t maxfit, Int_t rebin);
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Bool_t inoutanis, TGraphAsymmErrors *gRelSystEff);
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3);
Bool_t DefinePtBins(AliRDHFCuts *cutobj);
Int_t LoadRefGraphs(TString reffilename, TH1F **hRawYieldRef, TH1F **hRawYieldfsRef, TH1F **hRawYieldBC1Ref, TH1F **hRawYieldBC2Ref, TGraphAsymmErrors *&gv2Ref, TGraphAsymmErrors *&gv2fsRef, TGraphAsymmErrors *&gv2BC1Ref, TGraphAsymmErrors *&gv2BC2Ref, Bool_t inoutanis);
void GetMinMaxHisto(TH1F* histo, Double_t &xmin, Double_t &xmax);
void DivideCanvas(TCanvas* c, Int_t nPtBins);
void SetStyle(Int_t optfit=0);
Bool_t LoadD0toKpiMCHistos(TList *outlist);

//_________________________________________________________________
//METHODS IMPLEMENTATION
Int_t DmesonsFlowYieldSyst(Bool_t inoutanis){
  
  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}
  
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s",partname.Data(),suffix.Data());
  
  AliRDHFCuts *cutsobj=0x0;
  //Load input data from AliAnalysisTaskSEHFv2
  TFile *f = TFile::Open(filename.Data());
  if(!f){
    printf("file %s not found, please check file name\n",filename.Data());return 1;
  }
  TDirectoryFile* dir=(TDirectoryFile*)f->Get(dirname.Data());
  if(!dir){
    printf("Directory %s not found, please check dir name\n",dirname.Data());return 2;
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
    printf("list %s not found in file, please check list name\n",listname.Data());return 3;
  }
  if(!cutsobj){
    printf("cut object not found in file, please check keylist number\n");return 4;
  }
  //Define new pt bins
  if(!DefinePtBins(cutsobj)){
    printf("cut not define pt bins\n");return 5;
  }
  
  //Load mass histograms corresponding to the required centrality, pt range and phi binning
  printf("Load mass histos \n");
  TList *histlist=LoadMassHistos(list,inoutanis);
  
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
      return 6;
    }
    TH1F *histtofit = (TH1F*)hmasspt->ProjectionY("_py",binMin,binMax);
    Int_t nMassBins=histtofit->GetNbinsX();
    Double_t hmin=histtofit->GetBinLowEdge(2); // need wide range for <pt>
    Double_t hmax=histtofit->GetBinLowEdge(nMassBins-2); // need wide range for <pt>
    //AliHFMassFitter fitter(histtofit,hmin,hmax,1);
    AliHFMassFitterVAR *fitter=new AliHFMassFitterVAR(histtofit,hmin,hmax,1,typeb[0],types);
    if(useTemplD0Refl){
      Printf("USE TEMPLATE FOR AVERAGE Pt");
      TH1F *hrflTempl=(TH1F*)(histlist->FindObject(Form("histRfl_%d",ipt)))->Clone(Form("histrfl_%d",ipt));
      if(!hrflTempl) {Printf("histRfl_%d not found",ipt); return 7;}
      TH1F *hsigMC=(TH1F*)(histlist->FindObject(Form("histSgn_%d",ipt)))->Clone(Form("histsgn_%d",ipt));
      if(!hsigMC) {Printf("histSgn_%d not found",ipt); return 8;}
      fitter->SetTemplateReflections(hrflTempl);
      Float_t sOverRef=(hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999)))/(hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999)));
      Printf("R OVER S = %f",sOverRef);
      fitter->SetFixReflOverS(sOverRef,kTRUE);
    }
    fitter->MassFitter(kFALSE);
    Double_t massFromFit=fitter->GetMean();
    Double_t sigmaFromFit=fitter->GetSigma();
    TF1* funcB2=fitter->GetBackgroundRecalcFunc();
    utils->AveragePt(averagePt[ipt],errorPt[ipt],ptbinsnew[ipt],ptbinsnew[ipt+1],hmasspt,massFromFit,sigmaFromFit,funcB2,2.5,4.5,0.,3.,1);
  }
  printf("Average pt\n");
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++) printf("%f +- %f\n",averagePt[ipt],errorPt[ipt]);
  
  printf("Fill TGraphs for signal \n");
  //Fill TGraphs for signal
  TGraphAsymmErrors *gSignal[nptbinsnew];
  TGraphAsymmErrors *gSignalfs[nptbinsnew];
  TGraphAsymmErrors *gSignalBC1[nptbinsnew];
  TGraphAsymmErrors *gSignalBC2[nptbinsnew];
  TGraphAsymmErrors *gSigmaFree[nptbinsnew];
  TGraphAsymmErrors *gSigmaFixed[nptbinsnew];
  TGraphAsymmErrors *gChiSquareFree[nptbinsnew];
  TGraphAsymmErrors *gChiSquareFixed[nptbinsnew];
  for(Int_t i=0;i<nptbinsnew;i++){
    gSignal[i]=new TGraphAsymmErrors(nphi);
    gSignal[i]->SetName(Form("gasigpt%d",i));
    gSignal[i]->SetTitle(Form("Signal %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignal[i]->SetMarkerStyle(25);
    gSignalfs[i]=new TGraphAsymmErrors(nphi);
    gSignalfs[i]->SetName(Form("gasigfspt%d",i));
    gSignalfs[i]->SetTitle(Form("Signal (fixed sigma) %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignalfs[i]->SetMarkerStyle(21);
    gSignalBC1[i]=new TGraphAsymmErrors(nphi);
    gSignalBC1[i]->SetName(Form("gasigBC1pt%d",i));
    gSignalBC1[i]->SetTitle(Form("Signal (BC1) %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSignalBC2[i]=new TGraphAsymmErrors(nphi);
    gSignalBC2[i]->SetName(Form("gasigBC2pt%d",i));
    gSignalBC2[i]->SetTitle(Form("Signal (BC2) %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    
    gSigmaFree[i]=new TGraphAsymmErrors(nphi);
    gSigmaFree[i]->SetName(Form("gasigmafree%d",i));
    gSigmaFree[i]->SetTitle(Form("Sigma free %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gSigmaFixed[i]=new TGraphAsymmErrors(nphi);
    gSigmaFixed[i]->SetName(Form("gasigmafixed%d",i));
    gSigmaFixed[i]->SetTitle(Form("Sigma fixed %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));

    gChiSquareFree[i]=new TGraphAsymmErrors(nphi);
    gChiSquareFree[i]->SetName(Form("gachifree%d",i));
    gChiSquareFree[i]->SetTitle(Form("ChiSquare free %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
    gChiSquareFixed[i]=new TGraphAsymmErrors(nphi);
    gChiSquareFixed[i]->SetName(Form("gachifixed%d",i));
    gChiSquareFixed[i]->SetTitle(Form("ChiSquare fixed %.1f < #it{p}_{T} < %.1f GeV/c;#Delta#phi bin;Counts",ptbinsnew[i],ptbinsnew[i+1]));
  }
  
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
    Int_t nSubRes=1;
    TH1F* htestversion=(TH1F*)resolhist->FindObject(Form("hEvPlane%s%s",namereso[0].Data(),suffixcentr.Data()));
    if(htestversion){
      printf("Old version of the task\n");
    }else{
      printf("New version of the task\n");
      namereso[0]="Reso1";
    }
    if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
       evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
    for(Int_t ires=0;ires<nSubRes;ires++){
      hevplresos[ires]=(TH1F*)resolhist->FindObject(Form("hEvPlane%s%s",namereso[ires].Data(),suffixcentr.Data()));
    }
    resol=GetEventPlaneResolution(errorres,hevplresos[0],hevplresos[1],hevplresos[2]);
  }

  SetStyle();

  printf("Event plane resolution %f\n",resol);
  printf("Compute v2 \n");
  //compute v2
  
  TString fine="";
  if(nptbinsnew>=15) fine="_fineptbin";
  TString aniss="";
  if(nphi==2) aniss+="anis";
  
  //load reference graphs
  TGraphAsymmErrors* gv2Ref=0x0;
  TGraphAsymmErrors* gv2fsRef=0x0;
  TGraphAsymmErrors* gv2BC1Ref=0x0;
  TGraphAsymmErrors* gv2BC2Ref=0x0;
  TH1F** hRawYieldRef = new TH1F*[nphi];
  TH1F** hRawYieldfsRef = new TH1F*[nphi];
  TH1F** hRawYieldBC1Ref = new TH1F*[nphi];
  TH1F** hRawYieldBC2Ref = new TH1F*[nphi];
  TString reffilename = Form("%s/v2Output_%d_%d_%s%s%s.root",outputdir.Data(),minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data());
  Int_t loadref=LoadRefGraphs(reffilename, hRawYieldRef, hRawYieldfsRef, hRawYieldBC1Ref, hRawYieldBC2Ref, gv2Ref, gv2fsRef, gv2BC1Ref, gv2BC2Ref, inoutanis);
  if(loadref>0) {return 7;}
  
  TCanvas *cv2 =new TCanvas("cv2","v2 - systematic on yield extraction",1920,1080);
  DivideCanvas(cv2,nptbinsnew);
  TCanvas *cv2VsTrial =new TCanvas("cv2VsTrial","v2 vs. Trial - systematic on yield extraction",1920,1080);
  DivideCanvas(cv2VsTrial,nptbinsnew);

  TH1F** hv2 = new TH1F*[nptbinsnew];
  TH1F** hv2fs = new TH1F*[nptbinsnew];
  TH1F** hv2BC1 = new TH1F*[nptbinsnew];
  TH1F** hv2BC2 = new TH1F*[nptbinsnew];
  TH1F** hv2VsTrial = new TH1F*[nptbinsnew];
  TH1F** hv2fsVsTrial = new TH1F*[nptbinsnew];
  TH1F** hv2BC1VsTrial = new TH1F*[nptbinsnew];
  TH1F** hv2BC2VsTrial = new TH1F*[nptbinsnew];
  
  TCanvas** cRawYield=new TCanvas*[nphi];
  TCanvas** cSigma=new TCanvas*[nphi];
  TCanvas** cChiSquare=new TCanvas*[nphi];
  
  TH1F*** hRawYield=new TH1F**[nphi];
  TH1F*** hRawYieldfs=new TH1F**[nphi];
  TH1F*** hRawYieldBC1=new TH1F**[nphi];
  TH1F*** hRawYieldBC2=new TH1F**[nphi];
  TH1F*** hSigmaFixed=new TH1F**[nphi];
  TH1F*** hSigmaFree=new TH1F**[nphi];
  TH1F*** hChiSquareFixed=new TH1F**[nphi];
  TH1F*** hChiSquareFree=new TH1F**[nphi];
  
  for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
    hRawYield[iPhi]=new TH1F*[nptbinsnew];
    hRawYieldfs[iPhi]=new TH1F*[nptbinsnew];
    hRawYieldBC1[iPhi]=new TH1F*[nptbinsnew];
    hRawYieldBC2[iPhi]=new TH1F*[nptbinsnew];

    hSigmaFree[iPhi]=new TH1F*[nptbinsnew];
    hSigmaFixed[iPhi]=new TH1F*[nptbinsnew];
    hChiSquareFree[iPhi]=new TH1F*[nptbinsnew];
    hChiSquareFixed[iPhi]=new TH1F*[nptbinsnew];
    
    cRawYield[iPhi]=new TCanvas(Form("cRawYield_phi%d",iPhi),"Y - systematic on yield extraction",1920,1080);
    cSigma[iPhi]=new TCanvas(Form("cSigma_phi%d",iPhi),"Y - systematic on yield extraction",1920,1080);
    cChiSquare[iPhi]=new TCanvas(Form("cChiSquare_phi%d",iPhi),"Y - systematic on yield extraction",1920,1080);
    DivideCanvas(cRawYield[iPhi],nptbinsnew);
    DivideCanvas(cSigma[iPhi],nptbinsnew);
    DivideCanvas(cChiSquare[iPhi],nptbinsnew);
  }
  
  TH1F** hv2Syst = new TH1F*[nptbinsnew];
  TH1F** hv2fsSyst = new TH1F*[nptbinsnew];
  TH1F** hv2BC1Syst = new TH1F*[nptbinsnew];
  TH1F** hv2BC2Syst = new TH1F*[nptbinsnew];
  TH1F** hv2Syst2 = new TH1F*[nptbinsnew];
  TH1F** hv2fsSyst2 = new TH1F*[nptbinsnew];
  TH1F** hv2BC1Syst2 = new TH1F*[nptbinsnew];
  TH1F** hv2BC2Syst2 = new TH1F*[nptbinsnew];
  
  TH1F* hSyst1 = new TH1F("hSyst1","",nptbinsnew,ptbinsnew);
  hSyst1->SetLineWidth(2);
  hSyst1->SetStats(kFALSE);
  hSyst1->SetLineColor(colors[1]);
  hSyst1->GetYaxis()->SetTitle("#sqrt{RMS^{2} + shift^{2}} (fixed sigma)");
  hSyst1->GetYaxis()->SetTitleOffset(1.4);
  hSyst1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  
  const Int_t nbins=50;
  
  const Int_t nTrials=nReb*nMins*nMaxs*nBkgFcn;
  
  for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
    
    Double_t v2staterr = gv2fsRef->GetErrorY(iPt);
    
    hv2[iPt] = new TH1F(Form("hv2_%d",iPt),"",nbins,-3*v2staterr,3*v2staterr);
    hv2fs[iPt] = new TH1F(Form("hv2fs_%d",iPt),"",nbins,-3*v2staterr,3*v2staterr);
    hv2BC1[iPt] = new TH1F(Form("hv2BC1_%d",iPt),"",nbins,-3*v2staterr,3*v2staterr);
    hv2BC2[iPt] = new TH1F(Form("hv2BC2_%d",iPt),"",nbins,-3*v2staterr,3*v2staterr);
    
    hv2VsTrial[iPt] = new TH1F(Form("hv2VsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2fsVsTrial[iPt] = new TH1F(Form("hv2fsVsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2BC1VsTrial[iPt] = new TH1F(Form("hv2BC1VsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2BC2VsTrial[iPt] = new TH1F(Form("hv2BC2VsTrial_%d",iPt),"",nTrials,-0.5,nTrials-0.5);
    hv2VsTrial[iPt]->SetStats(kFALSE);
    hv2fsVsTrial[iPt]->SetStats(kFALSE);
    hv2BC1VsTrial[iPt]->SetStats(kFALSE);
    hv2BC2VsTrial[iPt]->SetStats(kFALSE);

    hv2Syst[iPt] = new TH1F(Form("hv2Syst_%d",iPt),"",4,-0.5,3.5);
    hv2fsSyst[iPt] = new TH1F(Form("hv2fsSyst_%d",iPt),"",4,-0.5,3.5);
    hv2BC1Syst[iPt] = new TH1F(Form("hv2BC1Syst_%d",iPt),"",4,-0.5,3.5);
    hv2BC2Syst[iPt] = new TH1F(Form("hv2BC2Syst_%d",iPt),"",4,-0.5,3.5);
    hv2Syst[iPt]->SetStats(kFALSE);
    hv2fsSyst[iPt]->SetStats(kFALSE);
    hv2BC1Syst[iPt]->SetStats(kFALSE);
    hv2BC2Syst[iPt]->SetStats(kFALSE);
    
    hv2Syst2[iPt] = new TH1F(Form("hv2Syst2_%d",iPt),"",4,-0.4,3.6);
    hv2fsSyst2[iPt] = new TH1F(Form("hv2fsSyst2_%d",iPt),"",4,-0.4,3.6);
    hv2BC1Syst2[iPt] = new TH1F(Form("hv2BC1Syst2_%d",iPt),"",4,-0.4,3.6);
    hv2BC2Syst2[iPt] = new TH1F(Form("hv2BC2Syst2_%d",iPt),"",4,-0.4,3.6);
    hv2Syst2[iPt]->SetStats(kFALSE);
    hv2fsSyst2[iPt]->SetStats(kFALSE);
    hv2BC1Syst2[iPt]->SetStats(kFALSE);
    hv2BC2Syst2[iPt]->SetStats(kFALSE);
    
    hv2[iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
    hv2[iPt]->GetYaxis()->SetTitle("Entries");
    hv2[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2fs[iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
    hv2fs[iPt]->GetYaxis()->SetTitle("Entries");
    hv2fs[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC1[iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
    hv2BC1[iPt]->GetYaxis()->SetTitle("Entries");
    hv2BC1[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC2[iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
    hv2BC2[iPt]->GetYaxis()->SetTitle("Entries");
    hv2BC2[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2[iPt]->SetLineColor(colors[0]);
    hv2fs[iPt]->SetLineColor(colors[1]);
    hv2BC1[iPt]->SetLineColor(colors[2]);
    hv2BC2[iPt]->SetLineColor(colors[3]);
    
    hv2VsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2VsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2VsTrial[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2fsVsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2fsVsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2fsVsTrial[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC1VsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2BC1VsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2BC1VsTrial[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC2VsTrial[iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
    hv2BC2VsTrial[iPt]->GetXaxis()->SetTitle("Trial #");
    hv2BC2VsTrial[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2VsTrial[iPt]->SetLineColor(colors[0]);
    hv2fsVsTrial[iPt]->SetLineColor(colors[1]);
    hv2BC1VsTrial[iPt]->SetLineColor(colors[2]);
    hv2BC2VsTrial[iPt]->SetLineColor(colors[3]);
    hv2VsTrial[iPt]->SetMarkerColor(colors[0]);
    hv2fsVsTrial[iPt]->SetMarkerColor(colors[1]);
    hv2BC1VsTrial[iPt]->SetMarkerColor(colors[2]);
    hv2BC2VsTrial[iPt]->SetMarkerColor(colors[3]);
    hv2VsTrial[iPt]->SetMarkerSize(0.5);
    hv2fsVsTrial[iPt]->SetMarkerSize(0.5);
    hv2BC1VsTrial[iPt]->SetMarkerSize(0.5);
    hv2BC2VsTrial[iPt]->SetMarkerSize(0.5);
    
    hv2Syst[iPt]->GetYaxis()->SetTitle("MEAN v_{2}-v_{2}^{ref}");
    hv2Syst[iPt]->GetYaxis()->SetTitleOffset(1.5);
    hv2Syst[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2fsSyst[iPt]->GetYaxis()->SetTitle("v_{2} syst. uncertainty");
    hv2fsSyst[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC1Syst[iPt]->GetYaxis()->SetTitle("v_{2} syst. uncertainty");
    hv2BC1Syst[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC2Syst[iPt]->GetYaxis()->SetTitle("v_{2} syst. uncertainty");
    hv2BC2Syst[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2Syst[iPt]->SetLineColor(colors[0]);
    hv2fsSyst[iPt]->SetLineColor(colors[1]);
    hv2BC1Syst[iPt]->SetLineColor(colors[2]);
    hv2BC2Syst[iPt]->SetLineColor(colors[3]);
    hv2Syst[iPt]->SetMarkerColor(colors[0]);
    hv2fsSyst[iPt]->SetMarkerColor(colors[1]);
    hv2BC1Syst[iPt]->SetMarkerColor(colors[2]);
    hv2BC2Syst[iPt]->SetMarkerColor(colors[3]);
    hv2Syst[iPt]->SetMarkerStyle(markers[0]);
    hv2fsSyst[iPt]->SetMarkerStyle(markers[1]);
    hv2BC1Syst[iPt]->SetMarkerStyle(markers[2]);
    hv2BC2Syst[iPt]->SetMarkerStyle(markers[3]);
    hv2Syst[iPt]->GetXaxis()->SetBinLabel(1,"free sigma");
    hv2Syst[iPt]->GetXaxis()->SetBinLabel(2,"fixed sigma");
    hv2Syst[iPt]->GetXaxis()->SetBinLabel(3,"BC1");
    hv2Syst[iPt]->GetXaxis()->SetBinLabel(4,"BC2");
    
    hv2Syst2[iPt]->GetYaxis()->SetTitle("MEAN v_{2}-v_{2}^{ref}");
    hv2Syst2[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2fsSyst2[iPt]->GetYaxis()->SetTitle("v_{2} syst. uncertainty");
    hv2fsSyst2[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC1Syst2[iPt]->GetYaxis()->SetTitle("v_{2} syst. uncertainty");
    hv2BC1Syst2[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2BC2Syst2[iPt]->GetYaxis()->SetTitle("v_{2} syst. uncertainty");
    hv2BC2Syst2[iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[iPt],ptbinsnew[iPt+1]));
    hv2Syst2[iPt]->SetLineColor(colors[0]);
    hv2fsSyst2[iPt]->SetLineColor(colors[1]);
    hv2BC1Syst2[iPt]->SetLineColor(colors[2]);
    hv2BC2Syst2[iPt]->SetLineColor(colors[3]);
    hv2Syst2[iPt]->SetMarkerColor(colors[0]);
    hv2fsSyst2[iPt]->SetMarkerColor(colors[1]);
    hv2BC1Syst2[iPt]->SetMarkerColor(colors[2]);
    hv2BC2Syst2[iPt]->SetMarkerColor(colors[3]);
    hv2Syst2[iPt]->SetMarkerStyle(markers[4]);
    hv2fsSyst2[iPt]->SetMarkerStyle(markers[5]);
    hv2BC1Syst2[iPt]->SetMarkerStyle(markers[6]);
    hv2BC2Syst2[iPt]->SetMarkerStyle(markers[7]);
    
    for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
      
      Double_t minraw=hRawYieldRef[iPhi]->GetBinContent(iPt+1)*(1-0.5);
      Double_t maxraw=hRawYieldRef[iPhi]->GetBinContent(iPt+1)*(1+0.5);
      
      hRawYield[iPhi][iPt] = new TH1F(Form("hRawYield_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYieldfs[iPhi][iPt] = new TH1F(Form("hRawYieldfs_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYieldBC1[iPhi][iPt] = new TH1F(Form("hRawYieldBC1_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYieldBC2[iPhi][iPt] = new TH1F(Form("hRawYieldBC2_%d_%d",iPhi,iPt),"",nbins,minraw,maxraw);
      hRawYield[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYield[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYield[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYieldfs[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYieldfs[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYieldfs[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYieldBC1[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYieldBC1[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYieldBC1[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYieldBC2[iPhi][iPt]->GetXaxis()->SetTitle("raw yield");
      hRawYieldBC2[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hRawYieldBC2[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hRawYield[iPhi][iPt]->SetLineColor(colors[0]);
      hRawYieldfs[iPhi][iPt]->SetLineColor(colors[1]);
      hRawYieldBC1[iPhi][iPt]->SetLineColor(colors[2]);
      hRawYieldBC2[iPhi][iPt]->SetLineColor(colors[3]);
      
      hSigmaFree[iPhi][iPt] = new TH1F(Form("hSigmaFree_%d_%d",iPhi,iPt),"",nbins,0,0.025);
      hSigmaFree[iPhi][iPt]->GetXaxis()->SetTitle("width (GeV/c^{2})");
      hSigmaFree[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hSigmaFree[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hSigmaFree[iPhi][iPt]->SetLineColor(colors[0]);
      hSigmaFixed[iPhi][iPt] = new TH1F(Form("hSigmaFixed_%d_%d",iPhi,iPt),"",nbins,0,0.025);
      hSigmaFixed[iPhi][iPt]->GetXaxis()->SetTitle("width (GeV/c^{2})");
      hSigmaFixed[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hSigmaFixed[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hSigmaFixed[iPhi][iPt]->SetLineColor(colors[1]);
      
      hChiSquareFree[iPhi][iPt] = new TH1F(Form("hChiSquareFree_%d_%d",iPhi,iPt),"",30,0,3.);
      hChiSquareFree[iPhi][iPt]->GetXaxis()->SetTitle("#chi^{2}");
      hChiSquareFree[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hChiSquareFree[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hChiSquareFree[iPhi][iPt]->SetLineColor(colors[0]);
      hChiSquareFixed[iPhi][iPt] = new TH1F(Form("hChiSquareFixed_%d_%d",iPhi,iPt),"",30,0,3.);
      hChiSquareFixed[iPhi][iPt]->GetXaxis()->SetTitle("#chi^{2}");
      hChiSquareFixed[iPhi][iPt]->GetYaxis()->SetTitle("Entries");
      hChiSquareFixed[iPhi][iPt]->SetTitle(Form("#phi%d, %.1f < #it{p}_{T} < %.1f GeV/c",iPhi,ptbinsnew[iPt],ptbinsnew[iPt+1]));
      hChiSquareFixed[iPhi][iPt]->SetLineColor(colors[1]);
    }
  }
  
  Int_t iTrial=0;
  for(Int_t iReb=0; iReb<nReb; iReb++) {
    for(Int_t iMin=0; iMin<nMins; iMin++) {
      for(Int_t iMax=0; iMax<nMaxs; iMax++) {
        for(Int_t iBkgFcn=0; iBkgFcn<nBkgFcn; iBkgFcn++) {
         
          Double_t chisquare[nptbinsnew][4];
          for(Int_t ipt=0; ipt<nptbinsnew; ipt++) {
            for(Int_t ichi=0; ichi<4; ichi++)
              chisquare[ipt][ichi] = 10.;
          }
          
          FillSignalGraph(histlist,gSignal,gSignalfs,gSignalBC1,gSignalBC2,gSigmaFree,gSigmaFixed,gChiSquareFree,gChiSquareFixed,inoutanis,typeb[iBkgFcn],minMassForFit[iMin],maxMassForFit[iMax],rebin[iReb]);
          
          TGraphAsymmErrors *gv2=Computev2(gSignal,resol,averagePt,inoutanis,0x0);
          TGraphAsymmErrors *gv2fs=Computev2(gSignalfs,resol,averagePt,inoutanis,0x0);
          TGraphAsymmErrors *gv2BC1=Computev2(gSignalBC1,resol,averagePt,inoutanis,0x0);
          TGraphAsymmErrors *gv2BC2=Computev2(gSignalBC2,resol,averagePt,inoutanis,0x0);

          //Fill multi-trial histos
          for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
            Double_t chifree[nphi];
            Double_t chifix[nphi];
            Bool_t isChiSquareFreeOk=kTRUE;
            Bool_t isChiSquareFixOk=kTRUE;
            for(Int_t iPhi=0; iPhi<nphi; iPhi++) {

              Double_t Y,sigma,phi;

              gSignal[iPt]->GetPoint(iPhi,phi,Y);
              gSigmaFree[iPt]->GetPoint(iPhi,phi,sigma);
              gChiSquareFree[iPt]->GetPoint(iPhi,phi,chifree[iPhi]);
              hRawYield[iPhi][iPt]->Fill(Y);
              hSigmaFree[iPhi][iPt]->Fill(sigma);
              hChiSquareFree[iPhi][iPt]->Fill(chifree[iPhi]);
              
              gSignalfs[iPt]->GetPoint(iPhi,phi,Y);
              gSigmaFixed[iPt]->GetPoint(iPhi,phi,sigma);
              gChiSquareFixed[iPt]->GetPoint(iPhi,phi,chifix[iPhi]);
              hRawYieldfs[iPhi][iPt]->Fill(Y);
              hSigmaFixed[iPhi][iPt]->Fill(sigma);
              hChiSquareFixed[iPhi][iPt]->Fill(chifix[iPhi]);
              
              gSignalBC1[iPt]->GetPoint(iPhi,phi,Y);
              hRawYieldBC1[iPhi][iPt]->Fill(Y);
              gSignalBC2[iPt]->GetPoint(iPhi,phi,Y);
              hRawYieldBC2[iPhi][iPt]->Fill(Y);
              
              if(chifree[iPhi]>maxchi) {isChiSquareFreeOk=kFALSE;}
              if(chifix[iPhi]>maxchi) {isChiSquareFixOk=kFALSE;}
            }
            
            Double_t v2, v2err, pt, v2ref, phi;
            gv2fsRef->GetPoint(iPt,pt,v2ref);
            
            if(isChiSquareFreeOk) {//fill only if chisquare is good
              gv2->GetPoint(iPt,pt,v2);
              v2err=gv2->GetErrorY(iPt);
              hv2[iPt]->Fill(v2-v2ref);
              hv2VsTrial[iPt]->SetBinContent(iTrial+1,v2);
              hv2VsTrial[iPt]->SetBinError(iTrial+1,v2err);
              gv2BC1->GetPoint(iPt,pt,v2);
              v2err=gv2BC1->GetErrorY(iPt);
              hv2BC1[iPt]->Fill(v2-v2ref);
              hv2BC1VsTrial[iPt]->SetBinContent(iTrial+1,v2);
              hv2BC1VsTrial[iPt]->SetBinError(iTrial+1,v2err);
              gv2BC2->GetPoint(iPt,pt,v2);
              v2err=gv2BC2->GetErrorY(iPt);
              hv2BC2[iPt]->Fill(v2-v2ref);
              hv2BC2VsTrial[iPt]->SetBinContent(iTrial+1,v2);
              hv2BC2VsTrial[iPt]->SetBinError(iTrial+1,v2err);
            }
            if(isChiSquareFixOk) {//fill only if chisquare is good
              gv2fs->GetPoint(iPt,pt,v2);
              v2err=gv2fs->GetErrorY(iPt);
              hv2fs[iPt]->Fill(v2-v2ref);
              hv2fsVsTrial[iPt]->SetBinContent(iTrial+1,v2);
              hv2fsVsTrial[iPt]->SetBinError(iTrial+1,v2err);
            }
          }
          iTrial++;
        }
      }
    }
  }
  
  for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
    hv2Syst[iPt]->SetBinContent(1,hv2[iPt]->GetMean());
    hv2Syst[iPt]->SetBinError(1,hv2[iPt]->GetRMS());
    hv2fsSyst[iPt]->SetBinContent(2,hv2fs[iPt]->GetMean());
    hv2fsSyst[iPt]->SetBinError(2,hv2fs[iPt]->GetRMS());
    hv2BC1Syst[iPt]->SetBinContent(3,hv2BC1[iPt]->GetMean());
    hv2BC1Syst[iPt]->SetBinError(3,hv2BC1[iPt]->GetRMS());
    hv2BC2Syst[iPt]->SetBinContent(4,hv2BC2[iPt]->GetMean());
    hv2BC2Syst[iPt]->SetBinError(4,hv2BC2[iPt]->GetRMS());
    hv2Syst2[iPt]->SetBinContent(1,hv2[iPt]->GetMean());
    Double_t xmin, xmax;
    GetMinMaxHisto(hv2[iPt],xmin,xmax);
    hv2Syst2[iPt]->SetBinError(1,(xmax-xmin)/TMath::Sqrt(12));
    hv2fsSyst2[iPt]->SetBinContent(2,hv2fs[iPt]->GetMean());
    GetMinMaxHisto(hv2fs[iPt],xmin,xmax);
    hv2fsSyst2[iPt]->SetBinError(2,(xmax-xmin)/TMath::Sqrt(12));
    hv2BC1Syst2[iPt]->SetBinContent(3,hv2BC1[iPt]->GetMean());
    GetMinMaxHisto(hv2BC1[iPt],xmin,xmax);
    hv2BC1Syst2[iPt]->SetBinError(3,(xmax-xmin)/TMath::Sqrt(12));
    hv2BC2Syst2[iPt]->SetBinContent(4,hv2BC2[iPt]->GetMean());
    GetMinMaxHisto(hv2BC2[iPt],xmin,xmax);
    hv2BC2Syst2[iPt]->SetBinError(4,(xmax-xmin)/TMath::Sqrt(12));
    
    hSyst1->SetBinContent(iPt+1,TMath::Sqrt(hv2fs[iPt]->GetMean()*hv2fs[iPt]->GetMean()+hv2fs[iPt]->GetRMS()*hv2fs[iPt]->GetRMS()));
  }
  
  //Prepare output file
  TFile *fout=new TFile(Form("%s/v2RawYieldSyst_%d_%d_%s%s%s.root",outputdir.Data(),minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()),"RECREATE");
  
  TPaveStats **pv2 = new TPaveStats*[nptbinsnew];
  TPaveStats **pv2fs = new TPaveStats*[nptbinsnew];
  TPaveStats **pv2BC1 = new TPaveStats*[nptbinsnew];
  TPaveStats **pv2BC2 = new TPaveStats*[nptbinsnew];

  TPaveStats ***pRawYield = new TPaveStats**[nphi];
  TPaveStats ***pRawYieldfs = new TPaveStats**[nphi];
  TPaveStats ***pRawYieldBC1 = new TPaveStats**[nphi];
  TPaveStats ***pRawYieldBC2 = new TPaveStats**[nphi];

  TPaveStats ***pSigmaFree = new TPaveStats**[nphi];
  TPaveStats ***pSigmaFixed = new TPaveStats**[nphi];
  
  TPaveStats ***pChiSquareFree = new TPaveStats**[nphi];
  TPaveStats ***pChiSquareFixed = new TPaveStats**[nphi];
  
  for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
    pRawYield[iPhi] = new TPaveStats*[nptbinsnew];
    pRawYieldfs[iPhi] = new TPaveStats*[nptbinsnew];
    pRawYieldBC1[iPhi] = new TPaveStats*[nptbinsnew];
    pRawYieldBC2[iPhi] = new TPaveStats*[nptbinsnew];
    pSigmaFree[iPhi] = new TPaveStats*[nptbinsnew];
    pSigmaFixed[iPhi] = new TPaveStats*[nptbinsnew];
    pChiSquareFree[iPhi] = new TPaveStats*[nptbinsnew];
    pChiSquareFixed[iPhi] = new TPaveStats*[nptbinsnew];
  }
  
  for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
    cv2->cd(iPt+1);
    hv2[iPt]->GetYaxis()->SetRangeUser(0.,hv2fs[iPt]->GetMaximum()*2.);
    hv2[iPt]->Draw();
    hv2fs[iPt]->Draw("sames");
    hv2BC1[iPt]->Draw("sames");
    hv2BC2[iPt]->Draw("sames");
    cv2->cd(iPt+1)->Update();
    pv2[iPt] = (TPaveStats*)hv2[iPt]->FindObject("stats");
    pv2[iPt]->SetTextColor(colors[0]);
    pv2[iPt]->SetY1NDC(0.74);
    pv2[iPt]->SetY2NDC(0.89);
    pv2fs[iPt] = (TPaveStats*)hv2fs[iPt]->FindObject("stats");
    pv2fs[iPt]->SetTextColor(colors[1]);
    pv2fs[iPt]->SetY1NDC(0.59);
    pv2fs[iPt]->SetY2NDC(0.74);
    pv2BC1[iPt] = (TPaveStats*)hv2BC1[iPt]->FindObject("stats");
    pv2BC1[iPt]->SetTextColor(colors[2]);
    pv2BC1[iPt]->SetY1NDC(0.44);
    pv2BC1[iPt]->SetY2NDC(0.59);
    pv2BC2[iPt] = (TPaveStats*)hv2BC2[iPt]->FindObject("stats");
    pv2BC2[iPt]->SetTextColor(colors[3]);
    pv2BC2[iPt]->SetY1NDC(0.29);
    pv2BC2[iPt]->SetY2NDC(0.44);
    cv2->cd(iPt+1)->Modified();
    
    cv2VsTrial->cd(iPt+1);
    hv2VsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2fsVsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2BC1VsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2BC2VsTrial[iPt]->GetYaxis()->SetRangeUser(-0.3,0.7);
    hv2VsTrial[iPt]->Draw();
    hv2fsVsTrial[iPt]->Draw("same");
    hv2BC1VsTrial[iPt]->Draw("same");
    hv2BC2VsTrial[iPt]->Draw("same");
    
    fout->cd();
    hv2[iPt]->Write();
    hv2fs[iPt]->Write();
    hv2BC1[iPt]->Write();
    hv2BC2[iPt]->Write();
    hv2VsTrial[iPt]->Write();
    hv2fsVsTrial[iPt]->Write();
    hv2BC1VsTrial[iPt]->Write();
    hv2BC2VsTrial[iPt]->Write();
    
    for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
      cRawYield[iPhi]->cd(iPt+1);
      hRawYield[iPhi][iPt]->GetYaxis()->SetRangeUser(0.,hRawYield[iPhi][iPt]->GetMaximum()*2.0);
      hRawYield[iPhi][iPt]->Draw();
      hRawYieldfs[iPhi][iPt]->Draw("sames");
      hRawYieldBC1[iPhi][iPt]->Draw("sames");
      hRawYieldBC2[iPhi][iPt]->Draw("sames");
      cRawYield[iPhi]->cd(iPt+1)->Update();
      pRawYield[iPhi][iPt] = (TPaveStats*)hRawYield[iPhi][iPt]->FindObject("stats");
      pRawYield[iPhi][iPt]->SetTextColor(colors[0]);
      pRawYield[iPhi][iPt]->SetY1NDC(0.74);
      pRawYield[iPhi][iPt]->SetY2NDC(0.89);
      pRawYieldfs[iPhi][iPt] = (TPaveStats*)hRawYieldfs[iPhi][iPt]->FindObject("stats");
      pRawYieldfs[iPhi][iPt]->SetTextColor(colors[1]);
      pRawYieldfs[iPhi][iPt]->SetY1NDC(0.59);
      pRawYieldfs[iPhi][iPt]->SetY2NDC(0.74);
      pRawYieldBC1[iPhi][iPt] = (TPaveStats*)hRawYieldBC1[iPhi][iPt]->FindObject("stats");
      pRawYieldBC1[iPhi][iPt]->SetTextColor(colors[2]);
      pRawYieldBC1[iPhi][iPt]->SetY1NDC(0.44);
      pRawYieldBC1[iPhi][iPt]->SetY2NDC(0.59);
      pRawYieldBC2[iPhi][iPt] = (TPaveStats*)hRawYieldBC2[iPhi][iPt]->FindObject("stats");
      pRawYieldBC2[iPhi][iPt]->SetTextColor(colors[3]);
      pRawYieldBC2[iPhi][iPt]->SetY1NDC(0.29);
      pRawYieldBC2[iPhi][iPt]->SetY2NDC(0.44);
      cRawYield[iPhi]->cd(iPt+1)->Modified();
      
      cSigma[iPhi]->cd(iPt+1);
      hSigmaFree[iPhi][iPt]->GetYaxis()->SetRangeUser(0.,hSigmaFixed[iPhi][iPt]->GetMaximum()*2.0);
      hSigmaFree[iPhi][iPt]->Draw();
      hSigmaFixed[iPhi][iPt]->Draw("sames");
      cSigma[iPhi]->cd(iPt+1)->Update();
      pSigmaFree[iPhi][iPt] = (TPaveStats*)hSigmaFree[iPhi][iPt]->FindObject("stats");
      pSigmaFree[iPhi][iPt]->SetTextColor(colors[0]);
      pSigmaFree[iPhi][iPt]->SetY1NDC(0.74);
      pSigmaFree[iPhi][iPt]->SetY2NDC(0.89);
      pSigmaFixed[iPhi][iPt] = (TPaveStats*)hSigmaFixed[iPhi][iPt]->FindObject("stats");
      pSigmaFixed[iPhi][iPt]->SetTextColor(colors[1]);
      pSigmaFixed[iPhi][iPt]->SetY1NDC(0.59);
      pSigmaFixed[iPhi][iPt]->SetY2NDC(0.74);
      cSigma[iPhi]->cd(iPt+1)->Modified();

      cChiSquare[iPhi]->cd(iPt+1);
      hChiSquareFree[iPhi][iPt]->GetYaxis()->SetRangeUser(0.,hChiSquareFixed[iPhi][iPt]->GetMaximum()*2.0);
      hChiSquareFree[iPhi][iPt]->Draw();
      hChiSquareFixed[iPhi][iPt]->Draw("sames");
      cChiSquare[iPhi]->cd(iPt+1)->Update();
      pChiSquareFree[iPhi][iPt] = (TPaveStats*)hChiSquareFree[iPhi][iPt]->FindObject("stats");
      pChiSquareFree[iPhi][iPt]->SetTextColor(colors[0]);
      pChiSquareFree[iPhi][iPt]->SetY1NDC(0.74);
      pChiSquareFree[iPhi][iPt]->SetY2NDC(0.89);
      pChiSquareFixed[iPhi][iPt] = (TPaveStats*)hChiSquareFixed[iPhi][iPt]->FindObject("stats");
      pChiSquareFixed[iPhi][iPt]->SetTextColor(colors[1]);
      pChiSquareFixed[iPhi][iPt]->SetY1NDC(0.59);
      pChiSquareFixed[iPhi][iPt]->SetY2NDC(0.74);
      cChiSquare[iPhi]->cd(iPt+1)->Modified();

      fout->cd();
      hRawYield[iPhi][iPt]->Write();
      hRawYieldfs[iPhi][iPt]->Write();
      hRawYieldBC1[iPhi][iPt]->Write();
      hRawYieldBC2[iPhi][iPt]->Write();
      hSigmaFree[iPhi][iPt]->Write();
      hSigmaFixed[iPhi][iPt]->Write();
      hChiSquareFree[iPhi][iPt]->Write();
      hChiSquareFixed[iPhi][iPt]->Write();
    }
  }
  
  TCanvas** cv2Syst = new TCanvas*[nptbinsnew];
  TLine** v2refline = new TLine*[nptbinsnew];
  TLine** v2reflinelow = new TLine*[nptbinsnew];
  TLine** v2reflinehigh = new TLine*[nptbinsnew];
  TLegend** leg = new TLegend*[nptbinsnew];
  
  for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
    cv2Syst[iPt] = new TCanvas(Form("cv2Syst_%d",iPt),"",800,800);
    cv2Syst[iPt]->SetLeftMargin(0.16);
    Double_t v2staterr = gv2fsRef->GetErrorY(iPt);
    cv2Syst[iPt]->cd();
    v2refline[iPt] = new TLine(-0.5,0,3.5,0);
    v2refline[iPt]->SetLineWidth(2);
    v2refline[iPt]->SetLineColor(colors[1]);
    v2reflinelow[iPt] = new TLine(-0.5,-v2staterr,3.5,-v2staterr);
    v2reflinelow[iPt]->SetLineWidth(2);
    v2reflinelow[iPt]->SetLineColor(colors[1]);
    v2reflinelow[iPt]->SetLineStyle(7);
    v2reflinehigh[iPt] = new TLine(-0.5,v2staterr,3.5,v2staterr);
    v2reflinehigh[iPt]->SetLineWidth(2);
    v2reflinehigh[iPt]->SetLineColor(colors[1]);
    v2reflinehigh[iPt]->SetLineStyle(7);
    leg[iPt] = new TLegend(0.5,0.7,0.89,0.89);
    leg[iPt]->SetFillStyle(0);
    leg[iPt]->AddEntry(v2reflinelow[iPt],"statistical uncertainty","l");
    leg[iPt]->AddEntry(hv2Syst[iPt],"RMS","p");
    leg[iPt]->AddEntry(hv2Syst2[iPt],"(x_{max}-x_{min})/#sqrt{12}","p");
    hv2Syst[iPt]->GetYaxis()->SetRangeUser(-v2staterr*2,v2staterr*3);
    hv2Syst[iPt]->Draw("E1X0");
    hv2fsSyst[iPt]->Draw("E1sameX0");
    hv2BC1Syst[iPt]->Draw("E1sameX0");
    hv2BC2Syst[iPt]->Draw("E1sameX0");
    hv2Syst2[iPt]->Draw("E1sameX0");
    hv2fsSyst2[iPt]->Draw("E1sameX0");
    hv2BC1Syst2[iPt]->Draw("E1sameX0");
    hv2BC2Syst2[iPt]->Draw("E1sameX0");
    v2refline[iPt]->Draw("same");
    v2reflinelow[iPt]->Draw("same");
    v2reflinehigh[iPt]->Draw("same");
    leg[iPt]->Draw("same");
    fout->cd();
    cv2Syst[iPt]->Write();
    cv2Syst[iPt]->SaveAs(Form("%s/v2RawYieldSyst_pt%.1f-%.1f_%d_%d_%s%s%s.pdf",outputdir.Data(),ptbinsnew[iPt],ptbinsnew[iPt+1],minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()));
  }
  fout->Close();

  TCanvas* cv2SystVsPt = new TCanvas("cv2SystVsPt","",1920,1080);
  cv2SystVsPt->SetLeftMargin(0.15);
  hSyst1->Draw("hist");
  
  cv2->SaveAs(Form("%s/v2RawYieldSyst_%d_%d_%s%s%s.pdf",outputdir.Data(),minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()));
  cv2SystVsPt->SaveAs(Form("%s/v2RawYieldSyst_RMSshift_%d_%d_%s%s%s.pdf",outputdir.Data(),minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()));
  for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
    cRawYield[iPhi]->SaveAs(Form("%s/RawYieldSyst_phi%d_%d_%d_%s%s%s.pdf",outputdir.Data(),iPhi,minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()));
    cSigma[iPhi]->SaveAs(Form("%s/RawYieldSyst_phi%d_%d_%d_%s%s%s.pdf",outputdir.Data(),iPhi,minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()));
    cChiSquare[iPhi]->SaveAs(Form("%s/RawYieldSyst_phi%d_%d_%d_%s%s%s.pdf",outputdir.Data(),iPhi,minCent,maxCent,aniss.Data(),suffix.Data(),fine.Data()));
  }

  return 0;
}

//______________________________________________________________
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
    for(Int_t ires=0;ires<3;ires++){
      resolSub[ires]=hevplresos[ires]->GetMean();
      errors[ires]=hevplresos[ires]->GetMeanError();
    }
    Double_t lowlim[3];for(Int_t ie=0;ie<3;ie++)lowlim[ie]=TMath::Abs(resolSub[ie]-errors[ie]);
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
      delete hMass;
      hMass=0x0;
    }
  }


  if(useTemplD0Refl){
    Bool_t retCode=LoadD0toKpiMCHistos(outlist);
    if(!retCode)Printf("ERROR: MC histograms loading failed");
    else Printf("******************************************\n MC HISTOGRAMS LOADED\n\n**********************************");
  }
  return outlist;
}
//______________________________________________________________
Bool_t DefinePtBins(AliRDHFCuts *cutobj){
  Int_t nPtBinsCuts=cutobj->GetNPtBins();
  Float_t *ptlimsCuts=(Float_t*)cutobj->GetPtBinLimits();
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
void FillSignalGraph(TList *histlist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2,TGraphAsymmErrors **gSigmaFree,TGraphAsymmErrors **gSigmaFixed, TGraphAsymmErrors **gChiSquareFree, TGraphAsymmErrors **gChiSquareFixed, Bool_t inoutanis, Int_t bkgfunc, Double_t minfit, Double_t maxfit, Int_t rebin){

  Int_t nphi=nphibins;
  if(inoutanis)nphi=2;

  TH1F *hrflTempl=0x0;
  TH1F *hsigMC=0x0;
  Float_t sOverRef=0.;
  
  Int_t nMassBins;
  for(Int_t ipt=0;ipt<nptbinsnew;ipt++){
    TH1F *histtofitfullsigma=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi0",ipt))->Clone();
    if(useTemplD0Refl){
      hrflTempl=(TH1F*)histlist->FindObject(Form("histRfl_%d",ipt));
      hsigMC=(TH1F*)histlist->FindObject(Form("histSgn_%d",ipt));
    }
    for(Int_t iphi=0;iphi<nphi;iphi++){
      Double_t signal=0,esignal=0;
      Double_t sigma=0, esigma=0;
      Double_t chisquare=0;
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      if(iphi>0)histtofitfullsigma->Add((TH1F*)histtofit->Clone());
      if(!histtofit){
        gSignal[ipt]->SetPoint(iphi,iphi,signal);
        gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
        return;
      }
      histtofit->SetTitle(Form("%.1f < #it{p}_{T} < %.1f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      nMassBins=histtofit->GetNbinsX();
      histtofit->Rebin(rebin);
      //AliHFMassFitter fitter(histtofit,minfit,maxfit,1,bkgfunc);
       AliHFMassFitterVAR *fitter=new AliHFMassFitterVAR(histtofit,minfit,maxfit,1,bkgfunc,types);
      if(useTemplD0Refl){
	fitter->SetTemplateReflections(hrflTempl);
	sOverRef=(hrflTempl->Integral(hrflTempl->FindBin(minfit*1.0001),hrflTempl->FindBin(maxfit*0.999)))/(hsigMC->Integral(hsigMC->FindBin(minfit*1.0001),hsigMC->FindBin(maxfit*0.999)));
	fitter->SetFixReflOverS(sOverRef,kTRUE);
      }
      fitter->SetInitialGaussianMean(massD);
      fitter->SetInitialGaussianSigma(0.012);
      Bool_t ok=fitter->MassFitter(kFALSE);
      Double_t sigmaforcounting=0;
      Double_t meanforcounting=0;
      if(ok){
        signal=fitter->GetRawYield();
        esignal=fitter->GetRawYieldError();
        sigma=fitter->GetSigma();
        esigma=fitter->GetSigmaUncertainty();
        sigmaforcounting=sigma;
        meanforcounting=fitter->GetMean();
        chisquare=fitter->GetReducedChiSquare();
      }
      gSignal[ipt]->SetPoint(iphi,iphi,signal);
      gSignal[ipt]->SetPointError(iphi,0,0,esignal,esignal);
      gSigmaFree[ipt]->SetPoint(iphi,iphi,sigma);
      gSigmaFree[ipt]->SetPointError(iphi,0,0,esigma,esigma);
      gChiSquareFree[ipt]->SetPoint(iphi,iphi,chisquare);
      gChiSquareFree[ipt]->SetPointError(iphi,0,0,0,0);
      TF1* fB1=fitter->GetBackgroundFullRangeFunc();
      TF1* fB2=fitter->GetBackgroundRecalcFunc();
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
    histtofitfullsigma->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",ptbinsnew[ipt],ptbinsnew[ipt+1]));
    histtofitfullsigma->GetXaxis()->SetTitle("M_{K#pi#pi} (GeV/c^{2})");
    histtofitfullsigma->GetXaxis()->SetTitleSize(0.05);
    nMassBins=histtofitfullsigma->GetNbinsX();
    histtofitfullsigma->Rebin(rebin);
    //AliHFMassFitter fitter(histtofitfullsigma,minfit,maxfit,1,bkgfunc);
    AliHFMassFitterVAR *fitter=new AliHFMassFitterVAR(histtofitfullsigma,minfit,maxfit,1,bkgfunc,types);
    if(useTemplD0Refl){
      fitter->SetTemplateReflections(hrflTempl);
      sOverRef=hrflTempl->Integral(hrflTempl->FindBin(minfit*1.0001),hrflTempl->FindBin(maxfit*0.999))/hsigMC->Integral(hsigMC->FindBin(minfit*1.0001),hsigMC->FindBin(maxfit*0.999));
      fitter->SetFixReflOverS(sOverRef,kTRUE);
    }
    fitter->SetInitialGaussianMean(massD);
    Bool_t ok=fitter->MassFitter(kFALSE);
    Double_t sigmatot=fitter->GetSigma();
    Double_t massFromFit=fitter->GetMean();
    for(Int_t iphi=0;iphi<nphi;iphi++){
      TH1F *histtofit=(TH1F*)histlist->FindObject(Form("hMass_pt%d_phi%d",ipt,iphi))->Clone();
      histtofit->SetTitle(Form("%.1f < #it{p}_{T} < %.1f, #phi%d",ptbinsnew[ipt],ptbinsnew[ipt+1],iphi));
      nMassBins=histtofit->GetNbinsX();
      histtofit->Rebin(rebin);
      //AliHFMassFitter fitter2(histtofit,minfit,maxfit,1,bkgfunc);
      AliHFMassFitterVAR *fitter2=new AliHFMassFitterVAR(histtofit,minfit,maxfit,1,bkgfunc,types);
      if(useTemplD0Refl){
	fitter2->SetTemplateReflections(hrflTempl);
	sOverRef=hrflTempl->Integral(hrflTempl->FindBin(minfit*1.0001),hrflTempl->FindBin(maxfit*0.999))/hsigMC->Integral(hsigMC->FindBin(minfit*1.0001),hsigMC->FindBin(maxfit*0.999));
	fitter2->SetFixReflOverS(sOverRef,kTRUE);
      }
      fitter2->SetInitialGaussianMean(massD);
      fitter2->SetFixGaussianSigma(sigmatot);
      if(fixAlsoMass) fitter2->SetFixGaussianMean(massFromFit);
      Bool_t ok2=fitter2->MassFitter(kFALSE);
      Double_t signal=0,esignal=0;
      Double_t sigma=0, esigma=0;
      Double_t chisquare=0;
      if(ok2){
        signal=fitter2->GetRawYield();
        esignal=fitter2->GetRawYieldError();
        sigma=fitter->GetSigma();
        esigma=fitter->GetSigmaUncertainty();
        chisquare=fitter2->GetReducedChiSquare();
      }
      gSignalfs[ipt]->SetPoint(iphi,iphi,signal);
      gSignalfs[ipt]->SetPointError(iphi,0,0,esignal,esignal);
      gSigmaFixed[ipt]->SetPoint(iphi,iphi,sigma);
      gSigmaFixed[ipt]->SetPointError(iphi,0,0,esigma,esigma);
      gChiSquareFixed[ipt]->SetPoint(iphi,iphi,chisquare);
      gChiSquareFixed[ipt]->SetPointError(iphi,0,0,0,0);
    }
  }//end loop on pt bin
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
      Double_t anis=0;
      Double_t eAnis=0;
      Double_t v2=0.;
      Double_t ev2=0;
      if((nIn+nOut)!=0) {
        anis=(nIn-nOut)/(nIn+nOut);
        eAnis=2./((nIn+nOut)*(nIn+nOut))*TMath::Sqrt(nIn*nIn*enOut*enOut+nOut*nOut*enIn*enIn);
        v2=anis*TMath::Pi()/4./resol;
        ev2=eAnis*TMath::Pi()/4./resol;
      }
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

//___________________________________________________________
Int_t LoadRefGraphs(TString reffilename, TH1F **hRawYieldRef, TH1F **hRawYieldfsRef, TH1F **hRawYieldBC1Ref, TH1F **hRawYieldBC2Ref, TGraphAsymmErrors *&gv2Ref, TGraphAsymmErrors *&gv2fsRef, TGraphAsymmErrors *&gv2BC1Ref, TGraphAsymmErrors *&gv2BC2Ref, Bool_t inoutanis) {
  
  Int_t nphi=nphibins;
  if(inoutanis) nphi=2;
  
  TFile* reffile = TFile::Open(reffilename.Data(),"READ");
  if(reffile) {
    gv2Ref=(TGraphAsymmErrors*)reffile->Get("gav2");
    gv2fsRef=(TGraphAsymmErrors*)reffile->Get("gav2fs");
    gv2BC1Ref=(TGraphAsymmErrors*)reffile->Get("gav2BC1");
    gv2BC2Ref=(TGraphAsymmErrors*)reffile->Get("gav2BC2");
  
    if(!gv2Ref) {return 2;}
    if(!gv2fsRef) {return 3;}
    if(!gv2BC1Ref) {return 4;}
    if(!gv2BC2Ref) {return 5;}
    
    for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
      hRawYieldRef[iPhi] = new TH1F(Form("hRawYieldRef_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldfsRef[iPhi]= new TH1F(Form("hRawYieldfsRef_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldBC1Ref[iPhi]= new TH1F(Form("hRawYieldBC1Ref_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldBC2Ref[iPhi]= new TH1F(Form("hRawYieldBC2Ref_phi%d",iPhi),"",nptbinsnew,ptbinsnew);
      hRawYieldRef[iPhi]->SetDirectory(0);
      hRawYieldfsRef[iPhi]->SetDirectory(0);
      hRawYieldBC1Ref[iPhi]->SetDirectory(0);
      hRawYieldBC2Ref[iPhi]->SetDirectory(0);
    }
    
    for(Int_t iPt=0; iPt<nptbinsnew; iPt++) {
      TGraphAsymmErrors* gtmp = (TGraphAsymmErrors*)reffile->Get(Form("gasigpt%d",iPt));
      TGraphAsymmErrors* gtmpfs = (TGraphAsymmErrors*)reffile->Get(Form("gasigfspt%d",iPt));
      TGraphAsymmErrors* gtmpBC1 = (TGraphAsymmErrors*)reffile->Get(Form("gasigBC1pt%d",iPt));
      TGraphAsymmErrors* gtmpBC2 = (TGraphAsymmErrors*)reffile->Get(Form("gasigBC2pt%d",iPt));
      
      if(gtmp && gtmpfs && gtmpBC1 && gtmpBC2) {
        for(Int_t iPhi=0; iPhi<nphi; iPhi++) {
          Double_t phi,Y, errY;
          gtmp->GetPoint(iPhi,phi,Y);
          hRawYieldRef[iPhi]->SetBinContent(iPt+1,Y);
          gtmpfs->GetPoint(iPhi,phi,Y);
          hRawYieldfsRef[iPhi]->SetBinContent(iPt+1,Y);
          gtmpBC1->GetPoint(iPhi,phi,Y);
          hRawYieldBC1Ref[iPhi]->SetBinContent(iPt+1,Y);
          gtmpBC2->GetPoint(iPhi,phi,Y);
          hRawYieldBC2Ref[iPhi]->SetBinContent(iPt+1,Y);
        }
      }
    }
    reffile->Close();
  }
  else {return 1;}
  
  return 0;
}

//___________________________________________________________
void GetMinMaxHisto(TH1F* histo, Double_t &xmin, Double_t &xmax) {

  Int_t nbins=histo->GetNbinsX();
  xmin=histo->GetBinLowEdge(nbins)+histo->GetBinWidth(nbins);
  Int_t content=0;
  for(Int_t iBin=0; iBin<nbins; iBin++) {
    content=histo->GetBinContent(iBin+1);
    if(content>0) {
      xmin=histo->GetBinLowEdge(iBin+1);
      break;
    }
  }
  content=0;
  xmax=histo->GetBinLowEdge(1);
  for(Int_t iBin=nbins; iBin>0; iBin--) {
    content=histo->GetBinContent(iBin);
    if(content>0) {
      xmax=histo->GetBinLowEdge(iBin)+histo->GetBinWidth(iBin);
      break;
    }
  }
}

//__________________________________________________________________________________________________________________
void DivideCanvas(TCanvas* c, const Int_t nPtBins) {
  if(nPtBins<2)
    c->Divide(1,1);
  if(nPtBins==2 || nPtBins==3)
    c->Divide(nPtBins,1);
  else if(nPtBins==4 || nPtBins==6 || nPtBins==8)
    c->Divide(nPtBins/2,2);
  else if(nPtBins==5 || nPtBins==7)
    c->Divide((nPtBins+1)/2,2);
  else if(nPtBins==9 || nPtBins==12 || nPtBins==15)
    c->Divide(nPtBins/3,3);
  else if(nPtBins==10 || nPtBins==11)
    c->Divide(4,3);
  else if(nPtBins==13 || nPtBins==14)
    c->Divide(5,3);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4==0)
    c->Divide(nPtBins/4,4);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4!=0)
    c->Divide(5,4);
  else if(nPtBins==21)
    c->Divide(7,3);
  else if(nPtBins>21 && nPtBins<=25)
    c->Divide(5,5);
  else if(nPtBins>25 && nPtBins%2==0)
    c->Divide(nPtBins/2,2);
  else
    c->Divide((nPtBins+1)/2,2);
}

//___________________________________________________________
void SetStyle(Int_t optfit) {
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(optfit);
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetStatFont(42);
  gStyle->SetStatY(0.89);
  gStyle->SetStatX(0.89);
  gStyle->SetTitleFont(42,"xyzg");
  gStyle->SetHistLineWidth(1);
  gStyle->SetLegendBorderSize(0);
}
Bool_t LoadD0toKpiMCHistos(TList *outlist){
  
  TFile *f=new TFile(fileNameMCD0refl.Data(),"READ");
  if(!f){
    printf("ERROR: file %s does not exist\n",fileNameMCD0refl.Data());
    return kFALSE;
  }
  f->ls();
  TH1F** hsig=new TH1F*[nptbinsnew];
  TH1F** hrfl=new TH1F*[nptbinsnew];
  for(Int_t j=0;j<nptbinsnew;j++){
    hsig[j]=(TH1F*)f->Get(Form("histSgn_%d",j));
    if(!hsig[j]) {Printf("histSgn_%d NOT FOUND",j); return kFALSE;}
    hrfl[j]=(TH1F*)f->Get(Form("histRflFitted%s_ptBin%d",rflFitType.Data(),j));
    if(!hrfl[j]) {Printf("histRflFitted%s_ptBin%d",rflFitType.Data(),j); return kFALSE;}
  }
  for(Int_t k=0;k<nptbinsnew;k++){
    outlist->Add(hsig[k]->Clone(Form("histSgn_%d",k)));
    outlist->Add(hrfl[k]->Clone(Form("histRfl_%d",k)));
  }
  outlist->ls();
  //for(Int_t p=0;p<nptbinsnew;p++){
  //delete hsig[p];
  //delete hrfl[p];
  //}
  //f->Close();
  return kTRUE;
}
