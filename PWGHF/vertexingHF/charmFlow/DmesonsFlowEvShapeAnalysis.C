#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TDatabasePDG.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <vector>

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

//methods for the analysis of AliAnalysisTaskSEHFv2 output in case of kEvShape method
//Author: Fabrizio Grosa, INFN Turin grosa@to.infn.it

//*******************************************************//
//                                                       //
//      Main Function: DmesonsFlowEvShapeAnalysis()      //
//                                                       //
//*******************************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIABLES
//to be set
//input file name
const TString infilename="$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/AnalysisResults_v2_EvShape_TwoEtaHalves.root";
const TString suffix="_Topod0Cut_QoverM_TwoEtaHalves_VZERO_EvShape";
//const TString infilename="$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/AnalysisResults_v2_EvShape_NonFlowTests_EP_ResoNonFlow.root";
//const TString suffix="_Topod0Cut_QoverM_DPosEta_q2NegTPC_VZERO_EvShape";
//const TString infilename = "$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/AnalysisResults_v2_3050_step2_EvShape_VZERO_CentAxis.root";
//const TString suffix = "_Topod0Cut_QoverM_VZERO_EvShape";
//const TString infilename = "$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/AnalysisResults_v2_3050_EvShape_q2TPCSmearing_2.root";
//const TString suffix = "_Topod0Cut_QoverM_q2TPC_VZERO_EvShape";
const TString partname="Dplus";
const Int_t minCent=30;
const Int_t maxCent=50;

const TString outputdir="Cent3050/v2/EvShape/q2Smearing";

//ptbins of the analysis
const Int_t nPtBins=5;
const Int_t nPtLims=nPtBins+1;
const Double_t PtLims[nPtLims] = {3.,4.,6.,8.,12.,16.};

//phi bins
const Int_t nPhiBins=4;
const Int_t nPhiLims=nPhiBins+1;
const Double_t PhiLims[nPhiLims]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};

//q2 cut values (absolute cut)
const Double_t q2smalllimit=2.2;
const Double_t q2largelimit=3.2;

//percentage of events with smaller/larger q2 (both for centrality integrated and centrality dependent cut)
const Double_t q2smallpercevents=0.60;
const Double_t q2largepercevents=0.20;

//EP resolution
const Int_t evPlane=AliEventPlaneResolutionHandler::kVZERO;
Int_t evPlaneRes=AliEventPlaneResolutionHandler::kThreeSub;
const Bool_t useNcollWeight=kFALSE;
const Bool_t useAliHandlerForRes=kFALSE;

// mass fit configuration
const Int_t rebin[]={3,4,4,6,6};
enum {kGaus=0, kDoubleGaus, kReflTempl};
const Int_t types=kGaus;//kReflTempl;
const Int_t typeb=AliHFMassFitter::kExpo; //Background: 0=expo, 1=linear, 2=pol2
Bool_t useTemplD0Refl=kFALSE;
TString rflFitType="DoubleGaus";
TString fileNameMCD0refl="./reflections/reflections_fitted_DoubleGaus.root";
const Bool_t fixAlsoMass=kFALSE;
const Double_t minMassForFit[]={1.72,1.70,1.70,1.70,1.70};
const Double_t maxMassForFit[]={2.02,2.05,2.05,2.15,2.15};
const Double_t nSigmaForCounting=3.5;

//not to be set
enum CutMethod{kAbsCut,kPercCut,kPercCutVsCent}; //kAbsCut->absolute cut values, kPercCut->cut according to the % of events with smaller/larger q2, kPercCutVsCent->cut according to the % of events with smaller/larger q2 in finer centrality bins
enum SmallOrLarge{kSmall,kLarge,kIntegrated};
enum AnalysisMethod{kEventPlane,kEventPlaneInOut,kScalarProd};

//in-out efficiency
const Double_t effInOverEffOut=1.03;

const Int_t colors[]={kRed,kBlue,kGray};
const Int_t markers[]={kFullSquare,kFullCircle,kFullDiamond,kOpenSquare,kOpenCircle,kOpenDiamond};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
Int_t DmesonsFlowEvShapeAnalysis(Int_t cutmeth=kPercCutVsCent, Int_t analysismeth=kEventPlaneInOut);
void Drawq2VsCent(Int_t cutmeth=kPercCutVsCent);
void DrawEventPlaneResolutionAndDistribution(Int_t cutmeth=kPercCutVsCent);
TList* LoadTList(Int_t &ncentbins);
THnSparseF* LoadSparseFromList(TList* inputlist);
TH2F* GetHistoq2VsCentr(TList* inputlist);
TList* LoadMassHistos(THnSparseF* sparse, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge, Int_t analysismeth);
TList* LoadResolutionHistos(TList *inputlist, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge);
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Int_t analysismeth, TGraphAsymmErrors *gRelSystEff);
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, Int_t smallorlarge, Int_t analysismeth, Int_t cutmeth);
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3);
void Applyq2Cut(THnSparseF* sparse, Double_t smallcutvalues, Double_t largecutvalues, Int_t smallorlarge);
void ApplyCut(THnSparseF* sparse, Double_t min, Double_t max, UInt_t axnum);
Bool_t Defineq2Cuts(TH2F* hq2VsCentr,vector<Double_t> &smallcutvalues,vector<Double_t> &largecutvalues, Int_t cutmeth, Int_t fSparseVers);
Bool_t Getq2CutValuePercEvents(TH1F* hq2, Double_t cutvalues[2]);
void ResetAxes(THnSparseF* sparse, Int_t axnum=-1);
void SetStyle();
Bool_t LoadD0toKpiMCHistos(TList *outlist);

//_____________________________________________________________________________________________
//ANALYSIS FUNCTION
Int_t DmesonsFlowEvShapeAnalysis(Int_t cutmeth, Int_t analysismeth) {

  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}

  Int_t ncentbins=0;
  TList* datalist=(TList*)LoadTList(ncentbins);
  if(!datalist){return 1;}
  THnSparseF* hMassPtPhiq2Centr=(THnSparseF*)LoadSparseFromList(datalist);
  if(!hMassPtPhiq2Centr) {return 2;}
  Int_t fSparseVers = 0; //check version of the task output
  if(hMassPtPhiq2Centr->GetAxis(4)) {
    TString centtitle=hMassPtPhiq2Centr->GetAxis(4)->GetTitle();
    if(centtitle.Contains("Centrality") || centtitle.Contains("centrality")) {
      cout << "Version of the sparse with centrality axis!" << endl;
      fSparseVers=1;
    }
    else {
      cout << "Version of the sparse without centrality axis!" << endl;
    }
  }
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) {return 3;}

  vector<Double_t> smallcutvalues;
  vector<Double_t> largecutvalues;
  Bool_t defq2cut=Defineq2Cuts(hq2VsCentr,smallcutvalues,largecutvalues,cutmeth,fSparseVers);
  if(!defq2cut) {return 4;}

  TString analysismethname="";
  if(analysismeth==kEventPlane) {analysismethname="EP";}
  else if(analysismeth==kEventPlaneInOut) {analysismethname="InOut";}
  else if(analysismeth==kScalarProd) {analysismethname="ScalarProduct";}

  Int_t q2region[3] = {kSmall,kLarge,kIntegrated};
  TString q2regionname[3] = {"q2Small","q2Large","q2Int"};

  SetStyle();

  //graphs for v2 in the two q2 regions and q2-integrated
  TCanvas *cv2fs =new TCanvas("cv2fs","v2 Fit with fixed sigma",1920,1080);
  TCanvas **cv2 =new TCanvas*[3];
  TGraphAsymmErrors **gv2=new TGraphAsymmErrors*[3];
  TGraphAsymmErrors **gv2fs=new TGraphAsymmErrors*[3];
  TGraphAsymmErrors **gv2BC1=new TGraphAsymmErrors*[3];
  TGraphAsymmErrors **gv2BC2=new TGraphAsymmErrors*[3];
  TGraphAsymmErrors **gRelSystEff=new TGraphAsymmErrors*[3];
  Double_t ymin[3]={1.,1.,1.};
  Double_t ymax[3]={-1.,-1.,-1.};

  TLegend* leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);

  TLegend** legsyst = new TLegend*[3];

  TString outname=Form("%s/v2Output_%d_%d_%s%s_q2Small%.2f_q2Large%.2f.root",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2smalllimit,q2largelimit);
  TString percsuffix="perc";
  if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/v2Output_%d_%d_%s%s_q2Small%0.f%s_q2Large%0.f%s.root",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
  TFile outfile(outname.Data(),"RECREATE"); //outputfile

  Double_t q2cut[3] = {q2smalllimit,q2largelimit,0};
  Double_t q2perccut[3] = {q2smallpercevents*100,q2largepercevents*100,100};

  for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
    if(analysismeth!=kScalarProd) {
      //load resolution and mass lists
      TList* masslist=(TList*)LoadMassHistos(hMassPtPhiq2Centr,smallcutvalues,largecutvalues,q2region[iq2],analysismeth);

      outname=Form("%s/v2Histograms_%d_%d_%s%s_%s%.2f.root",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[iq2].Data(),q2cut[iq2]);
      if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/v2Histograms_%d_%d_%s%s_%s%0.f%s.root",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[iq2].Data(),q2perccut[iq2],percsuffix.Data());}
      masslist->SaveAs(outname.Data(),"RECREATE");

      cout <<"Average pt for pt bin \n"<<endl;
      //average pt for pt bin
      AliVertexingHFUtils *utils=new AliVertexingHFUtils();
      Int_t minCentTimesTen=minCent*10;
      Int_t maxCentTimesTen=maxCent*10;
      TH2F* hmasspt=0x0;
      for(Int_t iCent=0; iCent<ncentbins; iCent++) {
        if(iq2!=kIntegrated) {Applyq2Cut(hMassPtPhiq2Centr,smallcutvalues[iCent],largecutvalues[iCent],q2region[iq2]);}
        TH2F* htmp=(TH2F*)hMassPtPhiq2Centr->Projection(0,1);
        if(iCent==0) {hmasspt=(TH2F*)htmp->Clone();}
        else {hmasspt->Add(htmp);}
        ResetAxes(hMassPtPhiq2Centr,-1.);
      }
      Float_t averagePt[nPtBins];
      Float_t errorPt[nPtBins];
      for(Int_t iPt=0;iPt<nPtBins;iPt++){
        Int_t binMin=hmasspt->FindBin(PtLims[iPt]);
        Int_t binMax=hmasspt->FindBin(PtLims[iPt+1]-0.001);
        if(TMath::Abs(hmasspt->GetXaxis()->GetBinLowEdge(binMin)-PtLims[iPt])>0.001 ||
           TMath::Abs(hmasspt->GetXaxis()->GetBinUpEdge(binMax)-PtLims[iPt+1])>0.001){
          cout << "Error in pt bin limits for projection!\n" << endl;
          return 4;
        }
        TH1F *histtofit = (TH1F*)hmasspt->ProjectionY("_py",binMin,binMax);
        Int_t nMassBins=histtofit->GetNbinsX();
        Double_t hmin=histtofit->GetBinLowEdge(2); // need wide range for <pt>
        Double_t hmax=histtofit->GetBinLowEdge(nMassBins-2); // need wide range for <pt>
	if (partname.Contains("Dstar")) {
	  if (hmin < 0.140) hmin=0.140;
	  if (hmax > 0.175) hmax=0.175;
	}
        AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(histtofit,hmin,hmax,1,typeb,types);
        if(useTemplD0Refl){
          Printf("USE TEMPLATE FOR AVERAGE Pt");
          TH1F *hrflTempl=(TH1F*)(masslist->FindObject(Form("histRfl_%d",iPt)))->Clone(Form("histrfl_%d",iPt));
          if(!hrflTempl) {Printf("histRfl_%d not found",iPt); return 15;}
          TH1F *hsigMC=(TH1F*)(masslist->FindObject(Form("histSgn_%d",iPt)))->Clone(Form("histsgn_%d",iPt));
          if(!hsigMC) {Printf("histSgn_%d not found",iPt); return 15;}
          fitter->SetTemplateReflections(hrflTempl);
          Float_t sOverRef=(hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999)))/(hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999)));
          Printf("R OVER S = %f",sOverRef);
          fitter->SetFixReflOverS(sOverRef,kTRUE);
        }
        fitter->SetUseLikelihoodFit();
        fitter->MassFitter(kFALSE);
        Double_t massFromFit=fitter->GetMean();
        Double_t sigmaFromFit=fitter->GetSigma();
        TF1* funcB2=fitter->GetBackgroundRecalcFunc();
        utils->AveragePt(averagePt[iPt],errorPt[iPt],PtLims[iPt],PtLims[iPt+1],hmasspt,massFromFit,sigmaFromFit,funcB2,2.5,4.5,0.,3.,1);
        if(averagePt[iPt]==0) averagePt[iPt] = (PtLims[iPt+1]+PtLims[iPt])/2;
      }
      cout << Form("Average pt %s region \n",q2regionname[iq2].Data()) << endl;
      for(Int_t iPt=0;iPt<nPtBins;iPt++) cout <<Form("%f +- %f\n",averagePt[iPt],errorPt[iPt])<<endl;

      Int_t nPhi=nPhiBins;
      if(analysismeth==kEventPlaneInOut) nPhi=2;

      printf("Fill TGraphs for signal \n");
      //Fill TGraphs for signal
      TGraphAsymmErrors *gSignal[nPtBins];
      TGraphAsymmErrors *gSignalfs[nPtBins];
      TGraphAsymmErrors *gSignalBC1[nPtBins];
      TGraphAsymmErrors *gSignalBC2[nPtBins];
      for(Int_t iPt=0;iPt<nPtBins;iPt++){
        gSignal[iPt]=new TGraphAsymmErrors(nPhi);
        gSignal[iPt]->SetName(Form("gasigpt%d_%s",iPt,q2regionname[iq2].Data()));
        gSignal[iPt]->SetTitle(Form("Signal %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",PtLims[iPt],PtLims[iPt+1]));
        gSignal[iPt]->SetMarkerStyle(25);
        gSignalfs[iPt]=new TGraphAsymmErrors(nPhi);
        gSignalfs[iPt]->SetName(Form("gasigfspt%d_%s",iPt,q2regionname[iq2].Data()));
        gSignalfs[iPt]->SetTitle(Form("Signal (fixed sigma) %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",PtLims[iPt],PtLims[iPt+1]));
        gSignalfs[iPt]->SetMarkerStyle(21);
        gSignalBC1[iPt]=new TGraphAsymmErrors(nPhi);
        gSignalBC1[iPt]->SetName(Form("gasigBC1pt%d_%s",iPt,q2regionname[iq2].Data()));
        gSignalBC1[iPt]->SetTitle(Form("Signal (BC1) %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",PtLims[iPt],PtLims[iPt+1]));
        gSignalBC2[iPt]=new TGraphAsymmErrors(nPhi);
        gSignalBC2[iPt]->SetName(Form("gasigBC2pt%d_%s",iPt,q2regionname[iq2].Data()));
        gSignalBC2[iPt]->SetTitle(Form("Signal (BC2) %.1f<#it{p}_{T}<%.1f GeV/c;#Delta#phi bin;Counts",PtLims[iPt],PtLims[iPt+1]));
      }
      FillSignalGraph(masslist,gSignal,gSignalfs,gSignalBC1,gSignalBC2,q2region[iq2],analysismeth,cutmeth);
      outfile.cd();
      for(Int_t iPt=0;iPt<nPtBins;iPt++){
        gSignal[iPt]->Write();
        gSignalfs[iPt]->Write();
        gSignalBC1[iPt]->Write();
        gSignalBC2[iPt]->Write();
      }

      //EP resolution
      Double_t resol=-1.;
      Double_t errorres=-1.;

      if(useAliHandlerForRes) {
        cerr << "Error: AliHandler for resolution not yet implemented. Exit" << endl;
        return 6;
      }
      else {
        TList* resolist=(TList*)LoadResolutionHistos(datalist,smallcutvalues,largecutvalues,q2region[iq2]);
        TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
        TH1F* hevplresos[3];
        TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
        Int_t nSubRes=1;
        TH2F* htestversion=(TH2F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[0].Data(),suffixcentr.Data()));
        if(htestversion){
          printf("Old version of the task\n");
        }else{
          printf("New version of the task\n");
          namereso[0]="Reso1Vsq2";
        }
        if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
           evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
        for(Int_t iRes=0;iRes<nSubRes;iRes++){
          hevplresos[iRes]=(TH1F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
        }
        resol=GetEventPlaneResolution(errorres,hevplresos[0],hevplresos[1],hevplresos[2]);
      }

      cout << "Event plane resolution: " << resol << " +/- " << errorres << endl;
      cout << "Compute v2" << endl;
      gRelSystEff[iq2] = new TGraphAsymmErrors(nPtBins);

      gv2[iq2] = (TGraphAsymmErrors*)Computev2(gSignal,resol,averagePt,analysismeth,0x0);
      gv2fs[iq2] = (TGraphAsymmErrors*)Computev2(gSignalfs,resol,averagePt,analysismeth,gRelSystEff[iq2]);
      gv2BC1[iq2] = (TGraphAsymmErrors*)Computev2(gSignalBC1,resol,averagePt,analysismeth,0x0);
      gv2BC2[iq2] = (TGraphAsymmErrors*)Computev2(gSignalBC2,resol,averagePt,analysismeth,0x0);

      gv2[iq2]->SetName(Form("gav2_%s",q2regionname[iq2].Data()));
      gv2[iq2]->SetTitle("");
      gv2[iq2]->GetYaxis()->SetTitle("v_{2} {EP}");
      gv2[iq2]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      gv2fs[iq2]->SetName(Form("gav2fs_%s",q2regionname[iq2].Data()));
      gv2fs[iq2]->SetTitle("");
      gv2fs[iq2]->GetYaxis()->SetTitle("v_{2} {EP}");
      gv2fs[iq2]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      gv2BC1[iq2]->SetName(Form("gav2BC1_%s",q2regionname[iq2].Data()));
      gv2BC1[iq2]->SetTitle("");
      gv2BC1[iq2]->GetYaxis()->SetTitle("v_{2} {EP}");
      gv2BC1[iq2]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      gv2BC2[iq2]->SetName(Form("gav2BC2_%s",q2regionname[iq2].Data()));
      gv2BC2[iq2]->SetTitle("");
      gv2BC2[iq2]->GetYaxis()->SetTitle("v_{2} {EP}");
      gv2BC2[iq2]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      gRelSystEff[iq2]->SetName(Form("gRelSystEff_%s",q2regionname[iq2].Data()));
      gRelSystEff[iq2]->SetTitle("SystErrEff");
      gv2[iq2]->SetMarkerStyle(markers[0]);
      gv2[iq2]->SetMarkerSize(1.5);
      gv2[iq2]->SetMarkerColor(colors[iq2]);
      gv2[iq2]->SetLineColor(colors[iq2]);
      gv2[iq2]->SetLineWidth(2);
      gv2fs[iq2]->SetMarkerStyle(markers[1]);
      gv2fs[iq2]->SetMarkerSize(1.5);
      gv2fs[iq2]->SetMarkerColor(colors[iq2]+1);
      gv2fs[iq2]->SetLineColor(colors[iq2]+1);
      gv2fs[iq2]->SetLineWidth(2);
      gv2BC1[iq2]->SetMarkerStyle(markers[3]);
      gv2BC1[iq2]->SetMarkerSize(1.5);
      gv2BC1[iq2]->SetMarkerColor(colors[iq2]+2);
      gv2BC1[iq2]->SetLineColor(colors[iq2]+2);
      gv2BC1[iq2]->SetLineWidth(2);
      gv2BC2[iq2]->SetMarkerStyle(markers[4]);
      gv2BC2[iq2]->SetMarkerSize(1.5);
      gv2BC2[iq2]->SetMarkerColor(colors[iq2]+3);
      gv2BC2[iq2]->SetLineColor(colors[iq2]+3);
      gv2BC2[iq2]->SetLineWidth(2);
      gRelSystEff[iq2]->SetMarkerStyle(markers[0]);
      gRelSystEff[iq2]->SetMarkerSize(1.5);
      gRelSystEff[iq2]->SetMarkerColor(colors[0]);
      gRelSystEff[iq2]->SetLineColor(colors[iq2]+4);
      gRelSystEff[iq2]->SetLineWidth(2);

      TString title="";
      if(iq2!=kIntegrated) {
        if(cutmeth==kAbsCut) {
          TString sign[2] = {"<",">"};
          title=Form("%s %s %0.1f",hq2VsCentr->GetYaxis()->GetTitle(),sign[iq2].Data(),q2cut[iq2]);
        }
        else {
          Double_t perc[2] = {q2smallpercevents*100,q2largepercevents*100};
          TString reg[2] = {"small-","large-"};
          title=Form("%0.f%% %s %s",perc[iq2],reg[iq2].Data(),hq2VsCentr->GetYaxis()->GetTitle());
        }
      }
      else {
        title=Form("%s-integrated",hq2VsCentr->GetYaxis()->GetTitle());
      }
      leg->AddEntry(gv2fs[iq2],title.Data(),"lpe");

      legsyst[iq2] = new TLegend(0.6,0.6,0.85,0.85);
      legsyst[iq2]->SetFillStyle(0);
      legsyst[iq2]->SetTextSize(0.05);
      legsyst[iq2]->AddEntry(gv2[iq2],"free sigma","lpe");
      legsyst[iq2]->AddEntry(gv2fs[iq2],"fixed sigma","lpe");
      legsyst[iq2]->AddEntry(gv2BC1[iq2],"bin counting 1","lpe");
      legsyst[iq2]->AddEntry(gv2BC2[iq2],"bin counting 2","lpe");

      cv2[iq2] = new TCanvas(Form("cv2_%s",q2regionname[iq2].Data()),"v2 - systematic on yield extraction",1920,1080);
      cv2[iq2]->cd();
      gv2[iq2]->SetTitle(title.Data());
      gv2[iq2]->Draw("AP");
      gv2fs[iq2]->Draw("P");
      gv2BC1[iq2]->Draw("P");
      gv2BC2[iq2]->Draw("P");
      legsyst[iq2]->Draw("same");

      Double_t tmp=TMath::MinElement(gv2[iq2]->GetN(),gv2[iq2]->GetY());
      if(ymin[iq2]>tmp-1.5*TMath::Abs(tmp)) ymin[iq2]=tmp-1.5*TMath::Abs(tmp);
      tmp=TMath::MinElement(gv2fs[iq2]->GetN(),gv2fs[iq2]->GetY());
      if(ymin[iq2]>tmp-1.5*TMath::Abs(tmp)) ymin[iq2]=tmp-1.5*TMath::Abs(tmp);
      tmp=TMath::MinElement(gv2BC1[iq2]->GetN(),gv2BC1[iq2]->GetY());
      if(ymin[iq2]>tmp-1.5*TMath::Abs(tmp)) ymin[iq2]=tmp-1.5*TMath::Abs(tmp);
      tmp=TMath::MinElement(gv2BC2[iq2]->GetN(),gv2BC2[iq2]->GetY());
      if(ymin[iq2]>tmp-1.5*TMath::Abs(tmp)) ymin[iq2]=tmp-1.5*TMath::Abs(tmp);

      tmp=TMath::MaxElement(gv2[iq2]->GetN(),gv2[iq2]->GetY());
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.2*TMath::Abs(tmp);
      tmp=TMath::MaxElement(gv2fs[iq2]->GetN(),gv2fs[iq2]->GetY());
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.2*TMath::Abs(tmp);
      tmp=TMath::MaxElement(gv2BC1[iq2]->GetN(),gv2BC1[iq2]->GetY());
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.2*TMath::Abs(tmp);
      tmp=TMath::MaxElement(gv2BC2[iq2]->GetN(),gv2BC2[iq2]->GetY());
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.2*TMath::Abs(tmp);

      gv2[iq2]->GetYaxis()->SetRangeUser(ymin[iq2],ymax[iq2]);

      outname=Form("%s/v2Output_%d_%d_%s%s_%s%.2f.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[iq2].Data(),q2cut[iq2]);
      if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/v2Output_%d_%d_%s%s_%s%.f%s.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[iq2].Data(),q2perccut[iq2],percsuffix.Data());}
      cv2[iq2]->SaveAs(outname.Data());

      cv2fs->cd();
      TString drawopt="P";
      if(iq2==kSmall) drawopt+="A";
      gv2fs[iq2]->Draw(drawopt.Data());
      if(iq2==kIntegrated) leg->Draw("same");
    }
    else {
      cerr << "Scalar product not yet implemented. Exit" << endl;
    }
    outfile.cd();
    gv2[iq2]->Write();
    gv2fs[iq2]->Write();
    gv2BC1[iq2]->Write();
    gv2BC2[iq2]->Write();
    gRelSystEff[iq2]->Write();
  }//loop over the two q2 regions

  cv2fs->Write();
  cv2[0]->Write();
  cv2[1]->Write();
  outfile.Close();

  Double_t min=ymin[kSmall];
  if(min>ymin[kLarge]) min=ymin[kLarge];
  if(min>ymin[kIntegrated]) min=ymin[kIntegrated];
  Double_t max=ymax[kLarge];
  if(max<ymax[kSmall]) max=ymax[kSmall];
  if(max<ymax[kIntegrated]) max=ymax[kIntegrated];
  gv2fs[kSmall]->GetYaxis()->SetRangeUser(min,max);


  outname=Form("%s/v2fsOutput_%d_%d_%s%s_q2Small%.2f_q2Large%.2f.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2smalllimit,q2largelimit);
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/v2fsOutput_%d_%d_%s%s_q2Small%.0f%s_q2Large%.0f%s.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
  cv2fs->SaveAs(outname.Data());

  return 0;
}

//_____________________________________________________________________________________________
//DRAW q_2 VS CENTRALITY
void Drawq2VsCent(Int_t cutmeth) {

  TGaxis::SetMaxDigits(3);

  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}

  SetStyle();

  Int_t ncentbins=0;
  TList* datalist=(TList*)LoadTList(ncentbins);
  if(!datalist){return;}
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) return;

  Double_t centwidth=(Double_t)(maxCent-minCent)/ncentbins;
  vector<Double_t> smallcutvalues;
  vector<Double_t> largecutvalues;
  Bool_t defq2cut=Defineq2Cuts(hq2VsCentr,smallcutvalues,largecutvalues,cutmeth,1);
  if(!defq2cut) return;

  if(ncentbins!=(Int_t)smallcutvalues.size() || ncentbins!=(Int_t)largecutvalues.size()) {
    cerr << "Number of centrality bins not consistent. Exit." << endl;
    return;
  }

  TH1F* hq2=(TH1F*)hq2VsCentr->ProjectionY();
  TString q2name = hq2->GetXaxis()->GetTitle();
  hq2->SetTitle("");
  hq2->GetYaxis()->SetTitle(Form("dN_{ev}/d%s",hq2->GetXaxis()->GetTitle()));
  hq2->SetLineColor(kBlack);

  TH1F** hq2centbin = new TH1F*[ncentbins];

  TH1F* hSmallLimitVsCent = new TH1F("hSmallLimitVsCent","",ncentbins,minCent,maxCent);
  TH1F* hLargeLimitVsCent = new TH1F("hLargeLimitVsCent","",ncentbins,minCent,maxCent);
  hSmallLimitVsCent->SetLineStyle(7);
  hLargeLimitVsCent->SetLineStyle(9);
  hSmallLimitVsCent->SetLineWidth(3);
  hLargeLimitVsCent->SetLineWidth(3);

  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    hq2centbin[iCent]=(TH1F*)hq2VsCentr->ProjectionY(Form("hq2centbin%d",iCent),iCent+1,iCent+1);
    hq2centbin[iCent]->GetYaxis()->SetTitle(Form("dN_{ev}/d%s",hq2->GetXaxis()->GetTitle()));
    hSmallLimitVsCent->SetBinContent(iCent+1,smallcutvalues[iCent]);
    hLargeLimitVsCent->SetBinContent(iCent+1,largecutvalues[iCent]);
  }

  TLegend* leg = new TLegend(0.55,0.6,0.89,0.89);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);

  TLegend* legVsCentSmall = new TLegend(0.25,0.75,0.55,0.89);
  legVsCentSmall->SetTextSize(0.045);
  legVsCentSmall->SetFillStyle(0);
  TLegend* legVsCentLarge = new TLegend(0.55,0.75,0.85,0.89);
  legVsCentLarge->SetTextSize(0.045);
  legVsCentSmall->SetFillStyle(0);

  TCanvas *cq2VsCent = new TCanvas("cq2VsCent","q_{2} vs. centrality",800,800);
  cq2VsCent->SetLogz();
  hSmallLimitVsCent->SetLineColor(kBlack);
  hLargeLimitVsCent->SetLineColor(kBlack);
  legVsCentSmall->AddEntry(hSmallLimitVsCent,Form("Small-%s",q2name.Data()),"l");
  legVsCentLarge->AddEntry(hLargeLimitVsCent,Form("Large-%s",q2name.Data()),"l");
  hq2VsCentr->Draw("colz");
  hSmallLimitVsCent->Draw("same");
  hLargeLimitVsCent->Draw("same");
  legVsCentSmall->Draw("same");
  legVsCentLarge->Draw("same");

  TCanvas *cq2 = new TCanvas("cq2","q_{2}",800,800);
  cq2->SetLogy();
  hq2->Draw();
  leg->AddEntry(hq2,Form("%d - %d %%",minCent,maxCent));
  Int_t selcentbins[3] = {0,ncentbins/2-1,ncentbins-1};
  Int_t colorline[3] = {kRed,kBlue,kGreen+2};
  for(Int_t iCent=0; iCent<3; iCent++) {
    hq2centbin[selcentbins[iCent]]->SetLineColor(colorline[iCent]);
    hq2centbin[selcentbins[iCent]]->Draw("same");
    Double_t min=(Double_t)minCent+selcentbins[iCent]*centwidth;
    Double_t max=(Double_t)minCent+(selcentbins[iCent]+1)*centwidth;
    leg->AddEntry(hq2centbin[selcentbins[iCent]],Form("%.1f - %.1f %%",min,max));
  }
  leg->Draw("same");

  TLatex* lat = new TLatex();
  lat->SetTextColor(kBlack);
  lat->SetTextFont(42);
  lat->SetTextSize(0.05);

  TCanvas* cq2Cut = 0x0;

  if(cutmeth==kAbsCut || cutmeth==kPercCut) {
    Double_t smallperc = -1.;
    Double_t largeperc = -1.;
    Int_t smallthresholdbin = hq2->GetXaxis()->FindBin(smallcutvalues[0]-0.0001);
    Int_t largethresholdbin = hq2->GetXaxis()->FindBin(largecutvalues[0]+0.0001);
    TH1F* hq2smallcut=(TH1F*)hq2->Clone();
    TH1F* hq2largecut=(TH1F*)hq2->Clone();
    hq2smallcut->SetLineWidth(0);
    hq2smallcut->SetFillColorAlpha(colors[0]+1,0.6);
    hq2largecut->SetLineWidth(0);
    hq2largecut->SetFillColorAlpha(colors[1]+1,0.6);

    Double_t totnum= 0.;
    for(Int_t iBin=1; iBin<=smallthresholdbin; iBin++) {
      totnum+=hq2->GetBinContent(iBin);
    }
    for(Int_t iBin=smallthresholdbin+1; iBin<=hq2->GetNbinsX(); iBin++) {
      hq2smallcut->SetBinContent(iBin,0);
    }
    smallperc=totnum/hq2->GetEntries()*100;
    cout <<"The percentage of events with q2 < "<< smallcutvalues[0] << " for the centrality class "<<minCent<<"-"<<maxCent<< " is equal to " << smallperc<<"%"<<endl;

    totnum=0;
    for(Int_t iBin=hq2->GetNbinsX(); iBin>=largethresholdbin; iBin--) {
      totnum+=hq2->GetBinContent(iBin);
    }
    for(Int_t iBin=1; iBin<largethresholdbin; iBin++) {
      hq2largecut->SetBinContent(iBin,0);
    }
    largeperc=totnum/hq2->GetEntries()*100;
    cout <<"The percentage of events with q2 > "<< largecutvalues[0] << " for the centrality class "<<minCent<<"-"<<maxCent<<" is equal to " << largeperc<<"%"<<endl;

    cq2Cut = new TCanvas("cq2Cut","",800,800);
    cq2Cut->SetLogy();
    hq2->Draw();
    hq2smallcut->Draw("same");
    hq2largecut->Draw("same");
    lat->DrawLatex(5.5,hq2->GetMaximum()*0.5,Form("Small-%s: %0.1f%%",q2name.Data(),smallperc));
    lat->DrawLatex(5.5,hq2->GetMaximum()*0.1,Form("Large-%s: %0.1f%%",q2name.Data(),largeperc));
  }

  TH1F** hq2smallcutcentbin = new TH1F*[ncentbins];
  TH1F** hq2largecutcentbin = new TH1F*[ncentbins];
  TH1F** hq2centbinclone = new TH1F*[ncentbins];
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    hq2smallcutcentbin[iCent]=(TH1F*)hq2centbin[iCent]->Clone();
    hq2largecutcentbin[iCent]=(TH1F*)hq2centbin[iCent]->Clone();
    hq2smallcutcentbin[iCent]->SetLineWidth(0);
    hq2smallcutcentbin[iCent]->SetLineColor(kBlack);
    hq2smallcutcentbin[iCent]->SetFillColorAlpha(colors[0]+1,0.6);
    hq2largecutcentbin[iCent]->SetLineWidth(0);
    hq2largecutcentbin[iCent]->SetLineColor(kBlack);
    hq2largecutcentbin[iCent]->SetFillColorAlpha(colors[1]+1,0.6);
  }

  TCanvas* cq2CutVsCent = new TCanvas("cq2CutVsCent","",1920,1080);
  if(ncentbins%2==0) {cq2CutVsCent->Divide(ncentbins/2,2);}
  else {cq2CutVsCent->Divide((ncentbins+1)/2,2);}

  cout << endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    Double_t totnum=0.;
    Double_t smallthresholdbin=hq2centbin[iCent]->GetXaxis()->FindBin(smallcutvalues[iCent]-0.0001);
    for(Int_t iBin=1; iBin<=smallthresholdbin; iBin++) {
      totnum+=hq2centbin[iCent]->GetBinContent(iBin);
    }
    for(Int_t iBin=smallthresholdbin+1; iBin<hq2centbin[iCent]->GetNbinsX(); iBin++) {
      hq2smallcutcentbin[iCent]->SetBinContent(iBin,0);
    }
    Double_t smallperc=totnum/hq2centbin[iCent]->GetEntries()*100;
    cout <<"The percentage of events with q2 < "<< smallcutvalues[iCent] << " for the centrality class "<<minCent+(centwidth*iCent)<<"-"<<minCent+(centwidth*(iCent+1))<< " is equal to " << smallperc<<"%"<<endl;
    cq2CutVsCent->cd(iCent+1)->SetLogy();
    hq2centbinclone[iCent] = (TH1F*)hq2centbin[iCent]->Clone();
    hq2centbinclone[iCent]->SetLineColor(kBlack);
    hq2centbinclone[iCent]->SetTitle(Form("%0.1f-%0.1f%%",minCent+(centwidth*iCent),minCent+(centwidth*(iCent+1))));
    hq2centbinclone[iCent]->Draw();
    hq2smallcutcentbin[iCent]->Draw("same");
    lat->DrawLatex(5.5,hq2centbin[iCent]->GetMaximum()*0.5,Form("Small-%s: %0.1f%%",q2name.Data(),smallperc));
  }
  cout << endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    Double_t totnum=0.;
    Double_t largethresholdbin=hq2centbin[iCent]->GetXaxis()->FindBin(largecutvalues[iCent]+0.0001);
    for(Int_t iBin=largethresholdbin; iBin<=hq2centbin[iCent]->GetNbinsX(); iBin++) {
      totnum+=hq2centbin[iCent]->GetBinContent(iBin);
    }
    for(Int_t iBin=1; iBin<largethresholdbin; iBin++) {
      hq2largecutcentbin[iCent]->SetBinContent(iBin,0);
    }
    Double_t largeperc=totnum/hq2centbin[iCent]->GetEntries()*100;
    cout <<"The percentage of events with q2 > "<< largecutvalues[iCent] << " for the centrality class "<<minCent+(centwidth*iCent)<<"-"<<minCent+(centwidth*(iCent+1))<< " is equal to " << largeperc<<"%"<<endl;
    cq2CutVsCent->cd(iCent+1);
    hq2largecutcentbin[iCent]->Draw("same");
    lat->DrawLatex(5.5,hq2centbin[iCent]->GetMaximum()*0.1,Form("Large-%s: %0.1f%%",q2name.Data(),largeperc));
  }
  cout << endl;

  cq2->SaveAs(Form("%s/Nev_vs_q2%s.pdf",outputdir.Data(),suffix.Data()));
  cq2->SaveAs(Form("%s/Nev_vs_q2%s.root",outputdir.Data(),suffix.Data()));

  TString outname=Form("%s/q2Selection%s_q2Small%0.2f_q2Large%0.2f.pdf",outputdir.Data(),suffix.Data(),q2smalllimit,q2largelimit);
  TString percsuffix="perc";
  if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/q2Selection%s_q2Small%0.f%s_q2Large%0.f%s.pdf",outputdir.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
  if(cutmeth!=kPercCutVsCent) {cq2Cut->SaveAs(outname.Data());}
  outname.ReplaceAll("q2Selection","q2SelectionVsCent");
  cq2CutVsCent->SaveAs(outname.Data());
  outname.ReplaceAll("q2SelectionVsCent","q2_vs_Centrality");
  cq2VsCent->SaveAs(outname.Data());
  outname.ReplaceAll(".pdf",".root");
  cq2VsCent->SaveAs(outname.Data());
}

//_____________________________________________________________________________________________
//DRAW RESOLUTION
void DrawEventPlaneResolutionAndDistribution(Int_t cutmeth) {

  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}

  SetStyle();

  Int_t ncentbins=0;
  TList* datalist=(TList*)LoadTList(ncentbins);
  if(!datalist){return;}
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) return;

  Double_t centwidth=(Double_t)(maxCent-minCent)/ncentbins;
  vector<Double_t> smallcutvalues;
  vector<Double_t> largecutvalues;
  vector<Double_t> intsmallcutvalues;
  vector<Double_t> intlargecutvalues;
  Bool_t defq2cut=Defineq2Cuts(hq2VsCentr,smallcutvalues,largecutvalues,cutmeth,1);
  if(!defq2cut) return;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    intsmallcutvalues.push_back(11.);
    intlargecutvalues.push_back(-1.);
  }

  Int_t q2region[3] = {kSmall,kLarge,kIntegrated};

  TGraphErrors** gRes = new TGraphErrors*[3];
  TF1** fRes = new TF1*[3];
  TList *resolist=0x0;
  Double_t error[3]={0.,0.,0.};
  Double_t resoFull[3]={0.,0.,0.};
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
  TString q2suffix[3] = {"q2Small","q2Large","q2Int"};

  for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
    resolist=(TList*)LoadResolutionHistos(datalist,smallcutvalues,largecutvalues,q2region[iq2]);
    gRes[iq2]=(TGraphErrors*)resolist->FindObject("gResolVsCent");
    gRes[iq2]->SetName(Form("gResolVsCent%s",q2suffix[iq2].Data()));
    gRes[iq2]->GetYaxis()->SetTitle("Event Plane resolution R_{2}");
    gRes[iq2]->GetXaxis()->SetTitle("Centrality (%)");
    gRes[iq2]->SetTitle("");
    gRes[iq2]->SetMarkerSize(1.5);
    gRes[iq2]->SetLineWidth(2);
    gRes[iq2]->SetMarkerStyle(markers[iq2]);
    gRes[iq2]->SetMarkerColor(colors[iq2]+1);
    gRes[iq2]->SetLineColor(colors[iq2]+1);
    TH1F* hevplresos[3];
    TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
    Int_t nSubRes=1;
    TH2F* htestversion=(TH2F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[0].Data(),suffixcentr.Data()));
    if(htestversion){
      printf("Old version of the task\n");
    }else{
      printf("New version of the task\n");
      namereso[0]="Reso1Vsq2";
    }
    if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
       evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
    for(Int_t iRes=0;iRes<nSubRes;iRes++){
      hevplresos[iRes]=(TH1F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
    }
    resoFull[iq2]=GetEventPlaneResolution(error[iq2],hevplresos[0],hevplresos[1],hevplresos[2]);
    fRes[iq2] = new TF1(Form("fResol%s",q2suffix[iq2].Data()),"pol0",minCent,maxCent);
    fRes[iq2]->SetParameter(0,resoFull[iq2]);
    fRes[iq2]->SetLineWidth(2);
    fRes[iq2]->SetLineColor(colors[iq2+2]);
  }

  TLegend* leg = new TLegend(0.2,0.2,0.8,0.5);
  leg->SetFillStyle(0);
  if(cutmeth==kAbsCut) {
    leg->AddEntry(gRes[0],Form("%s < %0.1f, <R_{2}> = %0.4f",hq2VsCentr->GetYaxis()->GetTitle(),q2smalllimit,resoFull[0]),"lpe");
    leg->AddEntry(gRes[1],Form("%s > %0.1f, <R_{2}> = %0.4f",hq2VsCentr->GetYaxis()->GetTitle(),q2largelimit,resoFull[1]),"lpe");
  }
  else {
    Double_t perc[2] = {q2smallpercevents*100,q2largepercevents*100};
    leg->AddEntry(gRes[0],Form("%0.f%% small %s, <R_{2}> = %0.4f",perc[0],hq2VsCentr->GetYaxis()->GetTitle(),resoFull[0]),"lpe");
    leg->AddEntry(gRes[1],Form("%0.f%% large %s, <R_{2}> = %0.4f",perc[1],hq2VsCentr->GetYaxis()->GetTitle(),resoFull[1]),"lpe");
  }
  leg->AddEntry(gRes[2],Form("%s-integrated, <R_{2}> = %0.4f",hq2VsCentr->GetYaxis()->GetTitle(),resoFull[2]),"lpe");

  TCanvas* cRes = new TCanvas("cRes","",800,800);
  Int_t npoints = gRes[0]->GetN();
  gRes[0]->GetYaxis()->SetRangeUser(0.,1.);
  gRes[0]->Draw("AP");
  gRes[1]->Draw("P");
  gRes[2]->Draw("P");
  leg->Draw("same");

  TString outname=Form("%s/EP_resolution%s_q2Small%.2f_q2Large%.2f.pdf",outputdir.Data(),suffix.Data(),q2smalllimit,q2largelimit);
  TString percsuffix="perc";
  if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/EP_resolution%s_q2Small%.0f%s_q2Large%.0f%s.pdf",outputdir.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}

  cRes->SaveAs(outname.Data());
  outname.ReplaceAll(".pdf",".root");
  TFile outfile(outname.Data(),"RECREATE");
  cRes->Write();
  gRes[0]->Write();
  gRes[1]->Write();
  gRes[2]->Write();
  fRes[0]->Write();
  fRes[1]->Write();
  fRes[2]->Write();
  outfile.Close();

  TH3F* hEvPlaneVsq2VsCent[6];
  TString DetNames[6] = {"TPC","TPCPosEta","TPCNegEta","VZERO","VZEROA","VZEROC"};
  hEvPlaneVsq2VsCent[0] = (TH3F*)datalist->FindObject(Form("hEvPlaneQncorr%sQoverMVsq2VsCent",DetNames[0].Data()));
  if(hEvPlaneVsq2VsCent[0]) {
    TH1F* hEvPlaneq2Int[6];
    TH1F* hEvPlaneLargeq2[6];
    TH1F* hEvPlaneSmallq2[6];
    TF1* fCosq2Int[6];
    TF1* fSinq2Int[6];
    TF1* fCosLargeq2[6];
    TF1* fSinLargeq2[6];
    TF1* fCosSmallq2[6];
    TF1* fSinSmallq2[6];
    TCanvas* cEP = new TCanvas("cEP","",1920,1080);
    cEP->Divide(3,2);
    TCanvas* cCosSin2Psi = new TCanvas("cCosSin2Psi","",1920,1080);
    cCosSin2Psi->Divide(3,2);
    TH1F* hCos2Psiq2Int = new TH1F("hCos2Psiq2Int","",6,-0.5,5.5);
    TH1F* hSin2Psiq2Int = new TH1F("hSin2Psiq2Int","",6,-0.5,5.5);
    hCos2Psiq2Int->SetLineColor(colors[2]);
    hSin2Psiq2Int->SetLineColor(colors[2]);
    hCos2Psiq2Int->SetMarkerColor(colors[2]);
    hSin2Psiq2Int->SetMarkerColor(colors[2]);
    hCos2Psiq2Int->SetMarkerStyle(kFullCircle);
    hSin2Psiq2Int->SetMarkerStyle(kOpenSquare);
    TH1F* hCos2PsiSmallq2 = new TH1F("hCos2PsiSmallq2","",6,-0.5,5.5);
    TH1F* hSin2PsiSmallq2 = new TH1F("hSin2PsiSmallq2","",6,-0.5,5.5);
    hCos2PsiSmallq2->SetLineColor(colors[0]+1);
    hSin2PsiSmallq2->SetLineColor(colors[0]+1);
    hCos2PsiSmallq2->SetMarkerColor(colors[0]+1);
    hSin2PsiSmallq2->SetMarkerColor(colors[0]+1);
    hCos2PsiSmallq2->SetMarkerStyle(kFullCircle);
    hSin2PsiSmallq2->SetMarkerStyle(kOpenSquare);
    TH1F* hCos2PsiLargeq2 = new TH1F("hCos2PsiLargeq2","",6,-0.5,5.5);
    TH1F* hSin2PsiLargeq2 = new TH1F("hSin2PsiLargeq2","",6,-0.5,5.5);
    hCos2PsiLargeq2->SetLineColor(colors[1]+1);
    hSin2PsiLargeq2->SetLineColor(colors[1]+1);
    hCos2PsiLargeq2->SetMarkerColor(colors[1]+1);
    hSin2PsiLargeq2->SetMarkerColor(colors[1]+1);
    hCos2PsiLargeq2->SetMarkerStyle(kFullCircle);
    hSin2PsiLargeq2->SetMarkerStyle(kOpenSquare);

    TH2F* hFrameBias = new TH2F("hFrameBias","",6,-0.5,5.5,1,-0.06,0.06);
    hFrameBias->GetYaxis()->SetTitle("< #it{f} (2#psi_{EP}) >");

    TLegend* leg2 = new TLegend(0.55,0.7,0.85,0.85);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(hCos2PsiSmallq2,"small q_{2}","l");
    leg2->AddEntry(hCos2PsiLargeq2,"large q_{2}","l");
    leg2->AddEntry(hCos2Psiq2Int,"q_{2}-integrated","l");

    TLegend* leg3 = new TLegend(0.55,0.2,0.85,0.3);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.04);
    leg3->AddEntry(hCos2Psiq2Int,"< cos(2#psi_{EP}) >","p");
    leg3->AddEntry(hSin2Psiq2Int,"< sin(2#psi_{EP}) >","p");

    for(Int_t iDet=0; iDet<6; iDet++) {
      if(iDet>0) {hEvPlaneVsq2VsCent[iDet]=(TH3F*)datalist->FindObject(Form("hEvPlaneQncorr%sQoverMVsq2VsCent",DetNames[iDet].Data()));}
      if(!hEvPlaneVsq2VsCent[iDet]) {continue;}
      hEvPlaneq2Int[iDet] = (TH1F*)hEvPlaneVsq2VsCent[iDet]->ProjectionZ(Form("hEvPlaneQncorr%sQoverMVsq2VsCent_q2Int",DetNames[iDet].Data()));
      hEvPlaneSmallq2[iDet] = (TH1F*)hEvPlaneq2Int[iDet]->Clone(Form("hEvPlaneQncorr%sQoverMVsq2VsCent_q2Int",DetNames[iDet].Data()));
      hEvPlaneLargeq2[iDet] = (TH1F*)hEvPlaneq2Int[iDet]->Clone(Form("hEvPlaneQncorr%sQoverMVsq2VsCent_q2Int",DetNames[iDet].Data()));
      TAxis* q2axis = (TAxis*)hEvPlaneVsq2VsCent[iDet]->GetYaxis();
      for(Int_t iEP=0; iEP<hEvPlaneq2Int[iDet]->GetNbinsX(); iEP++) {
        Double_t nentriesq2small=0;
        Double_t nentriesq2large=0;
        for(Int_t iCent=0; iCent<ncentbins; iCent++) {
          for(Int_t iq2=0; iq2<q2axis->GetNbins(); iq2++) {
            Double_t q2binlowerlimit = q2axis->GetBinLowEdge(iq2+1);
            Double_t q2binupperlimit = q2axis->GetBinLowEdge(iq2+1)+q2axis->GetBinWidth(iq2+1);
            if(q2binlowerlimit>=largecutvalues[iCent]) {nentriesq2large += hEvPlaneVsq2VsCent[iDet]->GetBinContent(iCent+1,iq2+1,iEP+1);}
            if(q2binupperlimit<=smallcutvalues[iCent]) {nentriesq2small += hEvPlaneVsq2VsCent[iDet]->GetBinContent(iCent+1,iq2+1,iEP+1);}
          }
        }
        hEvPlaneSmallq2[iDet]->SetBinContent(iEP+1,nentriesq2small);
        hEvPlaneLargeq2[iDet]->SetBinContent(iEP+1,nentriesq2large);
        hEvPlaneSmallq2[iDet]->SetBinError(iEP+1,TMath::Sqrt(nentriesq2small));
        hEvPlaneLargeq2[iDet]->SetBinError(iEP+1,TMath::Sqrt(nentriesq2large));
        hEvPlaneq2Int[iDet]->SetBinError(iEP+1,TMath::Sqrt(hEvPlaneq2Int[iDet]->GetBinContent(iEP+1)));
      }
      hEvPlaneSmallq2[iDet]->Sumw2();
      hEvPlaneLargeq2[iDet]->Sumw2();
      hEvPlaneq2Int[iDet]->Scale(1./hEvPlaneq2Int[iDet]->Integral());
      hEvPlaneSmallq2[iDet]->Scale(1./hEvPlaneSmallq2[iDet]->Integral());
      hEvPlaneLargeq2[iDet]->Scale(1./hEvPlaneLargeq2[iDet]->Integral());
      hEvPlaneq2Int[iDet]->SetLineColor(colors[2]);
      hEvPlaneSmallq2[iDet]->SetLineColor(colors[0]+1);
      hEvPlaneLargeq2[iDet]->SetLineColor(colors[1]+1);
      hEvPlaneq2Int[iDet]->GetYaxis()->SetTitle("Normalised entries");
      hEvPlaneSmallq2[iDet]->GetYaxis()->SetTitle("Normalised entries");
      hEvPlaneLargeq2[iDet]->GetYaxis()->SetTitle("Normalised entries");
      hEvPlaneq2Int[iDet]->GetYaxis()->SetTitleOffset(1.7);
      hEvPlaneSmallq2[iDet]->GetYaxis()->SetTitleOffset(1.7);
      hEvPlaneLargeq2[iDet]->GetYaxis()->SetTitleOffset(1.7);
      hEvPlaneq2Int[iDet]->GetXaxis()->SetTitle("#psi_{EP}");
      hEvPlaneSmallq2[iDet]->GetXaxis()->SetTitle("#psi_{EP}");
      hEvPlaneLargeq2[iDet]->GetXaxis()->SetTitle("#psi_{EP}");

      fCosq2Int[iDet] = new TF1(Form("fCosq2Int%s",DetNames[iDet].Data()),"[0]*(1+2*[1]*cos(2*x))",0.,TMath::Pi());
      fSinq2Int[iDet] = new TF1(Form("fSinq2Int%s",DetNames[iDet].Data()),"[0]*(1+2*[1]*sin(2*x))",0.,TMath::Pi());
      fCosq2Int[iDet]->SetLineColor(colors[2]);
      fSinq2Int[iDet]->SetLineColor(colors[2]);
      fCosq2Int[iDet]->SetLineStyle(5);
      fSinq2Int[iDet]->SetLineStyle(9);

      fCosLargeq2[iDet] = new TF1(Form("fCosLargeq2%s",DetNames[iDet].Data()),"[0]*(1+2*[1]*cos(2*x))",0.,TMath::Pi());
      fSinLargeq2[iDet] = new TF1(Form("fSinLargeq2%s",DetNames[iDet].Data()),"[0]*(1+2*[1]*sin(2*x))",0.,TMath::Pi());
      fCosLargeq2[iDet]->SetLineColor(colors[1]+1);
      fSinLargeq2[iDet]->SetLineColor(colors[1]+1);
      fCosLargeq2[iDet]->SetLineStyle(5);
      fSinLargeq2[iDet]->SetLineStyle(9);

      fCosSmallq2[iDet] = new TF1(Form("fCosSmallq2%s",DetNames[iDet].Data()),"[0]*(1+2*[1]*cos(2*x))",0.,TMath::Pi());
      fSinSmallq2[iDet] = new TF1(Form("fSinSmallq2%s",DetNames[iDet].Data()),"[0]*(1+2*[1]*sin(2*x))",0.,TMath::Pi());
      fCosSmallq2[iDet]->SetLineColor(colors[0]+1);
      fSinSmallq2[iDet]->SetLineColor(colors[0]+1);
      fCosSmallq2[iDet]->SetLineStyle(5);
      fSinSmallq2[iDet]->SetLineStyle(9);

      cEP->cd(iDet+1);
      hEvPlaneq2Int[iDet]->GetYaxis()->SetRangeUser(hEvPlaneLargeq2[iDet]->GetMinimum()*0.9,hEvPlaneLargeq2[iDet]->GetMaximum()*1.1);
      hEvPlaneq2Int[iDet]->Fit(Form("fCosq2Int%s",DetNames[iDet].Data()),"N");
      hEvPlaneSmallq2[iDet]->Fit(Form("fCosSmallq2%s",DetNames[iDet].Data()),"N");
      hEvPlaneLargeq2[iDet]->Fit(Form("fCosLargeq2%s",DetNames[iDet].Data()),"N");
      hEvPlaneq2Int[iDet]->Fit(Form("fSinq2Int%s",DetNames[iDet].Data()),"N");
      hEvPlaneSmallq2[iDet]->Fit(Form("fSinSmallq2%s",DetNames[iDet].Data()),"N");
      hEvPlaneLargeq2[iDet]->Fit(Form("fSinLargeq2%s",DetNames[iDet].Data()),"N");

      hEvPlaneq2Int[iDet]->Draw("E");
      hEvPlaneSmallq2[iDet]->Draw("Esame");
      hEvPlaneLargeq2[iDet]->Draw("Esame");
      if(iDet==0) {leg2->Draw("same");}

      hCos2Psiq2Int->SetBinContent(iDet+1,fCosq2Int[iDet]->GetParameter(1));
      hCos2Psiq2Int->SetBinError(iDet+1,fCosq2Int[iDet]->GetParError(1));
      hCos2PsiLargeq2->SetBinContent(iDet+1,fCosLargeq2[iDet]->GetParameter(1));
      hCos2PsiLargeq2->SetBinError(iDet+1,fCosLargeq2[iDet]->GetParError(1));
      hCos2PsiSmallq2->SetBinContent(iDet+1,fCosSmallq2[iDet]->GetParameter(1));
      hCos2PsiSmallq2->SetBinError(iDet+1,fCosSmallq2[iDet]->GetParError(1));

      hSin2Psiq2Int->SetBinContent(iDet+1,fSinq2Int[iDet]->GetParameter(1));
      hSin2Psiq2Int->SetBinError(iDet+1,fSinq2Int[iDet]->GetParError(1));
      hSin2PsiLargeq2->SetBinContent(iDet+1,fSinLargeq2[iDet]->GetParameter(1));
      hSin2PsiLargeq2->SetBinError(iDet+1,fSinLargeq2[iDet]->GetParError(1));
      hSin2PsiSmallq2->SetBinContent(iDet+1,fSinSmallq2[iDet]->GetParameter(1));
      hSin2PsiSmallq2->SetBinError(iDet+1,fSinSmallq2[iDet]->GetParError(1));

      hFrameBias->GetXaxis()->SetBinLabel(iDet+1,DetNames[iDet]);
    }
    cCosSin2Psi->cd();
    hCos2Psiq2Int->GetYaxis()->SetRangeUser(-0.05,0.05);
    hFrameBias->Draw();
    hCos2Psiq2Int->Draw("same");
    hCos2PsiLargeq2->Draw("same");
    hCos2PsiSmallq2->Draw("same");
    hSin2Psiq2Int->Draw("same");
    hSin2PsiLargeq2->Draw("same");
    hSin2PsiSmallq2->Draw("same");
    leg2->Draw("same");
    leg3->Draw("same");

    TString outnameEP=Form("%s/EP_distributions%s_q2Small%.2f_q2Large%.2f.pdf",outputdir.Data(),suffix.Data(),q2smalllimit,q2largelimit);
    TString percsuffix="perc";
    if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
    if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outnameEP=Form("%s/EP_resolution%s_q2Small%.0f%s_q2Large%.0f%s.pdf",outputdir.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}

    cEP->SaveAs(outnameEP.Data());

    outnameEP.ReplaceAll(".pdf",".root");
    TFile outfileEP(outnameEP.Data(),"RECREATE");
    for(Int_t iDet=0; iDet<6; iDet++) {
      hEvPlaneq2Int[iDet]->Write();
      hEvPlaneSmallq2[iDet]->Write();
      hEvPlaneLargeq2[iDet]->Write();
      fCosq2Int[iDet]->Write();
      fSinq2Int[iDet]->Write();
      fCosSmallq2[iDet]->Write();
      fSinSmallq2[iDet]->Write();
      fCosLargeq2[iDet]->Write();
      fSinLargeq2[iDet]->Write();
    }
    hCos2Psiq2Int->Write();
    hSin2Psiq2Int->Write();
    hCos2PsiSmallq2->Write();
    hSin2PsiSmallq2->Write();
    hCos2PsiLargeq2->Write();
    hSin2PsiLargeq2->Write();
    outfileEP.Close();
  }
}

//_____________________________________________________________________________________________
//LOAD DATA LIST
TList* LoadTList(Int_t &ncentbins) {

  cout << "Opening input file " << infilename << "..." <<endl;
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  TDirectoryFile* dir=0x0;
  TList* list=0x0;

  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s",partname.Data(),suffix.Data());

  if(infile) {dir=(TDirectoryFile*)infile->Get(dirname.Data()); cout << "File opened!" << endl;}
  else {cerr << "Error: File " << infilename << " not found. Exit." << endl; return 0x0;}
  if(dir) list=(TList*)dir->Get(listname.Data());
  else {cerr << "Error: Wrong TDirectoryFile name " << dirname << ". Exit." << endl; return 0x0;}
  if(!list) {cerr << "Error: Wrong TList name " << listname << ". Exit." << endl; return 0x0;}

  TH2F* hTest = (TH2F*)list->FindObject("hq2TPCFullEtaVsPosEta");
  TH2F* hTest2 = (TH2F*)list->FindObject("hMultVsCentFullTPC");
  TH2F* hTest3 = (TH2F*)list->FindObject("hEvPlaneQncorrTPCQoverMVsq2VsCent");
  if(hTest && hTest2 && hTest3) {
    ncentbins=(list->GetEntries()-10)/15;
  }
  else if(hTest && hTest2 && !hTest3) {
    ncentbins=(list->GetEntries()-10)/9;
  }
  else if(hTest && !hTest2 && !hTest3) {
    ncentbins=(list->GetEntries()-10)/6;
  }
  else {
    ncentbins=(list->GetEntries()-10)/3;
  }

  infile->Close();
  cout<<"Input file closed."<< endl;

  return list;
}

//_____________________________________________________________________________________________
//GET SPARSE
THnSparseF* LoadSparseFromList(TList* inputlist) {

  THnSparseF* sparse=(THnSparseF*)inputlist->FindObject("hMassPtPhiq2Centr");
  if(!sparse) {cerr << "Error: No THnSparse found. Check if the name hMassPtPhiq2Centr is right!" << endl; return 0x0;}
  else {cout << "THnSparse got!" << endl;}

  return sparse;
}

//_____________________________________________________________________________________________
//GET q2 VS. CENTRALITY HISTO
TH2F* GetHistoq2VsCentr(TList* inputlist) {

  Int_t ncentbins=0;
  TH2F* hTest = (TH2F*)inputlist->FindObject("hq2TPCFullEtaVsPosEta");
  TH2F* hTest2 = (TH2F*)inputlist->FindObject("hMultVsCentFullTPC");
  TH2F* hTest3 = (TH2F*)inputlist->FindObject("hEvPlaneQncorrTPCQoverMVsq2VsCent");
  if(hTest && hTest2 && hTest3) {
    ncentbins=(inputlist->GetEntries()-10)/15;
  }
  else if(hTest && hTest2 && !hTest3) {
    ncentbins=(inputlist->GetEntries()-10)/9;
  }
  else if(hTest && !hTest2 && !hTest3) {
    ncentbins=(inputlist->GetEntries()-10)/6;
  }
  else {
    ncentbins=(inputlist->GetEntries()-10)/3;
  }
  Double_t centwidth = (Double_t)(maxCent-minCent)/ncentbins;
  Double_t mincentpermil=minCent*10;
  Double_t maxcentpermil=(minCent+centwidth)*10;

  TH2F* hResVsq2=(TH2F*)inputlist->FindObject(Form("hEvPlaneReso1Vsq2centr%0.f_%0.f",mincentpermil,maxcentpermil));
  if(!hResVsq2) {return 0x0;}
  TH1F* hq2=(TH1F*)hResVsq2->ProjectionY();

  Int_t nq2bins = hq2->GetNbinsX();
  Int_t q2min = hq2->GetBinLowEdge(1);
  Int_t q2max = hq2->GetBinLowEdge(nq2bins)+hq2->GetBinWidth(nq2bins);

  TH2F* hq2VsCentr = new TH2F("hq2VsCentr","q_{2} vs. centrality",ncentbins,minCent,maxCent,nq2bins,q2min,q2max);
  hq2VsCentr->GetXaxis()->SetTitle("Centrality (%)");
  hq2VsCentr->GetYaxis()->SetTitle(hq2->GetXaxis()->GetTitle());

  for(Int_t iBin=0; iBin<nq2bins+1; iBin++) {
    hq2VsCentr->SetBinContent(1,iBin+1,hq2->GetBinContent(iBin+1)); //takes also overflow counts
  }

  for(Int_t iCent=1; iCent<ncentbins; iCent++) {
    mincentpermil=(minCent+iCent*centwidth)*10;
    maxcentpermil=(minCent+(iCent+1)*centwidth)*10;

    hResVsq2=(TH2F*)inputlist->FindObject(Form("hEvPlaneReso1Vsq2centr%0.f_%0.f",mincentpermil,maxcentpermil));
    if(!hResVsq2) return 0x0;
    hq2 = (TH1F*)hResVsq2->ProjectionY();
    for(Int_t iBin=0; iBin<hq2->GetNbinsX()+1; iBin++) {
      hq2VsCentr->SetBinContent(iCent+1,iBin+1,hq2->GetBinContent(iBin+1)); //takes also overflow counts
    }
    hq2VsCentr->GetYaxis()->SetTitle(hq2->GetXaxis()->GetTitle());
  }
  hq2VsCentr->SetStats(0);

  return hq2VsCentr;
}

//_____________________________________________________________________________________________
//GET MASS HISTOS IN SELECTED q2 RANGE
TList* LoadMassHistos(THnSparseF* sparse, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge, Int_t analysismeth) {

  Int_t ptaxis=1;
  Int_t phiaxis=2;
  Int_t q2axis=3;
  Int_t centaxis=4;

  const Int_t ncentbins = (const Int_t)smallcutvalues.size();
  if(ncentbins!=(const Int_t)largecutvalues.size()) {return 0x0;}
  Double_t centwidth = (Double_t)(maxCent-minCent)/ncentbins;

  TList *outlist = new TList();
  outlist->SetName("azimuthalhistoslist");

  TH1F* hMassTmp[ncentbins];
  TH1F* hMassTmp1[ncentbins];
  TH1F* hMassTmp2[ncentbins];

  Int_t fSparseVers = 0; //check version of the task output
  if(sparse->GetAxis(centaxis)) {
    TString centtitle=sparse->GetAxis(centaxis)->GetTitle();
    if(centtitle.Contains("Centrality") || centtitle.Contains("centrality")) {
      cout << "Version of the sparse with centrality axis!" << endl;
      fSparseVers=1;
    }
  }

  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    //get phi- and q2-integrated histos (for fixed sigma)
    ApplyCut(sparse,PtLims[iPt],PtLims[iPt+1],ptaxis); //apply pt cut
    hMassTmp[0]=(TH1F*)sparse->Projection(0);
    hMassTmp[0]->SetName(Form("hMass_pt%d",iPt));
    outlist->Add(hMassTmp[0]->Clone());
    delete hMassTmp[0];
    hMassTmp[0]=0x0;

    //get histos in deltaphi bins and with the q2 cut
    if(analysismeth==kEventPlane) {
      for(Int_t iPhi=0; iPhi<nPhiBins; iPhi++) {
        ApplyCut(sparse,PhiLims[iPhi],PhiLims[iPhi+1],phiaxis); //apply deltaphi cut
        if(fSparseVers==0) {
          Applyq2Cut(sparse,smallcutvalues[0],largecutvalues[0],smallorlarge); //apply q2 cut
          hMassTmp[0]=(TH1F*)sparse->Projection(0);
        }
        else if(fSparseVers==1) {
          for(Int_t iCent=0; iCent<ncentbins; iCent++) {
            ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
            if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
            hMassTmp[iCent]=(TH1F*)sparse->Projection(0);
            if(iCent>0) {hMassTmp[0]->Add(hMassTmp[iCent]);}
            ResetAxes(sparse,centaxis); //release centrality cut
            ResetAxes(sparse,q2axis); //release q2 cut
            if(iCent!=0) {
              delete hMassTmp[iCent];
              hMassTmp[iCent]=0x0;
            }
          }
        }
        hMassTmp[0]->SetName(Form("hMass_pt%d_phi%d",iPt,iPhi));
        outlist->Add(hMassTmp[0]->Clone());
        ResetAxes(sparse,q2axis); //release q2 cut
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        delete hMassTmp[0];
        hMassTmp[0]=0x0;
      }
    }
    else if(analysismeth==kEventPlaneInOut){
      if(fSparseVers==0) {
        if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[0],largecutvalues[0],smallorlarge);} //apply q2 cut
        ApplyCut(sparse,0.,TMath::Pi()/4,phiaxis); //apply deltaphi cut
        hMassTmp[0]=(TH1F*)sparse->Projection(0);
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        ApplyCut(sparse,3./4*TMath::Pi(),TMath::Pi(),phiaxis); //apply deltaphi cut
        hMassTmp1[0]=(TH1F*)sparse->Projection(0);
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        hMassTmp[0]->Add(hMassTmp1[0]);
        hMassTmp[0]->SetName(Form("hMass_pt%d_phi0",iPt));
        outlist->Add(hMassTmp[0]->Clone());
        ApplyCut(sparse,TMath::Pi()/4,3./4*TMath::Pi(),phiaxis); //apply deltaphi cut
        hMassTmp2[0]=(TH1F*)sparse->Projection(0);
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        ResetAxes(sparse,q2axis); //release q2 cut
        hMassTmp2[0]->SetName(Form("hMass_pt%d_phi1",iPt));
        outlist->Add(hMassTmp2[0]->Clone());
      }

      else if(fSparseVers==1) {
        ApplyCut(sparse,0.,TMath::Pi()/4,phiaxis); //apply deltaphi cut
        for(Int_t iCent=0; iCent<ncentbins; iCent++) {
          ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
          if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
          hMassTmp[iCent]=(TH1F*)sparse->Projection(0);
          if(iCent>0) {hMassTmp[0]->Add(hMassTmp[iCent]);}
          ResetAxes(sparse,centaxis); //release centrality cut
          ResetAxes(sparse,q2axis); //release q2 cut
          if(iCent!=0) {
            delete hMassTmp[iCent];
            hMassTmp[iCent]=0x0;
          }
        }
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        ApplyCut(sparse,3./4*TMath::Pi(),TMath::Pi(),phiaxis); //apply deltaphi cut
        for(Int_t iCent=0; iCent<ncentbins; iCent++) {
          ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
          if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
          hMassTmp1[iCent]=(TH1F*)sparse->Projection(0);
          if(iCent>0) {hMassTmp1[0]->Add(hMassTmp1[iCent]);}
          ResetAxes(sparse,centaxis); //release centrality cut
          ResetAxes(sparse,q2axis); //release q2 cut
          if(iCent!=0) {
            delete hMassTmp1[iCent];
            hMassTmp1[iCent]=0x0;
          }
        }
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        hMassTmp[0]->Add(hMassTmp1[0]);
        hMassTmp[0]->SetName(Form("hMass_pt%d_phi0",iPt));
        outlist->Add(hMassTmp[0]->Clone());
        ApplyCut(sparse,TMath::Pi()/4,3./4*TMath::Pi(),phiaxis); //apply deltaphi cut
        for(Int_t iCent=0; iCent<ncentbins; iCent++) {
          ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
          if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
          hMassTmp2[iCent]=(TH1F*)sparse->Projection(0);
          if(iCent>0) {hMassTmp2[0]->Add(hMassTmp2[iCent]);}
          ResetAxes(sparse,centaxis); //release centrality cut
          ResetAxes(sparse,q2axis); //release q2 cut
          if(iCent!=0) {
            delete hMassTmp2[iCent];
            hMassTmp2[iCent]=0x0;
          }
        }
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        hMassTmp2[0]->SetName(Form("hMass_pt%d_phi1",iPt));
        outlist->Add(hMassTmp2[0]->Clone());
      }
      delete hMassTmp[0];
      delete hMassTmp1[0];
      delete hMassTmp2[0];
      hMassTmp[0]=0x0;
      hMassTmp1[0]=0x0;
      hMassTmp2[0]=0x0;
      ResetAxes(sparse,ptaxis); //release pt cut
    }
  }
  ResetAxes(sparse,-1); //release all cuts

  if(useTemplD0Refl){
    Bool_t retCode=LoadD0toKpiMCHistos(outlist);
    if(!retCode)Printf("ERROR: MC histograms loading failed");
    else Printf("******************************************\n MC HISTOGRAMS LOADED\n\n**********************************");
  }

  return outlist;
}

//_____________________________________________________________________________________________
//GET RESOLUTION HISTOS IN SELECTED q2 RANGE
TList* LoadResolutionHistos(TList *inputlist, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge) {

  TList *outlist = new TList();
  outlist->SetName("eventplanemasslist");

  const Int_t nBins=20;
  Double_t ncoll[nBins]={1790.77,1578.44,1394.82,1236.17
    ,1095.08,969.86,859.571,759.959,669.648,589.588,516.039
    ,451.409,392.853,340.493,294.426,252.385,215.484,183.284
    ,155.101,130.963};

  Int_t ncentbins=0;
  TH2F* hTest = (TH2F*)inputlist->FindObject("hq2TPCFullEtaVsPosEta");
  TH2F* hTest2 = (TH2F*)inputlist->FindObject("hMultVsCentFullTPC");
  TH2F* hTest3 = (TH2F*)inputlist->FindObject("hEvPlaneQncorrTPCQoverMVsq2VsCent");
  if(hTest && hTest2 && hTest3) {
    ncentbins=(inputlist->GetEntries()-10)/15;
  }
  else if(hTest && hTest2 && !hTest3) {
    ncentbins=(inputlist->GetEntries()-10)/9;
  }
  else if(hTest && !hTest2 && !hTest3) {
    ncentbins=(inputlist->GetEntries()-10)/6;
  }
  else {
    ncentbins=(inputlist->GetEntries()-10)/3;
  }
  Double_t centwidth = (Double_t)(maxCent-minCent)/ncentbins;

  Double_t minCentTimesTen=minCent*10;
  Double_t maxCentTimesTen=maxCent*10;
  Double_t centwidthTimesTen=centwidth*10;

  TGraphErrors* gResolVsCent=new TGraphErrors(0);
  Int_t iCent=0;
  Int_t nSubRes=1;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
     evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
  TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
  TH2F* htestversion=(TH2F*)inputlist->FindObject(Form("hEvPlaneResocentr%0.f_%0.f",minCentTimesTen,minCentTimesTen+centwidthTimesTen));
  if(htestversion){
    printf("Old version of the task\n");
  }else{
    printf("New version of the task\n");
    namereso[0]="Reso1Vsq2";
  }
  TH1F* hevplresos[3];

  for(Int_t iCent=0; iCent<ncentbins; iCent++){
    TH2F* htmpresosvsq2[3];
    TH1F* htmpresos[3];

    for(Int_t iRes=0;iRes<nSubRes;iRes++){
      htmpresosvsq2[iRes]=(TH2F*)inputlist->FindObject(Form("hEvPlane%scentr%0.f_%0.f",namereso[iRes].Data(),minCentTimesTen+(iCent*centwidthTimesTen),minCentTimesTen+((iCent+1)*centwidthTimesTen)));
      Int_t binmin=-1;
      Int_t binmax=-1;

      //define bins for the q2 cut
      if(smallorlarge==kSmall) {
        binmin=1;
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(smallcutvalues[iCent]-0.0001);
      }
      else if(smallorlarge==kLarge) {
        binmin=htmpresosvsq2[iRes]->GetYaxis()->FindBin(largecutvalues[iCent]+0.0001);
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(htmpresosvsq2[iRes]->GetYaxis()->GetNbins()+1); //takes also overflow entries
        Int_t binmaxsmallregion=htmpresosvsq2[iRes]->GetYaxis()->FindBin(smallcutvalues[iCent]-0.0001);
        if(binmaxsmallregion==binmin) {
          binmin += 1; //avoid overlap between the two q2-regions
        }
      }
      else if(smallorlarge==kIntegrated) {
        binmin=1;
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(htmpresosvsq2[iRes]->GetYaxis()->GetNbins()+1); //takes also overflow entries
      }
      htmpresos[iRes]=(TH1F*)htmpresosvsq2[iRes]->ProjectionX(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()),binmin,binmax);
      if(htmpresos[iRes]){
        TString suffixcentrpart=Form("centr%0.1f_%0.1f",minCent+(iCent*centwidth),minCent+((iCent+1)*centwidth));
        htmpresos[iRes]->SetTitle(Form("Event Plane Resolution %s%s",namereso[iRes].Data(),suffixcentrpart.Data()));
        htmpresos[iRes]->SetName(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentrpart.Data()));
        if(useNcollWeight){
          cout << Form("Centr %0.f Bin %d  Ncoll %f\n",minCentTimesTen,ncentbins,ncoll[ncentbins]) << endl;
          htmpresos[iRes]->Scale(ncoll[ncentbins]);
        }
        outlist->Add(htmpresos[iRes]->Clone());
        if(iCent==0) {
          hevplresos[iRes]=(TH1F*)htmpresos[iRes]->Clone();
        }
        else {
          hevplresos[iRes]->Add(htmpresos[iRes]);
        }
      }
    }

    Double_t error;
    Double_t lowestRes=1;
    Double_t highestRes=0;
    Double_t resolBin=GetEventPlaneResolution(error,htmpresos[0],htmpresos[1],htmpresos[2]);

    if(resolBin<lowestRes) lowestRes=resolBin;
    if(resolBin>highestRes) highestRes=resolBin;

    Double_t binHalfWid=centwidth/2.;
    Double_t binCentr=(Double_t)minCent+(iCent*centwidth)+binHalfWid;
    gResolVsCent->SetPoint(iCent,binCentr,resolBin);
    gResolVsCent->SetPointError(iCent,binHalfWid,error);
  }

  for(Int_t iRes=0;iRes<nSubRes;iRes++){
    if(hevplresos[iRes]) {
      hevplresos[iRes]->SetTitle(Form("Event Plane Resolution %s%s",namereso[iRes].Data(),suffixcentr.Data()));
      hevplresos[iRes]->SetName(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
      outlist->Add(hevplresos[iRes]->Clone());
    }
  }
  gResolVsCent->SetName("gResolVsCent");
  gResolVsCent->SetTitle("Resolution vs. Centrality");
  outlist->Add(gResolVsCent->Clone());
  return outlist;
}

TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Int_t analysismeth, TGraphAsymmErrors *gRelSystEff) {

  TGraphAsymmErrors* gv2 = new TGraphAsymmErrors(nPtBins);

  if(analysismeth==kEventPlaneInOut) {
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
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
      gv2->SetPointError(iPt,averagePt[iPt]-PtLims[iPt],PtLims[iPt+1]-averagePt[iPt],ev2,ev2);
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
  else if(analysismeth==kEventPlane) {
    TF1 *flowFunc = new TF1("flow","[0]*(1.+2.*[1]*TMath::Cos(2.*x))");
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      //v2 from fit to Deltaphi distribution
      gSignal[iPt]->Fit(flowFunc);
      Double_t v2 = flowFunc->GetParameter(1)/resol;
      Double_t ev2=flowFunc->GetParError(1)/resol;
      gv2->SetPoint(iPt,averagePt[iPt],v2);
      gv2->SetPointError(iPt,averagePt[iPt]-PtLims[iPt],PtLims[iPt+1]-averagePt[iPt],ev2,ev2);
    }
    return gv2;
  }
  else {return 0x0;}
}

//_____________________________________________________________________________________________
//FILL SIGNAL GRAPHS
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, Int_t smallorlarge, Int_t analysismeth, Int_t cutmeth) {
  TString q2regionname="";
  Double_t q2cut=-1;
  Double_t q2perccut=-1.;
  if(smallorlarge==kSmall) {q2regionname="q2Small"; q2cut=q2smalllimit; q2perccut=q2smallpercevents*100;}
  else if(smallorlarge==kLarge) {q2regionname="q2Large"; q2cut=q2largelimit; q2perccut=q2largepercevents*100;}
  else if(smallorlarge==kIntegrated) {q2regionname="q2Int"; q2cut=0.; q2perccut=100;}

  TFile* infile = TFile::Open(infilename.Data(),"READ");
  TDirectoryFile* dir=0x0;
  TList* list=0x0;
  Double_t massD=-1;

  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());

  if(infile) {dir=(TDirectoryFile*)infile->Get(dirname.Data());}
  if(dir) {
    if(partname.Contains("Dzero")) {
      massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    }
    if(partname.Contains("Dplus")){
      massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
    }
    if(partname.Contains("Ds")) {
      massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
    }
    if(partname.Contains("Dstar")) {
      massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
    }
  }
  else {cerr << "Error: Wrong TDirectoryFile name " << dirname << ". Exit." << endl;}

  Int_t nPhi=nPhiBins;
  if(analysismeth==kEventPlaneInOut) nPhi=2;

  TH1F *hrflTempl=0x0;
  TH1F *hsigMC=0x0;
  Float_t sOverRef=0.;

  //Canvases for drawing histograms
  TCanvas *cDeltaPhi = new TCanvas(Form("cinvmassdeltaphi_%s",q2regionname.Data()),Form("Invariant mass distributions - %s",q2regionname.Data()),1920,1080);
  TCanvas *cDeltaPhifs = new TCanvas(Form("cinvmassdeltaphifs_%s",q2regionname.Data()),Form("Invariant mass distributions - fit with fixed sigma - %s",q2regionname.Data()),1920,1080);
  TCanvas *cPhiInteg = new TCanvas(Form("cinvmass_%s",q2regionname.Data()),Form("Invariant mass distributions - #phi and q_{2} integrated"),1920,1080);
  cDeltaPhi->Divide(nPtBins,nPhi);
  cDeltaPhifs->Divide(nPtBins,nPhi);
  Int_t nptpads = nPtBins/2;
  if(nPtBins%2!=0) nptpads = nPtBins/2+1;
  if(nPtBins>1) cPhiInteg->Divide(nptpads,2);
  Int_t nMassBins;
  Double_t hmin,hmax;
  for(Int_t iPt=0;iPt<nPtBins;iPt++){
    TH1F *histtofitfullsigma=(TH1F*)masslist->FindObject(Form("hMass_pt%d",iPt))->Clone();
    if(useTemplD0Refl){
      hrflTempl=(TH1F*)masslist->FindObject(Form("histRfl_%d",iPt));
      hsigMC=(TH1F*)masslist->FindObject(Form("histSgn_%d",iPt));
    }
    for(Int_t iPhi=0;iPhi<nPhi;iPhi++){
      Int_t ipad=(iPhi)*nPtBins+iPt+1;
      Double_t signal=0,esignal=0;
      TH1F *histtofit=(TH1F*)masslist->FindObject(Form("hMass_pt%d_phi%d",iPt,iPhi))->Clone();
      if(!histtofit){
        gSignal[iPt]->SetPoint(iPhi,iPhi,signal);
        gSignal[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
        return;
      }
      histtofit->SetTitle(Form("%.0f < #it{p}_{T} < %.0f, #phi%d",PtLims[iPt],PtLims[iPt+1],iPhi));
      nMassBins=histtofit->GetNbinsX();
      hmin=TMath::Max(minMassForFit[iPt],histtofit->GetBinLowEdge(2));
      hmax=TMath::Min(maxMassForFit[iPt],histtofit->GetBinLowEdge(nMassBins-2));
      histtofit->Rebin(rebin[iPt]);
      histtofit->GetYaxis()->SetTitle(Form("Entries / (%0.f MeV/c)",histtofit->GetBinWidth(5)*1000));
      AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(histtofit,hmin,hmax,1,typeb,types);
      if(useTemplD0Refl){
        fitter->SetTemplateReflections(hrflTempl);
        sOverRef=(hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999)))/(hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999)));
        fitter->SetFixReflOverS(sOverRef,kTRUE);
      }
      fitter->SetInitialGaussianMean(massD);
      fitter->SetInitialGaussianSigma(0.012);
      if (partname.Contains("Dstar")) {
	fitter->SetInitialGaussianSigma(0.0004);
      }
      fitter->SetUseLikelihoodFit();
      Bool_t ok=fitter->MassFitter(kFALSE);
      Double_t sigmaforcounting=0;
      Double_t meanforcounting=0;
      if(ok){
        fitter->DrawHere(cDeltaPhi->cd(ipad),3,1);
        signal = fitter->GetRawYield();
        esignal = fitter->GetRawYieldError();
        sigmaforcounting=fitter->GetSigma();
        meanforcounting=fitter->GetMean();
      }
      gSignal[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignal[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
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
      gSignalBC1[iPt]->SetPoint(iPhi,iPhi,cntSig1);
      gSignalBC1[iPt]->SetPointError(iPhi,0,0,cntErr,cntErr);
      gSignalBC2[iPt]->SetPoint(iPhi,iPhi,cntSig2);
      gSignalBC2[iPt]->SetPointError(iPhi,0,0,cntErr,cntErr);
    }
    //fit for fixed sigma
    histtofitfullsigma->SetTitle(Form("%.0f < #it{p}_{T} < %.0f GeV/c",PtLims[iPt],PtLims[iPt+1]));
    histtofitfullsigma->GetXaxis()->SetTitleSize(0.05);
    nMassBins=histtofitfullsigma->GetNbinsX();
    hmin=TMath::Max(minMassForFit[iPt],histtofitfullsigma->GetBinLowEdge(2));
    hmax=TMath::Min(maxMassForFit[iPt],histtofitfullsigma->GetBinLowEdge(nMassBins-2));
    histtofitfullsigma->Rebin(rebin[iPt]);
    histtofitfullsigma->GetYaxis()->SetTitle(Form("Entries / (%0.f MeV/c)",histtofitfullsigma->GetBinWidth(5)*1000));
    AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(histtofitfullsigma,hmin,hmax,1,typeb,types);
    if(useTemplD0Refl){
      fitter->SetTemplateReflections(hrflTempl);
      sOverRef=hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999))/hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999));
      fitter->SetFixReflOverS(sOverRef,kTRUE);
    }
    fitter->SetInitialGaussianMean(massD);
    if (partname.Contains("Dstar")) fitter->SetInitialGaussianSigma(0.0004);
    fitter->SetUseLikelihoodFit();
    Bool_t ok=fitter->MassFitter(kFALSE);
    if(ok){
      if(nPtBins==1) {fitter->DrawHere(cPhiInteg->cd(),3,1);}
      else {fitter->DrawHere(cPhiInteg->cd(iPt+1),3,1);}
    }
    Double_t sigma=fitter->GetSigma();
    Double_t massFromFit=fitter->GetMean();
    for(Int_t iPhi=0;iPhi<nPhi;iPhi++){
      Int_t ipad=(iPhi)*nPtBins+iPt+1;
      TH1F *histtofit=(TH1F*)masslist->FindObject(Form("hMass_pt%d_phi%d",iPt,iPhi))->Clone();
      histtofit->SetTitle(Form("%.f<#it{p}_{T}<%.f, #phi%d",PtLims[iPt],PtLims[iPt+1],iPhi));
      hmin=TMath::Max(minMassForFit[iPt],histtofit->GetBinLowEdge(2));
      hmax=TMath::Min(maxMassForFit[iPt],histtofit->GetBinLowEdge(nMassBins-2));
      nMassBins=histtofit->GetNbinsX();
      histtofit->Rebin(rebin[iPt]);
      histtofit->GetYaxis()->SetTitle(Form("Entries / (%0.f MeV/c)",histtofit->GetBinWidth(5)*1000));
      AliHFMassFitterVAR* fitter2=new AliHFMassFitterVAR(histtofit,hmin,hmax,1,typeb,types);
      if(useTemplD0Refl){
        fitter2->SetTemplateReflections(hrflTempl);
        sOverRef=hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999))/hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999));
        fitter2->SetFixReflOverS(sOverRef,kTRUE);
      }
      fitter2->SetInitialGaussianMean(massD);
      if (partname.Contains("Dstar")) fitter2->SetInitialGaussianSigma(0.0004);
      fitter2->SetFixGaussianSigma(sigma);
      if(fixAlsoMass) fitter2->SetFixGaussianMean(massFromFit);
      fitter2->SetUseLikelihoodFit();
      Bool_t ok2=fitter2->MassFitter(kFALSE);
      Double_t signal=0,esignal=0;
      if(ok2){
        fitter2->DrawHere(cDeltaPhifs->cd(ipad),3,1);
        signal = fitter2->GetRawYield();
        esignal = fitter2->GetRawYieldError();
      }
      gSignalfs[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignalfs[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
    }
  }//end loop on pt bin

  TString outnames[4] = {Form("%s/InvMassDeltaPhi%s_%s%.2f.pdf",outputdir.Data(),suffix.Data(),q2regionname.Data(),q2cut),
    Form("%s/InvMassDeltaPhi_fs%s_%s%.2f.pdf",outputdir.Data(),suffix.Data(),q2regionname.Data(),q2cut),
    Form("%s/InvMassDeltaPhi_fs%s_%s%.2f.root",outputdir.Data(),suffix.Data(),q2regionname.Data(),q2cut),
    Form("%s/InvMassfullphi%s.pdf",outputdir.Data(),suffix.Data())};

  if(smallorlarge==kIntegrated) {
    for(Int_t iFile=0; iFile<3; iFile++) {
      outnames[iFile].ReplaceAll("0.00","");
    }
  }

  TString percsuffix="perc";
  if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {
    outnames[0] = Form("%s/InvMassDeltaPhi%s_%s_%0.f%s.pdf",outputdir.Data(),suffix.Data(),q2regionname.Data(),q2perccut,percsuffix.Data());
    outnames[1] = Form("%s/InvMassDeltaPhi_fs%s_%s_%0.f%s.pdf",outputdir.Data(),suffix.Data(),q2regionname.Data(),q2perccut,percsuffix.Data());
    outnames[2] = Form("%s/InvMassDeltaPhi_fs%s_%s_%0.f%s.root",outputdir.Data(),suffix.Data(),q2regionname.Data(),q2perccut,percsuffix.Data());
  }

  cDeltaPhi->SaveAs(outnames[0].Data());
  cDeltaPhifs->SaveAs(outnames[1].Data());
  cDeltaPhifs->SaveAs(outnames[2].Data());
  cPhiInteg->SaveAs(outnames[3].Data());
}

//_____________________________________________________________________________________________
//GET EV PLANE RESOLUTION
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
    TH1F* hevplresos[3] = {hsubev1,hsubev2,hsubev3};
    for(Int_t iRes=0;iRes<3;iRes++){
      resolSub[iRes]=hevplresos[iRes]->GetMean();
      errors[iRes]=hevplresos[iRes]->GetMeanError();
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

//_____________________________________________________________________________________________
//q2 CUT APPLICATION
void Applyq2Cut(THnSparseF* sparse, Double_t smallcutvalues, Double_t largecutvalues, Int_t smallorlarge) {

  UInt_t q2axnum=3;

  TAxis* q2axis = (TAxis*)sparse->GetAxis(q2axnum);
  Int_t binmin = -1;
  Int_t binmax = -1;
  //define q2 cut bins
  if(smallorlarge==kSmall) {
    binmin=q2axis->FindBin(0.0001);
    binmax=q2axis->FindBin(smallcutvalues-0.0001);
  }
  else if(smallorlarge==kLarge) {
    binmin=q2axis->FindBin(largecutvalues+0.0001);
    binmax=q2axis->FindBin(q2axis->GetNbins()+1); //takes also overflow entries
    Int_t binmaxsmallregion=q2axis->FindBin(smallcutvalues-0.0001);
    if(binmaxsmallregion==binmin) {
      binmin += 1;
    }
  }
  q2axis->SetRange(binmin,binmax);
}

//_____________________________________________________________________________________________
//CUT APPLICATION
void ApplyCut(THnSparseF* sparse, Double_t min, Double_t max, UInt_t axnum) {

  TAxis* axis = (TAxis*)sparse->GetAxis(axnum);
  Int_t binmin = axis->FindBin(min+0.0001);
  Int_t binmax = axis->FindBin(max-0.0001);
  axis->SetRange(binmin,binmax);
}

//_____________________________________________________________________________________________
//DEFINE q2 CUTS
Bool_t Defineq2Cuts(TH2F* hq2VsCentr, vector<Double_t> &smallcutvalues, vector<Double_t> &largecutvalues, Int_t cutmeth, Int_t fSparseVers) {

  const Int_t ncentbins=hq2VsCentr->GetXaxis()->GetNbins();
  Double_t cutvalues[2] = {-1.,-1.};
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    if(cutmeth==kAbsCut) {smallcutvalues.push_back(q2smalllimit); largecutvalues.push_back(q2largelimit);}
    else if(cutmeth==kPercCut || (cutmeth==kPercCutVsCent && fSparseVers==0)) {
      TH1F* hq2=(TH1F*)hq2VsCentr->ProjectionY();
      Bool_t cuttaken = Getq2CutValuePercEvents(hq2,cutvalues);
      if(cuttaken) {
        smallcutvalues.push_back(cutvalues[0]);
        largecutvalues.push_back(cutvalues[1]);
      }
      else {return kFALSE;}
    }
    else if(cutmeth==kPercCutVsCent && fSparseVers==1) {
      TH1F* hq2centbin=(TH1F*)hq2VsCentr->ProjectionY(Form("hq2centbin%d",iCent+1),iCent+1,iCent+1);
      Bool_t cuttaken = Getq2CutValuePercEvents(hq2centbin,cutvalues);
      if(cuttaken) {
        smallcutvalues.push_back(cutvalues[0]);
        largecutvalues.push_back(cutvalues[1]);
      }
      else {return kFALSE;}
    }
  }

  cout << "Cut values for the small-q2 region:"<<endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    cout << "Centrality bin "<< iCent << " q2 < " << smallcutvalues[iCent] <<endl;
  }
  cout << "\nCut values for the large-q2 region:"<<endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    cout << "Centrality bin "<< iCent << " q2 > " << largecutvalues[iCent] <<endl;
  }
  cout << endl;

  return kTRUE;
}

//_____________________________________________________________________________________________
//q2 CUT VALUE SEARCH (ON A BASIS OF PERCENTAGE OF EVENTS WITH SMALLER/LARGER q2)
Bool_t Getq2CutValuePercEvents(TH1F* hq2, Double_t cutvalues[2]) {

  Int_t bincutvalues[2]={-1,-1};

  if(hq2) {
    Double_t SmallThresIntegral = q2smallpercevents*hq2->Integral(0.,1000.); //threshold integral value for small-q2 region (including overflow entries)
    Double_t LargeThresIntegral = q2largepercevents*hq2->Integral(0.,1000.); //threshold integral value for large-q2 region (including overflow entries)
    Double_t integral=0;
    Int_t q2bin=0;
    while(integral<SmallThresIntegral) {
      q2bin++;
      integral += hq2->GetBinContent(q2bin);
    }
    if(q2smallpercevents!=1.0) {bincutvalues[0]=q2bin-1;}//rounded down
    else {bincutvalues[0]=q2bin;}//exactly 100%
    if(q2smallpercevents==0.5 && q2largepercevents==0.5) bincutvalues[1]=bincutvalues[0]+1;
    else {
      integral=0;
      q2bin=hq2->GetNbinsX()+2;
      while(integral<LargeThresIntegral) {
        q2bin--;
        integral += hq2->GetBinContent(q2bin); //takes also overflow counts
      }
      if(q2largepercevents!=1.0) {bincutvalues[1]=q2bin+1;}//rounded down
      else {bincutvalues[1]=q2bin;}//exactly 100%
    }

    cutvalues[0]=hq2->GetBinLowEdge(bincutvalues[0])+hq2->GetBinWidth(bincutvalues[0]);
    cutvalues[1]=hq2->GetBinLowEdge(bincutvalues[1]);

    return kTRUE;
  }

  cout << "Warning: hq2 not set, impossible to define the cut values!" << endl;
  return kFALSE;
}

//_____________________________________________________________________________________________
//RESET AXES OF SPARSE
void ResetAxes(THnSparseF* sparse, Int_t axnum) {
  if(axnum<0) {
    for(Int_t iAxis=0; iAxis<sparse->GetNdimensions(); iAxis++) {
      sparse->GetAxis(iAxis)->SetRange(-1,-1);
    }
  }
  else {sparse->GetAxis(axnum)->SetRange(-1,-1);}
}

//_____________________________________________________________________________________________
//DRAWING STYLE
void SetStyle() {
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetStatFont(42);
  gStyle->SetStatY(0.89);
  gStyle->SetStatX(0.89);
  gStyle->SetTitleFont(42,"xyzg");
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPalette(53);
}

//_____________________________________________________________________________________________
//LOAD REFLECTION HISTOS FOR D0
Bool_t LoadD0toKpiMCHistos(TList *outlist){

  TFile *f=new TFile(fileNameMCD0refl.Data(),"READ");
  if(!f){
    printf("ERROR: file %s does not exist\n",fileNameMCD0refl.Data());
    return kFALSE;
  }
  f->ls();
  TH1F** hsig=new TH1F*[nPtBins];
  TH1F** hrfl=new TH1F*[nPtBins];
  for(Int_t j=0;j<nPtBins;j++){
    hsig[j]=(TH1F*)f->Get(Form("histSgn_%d",j));
    if(!hsig[j]) {Printf("histSgn_%d NOT FOUND",j); return kFALSE;}
    hrfl[j]=(TH1F*)f->Get(Form("histRflFitted%s_ptBin%d",rflFitType.Data(),j));
    if(!hrfl[j]) {Printf("histRflFitted%s_ptBin%d",rflFitType.Data(),j); return kFALSE;}
  }
  for(Int_t k=0;k<nPtBins;k++){
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
