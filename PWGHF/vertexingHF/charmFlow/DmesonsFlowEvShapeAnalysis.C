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
#include "AliHFVnVsMassFitter.h"

#endif

//methods for the analysis of AliAnalysisTaskSEHFv2/vn output in case of kEvShape method
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
const TString infilename = "$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/AnalysisResults_v2_1030_3050_1050_EvShape_VZERO_q2TPC_RemAllDau_NtrklDist.root";
const TString suffix="_3050_Topod0Cut_QoverM_q2TPC_q2RecalcRemAllDau_VZERO_EvShape";

const TString partname="Dplus";

//centrality
const Int_t ncentbins=20;
const Int_t minCent=30;
const Int_t maxCent=50;

const TString outputdir="Cent3050/EvShape/RemAllDau/q2Distr";

//ptbins of the analysis
const Int_t nPtBins=4;
const Int_t nPtLims=nPtBins+1;
const Double_t PtLims[nPtLims] = {3.,4.,6.,8.,12.};

//phi bins
const Int_t nPhiBins=4;
const Int_t nPhiLims=nPhiBins+1;
const Double_t PhiLims[nPhiLims]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};

//q2 cut values (absolute cut)
const Double_t q2smalllimit=2.2;
const Double_t q2largelimit=3.2;

//percentage of events with smaller/larger q2 (for kPercCut, kPercCutVsCent, kPercentileCut)
const Double_t q2smallpercevents=0.6;
const Double_t q2largepercevents=0.2;

//EP resolution
const Int_t evPlane=AliEventPlaneResolutionHandler::kVZERO;
Int_t evPlaneRes=AliEventPlaneResolutionHandler::kThreeSub;
const Bool_t useNcollWeight=kFALSE;
const Bool_t useAliHandlerForRes=kFALSE;

// mass fit configuration
//10-30
const Int_t rebin[]={3,3,4,4};
const Double_t minMassForFit[]={1.70,1.70,1.70,1.69};
const Double_t maxMassForFit[]={2.05,2.02,2.02,2.02};

//30-50
//const Int_t rebin[]={3,3,4,4,4,5,5,6,6,6};
//const Double_t minMassForFit[]={1.72,1.70,1.70,1.70,1.70};
//const Double_t maxMassForFit[]={2.02,2.05,2.05,2.05,2.05};

enum {kGaus=0, kDoubleGaus, kReflTempl};
const Int_t types=kGaus;//kReflTempl;
const Int_t typeb=AliHFMassFitter::kExpo; //Background: 0=expo, 1=linear, 2=pol2
Bool_t useTemplD0Refl=kFALSE;
TString rflFitType="DoubleGaus";
TString fileNameMCD0refl="./reflections/reflections_fitted_DoubleGaus.root";
const Bool_t fixAlsoMass=kFALSE;
const Double_t nSigmaForCounting=3.5;

//f(phi_D) vs. mass
const Int_t types_mass_simfit=AliHFVnVsMassFitter::kGaus;
const Int_t typeb_mass_simfit=AliHFVnVsMassFitter::kExpo;
const Int_t typeb_fphiD_simfit=AliHFVnVsMassFitter::kLin;

//not to be set
enum CutMethod{kAbsCut,kPercCut,kPercCutVsCent,kPercentileCut}; //kAbsCut->absolute cut values, kPercCut->cut according to the % of events with smaller/larger q2, kPercCutVsCent->cut according to the % of events with smaller/larger q2 in finer centrality bins, kPercentileCut->cut on percentile (if enabled in the task)
enum SmallOrLarge{kSmall,kLarge,kIntegrated};
enum AnalysisMethod{kEventPlane,kEventPlaneInOut,kScalarProd};

//in-out efficiency
const Double_t effInOverEffOut=1.03;

const Int_t colors[]={kRed,kBlue,kGray};
const Int_t markers[]={kFullSquare,kFullCircle,kFullDiamond,kOpenSquare,kOpenCircle,kOpenDiamond};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
Int_t DmesonsFlowEvShapeAnalysis(Int_t cutmeth=kPercentileCut, Int_t analysismeth=kEventPlaneInOut);
Int_t EvaluatePhiDModulations(Int_t cutmeth=kPercentileCut);
Int_t GetPhiDDistribution(Int_t cutmeth=kPercentileCut, Int_t nSigmaMinForSB=4, Int_t nSigmaMaxForSB=10);
Int_t Drawq2VsCent(Int_t cutmeth=kPercentileCut);
void DrawEventPlaneResolutionAndDistribution(Int_t cutmeth=kPercentileCut);
TList* LoadTList();
THnSparseF* LoadSparseFromList(TList* inputlist);
TH2F* GetHistoq2VsCentr(TList* inputlist);
TList* LoadMassHistos(THnSparseF* sparse, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge, Int_t analysismeth);
TList* LoadResolutionHistos(TList *inputlist, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge);
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Int_t analysismeth, TGraphAsymmErrors *gRelSystEff);
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, TH1F* hMean[], TH1F* hMeanfs[], TH1F* hSigmaFree[], TH1F* hSigmaFixed[],Int_t smallorlarge, Int_t analysismeth, Int_t cutmeth);
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

  TList* datalist=(TList*)LoadTList();
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
  else if(cutmeth==kPercentileCut) percsuffix="percentile";
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent || cutmeth==kPercentileCut) {outname=Form("%s/v2Output_%d_%d_%s%s_q2Small%0.f%s_q2Large%0.f%s.root",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
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
        delete htmp;
        htmp=0x0;
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
          if (hmax > 0.165) hmax=0.165;
        }
        histtofit->Rebin(rebin[iPt]);
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
        if (partname.Contains("Dstar")) {
          fitter->SetInitialGaussianMean(0.145);
          fitter->SetInitialGaussianSigma(0.0004);
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
      TH1F* hSigmaFixed[nPhi];
      TH1F* hSigmaFree[nPhi];
      TH1F* hMean[nPhi];
      TH1F* hMeanfs[nPhi];
      for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
        hSigmaFree[iPhi] = new TH1F(Form("hSigmaFree_phi%d_%s",iPhi,q2regionname[iq2].Data()),";#it{p}_{T} (GeV/c);width (GeV/c^{2})",nPtBins,PtLims);
        hSigmaFree[iPhi]->SetMarkerStyle(markers[3]);
        hSigmaFree[iPhi]->SetMarkerColor(colors[iq2]+2*iPhi);
        hSigmaFree[iPhi]->SetLineColor(colors[iq2]+2*iPhi);
        hSigmaFree[iPhi]->SetDirectory(0);
        hSigmaFixed[iPhi] = new TH1F(Form("hSigmaFixed_phi%d_%s",iPhi,q2regionname[iq2].Data()),";#it{p}_{T} (GeV/c);width (GeV/c^{2})",nPtBins,PtLims);
        hSigmaFixed[iPhi]->SetMarkerStyle(markers[1]);
        hSigmaFixed[iPhi]->SetMarkerColor(colors[iq2]+2*iPhi);
        hSigmaFixed[iPhi]->SetLineColor(colors[iq2]+2*iPhi);
        hSigmaFixed[iPhi]->SetDirectory(0);
        hMean[iPhi] = new TH1F(Form("hMean_phi%d_%s",iPhi,q2regionname[iq2].Data()),";#it{p}_{T} (GeV/c);mean (GeV/c^{2})",nPtBins,PtLims);
        hMean[iPhi]->SetMarkerStyle(markers[3]);
        hMean[iPhi]->SetMarkerColor(colors[iq2]+2*iPhi);
        hMean[iPhi]->SetLineColor(colors[iq2]+2*iPhi);
        hMean[iPhi]->SetDirectory(0);
        hMeanfs[iPhi] = new TH1F(Form("hMeanfs_phi%d_%s",iPhi,q2regionname[iq2].Data()),";#it{p}_{T} (GeV/c);mean (GeV/c^{2})",nPtBins,PtLims);
        hMeanfs[iPhi]->SetMarkerStyle(markers[1]);
        hMeanfs[iPhi]->SetMarkerColor(colors[iq2]+2*iPhi);
        hMeanfs[iPhi]->SetLineColor(colors[iq2]+2*iPhi);
        hMeanfs[iPhi]->SetDirectory(0);
      }
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
      FillSignalGraph(masslist,gSignal,gSignalfs,gSignalBC1,gSignalBC2,hMean,hMeanfs,hSigmaFree,hSigmaFixed,q2region[iq2],analysismeth,cutmeth);
      outfile.cd();
      for(Int_t iPt=0;iPt<nPtBins;iPt++){
        gSignal[iPt]->Write();
        gSignalfs[iPt]->Write();
        gSignalBC1[iPt]->Write();
        gSignalBC2[iPt]->Write();
      }
      for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
        hSigmaFixed[iPhi]->Write();
        hSigmaFree[iPhi]->Write();
        hMean[iPhi]->Write();
        hMeanfs[iPhi]->Write();
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
      TString q2name = hq2VsCentr->GetYaxis()->GetTitle();
      q2name.ReplaceAll("(%) ","");
      if(iq2!=kIntegrated) {
        if(cutmeth==kAbsCut) {
          TString sign[2] = {"<",">"};
          title=Form("%s %s %0.1f",q2name.Data(),sign[iq2].Data(),q2cut[iq2]);
        }
        else {
          Double_t perc[2] = {q2smallpercevents*100,q2largepercevents*100};
          TString reg[2] = {"small-","large-"};
          title=Form("%0.f%% %s %s",perc[iq2],reg[iq2].Data(),q2name.Data());
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
      if(cutmeth==kPercCut || cutmeth==kPercCutVsCent || cutmeth==kPercentileCut) {outname=Form("%s/v2Output_%d_%d_%s%s_%s%.f%s.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[iq2].Data(),q2perccut[iq2],percsuffix.Data());}
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
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent || cutmeth==kPercentileCut) {outname=Form("%s/v2fsOutput_%d_%d_%s%s_q2Small%.0f%s_q2Large%.0f%s.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
  cv2fs->SaveAs(outname.Data());

  return 0;
}

Int_t EvaluatePhiDModulations(Int_t cutmeth) {

  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}

  TList* datalist=(TList*)LoadTList();
  if(!datalist){return 1;}

  Double_t massD=-1;
  if(partname.Contains("Dzero")) {
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  else if(partname.Contains("Dplus")){
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
  else if(partname.Contains("Dstar")) {
    massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
  }
  else if(partname.Contains("Ds") && !partname.Contains("Dstar")) {
    massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
  }

  THnSparseF* hMassPtPhiq2Centr=(THnSparseF*)LoadSparseFromList(datalist);
  if(!hMassPtPhiq2Centr) {return 2;}
  Int_t fSparseVers = 0; //check version of the task output
  if(hMassPtPhiq2Centr->GetAxis(5)) {
    cout << "Version of the sparse with sin(2phiD) and cos(2phiD) axes!" << endl;
    fSparseVers=1;
  }
  else {
    cout << "Version of the sparse without sin(2phiD) and cos(2phiD) axes!" << endl;
    return 3;
  }
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) {return 4;}

  Double_t centwidth=(Double_t)(maxCent-minCent)/ncentbins;
  vector<Double_t> smallcutvalues;
  vector<Double_t> largecutvalues;
  Bool_t defq2cut=Defineq2Cuts(hq2VsCentr,smallcutvalues,largecutvalues,cutmeth,fSparseVers);
  if(!defq2cut) {return 5;}

  Int_t q2region[3] = {kSmall,kLarge,kIntegrated};
  TString q2regionname[3] = {"q2Small","q2Large","q2Int"};

  SetStyle();

  TGraphAsymmErrors* gSin2PhiD[3];
  TGraphAsymmErrors* gCos2PhiD[3];
  TCanvas* cSin2PhiDVsMass[3][nPtBins];
  TCanvas* cCos2PhiDVsMass[3][nPtBins];
  for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {

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

    Double_t sigmaFromFit[nPtBins];
    gSin2PhiD[iq2] = new TGraphAsymmErrors(0);
    gSin2PhiD[iq2]->SetName(Form("gSin2PhiD_%s",q2regionname[iq2].Data()));
    gSin2PhiD[iq2]->SetTitle("");
    gSin2PhiD[iq2]->SetLineColor(colors[iq2]+1);
    gSin2PhiD[iq2]->SetLineWidth(2);
    gSin2PhiD[iq2]->SetMarkerColor(colors[iq2]+1);
    gSin2PhiD[iq2]->SetMarkerStyle(kFullCircle);
    gCos2PhiD[iq2] = new TGraphAsymmErrors(0);
    gCos2PhiD[iq2]->SetName(Form("gCos2PhiD_%s",q2regionname[iq2].Data()));
    gCos2PhiD[iq2]->SetTitle("");
    gCos2PhiD[iq2]->SetLineColor(colors[iq2]+1);
    gCos2PhiD[iq2]->SetLineWidth(2);
    gCos2PhiD[iq2]->SetMarkerColor(colors[iq2]+1);
    gCos2PhiD[iq2]->SetMarkerStyle(kFullCircle);
    //load resolution and mass lists
    TList* masslist=(TList*)LoadMassHistos(hMassPtPhiq2Centr,smallcutvalues,largecutvalues,q2region[iq2],kEventPlaneInOut);

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
      delete htmp;
      htmp=0x0;
    }
    Float_t averagePt[nPtBins];
    Float_t errorPt[nPtBins];
    for(Int_t iPt=0;iPt<nPtBins;iPt++){
      Int_t ptbinMin=hmasspt->FindBin(PtLims[iPt]);
      Int_t ptbinMax=hmasspt->FindBin(PtLims[iPt+1]-0.001);
      if(TMath::Abs(hmasspt->GetXaxis()->GetBinLowEdge(ptbinMin)-PtLims[iPt])>0.001 ||
         TMath::Abs(hmasspt->GetXaxis()->GetBinUpEdge(ptbinMax)-PtLims[iPt+1])>0.001){
        cout << "Error in pt bin limits for projection!\n" << endl;
        return 7;
      }
      TH1F *hmass = (TH1F*)hmasspt->ProjectionY("_py",ptbinMin,ptbinMax);
      hmass->Rebin(rebin[iPt]);
      Int_t nMassBins=hmass->GetNbinsX();
      Double_t hmin=hmass->GetBinLowEdge(2); // need wide range for <pt>
      Double_t hmax=hmass->GetBinLowEdge(nMassBins-2); // need wide range for <pt>
      AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(hmass,hmin,hmax,1,typeb,types);
      if (partname.Contains("Dstar")) {
        if (hmin < 0.140) hmin=0.140;
        if (hmax > 0.175) hmax=0.175;
      }
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
      sigmaFromFit[iPt]=fitter->GetSigma();
      TF1* funcB2=fitter->GetBackgroundRecalcFunc();
      utils->AveragePt(averagePt[iPt],errorPt[iPt],PtLims[iPt],PtLims[iPt+1],hmasspt,massFromFit,sigmaFromFit[iPt],funcB2,2.5,4.5,0.,3.,1);
      if(averagePt[iPt]==0 || errorPt[iPt]>1) averagePt[iPt] = (PtLims[iPt+1]+PtLims[iPt])/2;

      //PhiD modulations (simultaneus fit)
      TH2F* hMassSin2PhiD=0x0;
      TH2F* hMassCos2PhiD=0x0;
      for(Int_t iCent=0; iCent<ncentbins; iCent++) {
        ApplyCut(hMassPtPhiq2Centr,PtLims[iPt],PtLims[iPt+1],1); //apply pt cut
        ApplyCut(hMassPtPhiq2Centr,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,4); //apply centrality cut (to have a q2 cut centrality dependent)
        if(iq2!=kIntegrated) {Applyq2Cut(hMassPtPhiq2Centr,smallcutvalues[iCent],largecutvalues[iCent],q2region[iq2]);}
        hMassPtPhiq2Centr->GetAxis(1)->SetRange(ptbinMin,ptbinMax);
        TH2F* htmpCos2PhiD=(TH2F*)hMassPtPhiq2Centr->Projection(0,5);
        TH2F* htmpSin2PhiD=(TH2F*)hMassPtPhiq2Centr->Projection(0,6);
        if(iCent==0) {
          hMassCos2PhiD=(TH2F*)htmpCos2PhiD->Clone();
          hMassSin2PhiD=(TH2F*)htmpSin2PhiD->Clone();
        }
        else {
          hMassCos2PhiD->Add(htmpCos2PhiD);
          hMassSin2PhiD->Add(htmpSin2PhiD);
        }
        ResetAxes(hMassPtPhiq2Centr,-1.);
        delete htmpCos2PhiD;
        htmpCos2PhiD=0x0;
        delete htmpSin2PhiD;
        htmpSin2PhiD=0x0;
      }
      TH1F* hMassForFit = (TH1F*)hMassCos2PhiD->ProjectionY(Form("hMass_Pt%d_%s",iPt,q2regionname[iq2].Data()));
      hMassForFit->Rebin(rebin[iPt]);
      TAxis* massaxis = (TAxis*)hMassPtPhiq2Centr->GetAxis(0);
      TAxis* massaxisreb = (TAxis*)hMassForFit->GetXaxis();
      Int_t nRebMassBins = massaxisreb->GetNbins();
      Double_t minmass = massaxisreb->GetBinLowEdge(1);
      Double_t maxmass = massaxisreb->GetBinLowEdge(nRebMassBins)+massaxisreb->GetBinWidth(nRebMassBins);
      TH1F* hSin2PhiDVsMass=new TH1F(Form("hSin2PhiDVsMass_Pt%d_%s",iPt,q2regionname[iq2].Data()),";M_{K#pi#pi} (GeV/c);< Sin(2#varphi) > / R_{2}",nRebMassBins,minmass,maxmass);
      TH1F* hCos2PhiDVsMass=new TH1F(Form("hCos2PhiDVsMass_Pt%d_%s",iPt,q2regionname[iq2].Data()),";M_{K#pi#pi} (GeV/c);< Cos(2#varphi) > / R_{2}",nRebMassBins,minmass,maxmass);
      for(Int_t iMass=0; iMass<nRebMassBins; iMass++){
        minmass = massaxisreb->GetBinLowEdge(iMass+1);
        maxmass = minmass+massaxisreb->GetBinWidth(iMass+1);
        Int_t binmassmin = massaxis->FindBin(minmass*1.0001);
        Int_t binmassmax = massaxis->FindBin(maxmass*0.9999);
        TH1F* hSin2PhiD = (TH1F*)hMassSin2PhiD->ProjectionX("hSin2PhiD",binmassmin,binmassmax);
        TH1F* hCos2PhiD = (TH1F*)hMassCos2PhiD->ProjectionX("hCos2PhiD",binmassmin,binmassmax);
        Double_t meansin2phiD = hSin2PhiD->GetMean();
        Double_t meansin2phiDerr = hSin2PhiD->GetMeanError();
        Double_t meancos2phiD = hCos2PhiD->GetMean();
        Double_t meancos2phiDerr = hCos2PhiD->GetMeanError();
        hSin2PhiDVsMass->SetBinContent(iMass+1,meansin2phiD/resol);
        hSin2PhiDVsMass->SetBinError(iMass+1,meansin2phiDerr/resol);
        hCos2PhiDVsMass->SetBinContent(iMass+1,meancos2phiD/resol);
        hCos2PhiDVsMass->SetBinError(iMass+1,meancos2phiDerr/resol);
      }
      AliHFVnVsMassFitter* fitterSin2PhiD = new AliHFVnVsMassFitter(hMassForFit,hSin2PhiDVsMass,minMassForFit[iPt],maxMassForFit[iPt],typeb_mass_simfit,types_mass_simfit,typeb_fphiD_simfit);
      fitterSin2PhiD->SetInitialGaussianMean(massD);
      fitterSin2PhiD->SetInitialGaussianSigma(sigmaFromFit[iPt]);
      fitterSin2PhiD->FixMeanFromMassFit();
      fitterSin2PhiD->FixSigmaFromMassFit();
      fitterSin2PhiD->SimultaneusFit(kFALSE);
      cSin2PhiDVsMass[iq2][iPt] = new TCanvas(Form("cSin2PhiDVsMass_Pt%d_%s",iPt,q2regionname[iq2].Data()),"",800,800);
      fitterSin2PhiD->DrawHere(cSin2PhiDVsMass[iq2][iPt]);
      Double_t sin2phiD = fitterSin2PhiD->GetVn();
      Double_t sin2phiDerr = fitterSin2PhiD->GetVnUncertainty();

      AliHFVnVsMassFitter* fitterCos2PhiD = new AliHFVnVsMassFitter(hMassForFit,hCos2PhiDVsMass,minMassForFit[iPt],maxMassForFit[iPt],typeb_mass_simfit,types_mass_simfit,typeb_fphiD_simfit);
      fitterCos2PhiD->SetInitialGaussianMean(massD);
      fitterCos2PhiD->SetInitialGaussianSigma(sigmaFromFit[iPt]);
      fitterCos2PhiD->FixMeanFromMassFit();
      fitterCos2PhiD->FixSigmaFromMassFit();
      fitterCos2PhiD->SimultaneusFit(kFALSE);
      cCos2PhiDVsMass[iq2][iPt] = new TCanvas(Form("cCos2PhiDVsMass_Pt%d_%s",iPt,q2regionname[iq2].Data()),"",800,800);
      fitterCos2PhiD->DrawHere(cCos2PhiDVsMass[iq2][iPt]);
      Double_t cos2phiD = fitterCos2PhiD->GetVn();
      Double_t cos2phiDerr = fitterCos2PhiD->GetVnUncertainty();

      gSin2PhiD[iq2]->SetPoint(iPt,averagePt[iPt],sin2phiD);
      gSin2PhiD[iq2]->SetPointError(iPt,averagePt[iPt]-PtLims[iPt],PtLims[iPt+1]-averagePt[iPt],sin2phiDerr,sin2phiDerr);
      gCos2PhiD[iq2]->SetPoint(iPt,averagePt[iPt],cos2phiD);
      gCos2PhiD[iq2]->SetPointError(iPt,averagePt[iPt]-PtLims[iPt],PtLims[iPt+1]-averagePt[iPt],cos2phiDerr,cos2phiDerr);
    }
    cout << Form("Average pt %s region \n",q2regionname[iq2].Data()) << endl;
    for(Int_t iPt=0;iPt<nPtBins;iPt++) cout <<Form("%f +- %f\n",averagePt[iPt],errorPt[iPt])<<endl;
  }

  TLine* lineatzero = new TLine(PtLims[0],0.,PtLims[nPtBins],0.);
  lineatzero->SetLineWidth(2);
  lineatzero->SetLineColor(kBlack);
  lineatzero->SetLineStyle(9);

  TLegend* leg = new TLegend(0.5,0.2,0.8,0.35);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gSin2PhiD[0],Form("%.0f%% Small-q_{2}",q2smallpercevents*100),"lpe");
  leg->AddEntry(gSin2PhiD[1],Form("%.0f%% Large-q_{2}",q2largepercevents*100),"lpe");
  leg->AddEntry(gSin2PhiD[2],"q_{2}-integrated","lpe");

  TCanvas* cCos2PhiD = new TCanvas("cCos2PhiD","",1920,1080);
  TString drawopt="AP";
  for(Int_t iq2=kIntegrated; iq2>=kSmall; iq2--) {
    gCos2PhiD[iq2]->GetYaxis()->SetRangeUser(-0.35,0.35);
    gCos2PhiD[iq2]->GetXaxis()->SetRangeUser(PtLims[0],PtLims[nPtBins]);
    gCos2PhiD[iq2]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gCos2PhiD[iq2]->GetYaxis()->SetTitle("< Cos(2#varphi_{D}) > / R_{2}");
    gCos2PhiD[iq2]->Draw(drawopt.Data());
    drawopt="P";
  }
  lineatzero->Draw("same");
  leg->Draw("same");
  TCanvas* cSin2PhiD = new TCanvas("cSin2PhiD","",1920,1080);
  drawopt="AP";
  for(Int_t iq2=kIntegrated; iq2>=kSmall; iq2--) {
    gSin2PhiD[iq2]->GetYaxis()->SetRangeUser(-0.35,0.35);
    gSin2PhiD[iq2]->GetXaxis()->SetRangeUser(PtLims[0],PtLims[nPtBins]);
    gSin2PhiD[iq2]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gSin2PhiD[iq2]->GetYaxis()->SetTitle("< Sin(2#varphi_{D}) > / R_{2}");
    gSin2PhiD[iq2]->Draw(drawopt.Data());
    drawopt="P";
  }
  lineatzero->Draw("same");
  leg->Draw("same");

  TString percsuffix="perc";
  if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
  else if(cutmeth==kPercentileCut) percsuffix="percentile";
  TString outname=Form("%s/Phi2D_modulations_%d_%d_%s%s_q2Small%.2f_q2Large%.2f.root",outputdir.Data(),minCent,maxCent,"InOut",suffix.Data(),q2smalllimit,q2largelimit);
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/Phi2D_modulations_%d_%d_%s%s_q2Small%.0f%s_q2Large%.0f%s.pdf",outputdir.Data(),minCent,maxCent,"InOut",suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
  TFile outfile(outname.Data(),"RECREATE");
  for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      cSin2PhiDVsMass[iq2][iPt]->Write();
      cSin2PhiDVsMass[iq2][iPt]->Write();
    }
    gCos2PhiD[iq2]->Write();
    gSin2PhiD[iq2]->Write();
  }
  cCos2PhiD->Write();
  cSin2PhiD->Write();
  outfile.Close();
  outname.ReplaceAll(".root",".pdf");
  outname.ReplaceAll("Phi2D","CosPhi2D");
  cCos2PhiD->SaveAs(outname.Data());
  outname.ReplaceAll("CosPhi2D","SinPhi2D");
  cSin2PhiD->SaveAs(outname.Data());

  return 0;
}
//_____________________________________________________________________________________________
Int_t GetPhiDDistribution(Int_t cutmeth, Int_t nSigmaMinForSB, Int_t nSigmaMaxForSB) {
  
  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}
  
  TList* datalist=(TList*)LoadTList();
  if(!datalist){return 1;}
  
  Double_t massD=-1;
  if(partname.Contains("Dzero")) {
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  else if(partname.Contains("Dplus")){
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
  else if(partname.Contains("Dstar")) {
    massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
  }
  else if(partname.Contains("Ds") && !partname.Contains("Dstar")) {
    massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
  }
  
  THnSparseF* hMassPtPhiq2Centr=(THnSparseF*)LoadSparseFromList(datalist);
  if(!hMassPtPhiq2Centr) {return 2;}
  Int_t fSparseVers = 0; //check version of the task output
  if(hMassPtPhiq2Centr->GetAxis(7)) {
    cout << "Version of the sparse with phiD axis!" << endl;
    fSparseVers=1;
  }
  else {
    cout << "Version of the sparse without phiD axis!" << endl;
    return 3;
  }
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) {return 4;}
  
  Double_t centwidth=(Double_t)(maxCent-minCent)/ncentbins;
  vector<Double_t> smallcutvalues;
  vector<Double_t> largecutvalues;
  Bool_t defq2cut=Defineq2Cuts(hq2VsCentr,smallcutvalues,largecutvalues,cutmeth,fSparseVers);
  if(!defq2cut) {return 5;}
  
  Int_t q2region[3] = {kSmall,kLarge,kIntegrated};
  TString q2regionname[3] = {"q2Small","q2Large","q2Int"};
  
  SetStyle();
  
  TH1F* hPhiDDistrAll[3];
  TH1F* hPhiDDistrAllNorm[3];
  TH1F* hPhiDDistrSBleft[3];
  TH1F* hPhiDDistrSBright[3];
  TH1F* hPhiDDistrBkg[3];
  TH1F* hPhiDDistrBkgNorm[3];
  TH1F* hPhiDDistrBkgScaled[3];
  TH1F* hPhiDDistrSignal[3];
  TH1F* hPhiDDistrSignalNorm[3];
  TH2F* hMassPhiD[3];
  TCanvas* cMassFit[3];
  
  for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
    //PhiD distributions in q2 classes
    ApplyCut(hMassPtPhiq2Centr,PtLims[0],PtLims[nPtBins],1); //apply pt cut -> Full pt range
    for(Int_t iCent=0; iCent<ncentbins; iCent++) {
      ApplyCut(hMassPtPhiq2Centr,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,4); //apply centrality cut (to have a q2 cut centrality dependent)
      if(iq2!=kIntegrated) {Applyq2Cut(hMassPtPhiq2Centr,smallcutvalues[iCent],largecutvalues[iCent],q2region[iq2]);}
      TH2F* htmpMassPhiD=(TH2F*)hMassPtPhiq2Centr->Projection(0,7);
      if(iCent==0) {
        hMassPhiD[iq2]=(TH2F*)htmpMassPhiD->Clone(Form("hMassPhiD_%s",q2regionname[iq2].Data()));
      }
      else {
        hMassPhiD[iq2]->Add(htmpMassPhiD);
      }
      delete htmpMassPhiD;
      htmpMassPhiD=0x0;
      ResetAxes(hMassPtPhiq2Centr,4);
    }
    TH1F* hMassForFit = (TH1F*)hMassPhiD[iq2]->ProjectionY(Form("hMass_%s",q2regionname[iq2].Data()));
    //hMassForFit->Rebin(rebin[0]);//--> first rebin
    Double_t hmin=hMassForFit->GetBinLowEdge(2); // need wide range for <pt>
    Double_t hmax=hMassForFit->GetBinLowEdge(hMassForFit->GetNbinsX()-2); // need wide range for <pt>
    if (partname.Contains("Dstar")) {
      if (hmin < 0.140) hmin=0.140;
      if (hmax > 0.175) hmax=0.175;
    }
    AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(hMassForFit,hmin,hmax,1,typeb,types);
    /* reflections for D0 to be added --> reflections for full pt range?
     if(useTemplD0Refl){
     hrflTempl=(TH1F*)masslist->FindObject(Form("histRfl_%d",iPt));
     hsigMC=(TH1F*)masslist->FindObject(Form("histSgn_%d",iPt));
     fitter->SetTemplateReflections(hrflTempl);
     sOverRef=(hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999)))/(hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999)));
     fitter->SetFixReflOverS(sOverRef,kTRUE);
     }
     */
    fitter->SetInitialGaussianMean(massD);
    fitter->SetInitialGaussianSigma(0.012);
    if (partname.Contains("Dstar")) {
      fitter->SetInitialGaussianSigma(0.0004);
    }
    fitter->SetUseLikelihoodFit();
    fitter->MassFitter(kFALSE);
    Double_t sigma = fitter->GetSigma();
    Double_t mean = fitter->GetMean();
    Double_t bkg, bkgerr;
    fitter->Background(3,bkg,bkgerr);
    cMassFit[iq2] = new TCanvas(Form("cMassFit_%s",q2regionname[iq2].Data()),"",800,800);
    fitter->DrawHere(gPad);
    
    Double_t minMassSBleft = mean-nSigmaMaxForSB*sigma;
    Double_t maxMassSBleft = mean-nSigmaMinForSB*sigma;
    Double_t minMassSBright = mean+nSigmaMinForSB*sigma;
    Double_t maxMassSBright = mean+nSigmaMaxForSB*sigma;
    Double_t minMassSignal = mean-3*sigma;
    Double_t maxMassSignal = mean+3*sigma;
    
    TBox* leftBox = new TBox(minMassSBleft,hMassForFit->GetMaximum()*0.4,maxMassSBleft,hMassForFit->GetMaximum()*0.8);
    leftBox->SetLineColor(kRed);
    leftBox->SetLineWidth(2);
    leftBox->SetFillStyle(0);
    leftBox->Draw("same");
    TBox* rightBox = new TBox(minMassSBright,hMassForFit->GetMaximum()*0.2,maxMassSBright,hMassForFit->GetMaximum()*0.6);
    rightBox->SetLineColor(kRed);
    rightBox->SetLineWidth(2);
    rightBox->SetFillStyle(0);
    rightBox->Draw("same");
    
    TAxis* massaxis = (TAxis*)hMassPhiD[iq2]->GetYaxis();
    Int_t binminMassSBleft = massaxis->FindBin(minMassSBleft*1.0001);
    Int_t binmaxMassSBleft = massaxis->FindBin(maxMassSBleft*0.9999);
    Int_t binminMassSBright = massaxis->FindBin(minMassSBright*1.0001);
    Int_t binmaxMassSBright = massaxis->FindBin(maxMassSBright*0.9999);
    Int_t binminMassSignal = massaxis->FindBin(minMassSignal*1.0001);
    Int_t binmaxMassSignal = massaxis->FindBin(maxMassSignal*0.9999);
    
    hPhiDDistrAll[iq2]=(TH1F*)hMassPhiD[iq2]->ProjectionX(Form("hPhiDDistrAll_%s",q2regionname[iq2].Data()),binminMassSignal,binmaxMassSignal);
    hPhiDDistrAll[iq2]->SetTitle("All candidates in mass range");
    hPhiDDistrAllNorm[iq2]=(TH1F*)hPhiDDistrAll[iq2]->Clone(Form("hPhiDDistrAllNorm_%s",q2regionname[iq2].Data()));
    hPhiDDistrAllNorm[iq2]->Sumw2();
    hPhiDDistrAllNorm[iq2]->Scale(1./hPhiDDistrAllNorm[iq2]->Integral());
    hPhiDDistrAllNorm[iq2]->SetLineColor(colors[iq2]);
    hPhiDDistrAllNorm[iq2]->GetYaxis()->SetTitle("Normalised entries");
    
    hPhiDDistrSBleft[iq2]=(TH1F*)hMassPhiD[iq2]->ProjectionX(Form("hPhiDDistrSBleft_%s",q2regionname[iq2].Data()),binminMassSBleft,binmaxMassSBleft);
    hPhiDDistrSBright[iq2]=(TH1F*)hMassPhiD[iq2]->ProjectionX(Form("hPhiDDistrSBright_%s",q2regionname[iq2].Data()),binminMassSBright,binmaxMassSBright);
    hPhiDDistrBkg[iq2]=(TH1F*)hPhiDDistrSBleft[iq2]->Clone(Form("hPhiDDistrBkg_%s",q2regionname[iq2].Data()));
    hPhiDDistrBkg[iq2]->SetTitle("Background");
    hPhiDDistrBkg[iq2]->Add(hPhiDDistrSBright[iq2]);
    hPhiDDistrBkgNorm[iq2]=(TH1F*)hPhiDDistrBkg[iq2]->Clone(Form("hPhiDDistrBkgNorm_%s",q2regionname[iq2].Data()));
    hPhiDDistrBkgNorm[iq2]->Sumw2();
    hPhiDDistrBkgNorm[iq2]->Scale(1./hPhiDDistrBkgNorm[iq2]->Integral());
    hPhiDDistrBkgNorm[iq2]->SetLineColor(colors[iq2]);
    hPhiDDistrBkgNorm[iq2]->GetYaxis()->SetTitle("Normalised entries");
    hPhiDDistrBkgScaled[iq2]=(TH1F*)hPhiDDistrBkgNorm[iq2]->Clone(Form("hPhiDDistrBkgScaled_%s",q2regionname[iq2].Data()));
    hPhiDDistrBkgScaled[iq2]->Scale(bkg);
    
    hPhiDDistrSignal[iq2]=(TH1F*)hPhiDDistrAll[iq2]->Clone(Form("hPhiDDistrSignal_%s",q2regionname[iq2].Data()));
    hPhiDDistrSignal[iq2]->SetTitle("Signal");
    hPhiDDistrSignal[iq2]->Add(hPhiDDistrBkgScaled[iq2],-1);
    hPhiDDistrSignalNorm[iq2]=(TH1F*)hPhiDDistrSignal[iq2]->Clone(Form("hPhiDDistrSignalNorm_%s",q2regionname[iq2].Data()));
    hPhiDDistrSignalNorm[iq2]->Sumw2();
    hPhiDDistrSignalNorm[iq2]->Scale(1./hPhiDDistrSignalNorm[iq2]->Integral());
    hPhiDDistrSignalNorm[iq2]->SetLineColor(colors[iq2]);
    hPhiDDistrSignalNorm[iq2]->GetYaxis()->SetTitle("Normalised entries");
    
    ResetAxes(hMassPtPhiq2Centr,-1);
  }
  
  TLegend* leg = new TLegend(0.5,0.75,0.8,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hPhiDDistrAllNorm[0],"Small-q_{2}","l");
  leg->AddEntry(hPhiDDistrAllNorm[1],"Large-q_{2}","l");
  leg->AddEntry(hPhiDDistrAllNorm[2],"q_{2}-integrated","l");
  
  TCanvas *cPhiDAll = new TCanvas("cPhiDAll","",800,800);
  hPhiDDistrAllNorm[0]->GetYaxis()->SetRangeUser(0.,0.15);
  hPhiDDistrAllNorm[0]->Draw();
  hPhiDDistrAllNorm[1]->Draw("same");
  hPhiDDistrAllNorm[2]->Draw("same");
  leg->Draw("same");
  
  TCanvas *cPhiDBkg = new TCanvas("cPhiDBkg","",800,800);
  hPhiDDistrBkgNorm[0]->GetYaxis()->SetRangeUser(0.,0.15);
  hPhiDDistrBkgNorm[0]->Draw();
  hPhiDDistrBkgNorm[1]->Draw("same");
  hPhiDDistrBkgNorm[2]->Draw("same");
  leg->Draw("same");
  
  TCanvas *cPhiDSignal = new TCanvas("cPhiDSignal","",800,800);
  hPhiDDistrSignalNorm[0]->GetYaxis()->SetRangeUser(0.,0.15);
  hPhiDDistrSignalNorm[0]->Draw();
  hPhiDDistrSignalNorm[1]->Draw("same");
  hPhiDDistrSignalNorm[2]->Draw("same");
  leg->Draw("same");
  
  TH1F* hPhiDDistrAllRatio[2];
  TH1F* hPhiDDistrBkgRatio[2];
  TH1F* hPhiDDistrSignalRatio[2];
  TF1* fCosAll[2];
  TF1* fCosBkg[2];
  TF1* fCosSignal[2];
  for(Int_t iq2=kSmall; iq2<=kLarge; iq2++) {
    hPhiDDistrAllRatio[iq2]=(TH1F*)hPhiDDistrAllNorm[iq2]->Clone(Form("hPhiDDistrAllRatio_%s",q2regionname[iq2].Data()));
    hPhiDDistrAllRatio[iq2]->Divide(hPhiDDistrAllNorm[iq2],hPhiDDistrAllNorm[2],1.,1.,"B");
    hPhiDDistrAllRatio[iq2]->GetYaxis()->SetTitle("Ratio");
    hPhiDDistrBkgRatio[iq2]=(TH1F*)hPhiDDistrBkgNorm[iq2]->Clone(Form("hPhiDDistrBkgRatio_%s",q2regionname[iq2].Data()));
    hPhiDDistrBkgRatio[iq2]->Divide(hPhiDDistrBkgNorm[iq2],hPhiDDistrBkgNorm[2],1.,1.,"B");
    hPhiDDistrBkgRatio[iq2]->GetYaxis()->SetTitle("Ratio");
    hPhiDDistrSignalRatio[iq2]=(TH1F*)hPhiDDistrSignalNorm[iq2]->Clone(Form("hPhiDDistrSignalRatio_%s",q2regionname[iq2].Data()));
    hPhiDDistrSignalRatio[iq2]->Divide(hPhiDDistrSignalNorm[iq2],hPhiDDistrSignalNorm[2],1.,1.,"B");
    hPhiDDistrSignalRatio[iq2]->GetYaxis()->SetTitle("Ratio");
    fCosAll[iq2] = new TF1(Form("fCosAll_%s",q2regionname[iq2].Data()),"[0]*(1+[1]*TMath::Cos([2]*x+[3]))",0.,2*TMath::Pi());
    fCosAll[iq2]->SetLineColor(colors[iq2]);
    fCosAll[iq2]->SetParLimits(2,1.,5.);
    fCosBkg[iq2] = new TF1(Form("fCosBkg_%s",q2regionname[iq2].Data()),"[0]*(1+[1]*TMath::Cos([2]*x+[3]))",0.,2*TMath::Pi());
    fCosBkg[iq2]->SetLineColor(colors[iq2]);
    fCosBkg[iq2]->SetParLimits(2,1.,5.);
    fCosSignal[iq2] = new TF1(Form("fCosSignal_%s",q2regionname[iq2].Data()),"[0]*(1+[1]*TMath::Cos([2]*x+[3]))",0.,2*TMath::Pi());
    fCosSignal[iq2]->SetParLimits(2,1.,5.);
    fCosSignal[iq2]->SetLineColor(colors[iq2]);
  }
  TLegend* legRatio = new TLegend(0.3,0.75,0.8,0.85);
  legRatio->SetFillStyle(0);
  legRatio->SetBorderSize(0);
  legRatio->SetTextSize(0.04);
  legRatio->AddEntry(hPhiDDistrAllRatio[0],"Small-q_{2} / q_{2}-integrated","l");
  legRatio->AddEntry(hPhiDDistrBkgRatio[1],"Large-q_{2} / q_{2}-integrated","l");
  
  TCanvas *cPhiDAllRatio = new TCanvas("cPhiDAllRatio","",800,800);
  hPhiDDistrAllRatio[0]->GetYaxis()->SetRangeUser(0.,2.);
  hPhiDDistrAllRatio[0]->Draw();
  hPhiDDistrAllRatio[0]->Fit(fCosAll[0]);
  hPhiDDistrAllRatio[1]->Draw("same");
  hPhiDDistrAllRatio[1]->Fit(fCosAll[1]);
  legRatio->Draw("same");
  
  TCanvas *cPhiDBkgRatio = new TCanvas("cPhiDBkgRatio","",800,800);
  hPhiDDistrBkgRatio[0]->GetYaxis()->SetRangeUser(0.,2.);
  hPhiDDistrBkgRatio[0]->Draw();
  hPhiDDistrBkgRatio[0]->Fit(fCosBkg[0]);
  hPhiDDistrBkgRatio[1]->Draw("same");
  hPhiDDistrBkgRatio[1]->Fit(fCosBkg[1]);
  legRatio->Draw("same");
  
  TCanvas *cPhiDSignalRatio = new TCanvas("cPhiDSignalRatio","",800,800);
  hPhiDDistrSignalRatio[0]->GetYaxis()->SetRangeUser(0.,2.);
  hPhiDDistrSignalRatio[0]->Draw();
  hPhiDDistrSignalRatio[0]->Fit(fCosSignal[0]);
  hPhiDDistrSignalRatio[1]->Draw("same");
  hPhiDDistrSignalRatio[1]->Fit(fCosSignal[1]);
  legRatio->Draw("same");
  
  TString percsuffix="perc";
  if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
  else if(cutmeth==kPercentileCut) percsuffix="percentile";
  TString outname=Form("%s/PhiD_distribution_%d_%d_%s%s_q2Small%.2f_q2Large%.2f.root",outputdir.Data(),minCent,maxCent,"InOut",suffix.Data(),q2smalllimit,q2largelimit);
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outname=Form("%s/PhiD_distribution_%d_%d_%s%s_q2Small%.0f%s_q2Large%.0f%s.root",outputdir.Data(),minCent,maxCent,"InOut",suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
  TFile outfile(outname.Data(),"RECREATE");
  for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
    hPhiDDistrAllNorm[iq2]->Write();
    hPhiDDistrBkgNorm[iq2]->Write();
    hPhiDDistrSignalNorm[iq2]->Write();
    if(iq2<kIntegrated) {
      hPhiDDistrAllRatio[iq2]->Write();
      hPhiDDistrBkgRatio[iq2]->Write();
      hPhiDDistrSignalRatio[iq2]->Write();
      fCosAll[iq2]->Write();
      fCosBkg[iq2]->Write();
      fCosSignal[iq2]->Write();
    }
  }
  cPhiDAll->Write();
  cPhiDBkg->Write();
  cPhiDSignal->Write();
  cPhiDAllRatio->Write();
  cPhiDBkgRatio->Write();
  cPhiDSignalRatio->Write();
  outfile.Close();
  
  outname.ReplaceAll(".root",".pdf");
  outname.ReplaceAll("PhiD_distribution","PhiD_distribution_Signal");
  cPhiDSignal->SaveAs(outname.Data());
  outname.ReplaceAll("PhiD_distribution_Signal","PhiD_distribution_Bkg");
  cPhiDBkg->SaveAs(outname.Data());
  outname.ReplaceAll("PhiD_distribution_Bkg","PhiD_distribution_All");
  cPhiDAll->SaveAs(outname.Data());
  outname.ReplaceAll("PhiD_distribution_All","PhiD_distribution_SignalRatio");
  cPhiDSignalRatio->SaveAs(outname.Data());
  outname.ReplaceAll("PhiD_distribution_SignalRatio","PhiD_distribution_BkgRatio");
  cPhiDBkgRatio->SaveAs(outname.Data());
  outname.ReplaceAll("PhiD_distribution_BkgRatio","PhiD_distribution_AllRatio");
  cPhiDAllRatio->SaveAs(outname.Data());
  
  return 0;
}

//_____________________________________________________________________________________________
//DRAW q_2 VS CENTRALITY
Int_t Drawq2VsCent(Int_t cutmeth) {
  
  TGaxis::SetMaxDigits(3);
  
  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}

  SetStyle();

  TList* datalist=(TList*)LoadTList();
  if(!datalist){return 1;}
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) {return 2;}

  Double_t centwidth=(Double_t)(maxCent-minCent)/ncentbins;
  vector<Double_t> smallcutvalues;
  vector<Double_t> largecutvalues;
  Bool_t defq2cut=Defineq2Cuts(hq2VsCentr,smallcutvalues,largecutvalues,cutmeth,1);
  if(!defq2cut) {return 3;}

  if(ncentbins!=(Int_t)smallcutvalues.size() || ncentbins!=(Int_t)largecutvalues.size()) {
    cerr << "Number of centrality bins not consistent. Exit." << endl;
    return 4;
  }

  TH1F* hq2=(TH1F*)hq2VsCentr->ProjectionY();
  TString q2name = hq2->GetXaxis()->GetTitle();
  TString q2namewoperc = q2name;
  q2namewoperc.ReplaceAll(" (%)","");
  hq2->SetTitle("");
  hq2->GetYaxis()->SetTitle("Entries");
  hq2->SetLineColor(kBlack);
  hq2->SetLineWidth(2);
  hq2->SetMarkerStyle(kFullCircle);

  TH1F** hq2centbin = new TH1F*[ncentbins];

  TH1F* hSmallLimitVsCent = new TH1F("hSmallLimitVsCent","",ncentbins,minCent,maxCent);
  TH1F* hLargeLimitVsCent = new TH1F("hLargeLimitVsCent","",ncentbins,minCent,maxCent);
  hSmallLimitVsCent->SetLineStyle(7);
  hLargeLimitVsCent->SetLineStyle(9);
  hSmallLimitVsCent->SetLineWidth(3);
  hLargeLimitVsCent->SetLineWidth(3);

  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    hq2centbin[iCent]=(TH1F*)hq2VsCentr->ProjectionY(Form("hq2centbin%d",iCent),iCent+1,iCent+1);
    hq2centbin[iCent]->GetYaxis()->SetTitle("Entries");
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
  legVsCentLarge->SetFillStyle(0);

  TCanvas *cq2VsCent = new TCanvas("cq2VsCent","q_{2} vs. centrality",800,800);
  cq2VsCent->SetLogz();
  hSmallLimitVsCent->SetLineColor(kBlack);
  hLargeLimitVsCent->SetLineColor(kBlack);
  legVsCentSmall->AddEntry(hSmallLimitVsCent,Form("Small-%s",q2namewoperc.Data()),"l");
  legVsCentLarge->AddEntry(hLargeLimitVsCent,Form("Large-%s",q2namewoperc.Data()),"l");
  hq2VsCentr->Draw("colz");
  hSmallLimitVsCent->Draw("same");
  hLargeLimitVsCent->Draw("same");
  legVsCentSmall->Draw("same");
  legVsCentLarge->Draw("same");

  TCanvas *cq2 = 0x0;
  if(cutmeth!=kPercentileCut) {
    cq2=new TCanvas("cq2","q_{2}",800,800);
    cq2->SetLogy();
    hq2->Draw("E");
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
  }
  
  TCanvas* cq2Cut = 0x0;

  if(cutmeth==kAbsCut || cutmeth==kPercCut || cutmeth==kPercentileCut) {
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
    if(cutmeth!=kPercentileCut) cq2Cut->SetLogy();
    else hq2->GetYaxis()->SetRangeUser(0.,hq2->GetMaximum()*1.5);
    hq2->Draw("E");
    hq2smallcut->Draw("same");
    hq2largecut->Draw("same");
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
    Int_t smallthresholdbin=hq2centbin[iCent]->GetXaxis()->FindBin(smallcutvalues[iCent]-0.0001);
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
  }
  cout << endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    Double_t totnum=0.;
    Int_t largethresholdbin=hq2centbin[iCent]->GetXaxis()->FindBin(largecutvalues[iCent]+0.0001);
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
  }
  cout << endl;

  TH2F* hq2SmallVsCentr = (TH2F*)hq2VsCentr->Clone("hq2SmallVsCentr");
  hq2SmallVsCentr->Clear();
  TH2F* hq2LargeVsCentr = (TH2F*)hq2VsCentr->Clone("hq2LargeVsCentr");
  hq2LargeVsCentr->Clear();

  TAxis* q2axis = (TAxis*)hq2VsCentr->GetYaxis();
  for(Int_t iCent=0; iCent<hq2VsCentr->GetXaxis()->GetNbins(); iCent++) {
    Int_t smallthresholdbin=q2axis->FindBin(smallcutvalues[iCent]-0.0001);
    Int_t largethresholdbin=q2axis->FindBin(largecutvalues[iCent]+0.0001);
    for(Int_t iq2=0; iq2<hq2VsCentr->GetYaxis()->GetNbins(); iq2++) {
      if(iq2+1<=smallthresholdbin) {
        hq2SmallVsCentr->SetBinContent(iCent+1,iq2+1,hq2VsCentr->GetBinContent(iCent+1,iq2+1));
        hq2LargeVsCentr->SetBinContent(iCent+1,iq2+1,0.);
      }
      else if(iq2+1>=largethresholdbin) {
        hq2LargeVsCentr->SetBinContent(iCent+1,iq2+1,hq2VsCentr->GetBinContent(iCent+1,iq2+1));
        hq2SmallVsCentr->SetBinContent(iCent+1,iq2+1,0.);
      }
      else {
        hq2LargeVsCentr->SetBinContent(iCent+1,iq2+1,0.);
        hq2SmallVsCentr->SetBinContent(iCent+1,iq2+1,0.);
      }
    }
  }

  TCanvas* cq2LargeVsCent = new TCanvas("cq2LargeVsCent","",800,800);
  hq2LargeVsCentr->Draw("colz");
  hLargeLimitVsCent->Draw("same");
  legVsCentLarge->Draw("same");

  TCanvas* cq2SmallVsCent = new TCanvas("cq2SmallVsCent","",800,800);
  hq2SmallVsCentr->Draw("colz");
  hSmallLimitVsCent->Draw("same");
  legVsCentSmall->Draw("same");

  TH1F* hCentLargeq2=(TH1F*)hq2LargeVsCentr->ProjectionX();
  hCentLargeq2->SetName("hCentLargeq2");
  hCentLargeq2->GetYaxis()->SetTitle("Normalised entries");
  hCentLargeq2->Sumw2();
  hCentLargeq2->Scale(1./hCentLargeq2->Integral());
  hCentLargeq2->SetLineColor(colors[1]+1);
  TH1F* hCentSmallq2=(TH1F*)hq2SmallVsCentr->ProjectionX();
  hCentSmallq2->SetName("hCentSmallq2");
  hCentSmallq2->GetYaxis()->SetTitle("Normalised entries");
  hCentSmallq2->Sumw2();
  hCentSmallq2->Scale(1./hCentSmallq2->Integral());
  hCentSmallq2->SetLineColor(colors[0]+1);

  TLegend* leg2 = new TLegend(0.55,0.7,0.89,0.89);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.045);
  leg2->AddEntry(hCentSmallq2,Form("Small-%s",q2namewoperc.Data()),"l");
  leg2->AddEntry(hCentLargeq2,Form("Large-%s",q2namewoperc.Data()),"l");

  TCanvas* cCent = new TCanvas("cCent","",800,800);
  hCentSmallq2->GetYaxis()->SetRangeUser(hCentSmallq2->GetMinimum()*0.9,hCentSmallq2->GetMaximum()*1.1);
  hCentSmallq2->Draw();
  hCentLargeq2->Draw("same");
  leg2->Draw("same");

  if(cutmeth!=kPercentileCut) {
    cq2->SaveAs(Form("%s/Nev_vs_q2%s.pdf",outputdir.Data(),suffix.Data()));
  }
  
  TString outname=Form("%s/q2Selection%s_q2Small%0.2f_q2Large%0.2f.pdf",outputdir.Data(),suffix.Data(),q2smalllimit,q2largelimit);
  TString percsuffix="perc";
  if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
  else if(cutmeth==kPercentileCut) percsuffix="percentile";

  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent || cutmeth==kPercentileCut) {outname=Form("%s/q2Selection%s_q2Small%0.f%s_q2Large%0.f%s.pdf",outputdir.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}
  if(cutmeth!=kPercCutVsCent) {
    cq2Cut->SaveAs(outname.Data());
  }
  outname.ReplaceAll("q2Selection","q2SelectionVsCent");
  cq2CutVsCent->SaveAs(outname.Data());
  outname.ReplaceAll("q2SelectionVsCent","q2_vs_Centrality");
  cq2VsCent->SaveAs(outname.Data());
  outname.ReplaceAll(".pdf",".root");

  TFile outfile(outname.Data(),"RECREATE");
  if(cutmeth!=kPercentileCut) {cq2->Write();}
  cq2CutVsCent->Write();
  cq2VsCent->Write();
  cq2LargeVsCent->Write();
  cq2SmallVsCent->Write();
  cCent->Write();
  outfile.Close();

  TCanvas* cNtrkl[3];
  TCanvas* cNtrklRatios[3];
  TCanvas* cNtrkVsq2[3];
  TCanvas* cMeanNtrkRatioVsq2[3];
  TH3F* hNtrklVsq2VsCent[3];
  TH2F* hNtrkVsq2[3];
  TH1F* hMeanNtrkl[3];
  TH1F* hRatioMeanNtrkl[3];
  TH1F* hRatioMeanNtrkl_Largeq2[3];
  TH1F* hRatioMeanNtrkl_Smallq2[3];

  hNtrklVsq2VsCent[0] = (TH3F*)datalist->FindObject("hNtrklVsq2VsCent");
  if(hNtrklVsq2VsCent[0]) {
    hNtrklVsq2VsCent[1] = (TH3F*)datalist->FindObject("hNtrklVsq2VsCentCand");
    hNtrklVsq2VsCent[2] = (TH3F*)datalist->FindObject("hNtrklVsq2VsCentCandInMass");
    TH1F* hNtrklUnbiased[3];
    TH1F* hNtrklq2Small[3];
    TH1F* hNtrklq2Large[3];
    TH1F* hNtrklRatioq2SmallUnbiased[3];
    TH1F* hNtrklRatioq2LargeUnbiased[3];
    TLegend* leg3 = new TLegend(0.55,0.65,0.85,0.85);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.04);

    TString Ntrklnames[3] = {"","Cand","CandInMass"};

    for(Int_t iHisto=0; iHisto<3; iHisto++) {
      hNtrklUnbiased[iHisto] = (TH1F*)hNtrklVsq2VsCent[iHisto]->ProjectionZ();

      for(Int_t iCent=0; iCent<hNtrklVsq2VsCent[iHisto]->GetXaxis()->GetNbins(); iCent++) {
        Int_t smallthresholdbin=q2axis->FindBin(smallcutvalues[iCent]-0.0001);
        Int_t largethresholdbin=q2axis->FindBin(largecutvalues[iCent]+0.0001);
        if(iCent==0) {
          hNtrklq2Small[iHisto] = (TH1F*)hNtrklVsq2VsCent[iHisto]->ProjectionZ(Form("hNtrklq2Small%s",Ntrklnames[iHisto].Data()),iCent+1,iCent+1,1,smallthresholdbin);
          hNtrklq2Large[iHisto] = (TH1F*)hNtrklVsq2VsCent[iHisto]->ProjectionZ(Form("hNtrklq2Large%s",Ntrklnames[iHisto].Data()),iCent+1,iCent+1,largethresholdbin,100000);
        }
        else {
          TH1F* hpart_q2Small = (TH1F*)hNtrklVsq2VsCent[iHisto]->ProjectionZ(Form("_partq2Small%s_Cent%d",Ntrklnames[iHisto].Data(),iCent),iCent+1,iCent+1,1,smallthresholdbin);
          TH1F* hpart_q2Large = (TH1F*)hNtrklVsq2VsCent[iHisto]->ProjectionZ(Form("_partq2Large%s_Cent%d",Ntrklnames[iHisto].Data(),iCent),iCent+1,iCent+1,largethresholdbin,100000);
          hNtrklq2Small[iHisto]->Add(hpart_q2Small);
          hNtrklq2Large[iHisto]->Add(hpart_q2Large);
        }
      }
      hNtrklUnbiased[iHisto]->SetName(Form("hNtrklUnbiased%s",Ntrklnames[iHisto].Data()));
      hNtrklUnbiased[iHisto]->SetTitle("");
      hNtrklUnbiased[iHisto]->GetYaxis()->SetTitle("Normalised entries");
      hNtrklUnbiased[iHisto]->GetXaxis()->SetTitle("N_{tracklets}");
      hNtrklUnbiased[iHisto]->GetXaxis()->SetTitleSize(0.05);
      hNtrklUnbiased[iHisto]->GetXaxis()->SetLabelSize(0.045);
      hNtrklUnbiased[iHisto]->Sumw2();
      hNtrklUnbiased[iHisto]->SetLineColor(colors[2]+1);
      hNtrklUnbiased[iHisto]->Scale(1./hNtrklUnbiased[iHisto]->Integral());
      hNtrklq2Small[iHisto]->SetTitle("");
      hNtrklq2Small[iHisto]->GetYaxis()->SetTitle("Normalised entries");
      hNtrklq2Small[iHisto]->GetXaxis()->SetTitle("N_{tracklets}");
      hNtrklq2Small[iHisto]->Sumw2();
      hNtrklq2Small[iHisto]->SetLineColor(colors[0]+1);
      hNtrklq2Small[iHisto]->Scale(1./hNtrklq2Small[iHisto]->Integral());
      hNtrklq2Small[iHisto]->GetXaxis()->SetTitleSize(0.05);
      hNtrklq2Small[iHisto]->GetXaxis()->SetLabelSize(0.045);
      hNtrklq2Large[iHisto]->SetTitle("");
      hNtrklq2Large[iHisto]->GetYaxis()->SetTitle("Normalised entries");
      hNtrklq2Large[iHisto]->GetXaxis()->SetTitle("N_{tracklets}");
      hNtrklq2Large[iHisto]->Sumw2();
      hNtrklq2Large[iHisto]->Scale(1./hNtrklq2Large[iHisto]->Integral());
      hNtrklq2Large[iHisto]->SetLineColor(colors[1]+1);
      hNtrklq2Large[iHisto]->GetXaxis()->SetTitleSize(0.05);
      hNtrklq2Large[iHisto]->GetXaxis()->SetLabelSize(0.045);

      if(iHisto==0) {
        leg3->AddEntry(hNtrklq2Small[iHisto],Form("Small-%s",q2namewoperc.Data()),"l");
        leg3->AddEntry(hNtrklq2Large[iHisto],Form("Large-%s",q2namewoperc.Data()),"l");
        leg3->AddEntry(hNtrklUnbiased[iHisto],Form("%s-integrated",q2namewoperc.Data()),"l");
      }

      hNtrklRatioq2SmallUnbiased[iHisto]=(TH1F*)hNtrklq2Small[iHisto]->Clone(Form("hNtrklRatioq2SmallUnbiased%s",Ntrklnames[iHisto].Data()));
      hNtrklRatioq2SmallUnbiased[iHisto]->GetYaxis()->SetTitle("Ratio w.r.t. unbiased");
      hNtrklRatioq2SmallUnbiased[iHisto]->Divide(hNtrklq2Small[iHisto],hNtrklUnbiased[iHisto],1.,1.,"B");
      hNtrklRatioq2SmallUnbiased[iHisto]->GetYaxis()->SetName("Weight");
      hNtrklRatioq2LargeUnbiased[iHisto]=(TH1F*)hNtrklq2Large[iHisto]->Clone(Form("hNtrklRatioq2LargeUnbiased%s",Ntrklnames[iHisto].Data()));
      hNtrklRatioq2LargeUnbiased[iHisto]->GetYaxis()->SetTitle("Ratio w.r.t. unbiased");
      hNtrklRatioq2LargeUnbiased[iHisto]->Divide(hNtrklq2Large[iHisto],hNtrklUnbiased[iHisto],1.,1.,"B");
      hNtrklRatioq2LargeUnbiased[iHisto]->GetYaxis()->SetName("Weight");

      cNtrkl[iHisto] = new TCanvas(Form("cNtrkl%s",Ntrklnames[iHisto].Data()),"",800,800);
      cNtrkl[iHisto]->SetLogy();
      hNtrklUnbiased[iHisto]->Draw();
      hNtrklq2Small[iHisto]->Draw("same");
      hNtrklq2Large[iHisto]->Draw("same");
      leg3->Draw("same");

      cNtrklRatios[iHisto] = new TCanvas(Form("cNtrklRatios%s",Ntrklnames[iHisto].Data()),"",800,800);
      hNtrklRatioq2SmallUnbiased[iHisto]->Draw();
      hNtrklRatioq2LargeUnbiased[iHisto]->Draw("same");
      leg2->Draw("same");

      TString outnameNtrkl = outname;
      outnameNtrkl.ReplaceAll(".root",".pdf");
      outnameNtrkl.ReplaceAll("q2_vs_Centrality",Form("NtrklDist%s",Ntrklnames[iHisto].Data()));
      cNtrkl[iHisto]->SaveAs(outnameNtrkl.Data());
      outnameNtrkl.ReplaceAll(Form("NtrklDist%s",Ntrklnames[iHisto].Data()),Form("NtrklWeights%s",Ntrklnames[iHisto].Data()));
      cNtrklRatios[iHisto]->SaveAs(outnameNtrkl.Data());
      
      if(cutmeth==kPercentileCut) {
        hNtrkVsq2[iHisto] = (TH2F*)hNtrklVsq2VsCent[iHisto]->Project3D("zy");
        hNtrkVsq2[iHisto]->SetTitle("");
        hNtrkVsq2[iHisto]->GetXaxis()->SetTitleSize(0.05);
        hNtrkVsq2[iHisto]->GetXaxis()->SetLabelSize(0.045);
        hNtrkVsq2[iHisto]->GetYaxis()->SetTitleSize(0.05);
        hNtrkVsq2[iHisto]->GetYaxis()->SetLabelSize(0.045);

        hMeanNtrkl[iHisto] = new TH1F(Form("hMeanNtrkl_%s",Ntrklnames[iHisto].Data()),"",100,0.,100.);
        hMeanNtrkl[iHisto]->GetYaxis()->SetTitle("<N_{tracklets}>");
        hMeanNtrkl[iHisto]->GetXaxis()->SetTitle(hNtrkVsq2[iHisto]->GetXaxis()->GetTitle());
        hMeanNtrkl[iHisto]->SetMarkerStyle(kFullCircle);
        hMeanNtrkl[iHisto]->SetMarkerColor(kBlack);
        hMeanNtrkl[iHisto]->SetLineColor(kBlack);
        hMeanNtrkl[iHisto]->SetLineWidth(2);
        
        hRatioMeanNtrkl[iHisto] = new TH1F(Form("hRatioMeanNtrkl_%s",Ntrklnames[iHisto].Data()),"",100,0.,100.);
        hRatioMeanNtrkl[iHisto]->GetYaxis()->SetTitle("<N_{tracklets}> ESE/unbiased");
        hRatioMeanNtrkl[iHisto]->GetXaxis()->SetTitle(hNtrkVsq2[iHisto]->GetXaxis()->GetTitle());
        hRatioMeanNtrkl[iHisto]->SetMarkerStyle(kFullCircle);
        hRatioMeanNtrkl[iHisto]->SetMarkerColor(kBlack);
        hRatioMeanNtrkl[iHisto]->SetLineColor(kBlack);
        hRatioMeanNtrkl[iHisto]->SetLineWidth(2);
        
        hRatioMeanNtrkl_Largeq2[iHisto] = new TH1F(Form("hRatioMeanNtrkl_Largeq2_%s",Ntrklnames[iHisto].Data()),"",1,(1-q2largepercevents)*100,100.);
        hRatioMeanNtrkl_Largeq2[iHisto]->GetYaxis()->SetTitle("<N_{tracklets}> ESE/unbiased");
        hRatioMeanNtrkl_Largeq2[iHisto]->GetXaxis()->SetTitle(hNtrkVsq2[iHisto]->GetXaxis()->GetTitle());
        hRatioMeanNtrkl_Largeq2[iHisto]->SetMarkerStyle(kFullDiamond);
        hRatioMeanNtrkl_Largeq2[iHisto]->SetMarkerColor(colors[1]+1);
        hRatioMeanNtrkl_Largeq2[iHisto]->SetLineColor(colors[1]+1);
        hRatioMeanNtrkl_Largeq2[iHisto]->SetLineWidth(2);
        hRatioMeanNtrkl_Largeq2[iHisto]->SetBinContent(1,hNtrklq2Large[iHisto]->GetMean()/hNtrklUnbiased[iHisto]->GetMean());
        hRatioMeanNtrkl_Largeq2[iHisto]->SetBinError(1,1.e-10);
        
        hRatioMeanNtrkl_Smallq2[iHisto] = new TH1F(Form("hRatioMeanNtrkl_Smallq2_%s",Ntrklnames[iHisto].Data()),"",1,0.,q2smallpercevents*100);
        hRatioMeanNtrkl_Largeq2[iHisto]->GetYaxis()->SetTitle("<N_{tracklets}> ESE/unbiased");
        hRatioMeanNtrkl_Largeq2[iHisto]->GetXaxis()->SetTitle(hNtrkVsq2[iHisto]->GetXaxis()->GetTitle());
        hRatioMeanNtrkl_Smallq2[iHisto]->SetMarkerStyle(kFullSquare);
        hRatioMeanNtrkl_Smallq2[iHisto]->SetMarkerColor(colors[0]+1);
        hRatioMeanNtrkl_Smallq2[iHisto]->SetLineColor(colors[0]+1);
        hRatioMeanNtrkl_Smallq2[iHisto]->SetLineWidth(2);
        hRatioMeanNtrkl_Smallq2[iHisto]->SetBinContent(1,hNtrklq2Small[iHisto]->GetMean()/hNtrklUnbiased[iHisto]->GetMean());
        hRatioMeanNtrkl_Smallq2[iHisto]->SetBinError(1,1.e-10);
        
        for(Int_t iq2=0; iq2<100; iq2++) {
          TH1F* hNtrklq2Bin = (TH1F*)hNtrklVsq2VsCent[iHisto]->ProjectionZ("hNtrklq2Bin",-1,-1,iq2+1,iq2+2);
          hMeanNtrkl[iHisto]->SetBinContent(iq2+1,hNtrklq2Bin->GetMean());
          hMeanNtrkl[iHisto]->SetBinError(iq2+1,hNtrklq2Bin->GetMeanError());
          hRatioMeanNtrkl[iHisto]->SetBinContent(iq2+1,hNtrklq2Bin->GetMean()/hNtrklUnbiased[iHisto]->GetMean());
          hRatioMeanNtrkl[iHisto]->SetBinError(iq2+1,1.e-10);
        }

        cNtrkVsq2[iHisto] = new TCanvas(Form("cNtrkVsq2_%s",Ntrklnames[iHisto].Data()),"",800,800);
        cNtrkVsq2[iHisto]->SetLogz();
        hNtrkVsq2[iHisto]->Draw("colz");
        hMeanNtrkl[iHisto]->Draw("same");
        
        cMeanNtrkRatioVsq2[iHisto] = new TCanvas(Form("cMeanNtrkRatioVsq2_%s",Ntrklnames[iHisto].Data()),"",800,800);
        hRatioMeanNtrkl[iHisto]->GetYaxis()->SetRangeUser(hRatioMeanNtrkl[iHisto]->GetMinimum()*0.95,hRatioMeanNtrkl[iHisto]->GetMaximum()*1.05);
        hRatioMeanNtrkl[iHisto]->Draw();
        hRatioMeanNtrkl_Largeq2[iHisto]->Draw("same");
        hRatioMeanNtrkl_Smallq2[iHisto]->Draw("same");
        leg2->Draw("same");
        
        outnameNtrkl.ReplaceAll(Form("NtrklWeights%s",Ntrklnames[iHisto].Data()),Form("NtrklVsq2%s",Ntrklnames[iHisto].Data()));
        cNtrkVsq2[iHisto]->SaveAs(outnameNtrkl.Data());
        outnameNtrkl.ReplaceAll(Form("NtrklVsq2%s",Ntrklnames[iHisto].Data()),Form("MeanNtrklRatioVsq2%s",Ntrklnames[iHisto].Data()));
        cMeanNtrkRatioVsq2[iHisto]->SaveAs(outnameNtrkl.Data());
      }
    }
    TString outnameNtrkl = outname;
    outnameNtrkl.ReplaceAll("q2_vs_Centrality","NtrklWeights");
    TFile outfile(outnameNtrkl.Data(),"RECREATE");
    for(Int_t iHisto=0; iHisto<3; iHisto++) {
      hNtrklRatioq2SmallUnbiased[iHisto]->Write();
      hNtrklRatioq2LargeUnbiased[iHisto]->Write();
    }
    outfile.Close();
  }

  return 0;
}

//_____________________________________________________________________________________________
//DRAW RESOLUTION
void DrawEventPlaneResolutionAndDistribution(Int_t cutmeth) {

  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}

  SetStyle();

  TList* datalist=(TList*)LoadTList();
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
    intsmallcutvalues.push_back(-1);
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
    gRes[iq2]->GetYaxis()->SetTitle("Event Plane Resolution R_{2}");
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
  else if(cutmeth==kPercentileCut) percsuffix="percentile";
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

    TH2F* hFrameBias = new TH2F("hFrameBias","",6,-0.5,5.5,1,-0.08,0.1);
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
    else if(cutmeth==kPercentileCut) percsuffix="percentile";
    if(cutmeth==kPercCut || cutmeth==kPercCutVsCent) {outnameEP=Form("%s/EP_distributions%s_q2Small%.0f%s_q2Large%.0f%s.pdf",outputdir.Data(),suffix.Data(),q2smallpercevents*100,percsuffix.Data(),q2largepercevents*100,percsuffix.Data());}

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
TList* LoadTList() {

  cout << "Opening input file " << infilename << "..." <<endl;
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  TDirectoryFile* dir=0x0;
  TList* list=0x0;

  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s",partname.Data(),suffix.Data());

  if(infile) {dir=(TDirectoryFile*)infile->Get(dirname.Data()); cout << "File opened!" << endl;}
  else {cerr << "Error: File " << infilename << " not found. Exit." << endl; return 0x0;}
  if(dir) list=(TList*)dir->Get(listname.Data());
  else {
    dirname.ReplaceAll("HFv2","HFvn");
    dir=(TDirectoryFile*)infile->Get(dirname.Data());
    if(!dir) {cerr << "Error: Wrong TDirectoryFile name " << dirname << ". Exit." << endl; return 0x0;}
    list=(TList*)dir->Get(listname.Data());
  }
  if(!list) {cerr << "Error: Wrong TList name " << listname << ". Exit." << endl; return 0x0;}

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

  Double_t centwidth = (Double_t)(maxCent-minCent)/ncentbins;
  Double_t mincentpermil=minCent*10;
  Double_t maxcentpermil=(minCent+centwidth)*10;

  TH2F* hResVsq2=(TH2F*)inputlist->FindObject(Form("hEvPlaneReso1Vsq2centr%0.f_%0.f",mincentpermil,maxcentpermil));
  if(!hResVsq2) {return 0x0;}
  TH1F* hq2=(TH1F*)hResVsq2->ProjectionY();

  Int_t nq2bins = hq2->GetNbinsX();
  Int_t q2min = hq2->GetBinLowEdge(1);
  Int_t q2max = hq2->GetBinLowEdge(nq2bins)+hq2->GetBinWidth(nq2bins);
  TString q2name = hq2->GetXaxis()->GetTitle();
  
  TH2F* hq2VsCentr = new TH2F("hq2VsCentr","q_{2} vs. centrality",ncentbins,minCent,maxCent,nq2bins,q2min,q2max);
  hq2VsCentr->GetXaxis()->SetTitle("Centrality (%)");
  hq2VsCentr->GetYaxis()->SetTitle(q2name.Data());

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
    hq2VsCentr->GetYaxis()->SetTitle(q2name.Data());
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
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, TH1F* hMean[], TH1F* hMeanfs[], TH1F* hSigmaFree[], TH1F* hSigmaFixed[], Int_t smallorlarge, Int_t analysismeth, Int_t cutmeth) {
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
  if(!dir) {
    dirname.ReplaceAll("HFv2","HFvn");
    dir=(TDirectoryFile*)infile->Get(dirname.Data());
  }
  if(dir) {
    if(partname.Contains("Dzero")) {
      massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    }
    else if(partname.Contains("Dplus")){
      massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
    }
    else if(partname.Contains("Dstar")) {
      massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
    }
    else if(partname.Contains("Ds") && !partname.Contains("Dstar")) {
      massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
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
      Double_t sigmaerr=0;
      Double_t meanerr=0;
      if(ok){
        fitter->DrawHere(cDeltaPhi->cd(ipad),3,1);
        signal = fitter->GetRawYield();
        esignal = fitter->GetRawYieldError();
        sigmaforcounting=fitter->GetSigma();
        meanforcounting=fitter->GetMean();
        sigmaerr=fitter->GetSigmaUncertainty();
        meanerr=fitter->GetMeanUncertainty();
      }
      gSignal[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignal[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
      hMean[iPhi]->SetBinContent(iPt+1,meanforcounting);
      hMean[iPhi]->SetBinError(iPt+1,meanerr);
      hSigmaFree[iPhi]->SetBinContent(iPt+1,sigmaforcounting);
      hSigmaFree[iPhi]->SetBinError(iPt+1,sigmaerr);
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
    Double_t sigmaerr=fitter->GetSigmaUncertainty();
    Double_t massFromFit=fitter->GetMean();
    Double_t massFromFiterr=fitter->GetMeanUncertainty();
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
        massFromFit = fitter2->GetMean();
        if(!fixAlsoMass) massFromFiterr = fitter2->GetMeanUncertainty();
      }
      gSignalfs[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignalfs[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
      hMeanfs[iPhi]->SetBinContent(iPt+1,massFromFit);
      hMeanfs[iPhi]->SetBinError(iPt+1,massFromFiterr);
      hSigmaFixed[iPhi]->SetBinContent(iPt+1,sigma);
      hSigmaFixed[iPhi]->SetBinError(iPt+1,sigmaerr);
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
  else if(cutmeth==kPercentileCut) percsuffix="percentile";
  if(cutmeth==kPercCut || cutmeth==kPercCutVsCent || cutmeth==kPercentileCut) {
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
    else if(cutmeth==kPercentileCut) {
      smallcutvalues.push_back(q2smallpercevents*100);
      largecutvalues.push_back(100-q2largepercevents*100);
    }
  }

  TString q2name = "q2";
  if(cutmeth==kPercentileCut) {q2name = "q2 percentile";}
  cout << "Cut values for the small-q2 region:"<<endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    cout << "Centrality bin "<< iCent << Form(" %s < ",q2name.Data()) << smallcutvalues[iCent] <<endl;
  }
  cout << "\nCut values for the large-q2 region:"<<endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    cout << "Centrality bin "<< iCent << Form(" %s > ",q2name.Data()) << largecutvalues[iCent] <<endl;
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
  gStyle->SetPadLeftMargin(0.15);
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
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
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
