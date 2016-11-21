#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
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

#include "AliHFMassFitter.h"
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
const TString infilename="$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/highIR_pass1/AnalysisResultsv2_EP_pass0_diffdet_topod0cut_3050.root";
const TString suffix="_Topod0Cut_pass0_QoverM_VZERO";//"_3050_CentralCuts_TPC";
const TString partname="Dplus";
const Int_t minCent=30;
const Int_t maxCent=50;

//ptbins of the analysis
const Int_t nPtBins=7;
const Int_t nPtLims=nPtBins+1;
const Double_t PtLims[nPtLims] = {2.,3.,4.,5.,6.,8.,12.,16.};

//phi bins
const Int_t nPhiBins=4;
const Int_t nPhiLims=nPhiBins+1;
const Double_t PhiLims[nPhiLims]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};

//q2 cut values (absolute cut)
const Double_t q2smalllimit=2.;
const Double_t q2largelimit=2.;

//percentage of events with smaller/larger q2
const Double_t q2percevents=0.4;

//EP resolution
const Int_t evPlane=AliEventPlaneResolutionHandler::kVZERO;
Int_t evPlaneRes=AliEventPlaneResolutionHandler::kThreeSub;
const Bool_t useNcollWeight=kFALSE;
const Bool_t useAliHandlerForRes=kFALSE;

// mass fit configuration
const Int_t rebin[]={5,5,5,5,5,5,5,5};
const Int_t typeb=AliHFMassFitter::kExpo;  // Background: 0=expo, 1=linear, 2=pol2
const Bool_t fixAlsoMass=kFALSE;
const Double_t minMassForFit=1.69;
const Double_t maxMassForFit=2.02;
const Double_t nSigmaForCounting=3.5;

//not to be set
enum CutMethod{kAbsCut,kPercCut}; //kAbsCut->absolute cut values, kPercCut->cut according to the % of events with smaller/larger q2
enum SmallOrLarge{kSmall,kLarge};
enum AnalysisMethod{kEventPlane,kEventPlaneInOut,kScalarProd};

//in-out efficiency
const Double_t effInOverEffOut=1.03;

Int_t colors[]={kBlack,kBlue+2,kRed+1,kGreen+2,kMagenta+3,kOrange+7,kCyan+2,kViolet+5,kYellow+2,kBlue-7,kGreen,kMagenta,kAzure,kRed+2};
Int_t markers[]={kFullCircle,kFullSquare,kFullDiamond,kFullTriangleUp,kFullTriangleDown,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenTriangleUp,kOpenTriangleDown};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
Int_t DmesonsFlowEvShapeAnalysis(Int_t cutmeth=kAbsCut, Int_t analysismeth=kEventPlaneInOut);
void Drawq2VsCent();
void DrawResolution(Int_t cutmeth=kAbsCut);
TList* LoadTList();
THnSparseF* LoadSparseFromList(TList* inputlist);
TH2F* GetHistoq2VsCentr(TList* inputlist);
TList* LoadMassHistos(THnSparseF* sparse, Double_t cutvalues[2], Int_t smallorlarge, Int_t analysismeth);
TList* LoadResolutionHistos(TList *inputlist, Double_t cutvalues[2], Int_t smallorlarge);
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Int_t analysismeth, TGraphAsymmErrors *gRelSystEff);
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, Int_t smallorlarge, Int_t analysismeth);
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3);
void Applyq2Cut(THnSparseF* sparse, Double_t cutvalues[2], Int_t smallorlarge);
void ApplyCut(THnSparseF* sparse, Double_t min, Double_t max, UInt_t axnum);
Bool_t Getq2CutValuePercEvents(TH1F* hq2, Double_t cutvalues[2]);
void ResetAxes(THnSparseF* sparse, Int_t axnum=-1);
void SetStyle();

//_____________________________________________________________________________________________
//ANALYSIS FUNCTION
Int_t DmesonsFlowEvShapeAnalysis(Int_t cutmeth, Int_t analysismeth) {
  
  TList* datalist=(TList*)LoadTList();
  if(!datalist){return 1;}
  THnSparseF* hMassPtPhiq2Centr=(THnSparseF*)LoadSparseFromList(datalist);
  if(!hMassPtPhiq2Centr) {return 2;}
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) return 3;
  TH1F* hq2=(TH1F*)hq2VsCentr->ProjectionY();
  
  Double_t cutvalues[2] = {-1.,-1.};
  if(cutmeth==kAbsCut) {cutvalues[0] = q2smalllimit; cutvalues[1] = q2largelimit;}
  else if(cutmeth) {Bool_t cuttaken = Getq2CutValuePercEvents(hq2,cutvalues);}

  TString analysismethname="";
  if(analysismeth==kEventPlane) {analysismethname="EP";}
  else if(analysismeth==kEventPlaneInOut) {analysismethname="InOut";}
  else if(analysismeth==kScalarProd) {analysismethname="ScalarProduct";}
  
  Int_t q2region[2] = {kSmall,kLarge};
  TString q2regionname[2] = {"q2Small","q2Large"};
  
  SetStyle();
  
  //graphs for v2 in the two q2 regions
  TCanvas *cv2fs =new TCanvas("cv2fs","v2 Fit with fixed sigma",1920,1080);
  TCanvas **cv2 =new TCanvas*[2];
  TGraphAsymmErrors **gv2=new TGraphAsymmErrors*[2];
  TGraphAsymmErrors **gv2fs=new TGraphAsymmErrors*[2];
  TGraphAsymmErrors **gv2BC1=new TGraphAsymmErrors*[2];
  TGraphAsymmErrors **gv2BC2=new TGraphAsymmErrors*[2];
  TGraphAsymmErrors **gRelSystEff=new TGraphAsymmErrors*[2];
  Double_t ymin[2]={1.,1.};
  Double_t ymax[2]={-1.,-1.};
  
  TLegend* leg = new TLegend(0.65,0.7,0.85,0.85);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  
  TLegend** legsyst = new TLegend*[2];
  
  TFile outfile(Form("v2Output_%d_%d_%s%s.root",minCent,maxCent,analysismethname.Data(),suffix.Data()),"RECREATE"); //outputfile
  
  for(Int_t iq2=0; iq2<2; iq2++) {
    if(analysismeth!=kScalarProd) {
      //load resolution and mass lists
      TList* masslist=(TList*)LoadMassHistos(hMassPtPhiq2Centr,cutvalues,q2region[iq2],analysismeth);
      masslist->SaveAs(Form("v2Histograms_%d_%d_%s%s_%s.root",minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[iq2].Data()),"RECREATE");

      cout <<"Average pt for pt bin \n"<<endl;
      //average pt for pt bin
      AliVertexingHFUtils *utils=new AliVertexingHFUtils();
      Int_t minCentTimesTen=minCent*10;
      Int_t maxCentTimesTen=maxCent*10;
      Applyq2Cut(hMassPtPhiq2Centr,cutvalues,q2region[iq2]);
      TH2F* hmasspt=(TH2F*)hMassPtPhiq2Centr->Projection(0,1);
      ResetAxes(hMassPtPhiq2Centr,-1.);
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
        AliHFMassFitter fitter(histtofit,hmin,hmax,1);
        fitter.MassFitter(kFALSE);
        Double_t massFromFit=fitter.GetMean();
        Double_t sigmaFromFit=fitter.GetSigma();
        TF1* funcB2=fitter.GetBackgroundRecalcFunc();
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
      FillSignalGraph(masslist,gSignal,gSignalfs,gSignalBC1,gSignalBC2,q2region[iq2],analysismeth);
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
        TList* resolist=(TList*)LoadResolutionHistos(datalist,cutvalues,q2region[iq2]);
        TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
        TH1F* hevplresos[3];
        TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
        Int_t nSubRes=1;
        if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
           evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
        for(Int_t iRes=0;iRes<nSubRes;iRes++){
          hevplresos[iRes]=(TH1F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
        }
        resol=GetEventPlaneResolution(errorres,hevplresos[0],hevplresos[1],hevplresos[2]);
      }
      
      cout << "Event plane resolution: " << resol <<endl;
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
      gv2[iq2]->SetMarkerStyle(markers[iq2]);
      gv2[iq2]->SetMarkerSize(1.5);
      gv2[iq2]->SetMarkerColor(colors[iq2]);
      gv2[iq2]->SetLineColor(colors[iq2]);
      gv2[iq2]->SetLineWidth(2);
      gv2fs[iq2]->SetMarkerStyle(markers[iq2+2]);
      gv2fs[iq2]->SetMarkerSize(1.5);
      gv2fs[iq2]->SetMarkerColor(colors[iq2+2]);
      gv2fs[iq2]->SetLineColor(colors[iq2+2]);
      gv2fs[iq2]->SetLineWidth(2);
      gv2BC1[iq2]->SetMarkerStyle(markers[iq2+5]);
      gv2BC1[iq2]->SetMarkerSize(1.5);
      gv2BC1[iq2]->SetMarkerColor(colors[iq2+5]);
      gv2BC1[iq2]->SetLineColor(colors[iq2+5]);
      gv2BC1[iq2]->SetLineWidth(2);
      gv2BC2[iq2]->SetMarkerStyle(markers[iq2+7]);
      gv2BC2[iq2]->SetMarkerSize(1.5);
      gv2BC2[iq2]->SetMarkerColor(colors[iq2+7]);
      gv2BC2[iq2]->SetLineColor(colors[iq2+7]);
      gv2BC2[iq2]->SetLineWidth(2);
      gRelSystEff[iq2]->SetMarkerStyle(markers[iq2]);
      gRelSystEff[iq2]->SetMarkerSize(1.5);
      gRelSystEff[iq2]->SetMarkerColor(colors[iq2]);
      gRelSystEff[iq2]->SetLineColor(colors[iq2]);
      gRelSystEff[iq2]->SetLineWidth(2);
      
      TString title="";
      if(cutmeth==kAbsCut) {
        TString sign[2] = {"<",">"};
        title=Form("%s %s %0.1f",hq2->GetXaxis()->GetTitle(),sign[iq2].Data(),cutvalues[iq2]);
      }
      else {
        Double_t perc = q2percevents*100;
        TString reg[2] = {"small","large"};
        title=Form("%0.f%% %s %s",perc,reg[iq2].Data(),hq2->GetXaxis()->GetTitle());
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
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.5*TMath::Abs(tmp);
      tmp=TMath::MaxElement(gv2fs[iq2]->GetN(),gv2fs[iq2]->GetY());
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.5*TMath::Abs(tmp);
      tmp=TMath::MaxElement(gv2BC1[iq2]->GetN(),gv2BC1[iq2]->GetY());
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.5*TMath::Abs(tmp);
      tmp=TMath::MaxElement(gv2BC2[iq2]->GetN(),gv2BC2[iq2]->GetY());
      if(ymax[iq2]<tmp+1.5*TMath::Abs(tmp)) ymax[iq2]=tmp+1.5*TMath::Abs(tmp);
   
      gv2[iq2]->GetYaxis()->SetRangeUser(ymin[iq2],ymax[iq2]);
      
      cv2[iq2]->SaveAs(Form("v2Output_%d_%d_%s%s_%s.pdf",minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[iq2].Data()));

      cv2fs->cd();
      TString drawopt="AP";
      if(iq2==1) drawopt="P";
      gv2fs[iq2]->Draw(drawopt.Data());
      if(iq2==1) leg->Draw("same");
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
  
  Double_t min=ymin[0];
  if(min>ymin[1]) min=ymin[1];
  Double_t max=ymax[0];
  if(max<ymax[1]) max=ymax[1];
  gv2fs[0]->GetYaxis()->SetRangeUser(min,max);
  
  cv2fs->SaveAs(Form("v2fsOutput_%d_%d_%s%s.pdf",minCent,maxCent,analysismethname.Data(),suffix.Data()));
  
  return 0;
}

//_____________________________________________________________________________________________
//DRAW q_2 VS CENTRALITY
void Drawq2VsCent() {

  SetStyle();
  
  TList* datalist=(TList*)LoadTList();
  if(!datalist){return;}
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) return;
  
  const Int_t ncentbins = hq2VsCentr->GetXaxis()->GetNbins();
  TH1F* hq2=(TH1F*)hq2VsCentr->ProjectionY();
  hq2->SetTitle("");
  hq2->GetYaxis()->SetTitle(Form("dN_{ev}/d%s",hq2->GetXaxis()->GetTitle()));
  hq2->SetLineColor(colors[0]);
  
  TH1F** hq2centbin = new TH1F*[ncentbins];
  
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    hq2centbin[iCent]=(TH1F*)hq2VsCentr->ProjectionY(Form("hq2centbin%d",iCent),iCent+1,iCent+1);
    hq2centbin[iCent]->GetYaxis()->SetTitle(Form("dN_{ev}/d%s",hq2->GetXaxis()->GetTitle()));
    hq2centbin[iCent]->SetLineColor(colors[iCent+1]);
  }
  
  TLegend* leg = new TLegend(0.65,0.45,0.89,0.89);
  leg->SetFillStyle(0);
  
  TCanvas *cq2VsCent = new TCanvas("cq2VsCent","q_{2} vs. centrality",800,800);
  cq2VsCent->SetLogz();
  hq2VsCentr->Draw("colz");
  TCanvas *cq2 = new TCanvas("cq2","q_{2}",800,800);
  cq2->SetLogy();
  hq2->Draw();
  leg->AddEntry(hq2,Form("%d - %d",minCent,maxCent));
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    hq2centbin[iCent]->Draw("same");
    Double_t min=minCent+iCent*2.5;
    Double_t max=minCent+(iCent+1)*2.5;
    leg->AddEntry(hq2centbin[iCent],Form("%.1f - %.1f",min,max));
  }
  leg->Draw("same");
  
  cq2VsCent->SaveAs(Form("q2_vs_Centrality%s.pdf",suffix.Data()));
  cq2->SaveAs(Form("Nev_vs_q2%s.pdf",suffix.Data()));
  
}

//_____________________________________________________________________________________________
//DRAW RESOLUTION
void DrawResolution(Int_t cutmeth) {

  SetStyle();
  
  TList* datalist=(TList*)LoadTList();
  if(!datalist){return;}
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) return;
  TH1F* hq2=(TH1F*)hq2VsCentr->ProjectionY();

  Double_t cutvalues[2] = {-1.,-1.};
  if(cutmeth==kAbsCut) {cutvalues[0] = q2smalllimit; cutvalues[1] = q2largelimit;}
  else if(cutmeth==kPercCut) {Bool_t cuttaken = Getq2CutValuePercEvents(hq2,cutvalues);}
  
  Int_t q2region[2] = {kSmall,kLarge};
  
  TGraphErrors** gRes = new TGraphErrors*[3];
  TList *resolist=0x0;
  Double_t error[3]={0.,0.,0.};
  Double_t resoFull[3]={0.,0.,0.};
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);

  for(Int_t iq2=0; iq2<3; iq2++) {
    resolist=(TList*)LoadResolutionHistos(datalist,cutvalues,q2region[iq2]);
    if(iq2==3) {
      Double_t q2cuts[2] = {11.,-1.};
      resolist=(TList*)LoadResolutionHistos(datalist,q2cuts,kSmall); //resolution integrated in q2
    }
    gRes[iq2]=(TGraphErrors*)resolist->FindObject("gResolVsCent");
    if(iq2<3) {gRes[iq2]->SetName(Form("gResolVsCent%d",q2region[iq2]));}
    else {gRes[iq2]->SetName("gResolVsCent");}
    gRes[iq2]->GetYaxis()->SetTitle("R_{2}");
    gRes[iq2]->GetXaxis()->SetTitle("Centrality (%)");
    gRes[iq2]->SetTitle("");
    gRes[iq2]->SetMarkerSize(1.5);
    gRes[iq2]->SetLineWidth(2);
    gRes[iq2]->SetMarkerStyle(markers[iq2]);
    gRes[iq2]->SetMarkerColor(colors[iq2]);
    gRes[iq2]->SetLineColor(colors[iq2]);
    TH1F* hevplresos[3];
    TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
    Int_t nSubRes=1;
    if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
       evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
    for(Int_t iRes=0;iRes<nSubRes;iRes++){
      hevplresos[iRes]=(TH1F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
    }
    resoFull[iq2]=GetEventPlaneResolution(error[iq2],hevplresos[0],hevplresos[1],hevplresos[2]);
  }
  
  TLegend* leg = new TLegend(0.35,0.2,0.8,0.5);
  leg->SetFillStyle(0);
  if(cutmeth==kAbsCut) {
    leg->AddEntry(gRes[0],Form("%s < %0.1f, <R_{2}> = %0.4f",hq2->GetXaxis()->GetTitle(),cutvalues[0],resoFull[0]),"lpe");
    leg->AddEntry(gRes[1],Form("%s > %0.1f, <R_{2}> = %0.4f",hq2->GetXaxis()->GetTitle(),cutvalues[1],resoFull[1]),"lpe");
  }
  else {
    Double_t perc = q2percevents*100;
    leg->AddEntry(gRes[0],Form("%0.f%% small %s, <R_{2}> = %0.4f",perc,hq2->GetXaxis()->GetTitle(),resoFull[0]),"lpe");
    leg->AddEntry(gRes[1],Form("%0.f%% large %s, <R_{2}> = %0.4f",perc,hq2->GetXaxis()->GetTitle(),resoFull[1]),"lpe");
  }
  leg->AddEntry(gRes[2],Form("%s-integrated, <R_{2}> = %0.4f",hq2->GetXaxis()->GetTitle(),resoFull[2]),"lpe");
  
  TCanvas* cRes = new TCanvas("cRes","",1920,1080);
  Int_t npoints = gRes[0]->GetN();
  Double_t min = TMath::MinElement(npoints,gRes[0]->GetY())*0.5;
  if(min>TMath::MinElement(npoints,gRes[1]->GetY())) {min = TMath::MinElement(npoints,gRes[0]->GetY())*0.9;}
  gRes[0]->GetYaxis()->SetRangeUser(min,1);
  gRes[0]->Draw("AP");
  gRes[1]->Draw("P");
  gRes[2]->Draw("P");
  leg->Draw("same");
  
  cRes->SaveAs(Form("EP_resolution%s.pdf",suffix.Data()));
  TFile outfile(Form("EP_resolution%s.root",suffix.Data()),"RECREATE");
  cRes->Write();
  gRes[0]->Write();
  gRes[1]->Write();
  gRes[2]->Write();
  outfile.Close();
}

//_____________________________________________________________________________________________
//LOAD DATA LIST
TList* LoadTList() {
  
  cout << "Opening input file " << infilename << "..." <<endl;
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  TDirectoryFile* dir=0x0;
  TList* list=0x0;
  
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s_EvShape",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s_EvShape",partname.Data(),suffix.Data());

  if(infile) {dir=(TDirectoryFile*)infile->Get(dirname.Data()); cout << "File opened!" << endl;}
  else {cerr << "Error: File " << infilename << " not found. Exit." << endl; return 0x0;}
  if(dir) list=(TList*)dir->Get(listname.Data());
  else {cerr << "Error: Wrong TDirectoryFile name " << dirname << ". Exit." << endl; return 0x0;}
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
  
  Int_t nq2bins = 100;
  Int_t q2min = 0.;
  Int_t q2max = 10.;
  Double_t step = 2.5;
  Int_t ncentbins = (maxCent-minCent)/step;
  
  TH2F* hq2VsCentr = new TH2F("hq2VsCentr","q_{2} vs. centrality",ncentbins,minCent,maxCent,nq2bins,q2min,q2max);
  hq2VsCentr->GetXaxis()->SetTitle("Centrality (%)");
  
  Double_t entries=0;
  
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    Double_t mincentpermil=(minCent+iCent*step)*10;
    Double_t maxcentpermil=(minCent+(iCent+1)*step)*10;
    
    TH2F* hResVsq2=(TH2F*)inputlist->FindObject(Form("hEvPlaneResoVsq2centr%0.f_%0.f",mincentpermil,maxcentpermil));
    if(!hResVsq2) return 0x0;
    TH1F* hq2 = (TH1F*)hResVsq2->ProjectionY();
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
TList* LoadMassHistos(THnSparseF* sparse, Double_t cutvalues[2], Int_t smallorlarge, Int_t analysismeth) {
  
  TList *outlist = new TList();
  outlist->SetName("azimuthalhistoslist");
  
  Applyq2Cut(sparse,cutvalues,smallorlarge);
  
  TH1F* hMassTmp=0x0;
  TH1F* hMassTmp1=0x0;
  TH1F* hMassTmp2=0x0;
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    ApplyCut(sparse,PtLims[iPt],PtLims[iPt+1],1); //apply pt cut
    hMassTmp=(TH1F*)sparse->Projection(0);
    hMassTmp->SetName(Form("hMass_pt%d",iPt));
    outlist->Add(hMassTmp->Clone());
    if(analysismeth==kEventPlane) {
      for(Int_t iPhi=0; iPhi<nPhiBins; iPhi++) {
        ApplyCut(sparse,PhiLims[iPhi],PhiLims[iPhi+1],2); //apply deltaphi cut
        hMassTmp=(TH1F*)sparse->Projection(0);
        hMassTmp->SetName(Form("hMass_pt%d_phi%d",iPt,iPhi));
        outlist->Add(hMassTmp->Clone());
        ResetAxes(sparse,2); //release deltaphi cut
        delete hMassTmp;
        hMassTmp=0x0;
      }
    }
    else if(analysismeth==kEventPlaneInOut){
      ApplyCut(sparse,0.,TMath::Pi()/4,2); //apply deltaphi cut
      hMassTmp=(TH1F*)sparse->Projection(0);
      ResetAxes(sparse,2); //release deltaphi cut
      ApplyCut(sparse,3./4*TMath::Pi(),TMath::Pi(),2); //apply deltaphi cut
      hMassTmp1=(TH1F*)sparse->Projection(0);
      ResetAxes(sparse,2); //release deltaphi cut
      hMassTmp->Add(hMassTmp1);
      hMassTmp->SetName(Form("hMass_pt%d_phi0",iPt));
      outlist->Add(hMassTmp->Clone());
      ApplyCut(sparse,TMath::Pi()/4,3./4*TMath::Pi(),2); //apply deltaphi cut
      hMassTmp2=(TH1F*)sparse->Projection(0);
      ResetAxes(sparse,2); //release deltaphi cut
      hMassTmp2->SetName(Form("hMass_pt%d_phi1",iPt));
      outlist->Add(hMassTmp2->Clone());
      delete hMassTmp;
      delete hMassTmp1;
      delete hMassTmp2;
      hMassTmp=0x0;
      hMassTmp1=0x0;
      hMassTmp2=0x0;
    }
    ResetAxes(sparse,1); //release pt cut
  }
  
  ResetAxes(sparse,-1); //release all cuts
  
  return outlist;
}

//_____________________________________________________________________________________________
//GET RESOLUTION HISTOS IN SELECTED q2 RANGE
TList* LoadResolutionHistos(TList *inputlist, Double_t cutvalues[2], Int_t smallorlarge) {

  TList *outlist = new TList();
  outlist->SetName("eventplanemasslist");
  
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
  TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
  TH2F* hevplresosvsq2[3];
  TH1F* hevplresos[3];
  Int_t ncBin=minCentTimesTen/25;
  
  for(Int_t iRes=0;iRes<nSubRes;iRes++){
    hevplresosvsq2[iRes]=(TH2F*)inputlist->FindObject(Form("hEvPlane%scentr%d_%d",namereso[iRes].Data(),minCentTimesTen,minCentTimesTen+25));
    Int_t binmin=-1;
    Int_t binmax=-1;
    //define q2 cut bins
    if(smallorlarge==kSmall) {
      binmin=1;
      binmax=hevplresosvsq2[iRes]->GetYaxis()->FindBin(cutvalues[0]*0.9999);
      cout << "\nCut values for the small q2 region (resolution histograms):\nbin min: " << binmin << "\nbin max: " << binmax <<"\n"<<endl;
    }
    else if(smallorlarge==kLarge) {
      binmin=hevplresosvsq2[iRes]->GetYaxis()->FindBin(cutvalues[1]*1.0001);
      binmax=hevplresosvsq2[iRes]->GetYaxis()->FindBin(hevplresosvsq2[iRes]->GetYaxis()->GetNbins()+1); //takes also overflow entries
      Int_t binmaxsmallregion=hevplresosvsq2[iRes]->GetYaxis()->FindBin(cutvalues[0]*0.9999);
      if(binmaxsmallregion==binmin) {
        binmin += 1;
        cout << "Warning: maximum bin for the small q2 range equal to the minimum bin for the large q2 range: shift the minimum bin for the large q2 range in order to avoid overlap" << endl;
      }
      cout << "\nCut values for the large q2 region (resolution histograms):\nbin min: " << binmin << "\nbin max: " << binmax <<"\n"<<endl;
    }
    hevplresos[iRes]=(TH1F*)hevplresosvsq2[iRes]->ProjectionX(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()),binmin,binmax);
    if(hevplresos[iRes]){
      hevplresos[iRes]->SetName(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
      hevplresos[iRes]->SetTitle(Form("Event Plane Resolution %s%s",namereso[iRes].Data(),suffixcentr.Data()));
      if(useNcollWeight){
        cout << Form("Centr %d Bin %d  Ncoll %f\n",minCentTimesTen,ncBin,ncoll[ncBin]) << endl;
        hevplresos[iRes]->Scale(ncoll[ncBin]);
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
    TH2F* htmpresosvsq2[3];
    TH1F* htmpresos[3];
    for(Int_t iRes=0;iRes<nSubRes;iRes++){
      htmpresosvsq2[iRes]=(TH2F*)inputlist->FindObject(Form("hEvPlane%scentr%d_%d",namereso[iRes].Data(),icentr,icentr+25));
      if(!htmpresosvsq2[iRes])cout<<"skipping ev pl reso "<<icentr<<"_"<<icentr+25<<endl;
      Int_t binmin=-1;
      Int_t binmax=-1;
      //define q2 cut bins
      if(smallorlarge==kSmall) {
        binmin=1;
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(cutvalues[0]*0.9999);
      }
      else if(smallorlarge==kLarge) {
        binmin=htmpresosvsq2[iRes]->GetYaxis()->FindBin(cutvalues[1]*1.0001);
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(htmpresosvsq2[iRes]->GetYaxis()->GetNbins()+1); //takes also overflow entries
        Int_t binmaxsmallregion=htmpresosvsq2[iRes]->GetYaxis()->FindBin(cutvalues[0]*0.9999);
        if(binmaxsmallregion==binmin) {
          binmin += 1;
        }
      }
      htmpresos[iRes]=(TH1F*)htmpresosvsq2[iRes]->ProjectionX(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()),binmin,binmax);
    }
    resolBin=GetEventPlaneResolution(error,htmpresos[0],htmpresos[1],htmpresos[2]);
    if(resolBin<lowestRes) lowestRes=resolBin;
    if(resolBin>highestRes) highestRes=resolBin;
    binCentr=(Double_t)icentr/10.+binHalfWid;
    gResolVsCent->SetPoint(iPt,binCentr,resolBin);
    gResolVsCent->SetPointError(iPt,binHalfWid,error);
    ++iPt;
    ncBin=icentr/25;
    for(Int_t iRes=0;iRes<nSubRes;iRes++){
      if(htmpresos[iRes]){
        if(useNcollWeight){
          cout << Form("Centr %d Bin %d  Ncoll %f\n",icentr,ncBin,ncoll[ncBin]) << endl;
          htmpresos[iRes]->Scale(ncoll[ncBin]);
        }
        hevplresos[iRes]->Add(htmpresos[iRes]);
      }
    }
  }
  for(Int_t iRes=0;iRes<nSubRes;iRes++){
    if(hevplresos[iRes]) outlist->Add(hevplresos[iRes]->Clone());
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
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2, Int_t smallorlarge, Int_t analysismeth) {
  
  TString q2regionname="";
  if(smallorlarge==kSmall) {q2regionname="q2Small";}
  if(smallorlarge==kLarge) {q2regionname="q2Large";}
  
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  TDirectoryFile* dir=0x0;
  TList* list=0x0;
  Double_t massD=-1;
  
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s_EvShape",partname.Data(),suffix.Data());
  
  if(infile) {dir=(TDirectoryFile*)infile->Get(dirname.Data());}
  if(dir) {
    if(partname.Contains("Dzero")) {
      massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    }
    if(partname.Contains("Dplus")){
      massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
    }
    if(partname.Contains("Dstar")) {
      massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
    }
    if(partname.Contains("Ds")) {
      massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
    }
  }
  else {cerr << "Error: Wrong TDirectoryFile name " << dirname << ". Exit." << endl;}
  
  Int_t nPhi=nPhiBins;
  if(analysismeth==kEventPlaneInOut) nPhi=2;
  
  //Canvases for drawing histograms
  TCanvas *cDeltaPhi = new TCanvas(Form("cinvmassdeltaphi_%s",q2regionname.Data()),Form("Invariant mass distributions - %s",q2regionname.Data()),1920,1080);
  TCanvas *cDeltaPhifs = new TCanvas(Form("cinvmassdeltaphifs_%s",q2regionname.Data()),Form("Invariant mass distributions - fit with fixed sigma - %s",q2regionname.Data()),1920,1080);
  TCanvas *cPhiInteg = new TCanvas(Form("cinvmass_%s",q2regionname.Data()),Form("Invariant mass distributions - #phi integrated - %s",q2regionname.Data()),1920,1080);
  cDeltaPhi->Divide(nPtBins,nPhi);
  cDeltaPhifs->Divide(nPtBins,nPhi);
  Int_t nptpads = nPtBins/2;
  if(nPtBins%2!=0) nptpads = nPtBins/2+1;
  if(nPtBins>1) cPhiInteg->Divide(nptpads,2);
  Int_t nMassBins;
  Double_t hmin,hmax;
  for(Int_t iPt=0;iPt<nPtBins;iPt++){
    TH1F *histtofitfullsigma=(TH1F*)masslist->FindObject(Form("hMass_pt%d",iPt))->Clone();
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
      hmin=TMath::Max(minMassForFit,histtofit->GetBinLowEdge(2));
      hmax=TMath::Min(maxMassForFit,histtofit->GetBinLowEdge(nMassBins-2));
      histtofit->Rebin(rebin[iPt]);
      histtofit->GetYaxis()->SetTitle(Form("Entries / (%0.f MeV/c)",histtofit->GetBinWidth(5)*1000));
      AliHFMassFitter fitter(histtofit,hmin,hmax,1,typeb);
      fitter.SetInitialGaussianMean(massD);
      fitter.SetInitialGaussianSigma(0.012);
      Bool_t ok=fitter.MassFitter(kFALSE);
      Double_t sigmaforcounting=0;
      Double_t meanforcounting=0;
      if(ok){
        fitter.DrawHere(cDeltaPhi->cd(ipad),3,1);
        fitter.Signal(3,signal,esignal);
        sigmaforcounting=fitter.GetSigma();
        meanforcounting=fitter.GetMean();
      }
      gSignal[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignal[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
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
      gSignalBC1[iPt]->SetPoint(iPhi,iPhi,cntSig1);
      gSignalBC1[iPt]->SetPointError(iPhi,0,0,cntErr,cntErr);
      gSignalBC2[iPt]->SetPoint(iPhi,iPhi,cntSig2);
      gSignalBC2[iPt]->SetPointError(iPhi,0,0,cntErr,cntErr);
    }
    //fit for fixed sigma
    histtofitfullsigma->SetTitle(Form("%.0f < #it{p}_{T} < %.0f GeV/c",PtLims[iPt],PtLims[iPt+1]));
    histtofitfullsigma->GetXaxis()->SetTitleSize(0.05);
    nMassBins=histtofitfullsigma->GetNbinsX();
    hmin=TMath::Max(minMassForFit,histtofitfullsigma->GetBinLowEdge(2));
    hmax=TMath::Min(maxMassForFit,histtofitfullsigma->GetBinLowEdge(nMassBins-2));
    histtofitfullsigma->Rebin(rebin[iPt]);
    histtofitfullsigma->GetYaxis()->SetTitle(Form("Entries / (%0.f MeV/c)",histtofitfullsigma->GetBinWidth(5)*1000));
    AliHFMassFitter fitter(histtofitfullsigma,hmin,hmax,1,typeb);
    fitter.SetInitialGaussianMean(massD);
    Bool_t ok=fitter.MassFitter(kFALSE);
    if(ok){
      if(nPtBins==1) {fitter.DrawHere(cPhiInteg->cd(),3,1);}
      else {fitter.DrawHere(cPhiInteg->cd(iPt+1),3,1);}
    }
    Double_t sigma=fitter.GetSigma();
    Double_t massFromFit=fitter.GetMean();
    for(Int_t iPhi=0;iPhi<nPhi;iPhi++){
      Int_t ipad=(iPhi)*nPtBins+iPt+1;
      TH1F *histtofit=(TH1F*)masslist->FindObject(Form("hMass_pt%d_phi%d",iPt,iPhi))->Clone();
      histtofit->SetTitle(Form("%.f<#it{p}_{T}<%.f, #phi%d",PtLims[iPt],PtLims[iPt+1],iPhi));
      nMassBins=histtofit->GetNbinsX();
      hmin=TMath::Max(minMassForFit,histtofit->GetBinLowEdge(2));
      hmax=TMath::Min(maxMassForFit,histtofit->GetBinLowEdge(nMassBins-2));
      histtofit->Rebin(rebin[iPt]);
      histtofit->GetYaxis()->SetTitle(Form("Entries / (%0.f MeV/c)",histtofit->GetBinWidth(5)*1000));
      AliHFMassFitter fitter2(histtofit,hmin,hmax,1,typeb);
      fitter2.SetInitialGaussianMean(massD);
      fitter2.SetFixGaussianSigma(sigma);
      if(fixAlsoMass) fitter2.SetFixGaussianMean(massFromFit);
      Bool_t ok2=fitter2.MassFitter(kFALSE);
      Double_t signal=0,esignal=0;
      if(ok2){
        fitter2.DrawHere(cDeltaPhifs->cd(ipad),3,1);
        fitter2.Signal(3,signal,esignal);
      }
      gSignalfs[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignalfs[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
    }
  }//end loop on pt bin
  
  cDeltaPhi->SaveAs(Form("InvMassDeltaPhi%s_%s.pdf",suffix.Data(),q2regionname.Data()));
  cDeltaPhifs->SaveAs(Form("InvMassDeltaPhi_fs%s_%s.pdf",suffix.Data(),q2regionname.Data()));
  cDeltaPhifs->SaveAs(Form("InvMassDeltaPhi_fs%s_%s.root",suffix.Data(),q2regionname.Data()));
  cPhiInteg->SaveAs(Form("InvMassfullphi%s_%s.pdf",suffix.Data(),q2regionname.Data()));
  
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
void Applyq2Cut(THnSparseF* sparse, Double_t cutvalues[2], Int_t smallorlarge) {
  
  UInt_t q2axnum=3;
  
  TAxis* q2axis = (TAxis*)sparse->GetAxis(q2axnum);
  Int_t binmin = -1;
  Int_t binmax = -1;
  //define q2 cut bins
  if(smallorlarge==kSmall) {
    binmin=q2axis->FindBin(0.0001);
    binmax=q2axis->FindBin(cutvalues[0]-0.0001);
    cout << "\nCut values for the small q2 region:\nbin min: " << binmin << "\nbin max: " << binmax <<"\n"<<endl;
  }
  else if(smallorlarge==kLarge) {
    binmin=q2axis->FindBin(cutvalues[1]+1.0001);
    binmax=q2axis->FindBin(q2axis->GetNbins()+1); //takes also overflow entries
    Int_t binmaxsmallregion=q2axis->FindBin(cutvalues[0]-0.0001);
    if(binmaxsmallregion==binmin) {
      binmin += 1;
      cout << "Warning: maximum bin for the small q2 range equal to the minimum bin for the large q2 range: shift the minimum bin for the large q2 range in order to avoid overlap" << endl;
    }
    cout << "\nCut values for the large q2 region:\nbin min: " << binmin << "\nbin max: " << binmax <<"\n"<<endl;
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
//q2 CUT VALUE SEARCH (ON A BASIS OF PERCENTAGE OF EVENTS WITH SMALLER/LARGER q2)
Bool_t Getq2CutValuePercEvents(TH1F* hq2, Double_t cutvalues[2]) {

  Int_t bincutvalues[2]={-1,-1};
  
  if(hq2) {
    Double_t ThresIntegral = q2percevents*hq2->Integral(0.,1000.); //threshold integral value (including overflow entries)
    Double_t integral=0;
    Int_t q2bin=0;
    while(integral<ThresIntegral) {
      q2bin++;
      integral += hq2->GetBinContent(q2bin);
    }
    bincutvalues[0]=q2bin-1;
    if(q2percevents==0.5) bincutvalues[1]=bincutvalues[0];
    else {
      integral=0;
      q2bin=hq2->GetNbinsX()+2;
      while(integral<ThresIntegral) {
        q2bin--;
        integral += hq2->GetBinContent(q2bin); //takes also overflow counts
      }
      bincutvalues[1]=q2bin+1;
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
