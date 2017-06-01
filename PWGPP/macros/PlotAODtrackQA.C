#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include "TTree.h"
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#endif

void ProjectDrawAndFitK0(TH2F* hInvMass2d, TF1* fmassk0, Double_t xmin, Double_t xmax);
void ProjectDrawAndFit(TH3F* hInvMass3d, TString partName, TF1* fmass, Double_t xmin, Double_t xmax);
void ProjectDrawAndFitMomP(TH3F* hInvMassLambdaGoodHyp3d,TH3F* hInvMassLambdaBadHyp3d,TH3F* hInvMassAntiLambdaGoodHyp3d,TH3F* hInvMassAntiLambdaBadHyp3d,TF1* fmass, Double_t xmin, Double_t xmax);
void InitFuncAndFit(TH1D* hm, TF1* fmass);
void DrawDistrib(TH1D* h1, TH1D* h2, TH1D* h3, Bool_t showStat);
TH1D* ComputeMatchEff(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle);
TH1D* ComputeRatio(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle);
void FillMeanAndRms(TH2F* h2d, TGraphErrors* gMean, TGraphErrors* gRms, TGraphErrors* gRel);

const Int_t totTrending=36;
Float_t vecForTrend[totTrending];

void PlotAODtrackQA(TString filename="AnalysisResults.root", TString suffix="QA", Int_t runNumber=111111){

  TString varForTrending[totTrending];
  for(Int_t jbit=0; jbit<12; jbit++){
    varForTrending[2*jbit]=Form("averPtFB%d",jbit);
    varForTrending[2*jbit+1]=Form("RatioFB%dFB0",jbit);
  }
  varForTrending[24]="MatchEffPt350EtaPos";
  varForTrending[25]="MatchEffPt1000EtaPos";
  varForTrending[26]="MatchEffPt4000EtaPos";
  varForTrending[27]="MatchEffPt350EtaNeg";
  varForTrending[28]="MatchEffPt1000EtaNeg";
  varForTrending[29]="MatchEffPt4000EtaNeg";
  varForTrending[30]="MatchEffSPDPt350EtaPos";
  varForTrending[31]="MatchEffSPDPt1000EtaPos";
  varForTrending[32]="MatchEffSPDPt4000EtaPos";
  varForTrending[33]="MatchEffSPDPt350EtaNeg";
  varForTrending[34]="MatchEffSPDPt1000EtaNeg";
  varForTrending[35]="MatchEffSPDPt4000EtaNeg";
 
  TTree* trtree=new TTree("trending","tree of trending variables");
  trtree->Branch("nrun",&runNumber,"nrun/I");
  for(Int_t j=0; j<totTrending; j++){
    trtree->Branch(varForTrending[j].Data(),&vecForTrend[j],Form("%s/F",varForTrending[j].Data()));
    vecForTrend[j]=-99.;
  }
  
  TFile* f=new TFile(filename.Data());
  TDirectoryFile* df=(TDirectoryFile*)f->Get("CheckAODTracks");
  TList* l=(TList*)df->Get(Form("clistCheckAODTracks%s",suffix.Data()));

  TH3F* hEtaPhiPtTPCsel=(TH3F*)l->FindObject("hEtaPhiPtTPCsel");
  TH3F* hEtaPhiPtTPCselITSref=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSref");
  TH3F* hEtaPhiPtTPCselSPDany=(TH3F*)l->FindObject("hEtaPhiPtTPCselSPDany");

  Int_t etamin=hEtaPhiPtTPCsel->GetXaxis()->FindBin(-0.799);
  Int_t etamax=hEtaPhiPtTPCsel->GetXaxis()->FindBin(0.799);
  Int_t eta0p=hEtaPhiPtTPCsel->GetXaxis()->FindBin(0.01);
  Int_t eta0m=hEtaPhiPtTPCsel->GetXaxis()->FindBin(-0.01);
  printf("Bin %d Range %f %f\n",etamin,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(etamin),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(etamin));  printf("Bin %d Range %f %f\n",eta0m,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(eta0m),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(eta0m));
  printf("Bin %d Range %f %f\n",eta0p,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(eta0p),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(eta0p));
  printf("Bin %d Range %f %f\n",etamax,hEtaPhiPtTPCsel->GetXaxis()->GetBinLowEdge(etamax),hEtaPhiPtTPCsel->GetXaxis()->GetBinUpEdge(etamax));


  TH1D* hPhiEtaNegTPCsel=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaNegTPCsel",etamin,eta0m);
  TH1D* hPhiEtaPosTPCsel=hEtaPhiPtTPCsel->ProjectionY("hPhiEtaPosTPCsel",eta0p,etamax);
  TH1D* hPtEtaNegTPCsel=hEtaPhiPtTPCsel->ProjectionZ("hPtEtaNegTPCsel",etamin,eta0m);
  TH1D* hPtEtaPosTPCsel=hEtaPhiPtTPCsel->ProjectionZ("hPtEtaPosTPCsel",eta0p,etamax);
  TH1D* hPhiEtaNegTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaNegTPCselITSref",etamin,eta0m);
  TH1D* hPhiEtaPosTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionY("hPhiEtaPosTPCselITSref",eta0p,etamax);
  TH1D* hPtEtaNegTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionZ("hPtEtaNegTPCselITSref",etamin,eta0m);
  TH1D* hPtEtaPosTPCselITSref=hEtaPhiPtTPCselITSref->ProjectionZ("hPtEtaPosTPCselITSref",eta0p,etamax);
  TH1D* hPhiEtaNegTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaNegTPCselSPDany",etamin,eta0m);
  TH1D* hPhiEtaPosTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionY("hPhiEtaPosTPCselSPDany",eta0p,etamax);
  TH1D* hPtEtaNegTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionZ("hPtEtaNegTPCselSPDany",etamin,eta0m);
  TH1D* hPtEtaPosTPCselSPDany=hEtaPhiPtTPCselSPDany->ProjectionZ("hPtEtaPosTPCselSPDany",eta0p,etamax);

  hPhiEtaNegTPCsel->SetMinimum(0);
  hPhiEtaPosTPCsel->SetMinimum(0);
  hPhiEtaNegTPCsel->SetTitle("#varphi tracks - #eta<0");
  hPhiEtaPosTPCsel->SetTitle("#varphi tracks - #eta>0");
  hPhiEtaNegTPCsel->GetXaxis()->SetTitle("#varphi (rad)");
  hPhiEtaPosTPCsel->GetXaxis()->SetTitle("#varphi (rad)");
  hPtEtaNegTPCsel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtEtaPosTPCsel->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtEtaNegTPCsel->SetTitle("p_{T} tracks - #eta<0");
  hPtEtaPosTPCsel->SetTitle("p_{T} tracks - #eta>0");

  hPtEtaNegTPCsel->Sumw2();
  hPtEtaNegTPCselITSref->Sumw2();
  hPtEtaNegTPCselSPDany->Sumw2();
  hPtEtaPosTPCsel->Sumw2();
  hPtEtaPosTPCselITSref->Sumw2();
  hPtEtaPosTPCselSPDany->Sumw2();
  hPhiEtaNegTPCsel->Sumw2();
  hPhiEtaNegTPCselITSref->Sumw2();
  hPhiEtaNegTPCselSPDany->Sumw2();
  hPhiEtaPosTPCsel->Sumw2();
  hPhiEtaPosTPCselITSref->Sumw2();
  hPhiEtaPosTPCselSPDany->Sumw2();

  TCanvas* cdist=new TCanvas("cdist","Distrib",900,900);
  cdist->Divide(2,2);
  cdist->cd(1);
  gPad->SetLogy();
  DrawDistrib(hPtEtaNegTPCsel,hPtEtaNegTPCselITSref,hPtEtaNegTPCselSPDany,kTRUE);
  cdist->cd(2);
  gPad->SetLogy();
  DrawDistrib(hPtEtaPosTPCsel,hPtEtaPosTPCselITSref,hPtEtaPosTPCselSPDany,kTRUE);
  cdist->cd(3);
  DrawDistrib(hPhiEtaNegTPCsel,hPhiEtaNegTPCselITSref,hPhiEtaNegTPCselSPDany,kFALSE);
  cdist->cd(4);
  DrawDistrib(hPhiEtaPosTPCsel,hPhiEtaPosTPCselITSref,hPhiEtaPosTPCselSPDany,kFALSE);

  TH1D* hMatchEffVsPtNegEta=ComputeMatchEff(hPtEtaNegTPCselITSref,hPtEtaNegTPCsel,"hMatchEffVsPtNegEta",1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEta=ComputeMatchEff(hPtEtaPosTPCselITSref,hPtEtaPosTPCsel,"hMatchEffVsPtPosEta",1,20,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtNegEtaSPDany=ComputeMatchEff(hPtEtaNegTPCselSPDany,hPtEtaNegTPCsel,"hMatchEffVsPtNegEtaSPDAny",kBlue-7,33,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPtPosEtaSPDany=ComputeMatchEff(hPtEtaPosTPCselSPDany,hPtEtaPosTPCsel,"hMatchEffVsPtPosEtaSPDAny",kBlue-7,33,"p_{T} (GeV/c)");
  TH1D* hMatchEffVsPhiNegEta=ComputeMatchEff(hPhiEtaNegTPCselITSref,hPhiEtaNegTPCsel,"hMatchEffVsPhiNegEta",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEta=ComputeMatchEff(hPhiEtaPosTPCselITSref,hPhiEtaPosTPCsel,"hMatchEffVsPhiPosEta",1,20,"#varphi (rad)");
  TH1D* hMatchEffVsPhiNegEtaSPDany=ComputeMatchEff(hPhiEtaNegTPCselSPDany,hPhiEtaNegTPCsel,"hMatchEffVsPhiNegEtaSPDAny",kBlue-7,33,"#varphi (rad)");
  TH1D* hMatchEffVsPhiPosEtaSPDany=ComputeMatchEff(hPhiEtaPosTPCselSPDany,hPhiEtaPosTPCsel,"hMatchEffVsPhiPosEtaSPDAny",kBlue-7,33,"#varphi (rad)");

  hMatchEffVsPtNegEta->SetTitle("#eta<0");
  hMatchEffVsPtPosEta->SetTitle("#eta>0");
  hMatchEffVsPhiNegEta->SetTitle("#eta<0");
  hMatchEffVsPhiPosEta->SetTitle("#eta>0");
  hMatchEffVsPtNegEtaSPDany->SetTitle("#eta<0");
  hMatchEffVsPtPosEtaSPDany->SetTitle("#eta>0");
  hMatchEffVsPhiNegEtaSPDany->SetTitle("#eta<0");
  hMatchEffVsPhiPosEtaSPDany->SetTitle("#eta>0");

  Int_t theBin350=hMatchEffVsPtNegEta->GetXaxis()->FindBin(0.35);
  Int_t theBin950=hMatchEffVsPtNegEta->GetXaxis()->FindBin(0.95);
  Int_t theBin3950=hMatchEffVsPtNegEta->GetXaxis()->FindBin(3.95);
  vecForTrend[24]=hMatchEffVsPtPosEta->GetBinContent(theBin350);
  vecForTrend[25]=hMatchEffVsPtPosEta->GetBinContent(theBin950);
  vecForTrend[26]=hMatchEffVsPtPosEta->GetBinContent(theBin3950);
  vecForTrend[27]=hMatchEffVsPtNegEta->GetBinContent(theBin350);
  vecForTrend[28]=hMatchEffVsPtNegEta->GetBinContent(theBin950);
  vecForTrend[29]=hMatchEffVsPtNegEta->GetBinContent(theBin3950);
  vecForTrend[30]=hMatchEffVsPtPosEtaSPDany->GetBinContent(theBin350);
  vecForTrend[31]=hMatchEffVsPtPosEtaSPDany->GetBinContent(theBin950);
  vecForTrend[32]=hMatchEffVsPtPosEtaSPDany->GetBinContent(theBin3950);
  vecForTrend[33]=hMatchEffVsPtNegEtaSPDany->GetBinContent(theBin350);
  vecForTrend[34]=hMatchEffVsPtNegEtaSPDany->GetBinContent(theBin950);
  vecForTrend[35]=hMatchEffVsPtNegEtaSPDany->GetBinContent(theBin3950);
 
  TCanvas* cme=new TCanvas("cme","MatchEff all",900,900);
  cme->Divide(2,2);
  cme->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtNegEta->Draw("PE");
  hMatchEffVsPtNegEtaSPDany->Draw("samepe");
  TLegend* leg=new TLegend(0.27,0.17,0.6,0.39);
  leg->AddEntry(hMatchEffVsPtNegEta,"ITSrefit","P");
  leg->AddEntry(hMatchEffVsPtNegEtaSPDany,"SPD any","P");
  leg->Draw();
  cme->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPtPosEta->Draw("PE");
  hMatchEffVsPtPosEtaSPDany->Draw("samepe");
  cme->cd(3);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPhiNegEta->Draw("PE");
  hMatchEffVsPhiNegEtaSPDany->Draw("samepe");
  cme->cd(4);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTickx();
  gPad->SetTicky();
  hMatchEffVsPhiPosEta->Draw("PE");
  hMatchEffVsPhiPosEtaSPDany->Draw("samepe");
  cme->SaveAs("MatchEff.png");

  TH3F* hEtaPhiPtTPCselITSrefGood=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSrefGood");
  TH3F* hEtaPhiPtTPCselITSrefFake=(TH3F*)l->FindObject("hEtaPhiPtTPCselITSrefFake");
  TH1D* hPtGood=hEtaPhiPtTPCselITSrefGood->ProjectionZ("hPtGood",etamin,eta0m);
  TH1D* hPtFake=hEtaPhiPtTPCselITSrefFake->ProjectionZ("hPtFake",etamin,eta0m);
  TH1D* hPtAll=(TH1D*)hPtGood->Clone("hPtAll");
  TH1F* hratiofake=(TH1F*)hPtFake->Clone("hratiofake");
  if(hPtFake->GetEntries()>0){
    hPtAll->Add(hPtFake);
    hPtAll->SetLineColor(1);
    hPtGood->Sumw2();
    hPtAll->Sumw2();
    hPtFake->Sumw2();
    hratiofake->Divide(hPtFake,hPtAll,1,1,"B");
    hratiofake->SetStats(0);
    hPtAll->SetMinimum(hPtFake->GetMinimum()*0.9);
  }

  TH3F* hImpParXYPtMulTPCselSPDanyGood=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanyGood");
  TH3F* hImpParXYPtMulTPCselSPDanyFake=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanyFake");
  TH1D* hImpParGood=hImpParXYPtMulTPCselSPDanyGood->ProjectionY("hImpParGood");
  TH1D* hImpParFake=hImpParXYPtMulTPCselSPDanyFake->ProjectionY("hImpParFake");
  hImpParGood->SetLineColor(kGreen+1);
  hImpParFake->SetLineColor(2);
  TH1D* hImpParAllGF=(TH1D*)hImpParGood->Clone("hImpParAllGF");
  hImpParAllGF->Add(hImpParFake);
  hImpParAllGF->SetLineColor(1);
  hImpParAllGF->Sumw2();
  hImpParGood->Sumw2();
  hImpParFake->Sumw2();
  TH1F* hratiofakeip=(TH1F*)hImpParFake->Clone("hratiofakeip");
  hratiofakeip->Divide(hImpParFake,hImpParAllGF,1,1,"B");
  hratiofakeip->SetLineColor(1);
  hratiofakeip->SetStats(0);
  if(hImpParFake->Integral()>0 && hImpParGood->Integral()>0 ){
    TCanvas* c1=new TCanvas("c1","FakeGood",1200,900);
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetLogy();
    hPtAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPtAll->Draw("histo");
    gPad->Update();
    TPaveStats* stp1=(TPaveStats*)hPtAll->GetListOfFunctions()->FindObject("stats");
    if(stp1){
      stp1->SetTextColor(1);
      stp1->SetY1NDC(0.73);
      stp1->SetY2NDC(0.92);
    }
    gPad->Modified();
    hPtGood->SetLineColor(kGreen+1);
    hPtGood->Draw("histosames");
    gPad->Update();
    TPaveStats* stp2=(TPaveStats*)hPtGood->GetListOfFunctions()->FindObject("stats");
    stp2->SetTextColor(kGreen+1);
    stp2->SetY1NDC(0.72);
    stp2->SetY2NDC(0.53);
    gPad->Modified();
    hPtFake->SetLineColor(2);
    hPtFake->Draw("histosames");
    gPad->Update();
    TPaveStats* stp3=(TPaveStats*)hPtFake->GetListOfFunctions()->FindObject("stats");
    stp3->SetTextColor(2);
    stp3->SetY1NDC(0.52);
    stp3->SetY2NDC(0.33);
    gPad->Modified();
    c1->cd(2);
    gPad->SetLogy();
    hImpParGood->Draw("histo");
    gPad->Update();
    TPaveStats* st1=(TPaveStats*)hImpParGood->GetListOfFunctions()->FindObject("stats");
    st1->SetTextColor(kGreen+1);
    st1->SetY1NDC(0.73);
    st1->SetY2NDC(0.92);
    gPad->Modified();
    hImpParFake->Draw("histosames"); 
    gPad->Update();
    TPaveStats* st2=(TPaveStats*)hImpParFake->GetListOfFunctions()->FindObject("stats");
    st2->SetTextColor(2);
    st2->SetY1NDC(0.72);
    st2->SetY2NDC(0.53);
    gPad->Modified();
    c1->cd(3);
    hratiofake->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratiofake->GetYaxis()->SetTitle("Fraction of fakes");
    hratiofake->GetYaxis()->SetTitleOffset(1.2);
    hratiofake->Draw();
    c1->cd(4);
    hratiofakeip->GetYaxis()->SetTitle("Fraction of fakes");
    hratiofakeip->GetYaxis()->SetTitleOffset(1.2);
    hratiofakeip->Draw();
    c1->SaveAs("GoodFakeTracks.png");
  }
  
  TH3F* hImpParXYPtMulTPCselSPDanyPrim=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanyPrim");
  TH3F* hImpParXYPtMulTPCselSPDanySecDec=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanySecDec");
  TH3F* hImpParXYPtMulTPCselSPDanySecMat=(TH3F*)l->FindObject("hImpParXYPtMulTPCselSPDanySecMat");
  TH1D* hPtPrim=hImpParXYPtMulTPCselSPDanyPrim->ProjectionX("hPtPrim");
  TH1D* hPtSecDec=hImpParXYPtMulTPCselSPDanySecDec->ProjectionX("hPtSecDec");
  TH1D* hPtSecMat=hImpParXYPtMulTPCselSPDanySecMat->ProjectionX("hPtSecMat");
  TH1D* hImpParPrim=hImpParXYPtMulTPCselSPDanyPrim->ProjectionY("hImpParPrim");
  TH1D* hImpParSecDec=hImpParXYPtMulTPCselSPDanySecDec->ProjectionY("hImpParSecDec");
  TH1D* hImpParSecMat=hImpParXYPtMulTPCselSPDanySecMat->ProjectionY("hImpParSecMat");
  TH1D* hImpParAll=(TH1D*)hImpParSecDec->Clone("hImpParAll");
  hImpParAll->Add(hImpParSecMat);
  hImpParAll->Add(hImpParPrim);
  TH1D* hPtAllPS=(TH1D*)hPtSecDec->Clone("hPtAllPS");
  hPtAllPS->Add(hPtSecMat);
  hPtAllPS->Add(hPtPrim);
  hImpParAll->SetLineColor(1);
  hImpParPrim->SetLineColor(4);
  hImpParSecDec->SetLineColor(kOrange+1);
  hImpParSecMat->SetLineColor(6);
  hImpParPrim->SetLineWidth(2);
  hImpParSecDec->SetLineWidth(2);
  hImpParSecMat->SetLineWidth(2);
  hPtPrim->SetLineColor(4);
  hPtSecDec->SetLineColor(kOrange+1);
  hPtSecMat->SetLineColor(6);
  hPtPrim->SetLineWidth(2);
  hPtSecDec->SetLineWidth(2);
  hPtSecMat->SetLineWidth(2);
  hPtAllPS->Sumw2();
  hPtSecDec->Sumw2();
  hPtSecMat->Sumw2();
  hPtPrim->Sumw2();
  hImpParAll->Sumw2();
  hImpParSecDec->Sumw2();
  hImpParSecMat->Sumw2();
  hImpParPrim->Sumw2();
  TH1F* hratiosecdec=(TH1F*)hPtSecDec->Clone("hratiosecdec");
  hratiosecdec->Divide(hPtSecDec,hPtAllPS,1,1,"B");
  hratiosecdec->SetStats(0);
  TH1F* hratiosecdecip=(TH1F*)hImpParSecDec->Clone("hratiosecdecpi");
  hratiosecdecip->Divide(hImpParSecDec,hImpParAll,1,1,"B");
  hratiosecdecip->SetStats(0);
  TH1F* hratiosecmat=(TH1F*)hPtSecMat->Clone("hratiosecmat");
  hratiosecmat->Divide(hPtSecMat,hPtAllPS,1,1,"B");
  hratiosecmat->SetStats(0);
  TH1F* hratiosecmatip=(TH1F*)hImpParSecMat->Clone("hratiosecmatpi");
  hratiosecmatip->Divide(hImpParSecMat,hImpParAll,1,1,"B");
  hratiosecmatip->SetStats(0);
  if(hImpParSecDec->Integral()>0 && hImpParPrim->Integral()>0 ){
    TCanvas* cps1=new TCanvas("cps1","SecPrim",1200,900);
    cps1->Divide(2,2);
    cps1->cd(1);
    gPad->SetLogy();
    hPtAllPS->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPtAllPS->Draw("histo");
    gPad->Update();
    TPaveStats* stp1=(TPaveStats*)hPtAllPS->GetListOfFunctions()->FindObject("stats");
    if(stp1){
      stp1->SetTextColor(1);
      stp1->SetY1NDC(0.73);
      stp1->SetY2NDC(0.92);
    }
    gPad->Modified();
    hPtPrim->Draw("histosames");
    gPad->Update();
    TPaveStats* stp2=(TPaveStats*)hPtPrim->GetListOfFunctions()->FindObject("stats");
    stp2->SetTextColor(hPtPrim->GetLineColor());
    stp2->SetY1NDC(0.72);
    stp2->SetY2NDC(0.53);
    gPad->Modified();
    hPtSecDec->Draw("histosames");
    gPad->Update();
    TPaveStats* stp3=(TPaveStats*)hPtSecDec->GetListOfFunctions()->FindObject("stats");
    stp3->SetTextColor(hPtSecDec->GetLineColor());
    stp3->SetY1NDC(0.52);
    stp3->SetY2NDC(0.33);
    gPad->Modified();
    hPtSecMat->Draw("histosames");
    gPad->Update();
    TPaveStats* stp4=(TPaveStats*)hPtSecMat->GetListOfFunctions()->FindObject("stats");
    stp4->SetTextColor(hPtSecMat->GetLineColor());
    stp4->SetY1NDC(0.32);
    stp4->SetY2NDC(0.13);
    gPad->Modified();
    cps1->cd(2);
    gPad->SetLogy();
    hImpParPrim->Draw("histo");
    gPad->Update();
    TPaveStats* sti1=(TPaveStats*)hImpParPrim->GetListOfFunctions()->FindObject("stats");
    sti1->SetTextColor(hImpParPrim->GetLineColor());
    sti1->SetY1NDC(0.73);
    sti1->SetY2NDC(0.92);
    gPad->Modified();
    hImpParSecDec->Draw("histosames"); 
    gPad->Update();
    TPaveStats* sti2=(TPaveStats*)hImpParSecDec->GetListOfFunctions()->FindObject("stats");
    sti2->SetTextColor(hImpParSecDec->GetLineColor());
    sti2->SetY1NDC(0.72);
    sti2->SetY2NDC(0.53);
    gPad->Modified();
    hImpParSecMat->Draw("histosames"); 
    gPad->Update();
    TPaveStats* sti3=(TPaveStats*)hImpParSecMat->GetListOfFunctions()->FindObject("stats");
    sti3->SetTextColor(hImpParSecMat->GetLineColor());
    sti3->SetY1NDC(0.52);
    sti3->SetY2NDC(0.33);
    gPad->Modified();
    cps1->cd(3);
    hratiosecdec->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratiosecdec->GetYaxis()->SetTitle("Fraction of secondaries");
    hratiosecdec->GetYaxis()->SetTitleOffset(1.2);
    hratiosecdec->Draw();
    hratiosecmat->Draw("same");
    cps1->cd(4);
    //    hratiosecip->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hratiosecdecip->GetYaxis()->SetTitle("Fraction of secondaries");
    hratiosecdecip->GetYaxis()->SetTitleOffset(1.2);
    hratiosecdecip->Draw();
    hratiosecmatip->Draw("same");
    cps1->SaveAs("PrimSecTracks.gif");
  }
  
  TH2F* hPtResidVsPtTPCselITSrefPion=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefpi");
  TH2F* hPtResidVsPtTPCselITSrefKaon=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefK");
  TH2F* hPtResidVsPtTPCselITSrefProton=(TH2F*)l->FindObject("hPtResidVsPtTPCselITSrefp");  
  if(hPtResidVsPtTPCselITSrefPion){
    hPtResidVsPtTPCselITSrefPion->SetStats(0);
    hPtResidVsPtTPCselITSrefKaon->SetStats(0);
    hPtResidVsPtTPCselITSrefProton->SetStats(0);
    hPtResidVsPtTPCselITSrefPion->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefKaon->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefProton->GetYaxis()->SetTitleOffset(1.5);
    hPtResidVsPtTPCselITSrefPion->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefKaon->GetXaxis()->SetTitleOffset(1.1);
    hPtResidVsPtTPCselITSrefProton->GetXaxis()->SetTitleOffset(1.1);
  }

  TGraphErrors* gMeanPi=new TGraphErrors(0);
  TGraphErrors* gMeanK=new TGraphErrors(0);
  TGraphErrors* gMeanProt=new TGraphErrors(0);
  TGraphErrors* gRmsPi=new TGraphErrors(0);
  TGraphErrors* gRmsK=new TGraphErrors(0);
  TGraphErrors* gRmsProt=new TGraphErrors(0);
  TGraphErrors* gRelPi=new TGraphErrors(0);
  TGraphErrors* gRelK=new TGraphErrors(0);
  TGraphErrors* gRelProt=new TGraphErrors(0);
  gMeanPi->SetName("gMeanProt");
  gMeanK->SetName("gMeanProt");
  gMeanProt->SetName("gMeanProt");
  gMeanPi->SetTitle("");
  gMeanK->SetTitle("");
  gMeanProt->SetTitle("");
  gRmsPi->SetName("gRmsProt");
  gRmsK->SetName("gRmsProt");
  gRmsProt->SetName("gRmsProt");
  gRmsPi->SetTitle("");
  gRmsK->SetTitle("");
  gRmsProt->SetTitle("");
  gRelPi->SetName("gRelProt");
  gRelK->SetName("gRelProt");
  gRelProt->SetName("gRelProt");
  gRelPi->SetTitle("");
  gRelK->SetTitle("");
  gRelProt->SetTitle("");

  if(hPtResidVsPtTPCselITSrefPion && hPtResidVsPtTPCselITSrefPion->Integral()>0){
    FillMeanAndRms(hPtResidVsPtTPCselITSrefPion,gMeanPi,gRmsPi,gRelPi);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefKaon,gMeanK,gRmsK,gRelK);
    FillMeanAndRms(hPtResidVsPtTPCselITSrefProton,gMeanProt,gRmsProt,gRelProt);
    TCanvas* c2=new TCanvas("c2","Pt resol",1500,900);
    c2->Divide(3,2);
    c2->cd(1);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefPion->Draw("colz");
    TLatex* tpi=new TLatex(0.45,0.93,"Pions");
    tpi->SetNDC();
    tpi->Draw();
    c2->cd(2);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefKaon->Draw("colz");
    TLatex* tk=new TLatex(0.45,0.93,"Kaons");
    tk->SetNDC();
    tk->Draw();
    c2->cd(3);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.11);
    hPtResidVsPtTPCselITSrefProton->Draw("colz");
    TLatex* tpr=new TLatex(0.43,0.93,"Protons");
    tpr->SetNDC();
    tpr->Draw();
    c2->cd(4);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gMeanPi->SetMarkerStyle(20);
    gMeanPi->SetMarkerColor(1);
    gMeanPi->SetLineColor(1);
    gMeanK->SetMarkerStyle(25);
    gMeanK->SetMarkerColor(2);
    gMeanK->SetLineColor(2);
    gMeanProt->SetMarkerStyle(27);
    gMeanProt->SetMarkerColor(4);
    gMeanProt->SetLineColor(4);
    gMeanPi->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanK->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanProt->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gMeanPi->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanK->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanProt->GetYaxis()->SetTitle("<p_{T,reco}-p_{T,gen}> (GeV/c)");
    gMeanProt->GetYaxis()->SetTitleOffset(1.6);
    gMeanProt->GetXaxis()->SetTitleOffset(1.2);
    gMeanProt->Draw("AP");
    gMeanPi->Draw("psame");
    gMeanK->Draw("psame");
    c2->cd(5);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gRmsPi->SetMarkerStyle(20);
    gRmsPi->SetMarkerColor(1);
    gRmsPi->SetLineColor(1);
    gRmsK->SetMarkerStyle(25);
    gRmsK->SetMarkerColor(2);
    gRmsK->SetLineColor(2);
    gRmsProt->SetMarkerStyle(27);
    gRmsProt->SetMarkerColor(4);
    gRmsProt->SetLineColor(4);
    gRmsPi->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsK->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsProt->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRmsPi->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsK->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsProt->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen}) (GeV/c)");
    gRmsProt->GetYaxis()->SetTitleOffset(1.6);
    gRmsProt->GetXaxis()->SetTitleOffset(1.2);
    gRmsProt->Draw("AP");
    gRmsProt->Draw("AP");
    gRmsPi->Draw("psame");
    gRmsK->Draw("psame");
    c2->cd(6);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetTickx();
    gPad->SetTicky();
    gRelPi->SetMarkerStyle(20);
    gRelPi->SetMarkerColor(1);
    gRelPi->SetLineColor(1);
    gRelK->SetMarkerStyle(25);
    gRelK->SetMarkerColor(2);
    gRelK->SetLineColor(2);
    gRelProt->SetMarkerStyle(27);
    gRelProt->SetMarkerColor(4);
    gRelProt->SetLineColor(4);
    gRelPi->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRelK->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRelProt->GetXaxis()->SetTitle("p_{T,gen} (GeV/c)");
    gRelPi->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen})/p_{T,gen}");
    gRelK->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen})/p_{T,gen}");
    gRelProt->GetYaxis()->SetTitle("#sigma(p_{T,reco}-p_{T,gen})/p_{T,gen}");
    gRelProt->GetYaxis()->SetTitleOffset(1.6);
    gRelProt->GetXaxis()->SetTitleOffset(1.2);
    gRelProt->Draw("AP");
    gRelPi->Draw("psame");
    gRelK->Draw("psame");
    c2->SaveAs("PtResisuals.png");
  }

  TH2F* hFilterBits=(TH2F*)l->FindObject("hFilterBits");
  hFilterBits->SetStats(0);
  Double_t cnt0=hFilterBits->GetBinContent(1,1)+hFilterBits->GetBinContent(1,2);
  if(cnt0>0){
    for(Int_t jbit=0; jbit<12; jbit++){
      Double_t cntj=hFilterBits->GetBinContent(jbit+1,1)+hFilterBits->GetBinContent(jbit+1,2);
      vecForTrend[2*jbit+1]=cntj/cnt0;
    }
  }

  TCanvas* cip=new TCanvas("cip","FiltBits",1100,900);
  cip->Divide(2,2);
  cip->cd(1);
  gPad->SetRightMargin(0.13);
  hFilterBits->Draw("colz");
  cip->cd(2);
  gPad->SetLogy();
  Int_t colors[12]={kRed+1,kRed-7,kOrange+1,kYellow+1,kGreen+1,kGreen,kCyan,kBlue+1,kMagenta,kMagenta+1,kGray+1,1};
  Int_t lstyl[12]={1,9,1,3,1,8,2,5,7,1,1,9};
  Int_t lwid[12]={2,2,2,3,2,2,3,3,3,2,2,2};
  TLegend* leg2=new TLegend(0.65,0.2,0.89,0.89);
  leg2->SetMargin(0.3);
  for(Int_t jbit=0; jbit<12; jbit++){
    TString hname=Form("hImpParXYPtMulPionFiltBit%d",jbit);
    TH3F* h=(TH3F*)l->FindObject(hname.Data());
    TH1D* h1=h->ProjectionY(Form("hImpParXYFiltBit%d",jbit));
    h1->SetLineColor(colors[jbit]);
    h1->SetLineStyle(lstyl[jbit]);
    h1->SetLineWidth(lwid[jbit]);
    h1->SetStats(0);
    if(jbit==0){ 
      h1->DrawNormalized();
    }
    else h1->DrawNormalized("same");
    leg2->AddEntry(h1,Form("Filt Bit %d",jbit),"L")->SetTextColor(colors[jbit]);
  }

  for(Int_t jbit=0; jbit<12; jbit++){
    TString hname=Form("hEtaPhiPtFiltBit%d",jbit);
    TH3F* h=(TH3F*)l->FindObject(hname.Data());
    TH1D* hphi=h->ProjectionY(Form("hPhiFiltBit%d",jbit));
    hphi->SetLineColor(colors[jbit]);
    hphi->SetLineStyle(lstyl[jbit]);
    hphi->SetLineWidth(lwid[jbit]);
    hphi->SetStats(0);
    TH1D* hpt1=h->ProjectionZ(Form("hPtFiltBit%d",jbit));
    vecForTrend[2*jbit]=hpt1->GetMean();
    hpt1->SetLineColor(colors[jbit]);
    hpt1->SetLineStyle(lstyl[jbit]);
    hpt1->SetLineWidth(lwid[jbit]);
    hpt1->SetStats(0);
    cip->cd(3);
    if(jbit==0){ 
      hphi->SetMinimum(0.);
      hphi->SetMaximum(1.2*hphi->GetMaximum());
      hphi->Draw();
    }else{ 
      hphi->Draw("same");
    }
    cip->cd(4);
    if(jbit==0){ 
      hpt1->SetMinimum(0.);
      hpt1->SetMaximum(1.4*hpt1->GetMaximum());
      hpt1->Draw();
    }else{ 
      hpt1->Draw("same");
    }
    leg2->Draw();
  }
  cip->SaveAs("FilterBits.png");


  for(Int_t jb=0; jb<12; jb++){
    TCanvas* ccc=new TCanvas(Form("cfb%d",jb),Form("cfb%d",jb),1500,800);
    ccc->Divide(3,2);
    ccc->cd(1);
    gPad->SetLogz();
    TH2F* hdist1=(TH2F*)l->FindObject(Form("hITScluPtFiltBit%d",jb));
    hdist1->SetStats(0);
    hdist1->SetTitle(Form("Filter bit %d",jb ));
    hdist1->Draw("colz");
    ccc->cd(2);
    gPad->SetLogz();
    TH2F* hdist2=(TH2F*)l->FindObject(Form("hSPDcluPtFiltBit%d",jb));
    hdist2->SetTitle(Form("Filter bit %d",jb ));
    hdist2->SetStats(0);
    hdist2->Draw("colz");
    ccc->cd(3);
    gPad->SetLogz();
    TH2F* hdist3=(TH2F*)l->FindObject(Form("hTPCcluPtFiltBit%d",jb));
    hdist3->SetTitle(Form("Filter bit %d",jb ));
    hdist3->SetStats(0);
    hdist3->Draw("colz");
    ccc->cd(4);
    gPad->SetLogz();
    TH2F* hdist4=(TH2F*)l->FindObject(Form("hTPCcrrowsPtFiltBit%d",jb));
    hdist4->SetTitle(Form("Filter bit %d",jb ));
    hdist4->SetStats(0);
    hdist4->Draw("colz");
    ccc->cd(5);
    gPad->SetLogz();
    TH2F* hdist5=(TH2F*)l->FindObject(Form("hTPCCrowOverFindPtFiltBit%d",jb));
    hdist5->SetTitle(Form("Filter bit %d",jb ));
    hdist5->GetYaxis()->SetTitle("Crossed rows /findable clusters");
    hdist5->SetStats(0);
    hdist5->Draw("colz");
    ccc->cd(6);
    gPad->SetLogz();
    TH2F* hdist6=(TH2F*)l->FindObject(Form("hTPCChi2ndfPtFiltBit%d",jb));
    hdist6->SetTitle(Form("Filter bit %d",jb ));
    hdist6->SetStats(0);
    hdist6->Draw("colz");
    gPad->SetLogz();
    ccc->SaveAs(Form("VarDistFiltBit%d.png",jb));
  }
  TH2F*	hInvMassK0s=(TH2F*)l->FindObject("hInvMassK0s");
  TH3F*	hInvMassLambda=(TH3F*)l->FindObject("hInvMassLambda");
  TH3F*	hInvMassAntiLambda=(TH3F*)l->FindObject("hInvMassAntiLambda");
  TF1* fmass=new TF1("fmass","[0]+[1]*x+[2]/sqrt(2.*TMath::Pi())/[4]*TMath::Exp(-0.5*(x-[3])*(x-[3])/[4]/[4])",1.11,1.122);
  fmass->SetLineWidth(2);
  fmass->SetLineColor(kRed+1);
  TF1* fmassk0=new TF1("fmassk0","[0]+[1]*x+[2]/sqrt(2.*TMath::Pi())/[4]*TMath::Exp(-0.5*(x-[3])*(x-[3])/[4]/[4])",0.46,0.52);
  fmassk0->SetLineWidth(2);
  fmassk0->SetLineColor(kMagenta+1);

  TCanvas* ck0=new TCanvas("ck0","K0s",1300,800);
  ck0->Divide(3,2);
  ck0->cd(1);
  ProjectDrawAndFitK0(hInvMassK0s,fmassk0,0.6,0.8);
  ck0->cd(2);
  ProjectDrawAndFitK0(hInvMassK0s,fmassk0,0.8,1.0);
  ck0->cd(3);
  ProjectDrawAndFitK0(hInvMassK0s,fmassk0,1.,2.);
  ck0->cd(4);
  ProjectDrawAndFitK0(hInvMassK0s,fmassk0,2.,3.);
  ck0->cd(5);
  ProjectDrawAndFitK0(hInvMassK0s,fmassk0,3.,5.);
  ck0->cd(6);
  ProjectDrawAndFitK0(hInvMassK0s,fmassk0,5.,10.);
  ck0->SaveAs("K0InvMass.png");

  TCanvas* clam=new TCanvas("clam","lambda",1300,800);
  clam->Divide(3,2);
  clam->cd(1);
  ProjectDrawAndFit(hInvMassLambda,"Lambda",fmass,0.6,0.8);
  clam->cd(2);
  ProjectDrawAndFit(hInvMassLambda,"Lambda",fmass,0.8,1.0);
  clam->cd(3);
  ProjectDrawAndFit(hInvMassLambda,"Lambda",fmass,1.0,2.0);
  clam->cd(4);
  ProjectDrawAndFit(hInvMassAntiLambda,"AntiLambda",fmass,0.6,0.8);
  clam->cd(5);
  ProjectDrawAndFit(hInvMassAntiLambda,"AntiLambda",fmass,0.8,1.0);
  clam->cd(6);
  ProjectDrawAndFit(hInvMassAntiLambda,"AntiLambda",fmass,1.,2.0);
  clam->SaveAs("LambdaInvMass.png");

  trtree->Fill();

  TFile* fouttree=new TFile("trending.root","recreate");
  trtree->Write();
  fouttree->Close();
  delete fouttree;

}

void ProjectDrawAndFitK0(TH2F* hInvMass2d, TF1* fmassk0, Double_t xmin, Double_t xmax){
 
  TH1D* hInvMass=0x0;

  Int_t b0=0;
  Int_t b1=9999;
  TString ptrange="";
  if(xmin<-0.1 || xmax<-0.1){
    hInvMass=hInvMass2d->ProjectionX("hInvMassK0allpt");
    ptrange="all p_{T}";
  }else{
    b0=hInvMass2d->GetYaxis()->FindBin(xmin+0.001);
    b1=hInvMass2d->GetYaxis()->FindBin(xmax-0.001);
    printf("Inv mass limits = %f-%f --> bins %d %d\n",xmin,xmax,b0,b1);
    hInvMass=hInvMass2d->ProjectionX(Form("hInvMassK0%d%d",b0,b1),b0,b1);
    ptrange=Form("%.1f<p_{T}^{K0}<%.1f GeV/c",xmin,xmax);
  }    
  hInvMass->GetXaxis()->SetTitle("Invariant mass #pi^{+}#pi^{-} (GeV/c^{2})");
  hInvMass->SetTitle(Form("K0s - %s",ptrange.Data()));
  hInvMass->GetYaxis()->SetTitle("Entries");

  hInvMass->Draw();
  if(hInvMass->GetEntries()>200000){
    fmassk0->SetParameter(0,hInvMass->GetBinContent(1));
    fmassk0->SetParameter(1,0.);
    fmassk0->SetParameter(2,10000.);
    fmassk0->SetParameter(3,0.5);
    fmassk0->SetParLimits(3,0.49,0.51);
    fmassk0->SetParameter(4,0.002);
    fmassk0->SetParLimits(4,0.0006,0.01);
    hInvMass->Fit("fmassk0","R");
    TLatex* t1=new TLatex(0.55,0.7,Form("Mean = %.3f+-%.3f GeV/c^{2}",fmassk0->GetParameter(3),fmassk0->GetParError(3)));
    t1->SetTextSize(0.04);
    t1->SetNDC();
    t1->Draw();
    TLatex* t2=new TLatex(0.55,0.65,Form("Sigma = %.2f+-%.2f MeV/c^{2}",fmassk0->GetParameter(4)*1000.,fmassk0->GetParError(4)*1000.));
    t2->SetNDC();
    t2->SetTextSize(0.04);
    t2->Draw();
    TLatex* t3=new TLatex(0.55,0.6,Form("Yield = %.0f+-%.0f",fmassk0->GetParameter(2)/hInvMass->GetBinWidth(1),fmassk0->GetParError(2)/hInvMass->GetBinWidth(1)));
    t3->SetNDC();
    t3->SetTextSize(0.04);
    t3->Draw();
  }


}

void ProjectDrawAndFit(TH3F* hInvMass3d, TString partName, TF1* fmass, Double_t xmin, Double_t xmax){
 
  TH1D* hInvMass=0x0;

  Int_t b0=0;
  Int_t b1=9999;
  TString ptrange="";
  if(xmin<-0.1 || xmax<-0.1){
    hInvMass=hInvMass3d->ProjectionX(Form("hInvMass%sallpt",partName.Data()));
    ptrange="all p_{T}";
  }else{
    b0=hInvMass3d->GetYaxis()->FindBin(xmin+0.001);
    b1=hInvMass3d->GetYaxis()->FindBin(xmax-0.001);
    printf("Inv mass limits = %f-%f --> bins %d %d\n",xmin,xmax,b0,b1);
    hInvMass=hInvMass3d->ProjectionX(Form("hInvMass%s%d%d",partName.Data(),b0,b1),b0,b1);
    ptrange=Form("%.1f<p_{T}^{#Lambda}<%.1f GeV/c",xmin,xmax);
  }    
  if(partName=="Lambda"){
    hInvMass->GetXaxis()->SetTitle("Invariant mass p#pi^{-} (GeV/c^{2})");
    hInvMass->SetTitle(Form("#Lambda - %s",ptrange.Data()));
  }else{
    hInvMass->GetXaxis()->SetTitle("Invariant mass #bar{p}#pi^{+} (GeV/c^{2})");
    hInvMass->SetTitle(Form("#bar{#Lambda} - %s",ptrange.Data()));
  }
  hInvMass->GetYaxis()->SetTitle("Entries");

  hInvMass->Draw();
  if(hInvMass->GetEntries()>200000){
    InitFuncAndFit(hInvMass,fmass);
  }
}

void ProjectDrawAndFitMomP(TH3F* hInvMassLambdaGoodHyp3d,TH3F* hInvMassLambdaBadHyp3d,TH3F* hInvMassAntiLambdaGoodHyp3d,TH3F* hInvMassAntiLambdaBadHyp3d,TF1* fmass, Double_t xmin, Double_t xmax){
 
  TH1D* hInvMassLambdaGoodHyp=0x0;
  TH1D* hInvMassLambdaBadHyp=0x0;
  TH1D* hInvMassAntiLambdaGoodHyp=0x0;
  TH1D* hInvMassAntiLambdaBadHyp=0x0;
  Int_t b0=0;
  Int_t b1=9999;
  TCanvas* clam=0x0;
  if(xmin<-0.1 || xmax<-0.1){
    clam=new TCanvas("clallprpt","lambda all pt",1100,900);
    hInvMassLambdaGoodHyp=hInvMassLambdaGoodHyp3d->ProjectionX("hInvMassLambdaGoodHypallpt");
    hInvMassLambdaBadHyp=hInvMassLambdaBadHyp3d->ProjectionX("hInvMassLambdaBadHypallpt");
    hInvMassAntiLambdaGoodHyp=hInvMassAntiLambdaGoodHyp3d->ProjectionX("hInvMassAntiLambdaGoodHypallpt");
    hInvMassAntiLambdaBadHyp=hInvMassAntiLambdaBadHyp3d->ProjectionX("hInvMassAntiLambdaBadHypallpt");
    hInvMassLambdaGoodHyp->SetTitle("#Lambda - good proton - all p_{T}");
    hInvMassLambdaBadHyp->SetTitle("#Lambda - bad proton - all p_{T}");
    hInvMassAntiLambdaGoodHyp->SetTitle("#bar{#Lambda} - good proton - all p_{T}");
    hInvMassAntiLambdaBadHyp->SetTitle("#bar{#Lambda} - bad proton - all p_{T}");
  }else{
    b0=hInvMassLambdaGoodHyp3d->GetYaxis()->FindBin(xmin+0.001);
    b1=hInvMassLambdaGoodHyp3d->GetYaxis()->FindBin(xmax-0.001);
    printf("Inv mass limits = %f-%f --> bins %d %d\n",xmin,xmax,b0,b1);

    clam=new TCanvas(Form("clprpt%d%d",b0,b1),Form("lambdavsppt%d%d",b0,b1),1100,900);
    hInvMassLambdaGoodHyp=hInvMassLambdaGoodHyp3d->ProjectionX(Form("hInvMassLambdaGoodHypPPt%d%d",b0,b1),0,-1,b0,b1);
    hInvMassLambdaBadHyp=hInvMassLambdaBadHyp3d->ProjectionX(Form("hInvMassLambdaBadHypPPt%d%d",b0,b1),0,-1,b0,b1);
    hInvMassAntiLambdaGoodHyp=hInvMassAntiLambdaGoodHyp3d->ProjectionX(Form("hInvMassAntiLambdaGoodHypPPt%d%d",b0,b1),0,-1,b0,b1);
    hInvMassAntiLambdaBadHyp=hInvMassAntiLambdaBadHyp3d->ProjectionX(Form("hInvMassAntiLambdaBadHypPPt%d%d",b0,b1),0,-1,b0,b1);
    hInvMassLambdaGoodHyp->SetTitle(Form("#Lambda - good proton - %.1f<p_{T}^{p,TPC}<%.1f GeV/c",xmin,xmax));
    hInvMassLambdaBadHyp->SetTitle(Form("#Lambda - bad proton - %.1f<p_{T}^{p,TPC}<%.1f GeV/c",xmin,xmax));
    hInvMassAntiLambdaGoodHyp->SetTitle(Form("#bar{#Lambda} - good proton - %.1f<p_{T}^{p,TPC}<%.1f GeV/c",xmin,xmax));
    hInvMassAntiLambdaBadHyp->SetTitle(Form("#bar{#Lambda} - bad proton - %.1f<p_{T}^{p,TPC}<%.1f GeV/c",xmin,xmax));
  }    
  hInvMassLambdaGoodHyp->GetXaxis()->SetTitle("Invariant mass p#pi^{-} (GeV/c^{2})");
  hInvMassLambdaBadHyp->GetXaxis()->SetTitle("Invariant mass p#pi^{-} (GeV/c^{2})");
  hInvMassAntiLambdaGoodHyp->GetXaxis()->SetTitle("Invariant mass #bar{p}#pi^{+} (GeV/c^{2})");
  hInvMassAntiLambdaBadHyp->GetXaxis()->SetTitle("Invariant mass #bar{p}#pi^{+} (GeV/c^{2})");
  hInvMassLambdaGoodHyp->GetYaxis()->SetTitle("Entries");
  hInvMassLambdaBadHyp->GetYaxis()->SetTitle("Entries");
  hInvMassAntiLambdaGoodHyp->GetYaxis()->SetTitle("Entries");
  hInvMassAntiLambdaBadHyp->GetYaxis()->SetTitle("Entries");

  clam->Divide(2,2);
  clam->cd(1);
  hInvMassLambdaGoodHyp->Draw();
  if(hInvMassLambdaGoodHyp->GetEntries()>200000){
    InitFuncAndFit(hInvMassLambdaGoodHyp,fmass);
  }
  clam->cd(2);
  hInvMassLambdaBadHyp->Draw();
  if(hInvMassLambdaBadHyp->GetEntries()>200000){
    InitFuncAndFit(hInvMassLambdaBadHyp,fmass);
    hInvMassLambdaBadHyp->Fit("fmass","R");
  }
  clam->cd(3);
  hInvMassAntiLambdaGoodHyp->Draw();
  if(hInvMassAntiLambdaGoodHyp->GetEntries()>200000){
    InitFuncAndFit(hInvMassAntiLambdaGoodHyp,fmass);
  }
  clam->cd(4);
  hInvMassAntiLambdaBadHyp->Draw();
  if(hInvMassAntiLambdaBadHyp->GetEntries()>200000){
    InitFuncAndFit(hInvMassAntiLambdaGoodHyp,fmass);
  }
  clam->SaveAs(Form("LambdaInvMass-ptp%d%d.png",b0,b1));
  hInvMassLambdaGoodHyp->Write();
  hInvMassLambdaBadHyp->Write();
  hInvMassAntiLambdaGoodHyp->Write();
  hInvMassAntiLambdaBadHyp->Write();
}

void InitFuncAndFit(TH1D* hm, TF1* fmass){
  fmass->SetParameter(0,hm->GetBinContent(hm->FindBin(1.10)));
  fmass->SetParameter(1,0.);
  //  fmass->SetParLimits(1,-99999999999,0.);
  fmass->SetParameter(2,100.);
  fmass->SetParLimits(2,0.,1.e+9);
  fmass->SetParameter(3,1.116);
  fmass->SetParLimits(3,1.11,1.12);
  fmass->SetParameter(4,0.002);
  fmass->SetParLimits(4,0.0006,0.003);
  hm->Fit("fmass","R");
  TLatex* t1=new TLatex(0.14,0.8,Form("Mean = %.3f+-%.3f GeV/c^{2}",fmass->GetParameter(3),fmass->GetParError(3)));
  t1->SetTextSize(0.04);
  t1->SetNDC();
  t1->Draw();
  TLatex* t2=new TLatex(0.14,0.75,Form("Sigma = %.2f+-%.2f MeV/c^{2}",fmass->GetParameter(4)*1000.,fmass->GetParError(4)*1000.));
  t2->SetNDC();
  t2->SetTextSize(0.04);
  t2->Draw();
  TLatex* t3=new TLatex(0.14,0.7,Form("Yield = %.0f+-%.0f",fmass->GetParameter(2)/hm->GetBinWidth(1),fmass->GetParError(2)/hm->GetBinWidth(1)));
  t3->SetNDC();
  t3->SetTextSize(0.04);
  t3->Draw();

}


void DrawDistrib(TH1D* h1, TH1D* h2, TH1D* h3, Bool_t showStat){
  if(!showStat){
    h1->SetStats(0);
  }
  h1->SetLineColor(1);
  h1->Draw();
  if(showStat){
    gPad->Update();
    TPaveStats* tp1=(TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    tp1->SetY1NDC(0.73);
    tp1->SetY2NDC(0.92);
    gPad->Modified();
  }
  h2->SetLineColor(2);
  if(showStat){
    h2->Draw("sames"); 
    gPad->Update();
    TPaveStats* tp2=(TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
    tp2->SetY1NDC(0.53);
    tp2->SetY2NDC(0.72);
    tp2->SetTextColor(2);
    gPad->Modified();
  }else{
    h2->Draw("same"); 
  }
  h3->SetLineColor(4);
  if(showStat){
    h3->Draw("sames");
    gPad->Update();
    TPaveStats* tp3=(TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
    tp3->SetY1NDC(0.33);
    tp3->SetY2NDC(0.52);
    tp3->SetTextColor(4);
    gPad->Modified();
  }else{
    h3->Draw("same");
  }
}

TH1D* ComputeMatchEff(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle){
  TH1D* hratio=(TH1D*)hnumer->Clone(name.Data());
  hratio->Divide(hnumer,hdenom,1.,1.,"B");
  hratio->SetLineColor(iCol);
  hratio->SetMarkerColor(iCol);
  hratio->SetMarkerStyle(iMarker);
  hratio->GetXaxis()->SetTitle(xtitle.Data());
  hratio->GetYaxis()->SetTitle("TPC-ITS matching efficiency");
  hratio->GetYaxis()->SetTitleOffset(1.25);
  hratio->SetMinimum(0);
  hratio->SetMaximum(1.05);
  hratio->SetStats(0);
  return hratio;
}

TH1D* ComputeRatio(TH1D* hnumer, TH1D* hdenom, TString name, Int_t iCol, Int_t iMarker, TString xtitle){
  TH1D* hratio=(TH1D*)hnumer->Clone(name.Data());
  TH1D* hsum=(TH1D*)hnumer->Clone("tmp");
  hsum->Add(hnumer,hdenom);
  hratio->Divide(hnumer,hsum,1.,1.,"B");
  hratio->SetLineColor(iCol);
  hratio->SetMarkerColor(iCol);
  hratio->SetMarkerStyle(iMarker);
  hratio->GetXaxis()->SetTitle(xtitle.Data());
  hratio->GetYaxis()->SetTitle("Fraction of protons tracked with wrong mass");
  hratio->GetYaxis()->SetTitleOffset(1.25);
  hratio->SetMinimum(0.0005);
  hratio->SetMaximum(1.05);
  hratio->SetStats(0);
  delete hsum;
  return hratio;
}


void FillMeanAndRms(TH2F* h2d, TGraphErrors* gMean, TGraphErrors* gRms, TGraphErrors* gRel){
  Int_t jpt=0;
  for(Int_t j=1; j<=h2d->GetNbinsX(); j++){
    TH1D* htmp=h2d->ProjectionY("htmp",j,j);
    Double_t pt=h2d->GetXaxis()->GetBinCenter(j);
    if(htmp->Integral()>20){
      Double_t m=htmp->GetMean();
      Double_t em=htmp->GetMeanError();
      Double_t r=htmp->GetRMS();
      Double_t er=htmp->GetRMSError();
      gMean->SetPoint(jpt,pt,m);
      gMean->SetPointError(jpt,0.,em);
      gRms->SetPoint(jpt,pt,r);
      gRms->SetPointError(jpt,0.,er);
      gRel->SetPoint(jpt,pt,r/pt);
      gRel->SetPointError(jpt,0.,er/pt);
      ++jpt;
      delete htmp;
    }
  }
}
 
