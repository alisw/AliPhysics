#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TLegendEntry.h>
#endif

enum EPtWei{kFONLL5overLHC13d3,kFONLL7overLHC10f6a,kNoWei};

const Int_t nPtBins=6;
Double_t binLims[nPtBins+1]={0.,1.,2.,3.,4.,6.,8.};

TString fileNameMC="AnalysisResultsMC_train434432.root";
TString suffix="c3SigPID_Pt400_EM1";
TString fileNameToy="Acceptance_Toy_D0Kpi_yfid08_etaDau09_ptDau100_FONLL7ptshape.root";
Int_t ptWeight=kFONLL5overLHC13d3;



void ComputeEfficiencyFromCombinHF(){


  // multiplicity weights
  TString fileWeightName="trackletsWeightsMultInt_LHC13d3_08092014.root";
  TString histoWeightName[3];
  histoWeightName[0]="hNtrUnCorrEvWithCand";
  histoWeightName[1]="hNtrUnCorrEvWithD";
  histoWeightName[2]="hNtrUnCorrEvSel";
  Int_t wcol[3]={kRed+1,kGreen+1,4};
  Int_t wmark[3]={22,23,26};


  // pt wieghts  
  TF1* funcPtWeight=0x0;
  if(ptWeight==kFONLL5overLHC13d3){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,30.);
    funcPtWeight->SetParameters(2.94999e+00,3.47032e+00,2.81278e+00,2.5,1.93370e-02,3.86865e+00,-1.54113e-01,8.86944e-02,2.56267e-02);
  }else if (ptWeight==kFONLL7overLHC10f6a){
    funcPtWeight=new TF1("funcPtWeight","([0]*x)/TMath::Power([2],(1+TMath::Power([3],x/[1])))+[4]*TMath::Exp([5]+[6]*x)+[7]*TMath::Exp([8]*x)",0.15,40.);
    funcPtWeight->SetParameters(2.41522e+01,4.92146e+00,6.72495e+00,2.5,6.15361e-03,4.78995e+00,-4.29135e-01,3.99421e-01,-1.57220e-02);
  }else{
    funcPtWeight=new TF1("funcPtWeight","[0]");
    funcPtWeight->SetParameter(0,1.);
  }

  Int_t ptcol[nPtBins]={1,kRed+1,kGreen+2,4,kOrange+2,kMagenta};

  
  TString dirName=Form("PWG3_D2H_InvMassDzeroLowPt%s",suffix.Data());
  TString lstName=Form("coutputDzero%s",suffix.Data());
  TFile* fil=new TFile(fileNameMC.Data());
  TDirectoryFile* df=(TDirectoryFile*)fil->Get(dirName.Data());
  TList* l=(TList*)df->Get(lstName.Data());
  //  l->ls();

  TH3F* hPtVsYVsMultReco=(TH3F*)l->FindObject("hPtVsYVsMultReco");
  TH3F* hPtVsYVsMultGenAcc=(TH3F*)l->FindObject("hPtVsYVsMultGenAcc");
  TH3F* hPtVsYVsMultGenLimAcc=(TH3F*)l->FindObject("hPtVsYVsMultGenLimAcc");

  TH2D* hypt=(TH2D*)hPtVsYVsMultGenAcc->Project3D("yx");
  hypt->SetTitle("Generated in acceptance, all multiplcities");
  TH2D* hptmult=(TH2D*)hPtVsYVsMultGenLimAcc->Project3D("xz");
  hptmult->SetTitle("Generated in |y|<0.5");
  hptmult->SetStats(0);
  hypt->SetStats(0);
  hypt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hypt->GetYaxis()->SetTitle("y");
  hptmult->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  hptmult->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  
  TProfile* hMeanPtMult=hptmult->ProfileX("hMeanPtMult");
  TCanvas* c0=new TCanvas("c0","2D plots",1200,600);
  c0->Divide(2,1);
  c0->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.12);
  hypt->Draw("colz");
  c0->cd(2);
  gPad->SetLogz();
  gPad->SetRightMargin(0.13);
  hptmult->Draw("colz");
  hMeanPtMult->Draw("same");

  
  TH1D* hPtReco=hPtVsYVsMultReco->ProjectionX("hPtReco");
  TH1D* hPtGenAcc=hPtVsYVsMultGenAcc->ProjectionX("hPtGenAcc");
  TH1D* hPtGenLimAcc=hPtVsYVsMultGenLimAcc->ProjectionX("hPtGenLimAcc");
  TH1D* hPtRecoR=(TH1D*)hPtReco->Rebin(nPtBins,"hPtRecoReb",binLims);
  TH1D* hPtGenAccR=(TH1D*)hPtGenAcc->Rebin(nPtBins,"hPtRecoReb",binLims);
  TH1D* hPtGenLimAccR=(TH1D*)hPtGenLimAcc->Rebin(nPtBins,"hPtRecoReb",binLims);
  TH1D* hEffVsPt=(TH1D*)hPtReco->Clone("hEff");
  hEffVsPt->Divide(hPtReco,hPtGenAcc,1,1,"B");
  TH1D* hEffVsPtR=(TH1D*)hPtRecoR->Clone("hEff");
  hEffVsPtR->Divide(hPtRecoR,hPtGenAccR,1,1,"B");
  hEffVsPt->SetStats(0);
  hEffVsPtR->SetStats(0);
  TH1D* hAccVsPt=(TH1D*)hPtGenAcc->Clone("hAcc");
  hAccVsPt->Divide(hPtGenAcc,hPtGenLimAcc,1,1,"B");
  TH1D* hAccVsPtR=(TH1D*)hPtGenAccR->Clone("hAcc");
  hAccVsPtR->Divide(hPtGenAccR,hPtGenLimAccR,1,1,"B");
  hAccVsPt->SetStats(0);
  hAccVsPtR->SetStats(0);
 
  TFile* fileAccToy=new TFile(fileNameToy.Data());
  TH1F* hPtGenAccToy=(TH1F*)fileAccToy->Get("hPtGenAcc");
  TH1F* hPtGenLimAccToy=(TH1F*)fileAccToy->Get("hPtGenLimAcc");
  TH1F* hAccToyFine=(TH1F*)fileAccToy->Get("hAccVsPt");
  TH1F* hPtGenAccToyR=(TH1F*)hPtGenAccToy->Rebin(nPtBins,"hPtGenAccToyReb",binLims);
  TH1F* hPtGenLimAccToyR=(TH1F*)hPtGenLimAccToy->Rebin(nPtBins,"hPtGenLimAccToyReb",binLims);
  TH1F* hAccToy=(TH1F*)hPtGenAccToyR->Clone("hAccToy");
  hPtGenAccToyR->Sumw2();
  hPtGenLimAccToyR->Sumw2();
  hAccToy->Divide(hPtGenAccToyR,hPtGenLimAccToyR,1,1,"B");
  hAccToy->SetLineColor(kGreen+2);
  hAccToy->SetLineWidth(3);
  hAccToy->SetMarkerStyle(25);
  hAccToy->SetMarkerColor(kGreen+2);
  hAccToy->SetStats(0);
  hAccToyFine->SetLineColor(kGreen+2);

  TCanvas* c1=new TCanvas("c1","EffVsPt",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy();
  hPtGenLimAcc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtGenLimAcc->GetYaxis()->SetTitle("Entries");
  hPtGenLimAcc->SetLineColor(1);
  hPtGenLimAcc->Draw();
  gPad->Update();
  TPaveStats *st1=(TPaveStats*)hPtGenLimAcc->GetListOfFunctions()->FindObject("stats");
  st1->SetY1NDC(0.7);
  st1->SetY2NDC(0.89);
  hPtGenAcc->SetLineColor(2);
  hPtGenAcc->Draw("sames");
  gPad->Update();
  TPaveStats *st2=(TPaveStats*)hPtGenAcc->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.5);
  st2->SetY2NDC(0.69);
  st2->SetTextColor(2);
  hPtReco->SetLineColor(4);
  hPtReco->Draw("sames");
  gPad->Update();
  TPaveStats *st3=(TPaveStats*)hPtReco->GetListOfFunctions()->FindObject("stats");
  st3->SetY1NDC(0.3);
  st3->SetY2NDC(0.49);
  st3->SetTextColor(4);
  gPad->Modified();
  c1->cd(2);
  hEffVsPt->SetLineColor(4);
  hEffVsPt->SetMinimum(0);
  hEffVsPt->SetMaximum(1.95);
  hEffVsPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hEffVsPt->GetYaxis()->SetTitle("Ratio");
  hEffVsPt->Draw();
  hAccVsPt->SetLineColor(2);
  hAccVsPt->Draw("same");
  hAccToyFine->Draw("same");
  hAccToy->Draw("same");
  TLatex* tacc=new TLatex(0.16,0.83,"Acceptance (CombinHF, 13d3)");
  tacc->SetNDC();
  tacc->SetTextColor(hAccVsPt->GetLineColor());
  tacc->Draw();
  TLatex* tacct=new TLatex(0.16,0.76,"Acceptance (ToyMC)");
  tacct->SetNDC();
  tacct->SetTextColor(hAccToyFine->GetLineColor());
  tacct->Draw();
  TLatex* te=new TLatex(0.16,0.69,"Efficiency");
  te->SetNDC();
  te->SetTextColor(hEffVsPt->GetLineColor());
  te->Draw();
  c1->SaveAs(Form("figures/EfficVsPt_%s.eps",suffix.Data()));

  TH1D* hMultRecoAllPt=hPtVsYVsMultReco->ProjectionZ("hMultReco");
  TH1D* hMultGenAccAllPt=hPtVsYVsMultGenAcc->ProjectionZ("hMultGenAcc");
  TH1D* hMultGenLimAccAllPt=hPtVsYVsMultGenLimAcc->ProjectionZ("hMultGenLimAcc");
  TH1D* hEffVsMultAllPt=(TH1D*)hMultRecoAllPt->Clone("hEff");
  hEffVsMultAllPt->Divide(hMultRecoAllPt,hMultGenAccAllPt,1,1,"B");
  TH1D* hAccEffVsMultAllPt=(TH1D*)hMultRecoAllPt->Clone("hAccEff");
  hAccEffVsMultAllPt->Divide(hMultRecoAllPt,hMultGenLimAccAllPt,1,1,"B");
  TH1D* hAccVsMultAllPt=(TH1D*)hMultGenAccAllPt->Clone("hAcc");
  hAccVsMultAllPt->Divide(hMultGenAccAllPt,hMultGenLimAccAllPt,1,1,"B");
  hEffVsMultAllPt->SetStats(0);
  hAccVsMultAllPt->SetStats(0);
  hAccEffVsMultAllPt->SetStats(0);

  TCanvas* c2=new TCanvas("c2","EffVsMult",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLogy();
  hMultGenLimAccAllPt->SetLineColor(1);
  hMultGenLimAccAllPt->Draw();
  hMultGenLimAccAllPt->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  hMultGenLimAccAllPt->GetYaxis()->SetTitle("Entries");
  gPad->Update();
  TPaveStats *st11=(TPaveStats*)hMultGenLimAccAllPt->GetListOfFunctions()->FindObject("stats");
  st11->SetY1NDC(0.7);
  st11->SetY2NDC(0.89);
  hMultGenAccAllPt->SetLineColor(2);
  hMultGenAccAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st12=(TPaveStats*)hMultGenAccAllPt->GetListOfFunctions()->FindObject("stats");
  st12->SetY1NDC(0.5);
  st12->SetY2NDC(0.69);
  st12->SetTextColor(2);
  hMultRecoAllPt->SetLineColor(4);
  hMultRecoAllPt->Draw("sames");
  gPad->Update();
  TPaveStats *st13=(TPaveStats*)hMultRecoAllPt->GetListOfFunctions()->FindObject("stats");
  st13->SetY1NDC(0.3);
  st13->SetY2NDC(0.49);
  st13->SetTextColor(4);
  gPad->Modified();
  c2->cd(2);
  hEffVsMultAllPt->SetLineColor(4);
  hEffVsMultAllPt->SetMinimum(0);
  hEffVsMultAllPt->SetMaximum(1.6);
  hEffVsMultAllPt->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  hEffVsMultAllPt->GetYaxis()->SetTitle("Ratio");
  hEffVsMultAllPt->Draw();
  hAccVsMultAllPt->SetLineColor(2);
  hAccVsMultAllPt->Draw("same");
  //  hAccEffVsMultAllPt->SetLineColor(6);
  // hAccEffVsMultAllPt->Draw("same");
  TLatex* tacc2=new TLatex(0.16,0.8,"Acceptance (CombinHF, 13d3)");
  tacc2->SetNDC();
  tacc2->SetTextColor(hAccVsMultAllPt->GetLineColor());
  tacc2->Draw();
  TLatex* te2=new TLatex(0.16,0.72,"Efficiency");
  te2->SetNDC();
  te2->SetTextColor(hEffVsMultAllPt->GetLineColor());
  te2->Draw();
  c2->SaveAs(Form("figures/EfficVsMult_%s.eps",suffix.Data()));


  const Int_t nMultBins=6;
  Double_t mulLims[nMultBins+1]={0.,5.,12.,20.,40.,80.,200.};
  Int_t mulcol[nMultBins]={1,kRed+1,kGreen+2,4,kOrange+2,kMagenta};

  TH1D* hPtRecoM[nMultBins];
  TH1D* hPtGenAccM[nMultBins];
  //  TH1D* hPtGenLimAccM[nMultBins];
  TH1D* hEffVsPtM[nMultBins];
  for(Int_t j=0; j<nMultBins; j++){
    Int_t lowBin=hPtVsYVsMultReco->GetZaxis()->FindBin(mulLims[j]);
    Int_t hiBin=hPtVsYVsMultReco->GetZaxis()->FindBin(mulLims[j+1]-0.001);
    //    printf("%d (%f)  %d(%f)\n",lowBin,hPtVsYVsMultReco->GetZaxis()->GetBinLowEdge(lowBin),hiBin,hPtVsYVsMultReco->GetZaxis()->GetBinUpEdge(hiBin));
    
    hPtRecoM[j]=hPtVsYVsMultReco->ProjectionX(Form("hPtRecoM%d",j),0,-1,lowBin,hiBin);
    hPtGenAccM[j]=hPtVsYVsMultGenAcc->ProjectionX(Form("hPtGenAccM%d",j),0,-1,lowBin,hiBin);
    //    hPtGenLimAccM[j]=hPtVsYVsMultGenLimAcc->ProjectionX(Form("hPtGenLimAccM%d",j),0,-1,lowBin,hiBin);
    hEffVsPtM[j]=(TH1D*)hPtRecoM[j]->Clone(Form("hEffM%d",j));
    hEffVsPtM[j]->Divide(hPtRecoM[j],hPtGenAccM[j],1,1,"B");
  }

  hEffVsPtM[0]->SetStats(0);
  TCanvas* cwp=new TCanvas("cwp","Eff and Weight vs pt",1200,600);
  cwp->Divide(2,1);
  cwp->cd(1);
  gPad->SetLeftMargin(0.12);
  hEffVsPtM[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hEffVsPtM[0]->GetYaxis()->SetTitle("Efficiency");
  hEffVsPtM[0]->GetYaxis()->SetTitleOffset(1.4);
  hEffVsPtM[0]->SetLineColor(mulcol[0]);
  hEffVsPtM[0]->SetMinimum(0);
  hEffVsPtM[0]->SetMaximum(1.8);
  hEffVsPtM[0]->Draw();
  TLegend* legp=new TLegend(0.16,0.5,0.55,0.89);
  legp->SetFillStyle(0);
  legp->SetBorderSize(0);
  legp->AddEntry(hEffVsPtM[0],Form("%.0f<N_{tracklets}<%.0f",mulLims[0],mulLims[1]),"L")->SetTextColor(mulcol[0]);
  for(Int_t j=1; j<nMultBins; j++){
    hEffVsPtM[j]->SetLineColor(mulcol[j]);
    hEffVsPtM[j]->Draw("same");
    legp->AddEntry(hEffVsPtM[j],Form("%.0f<N_{tracklets}<%.0f",mulLims[j],mulLims[j+1]),"L")->SetTextColor(mulcol[j]);
  }
  legp->Draw();
  cwp->cd(2);
  gPad->SetLeftMargin(0.12);
  funcPtWeight->SetTitle("");
  funcPtWeight->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  funcPtWeight->GetYaxis()->SetTitle("Weight");
  funcPtWeight->Draw();
  cwp->SaveAs(Form("figures/EfficVsPtMultBins_%s.eps",suffix.Data()));


  TH1D* hMultReco[nPtBins];
  TH1D* hMultGenAcc[nPtBins];
  //  TH1D* hMultGenLimAcc[nPtBins];
  TH1D* hEffVsMult[nPtBins];
  for(Int_t j=0; j<nPtBins; j++){
    Int_t lowBin=hPtVsYVsMultReco->GetXaxis()->FindBin(binLims[j]);
    Int_t hiBin=hPtVsYVsMultReco->GetXaxis()->FindBin(binLims[j+1]-0.001);
    //printf("%d (%f)  %d(%f)\n",lowBin,hPtVsYVsMultReco->GetXaxis()->GetBinLowEdge(lowBin),hiBin,hPtVsYVsMultReco->GetXaxis()->GetBinUpEdge(hiBin));
    
    hMultReco[j]=hPtVsYVsMultReco->ProjectionZ(Form("hMultReco%d",j),lowBin,hiBin);
    hMultGenAcc[j]=hPtVsYVsMultGenAcc->ProjectionZ(Form("hMultGenAcc%d",j),lowBin,hiBin);
    //    hMultGenLimAcc[j]=hPtVsYVsMultGenLimAcc->ProjectionZ(Form("hMultGenLimAcc%d",j),lowBin,hiBin);
    hEffVsMult[j]=(TH1D*)hMultReco[j]->Clone("hEff");
    hEffVsMult[j]->Divide(hMultReco[j],hMultGenAcc[j],1,1,"B");
  }

  // TCanvas* c2e=new TCanvas("c2e");
  // hEffVsMult[0]->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  // hEffVsMult[0]->GetYaxis()->SetTitle("Efficiency");
  // hEffVsMult[0]->GetYaxis()->SetTitleOffset(1.4);
  // hEffVsMult[0]->SetLineColor(ptcol[0]);
  // hEffVsMult[0]->SetMinimum(0);
  // hEffVsMult[0]->SetMaximum(1.6);
  // hEffVsMult[0]->Draw();
  // TLegend* legm=new TLegend(0.16,0.5,0.4,0.89);
  // legm->SetFillStyle(0);
  // legm->SetBorderSize(0);
  // legm->AddEntry(hEffVsMult[0],Form("%.0f<p_{T}<%.0f GeV/c",binLims[0],binLims[1]),"L")->SetTextColor(ptcol[0]);
  // for(Int_t j=1; j<nPtBins; j++){
  //   hEffVsMult[j]->SetLineColor(ptcol[j]);
  //   hEffVsMult[j]->Draw("same");
  //   legm->AddEntry(hEffVsMult[j],Form("%.0f<p_{T}<%.0f GeV/c",binLims[j],binLims[j+1]),"L")->SetTextColor(ptcol[j]);
  // }
  // legm->Draw();


  TFile* filw = new TFile(fileWeightName.Data());
  TH1F** hWeight=new TH1F*[3];
  for(Int_t iw=0;iw<3; iw++){
    hWeight[iw]=(TH1F*)filw->Get(histoWeightName[iw].Data());
    hWeight[iw]->SetLineColor(wcol[iw]);
    hWeight[iw]->SetStats(0);
    hWeight[iw]->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
    hWeight[iw]->GetYaxis()->SetTitle("Weight");
    hWeight[iw]->GetYaxis()->SetTitleOffset(1.4);
    hWeight[iw]->SetMaximum(3.);
    hWeight[iw]->SetMinimum(0.);
  }

  hEffVsMult[0]->SetStats(0);
  TCanvas* cw=new TCanvas("cw","Eff and Weight vs mult",1200,600);
  cw->Divide(2,1);
  cw->cd(1);
  gPad->SetLeftMargin(0.12);
  hEffVsMult[0]->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  hEffVsMult[0]->GetYaxis()->SetTitle("Efficiency");
  hEffVsMult[0]->GetYaxis()->SetTitleOffset(1.4);
  hEffVsMult[0]->SetLineColor(ptcol[0]);
  hEffVsMult[0]->SetMinimum(0);
  hEffVsMult[0]->SetMaximum(1.6);
  hEffVsMult[0]->Draw();
  TLegend* legm=new TLegend(0.16,0.5,0.4,0.89);
  legm->SetFillStyle(0);
  legm->SetBorderSize(0);
  legm->AddEntry(hEffVsMult[0],Form("%.0f<p_{T}<%.0f GeV/c",binLims[0],binLims[1]),"L")->SetTextColor(ptcol[0]);
  for(Int_t j=1; j<nPtBins; j++){
    hEffVsMult[j]->SetLineColor(ptcol[j]);
    hEffVsMult[j]->Draw("same");
    legm->AddEntry(hEffVsMult[j],Form("%.0f<p_{T}<%.0f GeV/c",binLims[j],binLims[j+1]),"L")->SetTextColor(ptcol[j]);
  }
  legm->Draw();
  cw->cd(2);
  gPad->SetLeftMargin(0.12);
  hWeight[0]->Draw();
  hWeight[1]->Draw("same");
  hWeight[2]->Draw("same");
  TLegend* legw=new TLegend(0.55,0.65,0.89,0.89);
  legw->SetFillStyle(0);
  legw->AddEntry(hWeight[0],"Candidate Weight","L");
  legw->AddEntry(hWeight[1],"D Weight","L");
  legw->AddEntry(hWeight[2],"EvSel Weight","L");
  legw->Draw();
  cw->SaveAs(Form("figures/EfficVsMultPtBins_%s.eps",suffix.Data()));


  TH1D *hpteffNoWeight =new TH1D("hEffVsPtNoWeight","",nPtBins,binLims);
  TH1D **hpteffMultWeight =new TH1D*[3];
  hpteffMultWeight[0]=new TH1D("hEffVsPtCandWeight","",nPtBins,binLims);
  hpteffMultWeight[1]=new TH1D("hEffVsPtDWeight","",nPtBins,binLims);
  hpteffMultWeight[2]=new TH1D("hEffVsPtEvSelWeight","",nPtBins,binLims);
  TH1D **hpteffMultPtWeight =new TH1D*[3];
  hpteffMultPtWeight[0]=new TH1D("hEffVsPtCandFONLLWeight","",nPtBins,binLims);
  hpteffMultPtWeight[1]=new TH1D("hEffVsPtDFONLLWeight","",nPtBins,binLims);
  hpteffMultPtWeight[2]=new TH1D("hEffVsPtEvSelFONLLWeight","",nPtBins,binLims);

  TH1D* hMultRecoAllPtW=hPtVsYVsMultReco->ProjectionZ("hMultRecoW");
  TH1D* hMultGenAccAllPtW=hPtVsYVsMultGenAcc->ProjectionZ("hMultGenAccW");
  TH1D* hMultGenLimAccAllPtW=hPtVsYVsMultGenLimAcc->ProjectionZ("hMultGenLimAccW");
  hMultRecoAllPtW->Reset("ines");
  hMultGenAccAllPtW->Reset("ines");
  hMultGenLimAccAllPtW->Reset("ines");

  Double_t countNumer[nPtBins];
  Double_t countDenom[nPtBins];
  Double_t countNumerMultWei[nPtBins][3];
  Double_t countDenomMultWei[nPtBins][3];
  Double_t countNumerMultPtWei[nPtBins][3];
  Double_t countDenomMultPtWei[nPtBins][3];

  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    countNumer[iPtBin]=0.0;
    countDenom[iPtBin]=0.0;
    for(Int_t iw=0; iw<3; iw++){
      countNumerMultWei[iPtBin][iw]=0.0;
      countDenomMultWei[iPtBin][iw]=0.0;
      countNumerMultPtWei[iPtBin][iw]=0.0;
      countDenomMultPtWei[iPtBin][iw]=0.0;
    }
  }

  for(Int_t ibx=0; ibx<=hPtVsYVsMultGenAcc->GetNbinsX()+1; ibx++){
    for(Int_t iby=0; iby<=hPtVsYVsMultGenAcc->GetNbinsY()+1; iby++){
      for(Int_t ibz=0; ibz<=hPtVsYVsMultGenAcc->GetNbinsZ()+1; ibz++){
	Double_t pt=hPtVsYVsMultGenAcc->GetXaxis()->GetBinCenter(ibx);
	//	Double_t y=hPtVsYVsMultGenAcc->GetYaxis()->GetBinCenter(iby);
	Double_t mult=hPtVsYVsMultGenAcc->GetZaxis()->GetBinCenter(ibz);
	//	  printf("%f %f %f\n",pt,y,mult);
	Int_t jPtBin=TMath::BinarySearch(nPtBins+1,binLims,pt);
	if(jPtBin>=0 && jPtBin<nPtBins){
	  Double_t wpt=funcPtWeight->Eval(pt);
	  countNumer[jPtBin]+=hPtVsYVsMultReco->GetBinContent(ibx,iby,ibz);
	  countDenom[jPtBin]+=hPtVsYVsMultGenAcc->GetBinContent(ibx,iby,ibz);
	  for(Int_t iw=0; iw<3; iw++){
	    Int_t binw=hWeight[iw]->FindBin(mult);
	    Double_t w=0;
	    if(binw>=1 && binw<hWeight[iw]->GetNbinsX()+1){
	      w=hWeight[iw]->GetBinContent(binw);
	      //		printf("mult %.0f   bin %d   wei %f\n",mult,binw,w);
	    }else{
	      if(hPtVsYVsMultGenAcc->GetBinContent(ibx,iby,ibz)>0){
		printf("mult %.0f   bin %d   wei %f\n",mult,binw,w);
		getchar();
	      }
	    }
	    countNumerMultWei[jPtBin][iw]+=hPtVsYVsMultReco->GetBinContent(ibx,iby,ibz)*w;
	    countDenomMultWei[jPtBin][iw]+=hPtVsYVsMultGenAcc->GetBinContent(ibx,iby,ibz)*w;
	    countNumerMultPtWei[jPtBin][iw]+=hPtVsYVsMultReco->GetBinContent(ibx,iby,ibz)*w*wpt;
	    countDenomMultPtWei[jPtBin][iw]+=hPtVsYVsMultGenAcc->GetBinContent(ibx,iby,ibz)*w*wpt;
	    hMultRecoAllPtW->Fill(mult,hPtVsYVsMultReco->GetBinContent(ibx,iby,ibz)*w);
	    hMultGenAccAllPtW->Fill(mult,hPtVsYVsMultGenAcc->GetBinContent(ibx,iby,ibz)*w);
	    hMultGenLimAccAllPtW->Fill(mult,hPtVsYVsMultGenLimAcc->GetBinContent(ibx,iby,ibz)*w);
	  }
	}
      }
    }
  }

  for(Int_t iPtBin=0; iPtBin<nPtBins; iPtBin++){
    Double_t eff1=countNumer[iPtBin]/countDenom[iPtBin];
    Double_t erreff1=TMath::Sqrt(eff1*(1-eff1)/countDenom[iPtBin]);
    printf("---- Pt range %.0f - %.0f ----\n",binLims[iPtBin],binLims[iPtBin+1]);
    printf("Eff from Projection = %f/%f = %f+-%f\n",hPtRecoR->GetBinContent(iPtBin+1),hPtGenAccR->GetBinContent(iPtBin+1),hEffVsPtR->GetBinContent(iPtBin+1),hEffVsPtR->GetBinError(iPtBin+1));
    printf("Eff No weights      = %f/%f = %f+-%f\n",countNumer[iPtBin],countDenom[iPtBin],eff1,erreff1);
    hpteffNoWeight->SetBinContent(iPtBin+1,eff1);
    hpteffNoWeight->SetBinError(iPtBin+1,erreff1);
    for(Int_t iw=0;iw<3; iw++){
      Double_t eff2=countNumerMultWei[iPtBin][iw]/countDenomMultWei[iPtBin][iw];
      Double_t erreff2=TMath::Sqrt(eff2*(1-eff2)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
      printf("Eff With mult weights %d   = %f/%f = %f+-%f\n",iw,countNumerMultWei[iPtBin][iw],countDenomMultWei[iPtBin][iw],eff2,erreff2);
      hpteffMultWeight[iw]->SetBinContent(iPtBin+1,eff2);
      hpteffMultWeight[iw]->SetBinError(iPtBin+1,erreff2);
      Double_t eff3=countNumerMultPtWei[iPtBin][iw]/countDenomMultPtWei[iPtBin][iw];
      Double_t erreff3=TMath::Sqrt(eff3*(1-eff3)/countDenom[iPtBin]);// countDenom is NOT a typo, it has to be like this to get proper statistical errors from the no-weight case
      printf("Eff With mult+pt weights %d   = %f/%f = %f+-%f\n",iw,countNumerMultPtWei[iPtBin][iw],countDenomMultPtWei[iPtBin][iw],eff3,erreff3);
      hpteffMultPtWeight[iw]->SetBinContent(iPtBin+1,eff3);
      hpteffMultPtWeight[iw]->SetBinError(iPtBin+1,erreff3);
    }
  }

  TH1D* hEffVsMultAllPtW=(TH1D*)hMultRecoAllPtW->Clone("hEff");
  hEffVsMultAllPtW->Divide(hMultRecoAllPtW,hMultGenAccAllPtW,1,1,"B");
  TH1D* hAccVsMultAllPtW=(TH1D*)hMultGenAccAllPtW->Clone("hAcc");
  hAccVsMultAllPtW->Divide(hMultGenAccAllPtW,hMultGenLimAccAllPtW,1,1,"B");
  hEffVsMultAllPtW->SetStats(0);
  hAccVsMultAllPtW->SetStats(0);

  TCanvas* c2w=new TCanvas("c2w","EffVsMultW",1200,600);
  c2w->Divide(2,1);
  c2w->cd(1);
  gPad->SetLogy();
  hMultGenLimAccAllPtW->SetLineColor(1);
  hMultGenLimAccAllPtW->Draw();
  hMultGenLimAccAllPtW->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  hMultGenLimAccAllPtW->GetYaxis()->SetTitle("Entries");
  hMultGenAccAllPtW->SetLineColor(2);
  hMultGenAccAllPtW->Draw("same");
  hMultRecoAllPtW->SetLineColor(4);
  hMultRecoAllPtW->Draw("same");
  c2w->cd(2);
  hEffVsMultAllPtW->SetLineColor(4);
  hEffVsMultAllPtW->SetMinimum(0);
  hEffVsMultAllPtW->SetMaximum(1.6);
  hEffVsMultAllPtW->GetXaxis()->SetTitle("N_{tracklets} in |#eta|<1");
  hEffVsMultAllPtW->GetYaxis()->SetTitle("Ratio");
  hEffVsMultAllPtW->Draw();
  hAccVsMultAllPtW->SetLineColor(2);
  hAccVsMultAllPtW->Draw("same");
  tacc2->Draw();
  te2->Draw();




  hEffVsPtR->SetMarkerStyle(0);
  hEffVsPtR->SetMarkerColor(0);
  hEffVsPtR->SetMarkerSize(1.2);
  hEffVsPtR->SetLineColor(kGray);
  hEffVsPtR->SetLineWidth(4);
  hEffVsPtR->SetStats(0);
  hpteffNoWeight->SetMarkerStyle(20);

  TH1F** hRatio=new TH1F*[3];
   for(Int_t iw=0;iw<3; iw++){
    hRatio[iw]=(TH1F*)hpteffMultWeight[iw]->Clone("hRatio");
    hRatio[iw]->Divide(hpteffMultWeight[iw],hpteffNoWeight);
    hRatio[iw]->GetYaxis()->SetTitle("With Mult weight / Without weight");
    hRatio[iw]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatio[iw]->GetYaxis()->SetTitleOffset(1.4);
    hRatio[iw]->SetLineColor(wcol[iw]);
    hRatio[iw]->SetMarkerColor(wcol[iw]);
    hRatio[iw]->SetMarkerStyle(wmark[iw]);
    hRatio[iw]->SetStats(0);
    hpteffMultWeight[iw]->SetLineColor(wcol[iw]);
    hpteffMultWeight[iw]->SetMarkerColor(wcol[iw]);
    hpteffMultWeight[iw]->SetMarkerStyle(wmark[iw]);
  }
 TH1F** hRatioPtW=new TH1F*[3];
   for(Int_t iw=0;iw<3; iw++){
    hRatioPtW[iw]=(TH1F*)hpteffMultPtWeight[iw]->Clone("hRatio");
    hRatioPtW[iw]->Divide(hpteffMultPtWeight[iw],hpteffMultWeight[iw]);
    hRatioPtW[iw]->GetYaxis()->SetTitle("With Mult+p_{T} weight / With only Mult weight");
    hRatioPtW[iw]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hRatioPtW[iw]->GetYaxis()->SetTitleOffset(1.4);
    hRatioPtW[iw]->SetLineColor(wcol[iw]);
    hRatioPtW[iw]->SetMarkerColor(wcol[iw]);
    hRatioPtW[iw]->SetMarkerStyle(wmark[iw]);
    hRatioPtW[iw]->SetStats(0);
    hpteffMultPtWeight[iw]->SetLineColor(wcol[iw]);
    hpteffMultPtWeight[iw]->SetMarkerColor(wcol[iw]);
    hpteffMultPtWeight[iw]->SetMarkerStyle(wmark[iw]);
  }

  TCanvas* ceff=new TCanvas("ceff","Eff Mult Wei",1200,600);
  ceff->Divide(2,1);
  ceff->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetTickx();
  gPad->SetTicky();
  hEffVsPtR->SetMinimum(0);
  hEffVsPtR->SetMaximum(1);
  hEffVsPtR->Draw();
  hEffVsPtR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hEffVsPtR->GetYaxis()->SetTitle("Efficiency");
  hEffVsPtR->GetYaxis()->SetTitleOffset(1.4);
  hpteffNoWeight->Draw("same");
  for(Int_t iw=0;iw<3; iw++) hpteffMultWeight[iw]->Draw("same");
  TLegend* leg=new TLegend(0.55,0.15,0.89,0.45);
  leg->SetFillStyle(0);
  leg->AddEntry(hEffVsPtR,"TH3F::Project","L");
  leg->AddEntry(hpteffNoWeight,"Multiplcity slices - No Weight","PL");
  leg->AddEntry(hpteffMultWeight[0],"Multiplcity slices - Cand Weight","PL");
  leg->AddEntry(hpteffMultWeight[1],"Multiplcity slices - D Weight","PL");
  leg->AddEntry(hpteffMultWeight[2],"Multiplcity slices - EvSel Weight","PL");
  leg->Draw();
  ceff->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetTickx();
  gPad->SetTicky();
  hRatio[0]->Draw();
  hRatio[0]->SetMinimum(0.85);
  hRatio[0]->SetMaximum(1.15);
  for(Int_t iw=1;iw<3; iw++) hRatio[iw]->Draw("same");
  ceff->SaveAs(Form("figures/EfficWithMultWeights_%s.eps",suffix.Data()));


  TCanvas* ceff2=new TCanvas("ceff2","Eff Mult+Pt Wei",1200,600);
  ceff2->Divide(2,1);
  ceff2->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetTickx();
  gPad->SetTicky();
  hEffVsPtR->SetMinimum(0);
  hEffVsPtR->SetMaximum(1);
  hEffVsPtR->Draw();
  hEffVsPtR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hEffVsPtR->GetYaxis()->SetTitle("Efficiency");
  hEffVsPtR->GetYaxis()->SetTitleOffset(1.4);
  hpteffNoWeight->Draw("same");
  for(Int_t iw=0;iw<3; iw++) hpteffMultPtWeight[iw]->Draw("same");
  TLegend* leg2=new TLegend(0.55,0.15,0.89,0.45);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hEffVsPtR,"TH3F::Project","L");
  leg2->AddEntry(hpteffNoWeight,"Multiplcity slices - No Weight","PL");
  leg2->AddEntry(hpteffMultPtWeight[0],"Multiplcity slices - FONLL+Cand Weight","PL");
  leg2->AddEntry(hpteffMultPtWeight[1],"Multiplcity slices - FONLL+D Weight","PL");
  leg2->AddEntry(hpteffMultPtWeight[2],"Multiplcity slices - FONLL+EvSel Weight","PL");
  leg2->Draw();
  ceff2->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetTickx();
  gPad->SetTicky();
  hRatioPtW[0]->Draw();
  hRatioPtW[0]->SetMinimum(0.95);
  hRatioPtW[0]->SetMaximum(1.05);
  for(Int_t iw=1;iw<3; iw++) hRatioPtW[iw]->Draw("same");
  ceff2->SaveAs(Form("figures/EfficWithMultAndPtWeights_%s.eps",suffix.Data()));

  TFile* out=new TFile(Form("outputEff_%s.root",suffix.Data()),"recreate");
  hEffVsPtR->Write();
  hAccToy->Write();
  hpteffMultPtWeight[0]->Write();
  hpteffMultPtWeight[1]->Write();
  hpteffMultPtWeight[2]->Write();
  out->Close();

}
