#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THn.h"
#include "TList.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

//loadlibs();
//gSystem->Load("libANALYSIS");
//gSystem->Load("libANALYSISalice");
//#include "AliESDtrackCuts.h"

void MakeSensitivityPlots2();
//
TCanvas * GetSensitivityPlot(TString cutname = "Ncl",  
			     Int_t projectionAxis = 1,
			     Int_t particleType = 5,
			     Double_t lowCut = 80.,
			     Double_t highCut = 1e4, 
			     TString inFileNameData = "data_LHC10d.root",
			     TString inFileNameMC = "mc_LHC10d.root");
//
TH1D   * GetAcceptedFraction(TString cutname = "Ncl",  
			     Int_t projectionAxis = 1,
			     Int_t particleType = 5,
			     Double_t lowCut = 80.,
			     Double_t highCut = 1e4, 
			     TString inFileName = "data_LHC10d.root");
//
void PrintCutGallery(TString inFileNameData = "data_LHC10d.root",
		     TString inFileNameMC   = "mc_LHC10d.root",
		     TString outFileName    = "plots/10d.pdf");

//
TH1D *GetITSTPCMatchingEff(TString inFileNameData = "data_LHC10d.root",
                              Int_t particleType = 5,
                              Int_t projectionAxis = 1);

//
TH1D *GetITSTPCMatchingHisto(TString inFileNameData = "data_LHC10d.root",
                             Int_t particleType = 5,
                             Bool_t isMatched = kTRUE,
                             Int_t projectionAxis = 1);


//______________________________________________________________________________
void MakeSensitivityPlots2() {
  //
  // make all the senstivity plots
  //
  //PrintCutGallery("output/LHC10b_data.root", "output/LHC10b_MC.root", "plots/10b.pdf");
  PrintCutGallery("data_LHC10d.root", "mc_LHC10d.root", "plots/10d.pdf");
  //PrintCutGallery("output/LHC10e_data.root", "output/LHC10e_MC.root", "plots/10e.pdf");
  //PrintCutGallery("output/LHC10h_data.root", "output/LHC10h_MC.root", "plots/10h.pdf");

}

//______________________________________________________________________________
void PrintCutGallery(TString inFileNameData, TString inFileNameMC, TString outFileName) {
  //
  // print a gallery of cuts
  //
  TCanvas * canvMaster = new TCanvas("canvMaster","canvMaster");
  canvMaster->Divide(3,2);
  //TCanvas * canv = GetSensitivityPlot();
  //
  TCanvas * canvSens[6]; 
  for(Int_t i=0; i < 6; i++) {
    canvSens[i]= GetSensitivityPlot("Ncl", 1, 5, 50 + i*10, 1e4, inFileNameData.Data(),  inFileNameMC.Data());
    //canvSens[i]= GetSensitivityPlot("Chi2Tpc", 1, 5, -1e4, 2.5 + i*0.5, inFileNameData.Data(),  inFileNameMC.Data());
    canvSens[i]->SetName(Form("canvSens_%i",i));
    canvMaster->cd(i+1);
    canvSens[i]->DrawClonePad();
      
  }
  canvMaster->Print(outFileName.Data());

}


//______________________________________________________________________________
TCanvas * GetSensitivityPlot(TString cutname,  
			     Int_t projectionAxis,
			     Int_t particleType,
			     Double_t lowCut,
			     Double_t highCut, 
			     TString inFileNameData, 
			     TString inFileNameMC) {
  //
  // make a single plot
  //
  TH1D * nclAcceptedData = GetAcceptedFraction(cutname.Data(), projectionAxis, particleType, lowCut, highCut, inFileNameData);
  nclAcceptedData->SetNameTitle(Form("%s_DATA",nclAcceptedData->GetName()), Form("%s_DATA",nclAcceptedData->GetName()));
  //
  //
  nclAcceptedData->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  nclAcceptedData->GetYaxis()->SetTitle("accepted fraction");
  nclAcceptedData->GetYaxis()->SetLabelSize(0.04);
  nclAcceptedData->GetYaxis()->SetTitleSize(0.06);
  nclAcceptedData->GetYaxis()->SetTitleOffset(1);
  nclAcceptedData->SetMaximum(1.07);
  nclAcceptedData->SetMinimum(0.625);
  nclAcceptedData->SetLineColor(kRed -3);
  //
  TH1D * nclAcceptedMc   = GetAcceptedFraction(cutname.Data(), projectionAxis, particleType, lowCut, highCut, inFileNameMC);
  nclAcceptedMc->SetNameTitle(Form("%s_MC",nclAcceptedMc->GetName()), Form("%s_MC",nclAcceptedMc->GetName()));
  nclAcceptedMc->SetLineColor(kBlue -3);
  //
  // create the complicated splitted canvas
  //
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas * canvCut = new TCanvas("canvCut",Form("sensitivity to %s cut",cutname.Data()),600,600);
  canvCut->Divide(1,2);
  //
  TPad* canvas_up = (TPad*) canvCut->GetListOfPrimitives()->FindObject(Form("canvCut_1"));
  TPad* canvas_dw = (TPad*) canvCut->GetListOfPrimitives()->FindObject(Form("canvCut_2"));
  //
  canvas_up->SetLogx();
  canvas_dw->SetLogx();
  //
  // define the size
  double up_height     = 0.75; // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.30; // please tune so that the smaller canvas size will work in your environment
  //
  //double font_size_dw  = 0.1; // please tune the font size parameter for bottom figure
  double dw_height    = (1. - up_height) * dw_correction;
  // set pad size
  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
  canvas_dw->SetPad(0., 0., 1., dw_height);
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  //
  // set top margin 0 for bottom figure
  canvas_dw->SetTopMargin(0);
  canvas_dw->SetLeftMargin(0.12);
  canvas_dw->SetBottomMargin(0.3);
  //
  //set margins for top figure
  canvas_up->SetLeftMargin(0.12);
  //
  canvas_up->cd();
  nclAcceptedData->GetYaxis()->SetLabelFont(62);
  nclAcceptedData->SetObjectStat(0);
  nclAcceptedMc->SetObjectStat(0);
  //
  nclAcceptedData->DrawCopy();
  nclAcceptedMc->DrawCopy("SAME");
  //
  TLegend *leg = new TLegend(0.65,0.2,.85,0.4);
  leg->AddEntry(nclAcceptedData,"Data","f");
  leg->AddEntry(nclAcceptedMc,"MC","f");
  leg->SetBorderSize(0);
  gStyle->SetFillColor(0);
  leg->Draw();
  Float_t lowCutTitle = lowCut < -999. ? 0 : lowCut;
  TLatex *   tex = new TLatex(0.9654768,1.032485,Form("%4.2f < %s < %4.2f", lowCutTitle, cutname.Data(), highCut));
  tex->Draw();
  //
  // Draw the ratio
  //    
  TH1D * nclAcceptedMcDataRatio = (TH1D*)nclAcceptedData->Clone();
  //
  canvCut->cd(2)->SetLogx();
  // gPad->SetTicky(2);
  gPad->SetFillStyle(0);
  //
  nclAcceptedMcDataRatio->Divide(nclAcceptedMc);
  nclAcceptedMcDataRatio->GetYaxis()->SetTitle("ratio");
  nclAcceptedMcDataRatio->GetXaxis()->SetLabelSize(0.1);
  nclAcceptedMcDataRatio->GetXaxis()->SetTitleSize(0.11);
  nclAcceptedMcDataRatio->GetYaxis()->SetLabelSize(0.07);
  nclAcceptedMcDataRatio->GetYaxis()->SetTitleSize(0.11);
  nclAcceptedMcDataRatio->GetXaxis()->SetTitleOffset(1.23);
  nclAcceptedMcDataRatio->GetYaxis()->SetTitleOffset(0.5);
  nclAcceptedMcDataRatio->GetYaxis()->CenterTitle();
  nclAcceptedMcDataRatio->SetMaximum(1.18);
  nclAcceptedMcDataRatio->SetMinimum(0.8);
  nclAcceptedMcDataRatio->GetYaxis()->SetLabelFont(62);
  nclAcceptedMcDataRatio->SetLineColor(kGreen+2);
  nclAcceptedMcDataRatio->Draw();
  canvas_dw->SetFillStyle(1001);
  //
  //
  //
  return canvCut;
  
}


//______________________________________________________________________________
TH1D * GetAcceptedFraction(TString cutname,  
			   Int_t projectionAxis,
			   Int_t particleType,
			   Double_t lowCut,
			   Double_t highCut, 
			   TString inFileName) {
  //
  // accepted fraction of tracks for ncl cut vs. pT
  //
  TFile * inFileData = TFile::Open(inFileName);
  TList * l = (TList * ) inFileData->Get("akalweit_TrackingUncert");
  THnF * histNcl = (THnF *) l->FindObject(Form("hist%s",cutname.Data()));
  //
  // select particleType
  //
  histNcl->GetAxis(4)->SetRangeUser(particleType, particleType);
  //
  // restrict eta-range
  //
  histNcl->GetAxis(2)->SetRangeUser(-0.75, 0.75);
  //
  // determine sensitivities
  //
  TH1D * hAll = histNcl->Projection(projectionAxis);
  hAll->SetNameTitle("hAll","hAll");  
  //
  const Int_t kVeryBig = 10000;
  histNcl->GetAxis(0)->SetRangeUser(lowCut, highCut);
  TH1D * hAccepted = histNcl->Projection(1);
  hAccepted->SetNameTitle(Form("hAccepted_%s_%f_%f",cutname.Data(),lowCut,highCut),Form("hAccepted_%s_%f_%f",cutname.Data(),lowCut,highCut));
  //
  hAccepted->Divide(hAll);
  //
  // some cosmetics
  //
  hAccepted->SetLineWidth(3);
  hAccepted->SetDirectory(0);
  //
  delete l;
  inFileData->Close();
  //
  return hAccepted;

  
}

//______________________________________________________________________________
TH1D *GetITSTPCMatchingEff(TString inFileNameData,
                              Int_t particleType,
                              Int_t projectionAxis){
    
    //
    // make a single plot
    //
    TH1D * hMatching = GetITSTPCMatchingHisto(inFileNameData, particleType, kTRUE, projectionAxis);
    hMatching->SetNameTitle(Form("PID: %d\t proj: %d",particleType,projectionAxis), Form("PID: %d\t proj: %d",particleType,projectionAxis));
    hMatching->Sumw2();
    //
    TH1D *hNoMatching = GetITSTPCMatchingHisto(inFileNameData, particleType, kFALSE, projectionAxis);
    hNoMatching->SetNameTitle(Form("PID: %d\t proj: %d",particleType,projectionAxis), Form("PID: %d\t proj: %d",particleType,projectionAxis));
    hNoMatching->Sumw2();
    //
    hMatching->Divide(hMatching,hNoMatching,1,1,"B");
    
    //TCanvas * canvEff = new TCanvas("canvEff","matching efficiency",600,600);
    //hMatching->Draw();
    
    return hMatching;
    
    
}



//______________________________________________________________________________
TH1D *GetITSTPCMatchingHisto(TString inFileNameData,
                             Int_t particleType,
                             Bool_t isMatched,
                             Int_t projectionAxis){
    
  //
  // ITS-TPC matching histograms as a funct. of pT, eta, phi, for each species
  //    
  TFile * inFileData = TFile::Open(inFileNameData);
  TList * l = (TList * ) inFileData->Get("akalweit_TrackingUncert");
  THnF * histITSTPC = (THnF *) l->FindObject("histTpcItsMatch");

  //
  // select particleType
  //
  histITSTPC->GetAxis(4)->SetRangeUser(particleType, particleType);

  //
  // extract isMatched = kFALSE or kTRUE
  //
  if(isMatched)
    histITSTPC->GetAxis(0)->SetRangeUser(1,1);
  else
    histITSTPC->GetAxis(0)->SetRangeUser(0,0);
  
  TH1D * hProj  = histITSTPC->Projection(projectionAxis);
  hProj->SetDirectory(0);
  hProj->SetNameTitle(Form("h%d",projectionAxis),Form("h%d",projectionAxis));
    
  //
  delete l;
  inFileData->Close();
  //
    
    
  return hProj;
  
    
}
