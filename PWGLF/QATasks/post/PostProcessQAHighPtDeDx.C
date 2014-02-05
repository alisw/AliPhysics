/*

  Use AliRoot because of AliXRDPROOFtoolkit:



  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/Base")
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGLF/SPECTRA/IdentifiedHighPt/macros")
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGLF/SPECTRA/IdentifiedHighPt/grid")
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGLF/SPECTRA/IdentifiedHighPt/lib")
  gROOT->SetMacroPath(".:$ALICE_ROOT/PWGLF/SPECTRA/IdentifiedHighPt/macros:$ALICE_ROOT/PWGLF/SPECTRA/IdentifiedHighPt/grid:$ALICE_ROOT/PWGLF/SPECTRA/IdentifiedHighPt/lib/")
  .L my_functions.C+
  .L my_tools.C+
  .L PostProcessQAHighPtDeDx.C+

   PlotQA("FileRoot/AnalysisResults.root")
   MakeFitsExternalData("FileRoot/AnalysisResults.root", "HistosForBB")

  MakeFitsV0s("FileRoot/AnalysisResults.root", "HistosForBB/PrimaryElectrons.root", "HistosForBB",0)
  MakeFitsV0s("FileRoot/AnalysisResults.root", "HistosForBB/PrimaryElectrons.root", "HistosForBB",1)
  MakeFitsV0s("FileRoot/AnalysisResults.root", "HistosForBB/PrimaryElectrons.root", "HistosForBB",2)
  MakeFitsV0s("FileRoot/AnalysisResults.root", "HistosForBB/PrimaryElectrons.root", "HistosForBB",3)


  PlotParametrizations("HistosForBB")




  FitDeDxVsP("FileRoot/AnalysisResults.root", 3.0, 10.0, 0, 6, 13, kTRUE,  0,  2, 0,1, 0, "HistosForBB/hres_0_5_02.root",27);
  FitDeDxVsP("FileRoot/AnalysisResults.root", 3.0, 10.0, 0, 6, 13, kTRUE,  2,  4, 0,1, 0, "HistosForBB/hres_0_5_24.root",27);
  FitDeDxVsP("FileRoot/AnalysisResults.root", 3.0, 10.0, 0, 6, 13, kTRUE,  4,  6, 0,1, 0, "HistosForBB/hres_0_5_46.root",27);
  FitDeDxVsP("FileRoot/AnalysisResults.root", 3.0, 10.0, 0, 6, 13, kTRUE,  6,  8, 0,1, 0, "HistosForBB/hres_0_5_68.root",27);


  MakeNSigmaPlot("FileRoot/AnalysisResults.root","fitparameters/MB/02_dataPbPb.root",2,50, 0,  "02_dataPbPb.root");
  MakeNSigmaPlot("FileRoot/AnalysisResults.root","fitparameters/MB/24_dataPbPb.root",2,50, 1,  "24_dataPbPb.root");
  MakeNSigmaPlot("FileRoot/AnalysisResults.root","fitparameters/MB/46_dataPbPb.root",2,50, 2,  "46_dataPbPb.root");
  MakeNSigmaPlot("FileRoot/AnalysisResults.root","fitparameters/MB/68_dataPbPb.root",2,50, 3,  "68_dataPbPb.root");


  PlotNSigma("nsigma_results")



  //add particle fractions vs p


  ExtractUncPartFractvsP("FileRoot/AnalysisResults.root", 2, 50,0,0, "results/eta02","fractions");
  ExtractUncPartFractvsP("FileRoot/AnalysisResults.root", 2, 50,0,1, "results/eta24","fractions");
  ExtractUncPartFractvsP("FileRoot/AnalysisResults.root", 2, 50,0,2, "results/eta46","fractions");
  ExtractUncPartFractvsP("FileRoot/AnalysisResults.root", 2, 50,0,3, "results/eta68","fractions");


*/
#include <TLegend.h>
#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include "my_tools.C"
#include "my_functions.C"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;





//const Char_t* ending[4] = {"02", "24", "46", "68"};

const Char_t * legEtaCut[4]=
  {"|#eta_{lab}|<0.2",  "0.2<|#eta_{lab}|<0.4",   "0.4<|#eta_{lab}|<0.6",   "0.6<|#eta_{lab}|<0.8"};


const Char_t* endingCent[1] = 
  {"MB"};
const Char_t* CentLatex[1] =
  {"MB"};

  const Char_t *Pid[4]={"Pion", "Kaon", "Proton", "Electron"};
//const Char_t *Pid[3]={"Pion","Kaon","Proton"};
const Char_t *PidLatex[3]={"#pi^{+} + #pi^{-}","K^{+} + K^{-}","p + #bar{p}"};


const Color_t PidColor[3]={2,kGreen,4};

TF1* piFunc = 0;
TF1* kFunc  = 0;
TF1* pFunc = 0;
TF1* eFunc = 0;
TF1* sigmaFunc = 0;

const Int_t nPtBins = 68;

Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
			     0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			     1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			     2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			     4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			     11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			     26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };




const Int_t nPtBinsV0s = 25;
Double_t ptBinsV0s[nPtBinsV0s+1] = { 0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
				     1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
				     9.0 , 12.0, 15.0, 20.0 };


TH2D * AddTwoSameBinningTH2D(TH2D *hPos, TH2D *hNeg, const Char_t *nameHist);


void PlotNSigma(const Char_t *path){

  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);



  TFile *fin[2]={0,0};
  
  fin[0] = FindFileFresh(Form("%s/02_dataPbPb.root",path));
  if(!fin[0])
    return;

  fin[1] = FindFileFresh(Form("%s/68_dataPbPb.root",path));
  if(!fin[1])
    return;

  TH1D *hNPiP[2]={0,0}; 
  TH1D *hNPiK[2]={0,0}; 
  TH1D *hNPK[2]={0,0}; 
  Int_t colorestmp[2]={2,4};
  for(Int_t i=0; i<2; ++i){

    hNPiP[i]=(TH1D *)fin[i]->Get("hPionProtonSeparation");
    hNPiP[i]->SetLineColor(colorestmp[i]);
    hNPiP[i]->SetLineWidth(2);
    hNPiK[i]=(TH1D *)fin[i]->Get("hPionKaonSeparation");
    hNPiK[i]->SetLineColor(colorestmp[i]);
    hNPiK[i]->SetLineWidth(2);
    hNPK[i]=(TH1D *)fin[i]->Get("hKaonProtonSeparation");
    hNPK[i]->SetLineColor(colorestmp[i]);
    hNPK[i]->SetLineWidth(2);

  }


  TH1D *hframe[3]={0,0,0};
  for(Int_t i=0; i<3; ++i){

    hframe[i] = new TH1D(Form("hframe%d",i),";#it{p} (GeV/#it{c});",10,3,15);
    hframe[i]->GetYaxis()->SetRangeUser(-0.01,7.01);
    hframe[i]->GetXaxis()->SetTitleSize(0.05);
    hframe[i]->GetYaxis()->SetTitleSize(0.05);
  }


  TCanvas* vC1 = new TCanvas("vC1","vC1",200,10,1050,400);
  TPad *pad[3]={0,0,0};
  for(Int_t i=0; i<3; ++i){

    pad[i]=new TPad(Form("pad%d",i),"pad1",0.01+0.33*i,0.01,0.33+0.33*i,0.99,0);
    pad[i]->SetLeftMargin(0.12);
    pad[i]->SetRightMargin(0.01);
    pad[i]->SetTopMargin(0.01);
    pad[i]->SetBottomMargin(0.1);
    pad[i]->SetBorderSize(0);
    pad[i]->SetBorderMode(0);

  }

 
  
  
  for(Int_t i =0 ; i<3; ++i){
    
    


    vC1->cd();
    pad[i]->Draw();
    
    
    pad[i]->cd();
    
    
    pad[i]->SetTickx(1);
    pad[i]->SetTicky(1);
    
    if(i==0)
      hframe[i]->GetYaxis()->SetTitle("Separation #pi/p, N_{#sigma}");
    if(i==1)
      hframe[i]->GetYaxis()->SetTitle("Separation #pi/K, N_{#sigma}");
    if(i==2)
      hframe[i]->GetYaxis()->SetTitle("Separation p/K, N_{#sigma}");


    hframe[i]->Draw();

    if(i==0){
      for(Int_t j =0 ; j<2; ++j){

	hNPiP[j]->Draw("samel");
	
      }
      TF1 *fpip1=new TF1("fpip1","5.0+pol0",5,15);
      fpip1->SetLineColor(1);
      fpip1->SetLineWidth(2);
      fpip1->SetLineStyle(2);
      fpip1->Draw("same");

      TF1 *fpip2=new TF1("fpip2","3.5+pol0",5,15);
      fpip2->SetLineColor(1);
      fpip2->SetLineWidth(2);
      fpip2->SetLineStyle(4);
      fpip2->Draw("same");

      TLatex* latex0 = new TLatex();
      latex0->SetNDC();
      latex0->SetTextAlign(22);
      latex0->SetTextFont(42);
      latex0->SetTextSize(0.05);
      latex0->DrawLatex(0.5,0.77,"pp @ 2.76 TeV, 0.6<|#eta|<0.8 (final)");

      latex0->DrawLatex(0.5,0.57,"Pb-Pb, 0-5%, 0.6<|#eta|<0.8 (final)");

    }

    if(i==1){
      for(Int_t j =0 ; j<2; ++j){
	hNPiK[j]->Draw("samel");
	
      }


      TF1 *fpik1=new TF1("fpik1","3.5+pol0",5,15);
      fpik1->SetLineColor(1);
      fpik1->SetLineWidth(2);
      fpik1->SetLineStyle(2);
      fpik1->Draw("same");

      TF1 *fpik2=new TF1("fpik2","2.2+pol0",5,15);
      fpik2->SetLineColor(1);
      fpik2->SetLineWidth(2);
      fpik2->SetLineStyle(4);
      fpik2->Draw("same");

      TLegend* legend = new TLegend(0.51, 0.65, 0.85, 0.85);    
      legend->SetBorderSize(0);
      legend->SetFillColor(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.05);
      legend->SetLineColor(kWhite);
      legend->SetLineStyle(3);
      legend->SetShadowColor(kWhite);
      legend->SetFillColor(kWhite);
      legend->SetFillStyle(0);

      legend->AddEntry(hNPiK[0], "|#eta|<0.2", "L");
      legend->AddEntry(hNPiK[1], "0.6<|#eta|<0.8", "L");
      legend->Draw();


    }

    if(i==2){
      for(Int_t j =0 ; j<2; ++j){
	hNPK[j]->Draw("samel");
	
      }

      TF1 *fpk1=new TF1("fpk1","1.8+pol0",5,15);
      fpk1->SetLineColor(1);
      fpk1->SetLineWidth(2);
      fpk1->SetLineStyle(2);
      fpk1->Draw("same");

      TF1 *fpk2=new TF1("fpk2","1.2+pol0",5,15);
      fpk2->SetLineColor(1);
      fpk2->SetLineWidth(2);
      fpk2->SetLineStyle(4);
      fpk2->Draw("same");

    }

  }


    
  vC1->SaveAs("NSigma.png");
  vC1->SaveAs("NSigma.eps");
    






}



//____________________________________________________________________________
void MakeNSigmaPlot(const Char_t* inFile, const Char_t* fitFileName,
		    Double_t ptStart, Double_t ptStop, Int_t i_eta,
		    const Char_t* outFileName=0)
{

  const Char_t* ending[4] = {"02", "24", "46", "68"};
  gStyle->SetOptStat(0);
  
  if(fitFileName) {
    
    TFile* fitFile = FindFileFresh(fitFileName);
    if(!fitFile)
      return;
    DeDxFitInfo* fitPar = (DeDxFitInfo*)fitFile->Get("fitInfo");
    fitPar->Print();
    
    fixMIP      = fitPar->MIP;
    fixPlateau  = fitPar->plateau;

    Double_t dedxPar[6]  = {0, 0, 0, 0, 0, 0};
    Double_t sigmaPar[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    
    dedxPar[0] = fitPar->optionDeDx;
    for(Int_t i = 0; i < fitPar->nDeDxPar; i++) {
      dedxPar[i+1] = fitPar->parDeDx[i];
    }

    sigmaPar[0] = fitPar->optionSigma;
    for(Int_t i = 0; i < fitPar->nSigmaPar; i++) {
      sigmaPar[i+1] = fitPar->parSigma[i];
    }

    piFunc = new TF1("piFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    piFunc->SetParameters(dedxPar);

 
    kFunc = new TF1("kFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    kFunc->SetParameters(dedxPar);
    kFunc->SetParameter(0, kFunc->GetParameter(0)+10);

    pFunc = new TF1("pFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    pFunc->SetParameters(dedxPar);
    pFunc->SetParameter(0, pFunc->GetParameter(0)+20);

    eFunc = new TF1("eFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
    eFunc->SetParameters(dedxPar);
    eFunc->SetParameter(0, eFunc->GetParameter(0)+30);

    sigmaFunc = new TF1("sigmaFunc", SigmaFunc, 0, 100, fitPar->nSigmaPar+1); 
    sigmaFunc->SetParameters(sigmaPar);
  }
 
  kFunc->SetLineColor(kGreen+2);
  pFunc->SetLineColor(4);
  eFunc->SetLineColor(kMagenta);
 

  /*
  TString baseName(gSystem->BaseName(calibFileName));
  baseName.ReplaceAll(".root", "");
  
  TFile* calibFile = FindFileFresh(calibFileName);
  if(!calibFile)
    return;
  AliHighPtDeDxCalib* calib = (AliHighPtDeDxCalib*)GetObject(calibFile, filter, usePhiCut, run, etaAbs, etaLow, etaHigh);
  if(calib) {    
    fixMIP      = calib->GetHistDeDx(kTRUE, 0)->GetMean();
    fixPlateau  = calib->GetHistDeDx(kFALSE, 0)->GetMean();
    cout << "Plateau: " << fixPlateau << endl;
  } else {
    calib = (AliHighPtDeDxCalib*)calibFile->Get("v0phicut");
    fixMIP = 50.0;
    fixPlateau  = 83.461;
  }
  calib->Print();

  hDeDxVsP = calib->GetHistDeDxVsP(0);
  */

  TFile* dedxFile = FindFileFresh(inFile);
  if(!dedxFile)
    return;

  TList * list = (TList *)dedxFile->Get("outputdedx");
  TH2D *hDeDxVsPPlus = 0;
  hDeDxVsPPlus = (TH2D *)list->FindObject(Form("histAllCh%s",ending[i_eta]));
  hDeDxVsPPlus->Sumw2();

  TH2D *hDeDxVsPMinus = 0;
  hDeDxVsPMinus = (TH2D *)list->FindObject(Form("histAllCh%s",ending[i_eta]));
  hDeDxVsPMinus->Sumw2();
  
  gSystem->Exec("rm *.root");
  TFile *fplus = new TFile("hist2Dplus.root","recreate");
  fplus->cd();
  hDeDxVsPPlus->SetName(Form("histAllCh%s", ending[i_eta]));
  hDeDxVsPPlus->Write();
  fplus->Close();
  
  TFile *fminus = new TFile("hist2Dminus.root","recreate");
  fminus->cd();
  hDeDxVsPMinus->SetName(Form("histAllCh%s", ending[i_eta]));
  hDeDxVsPMinus->Write();
  fminus->Close();
  
  gSystem->Exec("hadd hist2Ddedx.root hist2Dplus.root hist2Dminus.root");
  
  TFile *ffinal = TFile::Open("hist2Ddedx.root"); 
  hDeDxVsP = (TH2D *)ffinal->Get(Form("histAllCh%s", ending[i_eta]));




   //in this case, sigma is not the relative resolution

  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeDxVsP->GetXaxis()->FindBin(ptStart+0.01);
  ptStart = hDeDxVsP->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeDxVsP->GetXaxis()->FindBin(ptStop-0.01);
  ptStop  = hDeDxVsP->GetXaxis()->GetBinUpEdge(binStop);
  //  const Int_t nBins    = binStop - binStart + 1;

  cout << "Doing fits from pTlow = " << ptStart << " (bin: " << binStart
       << ") to pThigh = " << ptStop << " (bin: " << binStop << ")" << endl;


  TH1D* hPionRatio =(TH1D*)hDeDxVsP->ProjectionX("hPionRatio", 1, 1);
  hPionRatio->Reset();

  TH1D* hPionKaonSeparation = (TH1D*)hPionRatio->Clone("hPionKaonSeparation");
  hPionKaonSeparation->SetTitle(Form("#pi/K separation, %s; #it{p} (GeV/#it{c}); N_{#sigma}",legEtaCut[i_eta]));
  TH1D* hPionProtonSeparation = (TH1D*)hPionRatio->Clone("hPionProtonSeparation");
  hPionProtonSeparation->SetTitle(Form("#pi/p separation, %s; #it{p} (GeV/#it{c}); N_{#sigma}",legEtaCut[i_eta]));
  TH1D* hKaonProtonSeparation = (TH1D*)hPionRatio->Clone("hKaonProtonSeparation");
  hKaonProtonSeparation->SetTitle(Form("K/p separation, %s; #it{p} (GeV/#it{c}); N_{#sigma}",legEtaCut[i_eta]));




  for(Int_t bin = binStart; bin <= binStop; bin++){
    
    Double_t pt=hPionRatio->GetBinCenter(bin);

    Double_t dEdx_pi = piFunc->Eval(pt);
    Double_t Sigma_pi = sigmaFunc->Eval(piFunc->Eval(pt));

    Double_t dEdx_k = kFunc->Eval(pt);
    Double_t Sigma_k = sigmaFunc->Eval(kFunc->Eval(pt));

    Double_t dEdx_p = pFunc->Eval(pt);
    Double_t Sigma_p = sigmaFunc->Eval(pFunc->Eval(pt));

    Double_t N1 = (dEdx_pi - dEdx_k)/( (Sigma_pi + Sigma_k)/2.0 );
    Double_t N2 = (dEdx_pi - dEdx_p)/( (Sigma_pi + Sigma_p)/2.0 );
    Double_t N3 = (dEdx_k - dEdx_p)/( (Sigma_k + Sigma_p)/2.0 );

    hPionKaonSeparation->SetBinContent(bin,N1);
    hPionProtonSeparation->SetBinContent(bin,N2);
    hKaonProtonSeparation->SetBinContent(bin,N3);

  }



  if(outFileName) {

    CreateDir("nsigma_results");
    TFile* fileOut = new TFile(Form("nsigma_results/%s", outFileName), "RECREATE");

    hPionKaonSeparation->Write();
    hPionProtonSeparation->Write();
    hKaonProtonSeparation->Write();



    fileOut->Close();
  }



}



//________________________________________________________________________
void PlotParametrizations(const Char_t *path){


  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);

  TFile *fin[2]={0,0};


  fin[0] = FindFileFresh(Form("%s/hres_0_5_02.root",path));
  if(!fin[0])
    return;

  fin[1] = FindFileFresh(Form("%s/hres_0_5_68.root",path));
  if(!fin[1])
    return;


  TH1D *hframeRes = new TH1D("hframeRes","; #LT dE/dx #GT; #sigma/#LT dE/dx #GT",50,40,79);
  hframeRes->GetYaxis()->SetRangeUser(0,0.105);
  hframeRes->GetYaxis()->SetTitleSize(0.05);
  hframeRes->GetXaxis()->SetTitleSize(0.05);

  TH1D *hframeBB = new TH1D("hframeBB","; #beta#gamma;#LT dE/dx #GT",50,2,59);
  hframeBB->GetYaxis()->SetRangeUser(40,81);
  hframeBB->GetYaxis()->SetTitleSize(0.05);
  hframeBB->GetXaxis()->SetTitleSize(0.05);

  TGraphErrors *gRes[2] = {0,0}; 

  gRes[0] = (TGraphErrors *)fin[0]->Get("res_allpions");
  gRes[1] = (TGraphErrors *)fin[1]->Get("res_allpions");


  TGraphErrors *gBB[2] = {0,0}; 

  gBB[0] = (TGraphErrors *)fin[0]->Get("gBB");
  gBB[1] = (TGraphErrors *)fin[1]->Get("gBB");

  Int_t colorestmp[2]={2,4};



  TCanvas* vC1 = new TCanvas("vC1","vC1",200,10,700,400);
  TPad * pad1=new TPad("pad1","pad1",0.01,0.01,0.5,0.99,0);
  TPad * pad2=new TPad("pad2","pad2",0.5,0.01,0.99,0.99,0);
  
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.01);
  pad1->SetTopMargin(0.01);
  pad1->SetBottomMargin(0.1);
  pad1->SetBorderSize(0);
  pad1->SetBorderMode(0);
  
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.01);
  pad2->SetBottomMargin(0.1);
  pad2->SetTopMargin(0);
  pad2->SetBorderSize(0);
  vC1->cd();
  pad1->Draw();
  
  
  pad1->cd();


  pad1->SetTickx(1);
  pad1->SetTicky(1);






  hframeRes->Draw();

  for(Int_t i =0 ; i<2; ++i){
    gRes[i]->SetMarkerStyle(20);
    gRes[i]->SetLineWidth(2);
    gRes[i]->SetLineColor(colorestmp[i]);
    gRes[i]->SetMarkerColor(colorestmp[i]);
    gRes[i]->Draw("samep");
  }


  TF1 *f0005 = new TF1("f0005","pol2",50,80);
  //  f0005->SetParameter(0,1.20000000000000000e+01);
  f0005->SetParameter(0,6.08803411659195118e-02);
  f0005->SetParameter(1,2.34760597152924448e-04);
  f0005->SetParameter(2,-2.81841683363273653e-06);

  f0005->SetLineColor(1);
  f0005->SetLineWidth(2);
  f0005->SetLineStyle(2);
  f0005->Draw("samel");


  TF1 *fpp = new TF1("fpp","pol2",50,80);
  fpp->SetParameter(0,1.06773399595441187e-01);
  fpp->SetParameter(1,-1.55087701882609280e-03);
  fpp->SetParameter(2,1.01579672586848698e-05);

  fpp->SetLineColor(1);
  fpp->SetLineWidth(2);
  fpp->SetLineStyle(3);
  fpp->Draw("samel");


  TLegend* legendRes = new TLegend(0.21, 0.15, 0.65, 0.35);    
  legendRes->SetBorderSize(0);
  legendRes->SetFillColor(0);
  legendRes->SetTextFont(42);
  legendRes->SetTextSize(0.05);
  legendRes->SetLineColor(kWhite);
  legendRes->SetLineStyle(3);
  legendRes->SetShadowColor(kWhite);
  legendRes->SetFillColor(kWhite);
  legendRes->SetFillStyle(0);
  legendRes->SetHeader("Final, 0.6<|#eta|<0.8");
  legendRes->AddEntry(f0005, "Pb-Pb 0-5%", "L");
  legendRes->AddEntry(fpp, "pp @ 2.76 TeV", "L");
  legendRes->Draw();

  vC1->cd();
  pad2->Draw();
  
  
  pad2->cd();


  pad2->SetTickx(1);
  pad2->SetTicky(1);






  hframeBB->Draw();

  for(Int_t i =0 ; i<2; ++i){
    gBB[i]->SetMarkerStyle(20);
    gBB[i]->SetLineWidth(2);
    gBB[i]->SetLineColor(colorestmp[i]);
    gBB[i]->SetMarkerColor(colorestmp[i]);
    gBB[i]->Draw("samep");
  }

  TLegend* legend = new TLegend(0.51, 0.8, 0.85, 0.95);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.05);
  legend->SetLineColor(kWhite);
  legend->SetLineStyle(3);
  legend->SetShadowColor(kWhite);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  
  legend->AddEntry(gBB[0], "|#eta|<0.2", "P");
  legend->AddEntry(gBB[1], "0.6<|#eta|<0.8", "P");
  legend->Draw();




  vC1->SaveAs("Parametrizations.png");
  vC1->SaveAs("Parametrizations.eps");



}
//____________________________________________________________________________
void ExtractUncPartFractvsP(const Char_t * inFile,  Double_t ptStart, Double_t ptStop,
			    Int_t i_cent=1,
			    Int_t i_eta=1,
			    const Char_t* dirName="debugfits", const Char_t* outFileName=0)


{


  const Char_t* ending[4] = {"02", "24", "46", "68"};
  const Double_t etaHigh[4]={2,4,6,8};
  gStyle->SetOptStat(0);
  
  
  
  TFile* fitFile = FindFileFresh(Form("fitparameters/MB/%s_dataPbPb.root",ending[i_eta]));
  if(!fitFile)
    return;
  
 
  DeDxFitInfo* fitPar = (DeDxFitInfo*)fitFile->Get("fitInfo");
  fitPar->Print();
  
  fixMIP      = fitPar->MIP;
 
  

  //return;

  Double_t dedxPar[6]  = {0, 0, 0, 0, 0, 0};
  Double_t sigmaPar[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  
  dedxPar[0] = fitPar->optionDeDx;
  for(Int_t i = 0; i < fitPar->nDeDxPar; i++) {
    dedxPar[i+1] = fitPar->parDeDx[i];
    cout<<"idedx_par"<<i+1<<"="<<dedxPar[i+1]<<endl;

  }
  fixPlateau  = dedxPar[5];
  cout<<"fixPlateau="<<fixPlateau<<"  fixMIP="<<fixMIP<<endl;
  //return;
  
  sigmaPar[0] = fitPar->optionSigma;
  for(Int_t i = 0; i < fitPar->nSigmaPar; i++) {
    sigmaPar[i+1] = fitPar->parSigma[i];
  }
  
  piFunc = new TF1("piFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
  piFunc->SetParameters(dedxPar);
  
  
  kFunc = new TF1("kFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
  kFunc->SetParameters(dedxPar);
  kFunc->SetParameter(0, kFunc->GetParameter(0)+10);
  
  pFunc = new TF1("pFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
  pFunc->SetParameters(dedxPar);
  pFunc->SetParameter(0, pFunc->GetParameter(0)+20);
  
  eFunc = new TF1("eFunc", FitFunc, 0, 100, fitPar->nDeDxPar+1);
  eFunc->SetParameters(dedxPar);
  eFunc->SetParameter(0, eFunc->GetParameter(0)+30);
  
  sigmaFunc = new TF1("sigmaFunc", SigmaFunc, 0, 100, fitPar->nSigmaPar+1); 
  sigmaFunc->SetParameters(sigmaPar);
  
  
  kFunc->SetLineColor(kGreen+2);
  pFunc->SetLineColor(4);
  eFunc->SetLineColor(kMagenta);
 
  /*
    TFile* calibFile = FindFileFresh(inFile);
    if(!calibFile)
    return;
    hDeDxVsP = (TH2D *)calibFile->Get(Form("hDeDxVsP%s_%s", ending[i_eta], endingCent[i_cent]));
    hDeDxVsP->Sumw2();
  */

  TFile* dedxFile = FindFileFresh(inFile);
  if(!dedxFile)
    return;

  TList * list = (TList *)dedxFile->Get("outputdedx");
  TH2D *hDeDxVsPPlus = 0;
  hDeDxVsPPlus = (TH2D *)list->FindObject(Form("histAllCh%s",ending[i_eta]));
  hDeDxVsPPlus->Sumw2();

  TH2D *hDeDxVsPMinus = 0;
  hDeDxVsPMinus = (TH2D *)list->FindObject(Form("histAllCh%s",ending[i_eta]));
  hDeDxVsPMinus->Sumw2();
  
  gSystem->Exec("rm *.root");
  TFile *fplus = new TFile("hist2Dplus.root","recreate");
  fplus->cd();
  hDeDxVsPPlus->SetName(Form("histAllCh%s", ending[i_eta]));
  hDeDxVsPPlus->Write();
  fplus->Close();
  
  TFile *fminus = new TFile("hist2Dminus.root","recreate");
  fminus->cd();
  hDeDxVsPMinus->SetName(Form("histAllCh%s", ending[i_eta]));
  hDeDxVsPMinus->Write();
  fminus->Close();
  
  gSystem->Exec("hadd hist2Ddedx.root hist2Dplus.root hist2Dminus.root");
  
  TFile *ffinal = TFile::Open("hist2Ddedx.root"); 
  hDeDxVsP = (TH2D *)ffinal->Get(Form("histAllCh%s", ending[i_eta]));
  






  CreateDir(Form("fit_results/%s_%s",dirName,endingCent[i_cent]));
  
  TCanvas* cDeDxVsPLogX = new TCanvas("cDeDxVsPLogX", "; #it{p} (GeV/#it{c}) ;d#it{E}/d#it{x}", 500, 400);
  cDeDxVsPLogX->Clear();
  cDeDxVsPLogX->cd();
  cDeDxVsPLogX->SetLogz();
  
  TH2D *hDeDxVsP2=(TH2D *)hDeDxVsP->Clone("h2DP");
  hDeDxVsP2->SetTitle("; #it{p} (GeV/#it{c}) ;d#it{E}/d#it{x}");
  hDeDxVsP2->GetXaxis()->SetRangeUser(2,19.5);
  hDeDxVsP2->Draw("COLZ");
  
  piFunc->SetLineColor(1);
  piFunc->Draw("samel");
  
  eFunc->SetLineColor(1);
  eFunc->Draw("samel");
  
  kFunc->SetLineColor(1);
  kFunc->Draw("samel");
  
  pFunc->SetLineColor(1);
  pFunc->Draw("samel");
  
  TLatex* latex0 = new TLatex();
  latex0->SetNDC();
  latex0->SetTextAlign(22);
  latex0->SetTextSize(0.05);
  latex0->DrawLatex(0.7,0.94+0.035,Form("p-Pb, %s",legEtaCut[i_eta]));
  latex0->DrawLatex(0.71,0.89+0.035,Form("#sqrt{s_{NN}}=5.02 TeV, %s",CentLatex[i_cent]));

  
  
  cDeDxVsPLogX->SaveAs(Form("fit_results/%s_%s/%s.png",dirName,endingCent[i_cent],cDeDxVsPLogX->GetName()));
  cDeDxVsPLogX->SaveAs(Form("fit_results/%s_%s/%s.eps",dirName,endingCent[i_cent],cDeDxVsPLogX->GetName()));
  //in this case, sigma is not the relative resolution
  
  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeDxVsP->GetXaxis()->FindBin(ptStart+0.01);
  ptStart = hDeDxVsP->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeDxVsP->GetXaxis()->FindBin(ptStop-0.01);
  ptStop  = hDeDxVsP->GetXaxis()->GetBinUpEdge(binStop);
  //  const Int_t nBins    = binStop - binStart + 1;
  
  cout << "Doing fits from pTlow = " << ptStart << " (bin: " << binStart
       << ") to pThigh = " << ptStop << " (bin: " << binStop << ")" << endl;
  
  
  //cross check
  TCanvas* cFits = new TCanvas("cFits", "Fit comparison to data", 1200, 800);
  cFits->Clear();
  cFits->Divide(7, 4);
  
  
  TF1* pion = new TF1("pion", "gausn", 40, 100);
  
  TF1* kaon = new TF1("kaon", "gausn", 40, 100);
  
  TF1* proton = new TF1("proton", "gausn", 40, 100);
  //proton->SetLineWidth(2);
  //proton->SetLineColor(kBlue);
  TF1* electron = new TF1("electron", "gausn", 40, 100);
  //electron->SetLineWidth(2);
  //electron->SetLineColor(kMagenta);
  //  TF1* total = new TF1("total", "gausn(0)+gausn(3)+gausn(6)+gausn(9)", -30, 30);
  TF1* total = 0;
  total = new TF1("total", "([0]+[12])*gausn(1)+[4]*gausn(5)+[8]*gausn(9)+[12]*gausn(13)", 40, 100);
  
  total->SetLineColor(kBlack);
  total->SetLineWidth(2);
  total->SetLineStyle(2);
  
  TLegend* legend = new TLegend(0.11, 0.6, 0.35, 0.85);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(total, "4-Gauss fit", "L");
  legend->AddEntry(pion, "#pi", "L");
  legend->AddEntry(kaon, "K", "L");
  legend->AddEntry(proton, "p", "L");
  legend->AddEntry(electron, "e", "L");
  
  TCanvas* cSingleFit = new TCanvas("cSingleFit", "single fit");
  cSingleFit->SetLeftMargin(0.11);
  cSingleFit->SetRightMargin(0.03);
  cSingleFit->SetTopMargin(0.01);
  cSingleFit->SetBottomMargin(0.1);
  cSingleFit->SetBorderSize(0);
  cSingleFit->SetBorderMode(0);
  
  
  TH1D* hPionRatio =(TH1D*)hDeDxVsP->ProjectionX("hPionRatio", 1, 1);
  hPionRatio->SetTitle("particle fractions; p [GeV/c]; particle fraction");
  hPionRatio->Reset();
  TH1D* hKaonRatio   = (TH1D*)hPionRatio->Clone("hKaonRatio");
  TH1D* hProtonRatio = (TH1D*)hPionRatio->Clone("hProtonRatio");
  TH1D* hElectronRatio = (TH1D*)hPionRatio->Clone("hElectronRatio");
  
    
    
  TH1D* hPionYield =(TH1D*)hDeDxVsP->ProjectionX("hPionYield", 1, 1);
  hPionYield->SetTitle("particle fractions; p [GeV/c]; particle fraction");
  hPionYield->Reset();
  TH1D* hKaonYield   = (TH1D*)hPionYield->Clone("hKaonYield");
  TH1D* hProtonYield = (TH1D*)hPionYield->Clone("hProtonYield");
  TH1D* hElectronYield = (TH1D*)hPionYield->Clone("hElectronYield");
  
  
  TF1 *fmip=new TF1("fmip","pol0",0,1);
  fmip->SetParameter(0,fixMIP);
  
  
  
  TF1* fElectronFraction = 0;
  if(etaHigh[i_eta]==8||etaHigh[i_eta]==-6)
    fElectronFraction = new TF1("fElectronFraction", "[0]+(x<9.0)*[1]*(x-9.0)", 0, ptStop);
  if(etaHigh[i_eta]==6||etaHigh[i_eta]==-4)
    fElectronFraction = new TF1("fElectronFraction", "[0]+(x<8.0)*[1]*(x-8.0)", 0, ptStop);
  if(etaHigh[i_eta]==4||etaHigh[i_eta]==-2)
    fElectronFraction = new TF1("fElectronFraction", "[0]+(x<7.5)*[1]*(x-7.5)", 0, ptStop);
  if(etaHigh[i_eta]==2||etaHigh[i_eta]==0)
    fElectronFraction = new TF1("fElectronFraction", "[0]+(x<7.0)*[1]*(x-7.0)", 0, ptStop);
  
  //TF1* fElectronFraction = new TF1("fElectronFraction", "[0]+(x<7.0)*[1]*(x-7.0)", 0, ptStop);
  fElectronFraction->SetParameters(0.01, 0.0);
  
  
  TH1D *dEdxpoints[binStop-binStart];
  TF1 *fPion[binStop-binStart];
  TF1 *fKaon[binStop-binStart];
  TF1 *fProton[binStop-binStart];
  TF1 *fElectron[binStop-binStart];
  TF1 *fTotal[binStop-binStart];
  
  
  for(Int_t bin = binStart; bin <= binStop; bin++){
      
    Double_t pt=hPionRatio->GetBinCenter(bin);
    
    cout << "Making projection for bin: " << bin << endl;
    
    const Int_t j = bin-binStart;
    
    TH1D* hDeltaPiVsPtProj =(TH1D*)hDeDxVsP->ProjectionY(Form("hDeltaPiVsPtProj%d", bin), bin, bin);
    //    hDeltaPiVsPtProj->GetXaxis()->SetRangeUser(-25, 20);
    hDeltaPiVsPtProj->GetXaxis()->SetRangeUser(40, 90);
    hDeltaPiVsPtProj->SetTitle(Form("%.1f<p<%.1f GeV/c", 
				    hDeDxVsP->GetXaxis()->GetBinLowEdge(bin),
				    hDeDxVsP->GetXaxis()->GetBinUpEdge(bin)));
    
    dEdxpoints[j]=0;
    dEdxpoints[j]=hDeltaPiVsPtProj;
    fPion[j]=0;
    fPion[j]=new TF1(Form("fPiongauss_%d",bin), "gausn", 40, 100);
    
    fKaon[j]=0;
    fKaon[j]=new TF1(Form("fKaongauss_%d",bin), "gausn", 40, 100);
    
    fProton[j]=0;
    fProton[j]=new TF1(Form("fProtongauss_%d",bin), "gausn", 40, 100);
    
    fElectron[j]=0;
    fElectron[j]=new TF1(Form("fElectrongauss_%d",bin), "gausn", 40, 100);
    
    fTotal[j]=0;
    fTotal[j]=new TF1(Form("fTotalgauss_%d",bin), "([0]+[12])*gausn(1)+[4]*gausn(5)+[8]*gausn(9)+[12]*gausn(13)", 40, 100);
    
    
    Double_t all =  hDeltaPiVsPtProj->GetEntries();
    
    
    const Int_t nPar = 16;
    Double_t gausParams[nPar] = { 
      0.59,
      all,
      piFunc->Eval(pt), 
      //fpi->Eval(pt),
      sigmaFunc->Eval(piFunc->Eval(pt)),
      
      0.3,
      all,
      kFunc->Eval(pt), 
      //fp->Eval(pt),
      sigmaFunc->Eval(kFunc->Eval(pt)),
      
      0.1,
      all,
      pFunc->Eval(pt), 
      //fp->Eval(pt),
      sigmaFunc->Eval(pFunc->Eval(pt)),
      
      0.01,
      all,
      eFunc->Eval(pt),
      //fpi->Eval(pt),
      sigmaFunc->Eval(eFunc->Eval(pt)),
    };
    
   
    for(Int_t ipar=0;ipar<nPar;++ipar)
      cout<<gausParams[ipar]<<endl;
    
    
    
    cFits->cd();
    cFits->cd(j + 1);
    
    total->SetParameters(gausParams);
    
    
    for(Int_t i = 0; i < nPar; i++) {
      if(i!=0 && i!=4 && i!=8 && i!=12)
	total->FixParameter(i, gausParams[i]);
      else
	continue;
    }
    
    total->SetParLimits(7, gausParams[7]-0.05*gausParams[7],gausParams[7]+0.05*gausParams[7]);
    total->SetParLimits(8, 0.005,0.4);
    total->SetParLimits(4, 0.005,0.6);
    
    if(bin==48) {
      if(etaHigh[i_eta]==8||etaHigh[i_eta]==-6)
	hElectronRatio->Fit(fElectronFraction, "N", "", 4.0, 10.0);
      if(etaHigh[i_eta]==6||etaHigh[i_eta]==-4)
	hElectronRatio->Fit(fElectronFraction, "N", "", 3.66, 10.0);
      if(etaHigh[i_eta]==4||etaHigh[i_eta]==-2)
	hElectronRatio->Fit(fElectronFraction, "N", "", 3.33, 10.0);
      if(etaHigh[i_eta]==2||etaHigh[i_eta]==0)
	hElectronRatio->Fit(fElectronFraction, "N", "", 3.0, 10.0);
      fElectronFraction->SetRange(0, ptStop);
    }
    
    if(bin>=48) {
      total->FixParameter(12, fElectronFraction->Eval(hDeDxVsP->GetXaxis()->GetBinCenter(bin)));
    }
     
      
    hDeltaPiVsPtProj->Fit(total, "0L");
    
    
    hDeltaPiVsPtProj->Draw();
    total->DrawCopy("same");    
    
    Double_t parametersOut[nPar];
    total->GetParameters(parametersOut);
    const Double_t* parameterErrorsOut = total->GetParErrors();
    
    for(Int_t i = 0; i < nPar; i++) 
      cout << parametersOut[i] << ", ";
    cout << endl;
    fTotal[j]->SetParameters(parametersOut);
    
    

    pion->SetParameters(&parametersOut[1]);
    pion->SetParameter(0,all*(parametersOut[0]+parametersOut[12]));
    
    fPion[j]->SetParameters(&parametersOut[1]);;
    fPion[j]->SetParameter(0,all*(parametersOut[0]+parametersOut[12]));
    
    
    hPionRatio->SetBinContent(bin, parametersOut[0]);
    hPionRatio->SetBinError(bin, parameterErrorsOut[0]);
    hPionYield->SetBinContent(bin, parametersOut[0]*all);
    hPionYield->SetBinError(bin, parameterErrorsOut[0]*all);
      
    kaon->SetParameters(&parametersOut[5]);
    kaon->SetParameter(0,all*parametersOut[4]);
    
    fKaon[j]->SetParameters(&parametersOut[5]);
    fKaon[j]->SetParameter(0,all*parametersOut[4]);
    
    
    hKaonRatio->SetBinContent(bin, parametersOut[4]);
    hKaonRatio->SetBinError(bin, parameterErrorsOut[4]);
    hKaonYield->SetBinContent(bin, parametersOut[4]*all);
    hKaonYield->SetBinError(bin, parameterErrorsOut[4]*all);
    
    proton->SetParameters(&parametersOut[9]);
    proton->SetParameter(0,all*parametersOut[8]);
    //proton->DrawCopy("same");
    
    fProton[j]->SetParameters(&parametersOut[9]);
    fProton[j]->SetParameter(0,all*parametersOut[8]);
    
    
    hProtonRatio->SetBinContent(bin, parametersOut[8]);
    hProtonRatio->SetBinError(bin, parameterErrorsOut[8]);
    hProtonYield->SetBinContent(bin, parametersOut[8]*all);
    hProtonYield->SetBinError(bin, parameterErrorsOut[8]*all);
    
    electron->SetParameters(&parametersOut[13]);
    electron->SetParameter(0,all*parametersOut[12]);
    
    fElectron[j]->SetParameters(&parametersOut[13]);
    fElectron[j]->SetParameter(0,all*parametersOut[12]);
      
      
    //electron->DrawCopy("same");
    hElectronRatio->SetBinContent(bin, parametersOut[12]);
    hElectronRatio->SetBinError(bin, parameterErrorsOut[12]);
    hElectronYield->SetBinContent(bin, parametersOut[12]*all);
    hElectronYield->SetBinError(bin, parameterErrorsOut[12]*all);
    
    //DrawALICELogo(kFALSE, 0.71, 0.76, 0.82, 0.98);
    cSingleFit->cd();
    cSingleFit->Clear();
    cSingleFit->SetTickx(1);
    cSingleFit->SetTicky(1);
    cSingleFit->SetLogy(0);
    
    //    cSingleFit->SetLogy();
    hDeltaPiVsPtProj->GetXaxis()->SetRangeUser(40,90);
    hDeltaPiVsPtProj->SetMarkerStyle(25);
    hDeltaPiVsPtProj->SetMarkerColor(1);
    hDeltaPiVsPtProj->SetLineColor(1);
    hDeltaPiVsPtProj->SetLineWidth(2);
    hDeltaPiVsPtProj->GetYaxis()->SetLabelFont(42);
    hDeltaPiVsPtProj->GetXaxis()->SetLabelFont(42);
    hDeltaPiVsPtProj->GetZaxis()->SetLabelFont(42);
    hDeltaPiVsPtProj->GetYaxis()->SetTitleFont(42);
    hDeltaPiVsPtProj->GetXaxis()->SetTitleFont(42);
    hDeltaPiVsPtProj->GetZaxis()->SetTitleFont(42);
    hDeltaPiVsPtProj->GetYaxis()->SetTitle("Entries");
    hDeltaPiVsPtProj->GetXaxis()->SetTitle("d#it{E}/d#it{x} in TPC (arb. units)");
    hDeltaPiVsPtProj->GetXaxis()->SetTitleSize(0.05);
    hDeltaPiVsPtProj->GetXaxis()->SetTitleOffset(0.9);
    hDeltaPiVsPtProj->GetYaxis()->SetTitleSize(0.05);
    hDeltaPiVsPtProj->GetYaxis()->SetTitleOffset(0.9);
    hDeltaPiVsPtProj->Draw();
    
    total->SetLineColor(kGray+1);
    total->SetLineWidth(3);
    total->SetLineStyle(1);
    total->DrawCopy("same");
    
    pion->GetXaxis()->SetRangeUser(40,90);
    pion->SetLineColor(2);
    pion->SetLineWidth(2);
    //pion->SetLineStyle(2);
    pion->DrawCopy("same");
    
    kaon->SetLineColor(kGreen+2);
    kaon->SetLineWidth(2);
    //kaon->SetLineStyle(2);
    kaon->DrawCopy("same");
    
    
    proton->SetLineColor(4);
    proton->SetLineWidth(2);
    //proton->SetLineStyle(2);
    proton->DrawCopy("same");
    
    
    TLatex* latex = new TLatex();
    //  latex->SetTextFont(132);
    latex->SetNDC();
    latex->SetTextAlign(22);
    
    latex->SetTextSize(0.05);
    latex->SetTextFont(42);
    latex->SetTextColor(kRed+2);
    //latex->DrawLatex(0.8,0.3,"pp");
    //latex->DrawLatex(0.81,0.25,"#sqrt{s}=2.76 TeV");
    latex->DrawLatex(0.8,0.25+0.035,Form("p-Pb, %s",CentLatex[i_cent]));
    latex->DrawLatex(0.81,0.2+0.035,"#sqrt{s_{NN}}=5.02 TeV");
    
    latex->SetTextColor(1);
    latex->SetTextSize(0.04);
    //latex->DrawLatex(0.81,0.2,"25/07/2012");
    //latex->DrawLatex(0.80,0.3+0.036,"31/10/2012");
    
    latex->SetTextSize(0.05);
    latex->DrawLatex(0.8,0.71,legEtaCut[i_eta]);
    latex->DrawLatex(0.8,0.64,Form("%.1f<#it{p}<%.1f GeV/#it{c}", hDeDxVsP->GetXaxis()->GetBinLowEdge(bin), hDeDxVsP->GetXaxis()->GetBinUpEdge(bin)));
    
    
    
    //cSingleFit->SaveAs(Form("%s/ptspectrum_bin%d.gif", dirName, bin));
    cSingleFit->SaveAs(Form("fit_results/%s_%s/ptspectrum_bin%d.png", dirName, endingCent[i_cent],bin));
    cSingleFit->SaveAs(Form("fit_results/%s_%s/ptspectrum_bin%d.pdf", dirName, endingCent[i_cent],bin));
    cSingleFit->SaveAs(Form("fit_results/%s_%s/ptspectrum_bin%d.eps", dirName, endingCent[i_cent],bin));
    //    legend->Draw();
    
    
    
      
  }

    
  TCanvas* cRatio = new TCanvas("cRatio", "ratios/all vs p; pt", 600, 400);
  
  cRatio->Clear();
  cRatio->SetGridy();
  hElectronRatio->GetYaxis()->SetRangeUser(-0.05, 0.1);
  hElectronRatio->DrawCopy("P E");
  fElectronFraction->DrawCopy("SAME");
  
  gROOT->ProcessLine(".x drawText.C");
  cRatio->SaveAs(Form("fit_results/%s_%s/electron_ratio.gif", dirName,endingCent[i_cent]));
  cRatio->SaveAs(Form("fit_results/%s_%s/electron_ratio.pdf", dirName,endingCent[i_cent]));
  
  cRatio->Clear();
  cRatio->SetGridy(kFALSE);
  hPionRatio->SetMarkerStyle(20);
  hPionRatio->SetMarkerColor(2);
  hPionRatio->SetLineColor(2);
  hPionRatio->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hPionRatio->GetYaxis()->SetRangeUser(0.0, 1.0);
  hPionRatio->DrawCopy("P E");
  hPionYield->SetMarkerStyle(20);
  hPionYield->SetMarkerColor(2);
  hPionYield->SetLineColor(2);
  hPionYield->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  
  hKaonRatio->SetMarkerStyle(20);
  hKaonRatio->SetMarkerColor(kGreen);
  hKaonRatio->SetLineColor(kGreen);
  hKaonRatio->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hKaonYield->SetMarkerStyle(20);
  hKaonYield->SetMarkerColor(kGreen);
  hKaonYield->SetLineColor(kGreen);
  hKaonYield->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hKaonRatio->DrawCopy("SAME P E");
  
  hProtonRatio->SetMarkerStyle(20);
  hProtonRatio->SetMarkerColor(4);
  hProtonRatio->SetLineColor(4);
  hProtonRatio->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hProtonYield->SetMarkerStyle(20);
  hProtonYield->SetMarkerColor(4);
  hProtonYield->SetLineColor(4);
  hProtonYield->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hProtonRatio->DrawCopy("SAME P E");
  
  hElectronRatio->SetMarkerStyle(20);
  hElectronRatio->SetMarkerColor(6);
  hElectronRatio->SetLineColor(6);
  hElectronRatio->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  hElectronRatio->DrawCopy("SAME P E");
  hElectronYield->SetMarkerStyle(20);
  hElectronYield->SetMarkerColor(6);
  hElectronYield->SetLineColor(6);
  hElectronYield->GetXaxis()->SetRangeUser(0.0, ptStop-0.001);
  
  TH1D *hallRec=0;
  
  hallRec=(TH1D *)hPionRatio->Clone("hRecAll");
  hallRec->Add(hKaonRatio);
  hallRec->Add(hProtonRatio);
  hallRec->SetLineColor(kYellow-3);
  hallRec->SetMarkerColor(kYellow-3);
  hallRec->DrawCopy("SAME P E");
  
  
  gROOT->ProcessLine(".x drawText.C");
  cRatio->SaveAs(Form("fit_results/%s_%s/particle_ratios.gif", dirName,endingCent[i_cent]));
  cRatio->SaveAs(Form("fit_results/%s_%s/particle_ratios.pdf", dirName,endingCent[i_cent]));
  
  

  if(outFileName) {
    
    CreateDir(Form("fit_results/fit_yields_results_%s",endingCent[i_cent]));
    TFile* fileOut = new TFile(Form("fit_results/fit_yields_results_%s/%s%s.root", endingCent[i_cent], outFileName, ending[i_eta]), "RECREATE");
    
    
    hPionRatio->Write();
    hKaonRatio->Write();
    hProtonRatio->Write();
    hElectronRatio->Write();
    
    hPionYield->Write();
    hKaonYield->Write();
    hProtonYield->Write();
    hElectronYield->Write();
    
    
    
    fElectronFraction->Write();
    
    for(Int_t bin = binStart; bin <= binStop; bin++){
      
      Int_t j = bin-binStart;
      dEdxpoints[j]->Write();
      fPion[j]->Write();
      fKaon[j]->Write();
      fProton[j]->Write();
      fElectron[j]->Write();
      fTotal[j]->Write();
      
      
    }
    fmip->Write();
    fileOut->Close();
  }
  
  
  
  
}






//____________________________________________________________________________
void FitDeDxVsP(const Char_t * dedxFileName,
		Double_t pStart, Double_t pStop,  Int_t i_cent,
		Int_t optionDeDx, Int_t optionSigma,
		Bool_t performFit = kFALSE,
		Int_t etaLow=0, Int_t etaHigh=8, 
		Bool_t fixAllParBB=kFALSE,
		Bool_t fixAllParSigma=kFALSE,
		Bool_t fixKtoZero=kFALSE,
		const Char_t* initialFitFileName=0,
		Int_t fixParametersDedx=-1,
		Int_t fixParametersSigma=-1)
{
  //gStyle->SetOptStat(0);



  TFile* dedxFile = FindFileFresh(dedxFileName);
  if(!dedxFile)
    return;

  TList * list = (TList *)dedxFile->Get("outputdedx");
  TH2D *hDeDxVsPPlus = 0;
  hDeDxVsPPlus = (TH2D *)list->FindObject(Form("histAllCh%d%d", etaLow, etaHigh));
  hDeDxVsPPlus->Sumw2();

  TH2D *hDeDxVsPMinus = 0;
  hDeDxVsPMinus = (TH2D *)list->FindObject(Form("histAllCh%d%d", etaHigh, etaLow));
  hDeDxVsPMinus->Sumw2();

  gSystem->Exec("rm *.root");
  TFile *fplus = new TFile("hist2Dplus.root","recreate");
  fplus->cd();
  hDeDxVsPPlus->SetName(Form("histAllCh%d%d", etaLow, etaHigh));
  hDeDxVsPPlus->Write();
  fplus->Close();

  TFile *fminus = new TFile("hist2Dminus.root","recreate");
  fminus->cd();
  hDeDxVsPMinus->SetName(Form("histAllCh%d%d", etaLow, etaHigh));
  hDeDxVsPMinus->Write();
  fminus->Close();

  gSystem->Exec("hadd hist2Ddedx.root hist2Dplus.root hist2Dminus.root");

  TFile *ffinal = TFile::Open("hist2Ddedx.root"); 
  hDeDxVsP = (TH2D *)ffinal->Get(Form("histAllCh%d%d", etaLow, etaHigh));



  CreateDir("old");
  gSystem->Exec(Form("mv results/calibdedx_%d%d/* old/",etaLow, etaHigh));


  TCanvas* cDeDxVsP = new TCanvas("cDeDxVsP", "dE/dx vs p", 400, 300);
  cDeDxVsP->Clear();
  cDeDxVsP->cd();
  cDeDxVsP->SetLogz();
  hDeDxVsP->SetTitle(0);
  hDeDxVsP->Draw("COLZ");

  TCanvas* cDeDxVsPLogX = new TCanvas("cDeDxVsPLogX", "dE/dx vs p", 400, 300);
  cDeDxVsPLogX->Clear();
  cDeDxVsPLogX->cd();
  cDeDxVsPLogX->SetLogz();
  cDeDxVsPLogX->SetLogx();
  hDeDxVsP->Draw("COLZ");

  // Root is a bit stupid with finidng bins so we have to add and subtract a
  // little to be sure we get the right bin as we typically put edges as
  // limits
  const Int_t binStart = hDeDxVsP->GetXaxis()->FindBin(pStart+0.01);
  pStart = hDeDxVsP->GetXaxis()->GetBinLowEdge(binStart);
  const Int_t binStop  = hDeDxVsP->GetXaxis()->FindBin(pStop-0.01);
  pStop = hDeDxVsP->GetXaxis()->GetBinUpEdge(binStop);
  const Int_t nBins    = binStop - binStart + 1;

  cout << "Doing 2d fit from pTlow = " << pStart << " (bin: " << binStart
       << ") to pThigh = " << pStop << " (bin: " << binStop << ")" << endl;
  
  // Double_t binSize = (histdEdxvsP->GetXaxis()->GetXmax() - histdEdxvsP->GetXaxis()->GetXmin())/ histdEdxvsP->GetXaxis()->GetNbins();
  
  Double_t parDeDx[6]  = {0, 0, 0, 0, 0, 0};
  Double_t parSigma[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  
  const Int_t nLocalParams  = 3; // pi, K, p yields
  Int_t nDeDxPar      = 0;
  Int_t nSigmaPar     = 0;

  switch(optionDeDx) {
    
  case 1:
    nDeDxPar = 2;
    parDeDx[0] = 39.7;
    parDeDx[1] =  6.3;
    break;
  case 2:
    nDeDxPar = 1;
    parDeDx[0] =  7.3;
    break;
  case 3:
    nDeDxPar = 2;
    parDeDx[0] =  6.85097;
    parDeDx[1] =  29.5933;
    break;
  case 4:
    nDeDxPar = 3;
    parDeDx[0] =  31.2435;
    parDeDx[1] =  9.73037;
    parDeDx[2] =  1.95832;
    break;
  case 5:
    nDeDxPar = 4;
    parDeDx[0] =  31.384;
    parDeDx[1] =  9.5253;
    parDeDx[2] =  7.3542;
    parDeDx[3] =  1.4;
    break;
  case 6:
    nDeDxPar = 5;
    parDeDx[0] =  31.384;
    parDeDx[1] =  9.5253;
    parDeDx[2] =  7.3542;
    parDeDx[3] =  1.5;
    parDeDx[4] =  81;
    break;
  case 7:
    nDeDxPar = 3;
    parDeDx[0] =  31.3;
    parDeDx[1] =  9.5;
    parDeDx[2] =  1.4;
    break;
  default:

    cout << "optionDeDx does not support option: " << optionDeDx << endl;
    return;
    break;
  }

  switch(optionSigma) {
    
  case 1:
    nSigmaPar = 1;
    parSigma[0] = 3.0;
    break;
  case 2:
    nSigmaPar = 1;
    parSigma[0] = 0.0655688;
    break;
  case 3:
    nSigmaPar = 2;
    parSigma[0] = 0.06;
    parSigma[1] = -0.001;
  case 4:
    nSigmaPar = 2;
    parSigma[0] = 0.06;
    parSigma[1] = 1.0;
    break;
  case 5:
    nSigmaPar = 2;
    parSigma[0] = 0.06;
    parSigma[1] = 1.0;
    break;
  case 6:
    nSigmaPar = 2;
    parSigma[0] = 0.06;
    parSigma[1] = 0.0;
    break;
  case 7:
    nSigmaPar = 3;
    parSigma[0] = 0.06;
    parSigma[1] = 0.0;
    parSigma[2] = 2.0;
    break;
  case 8:
    nSigmaPar = 3;
    parSigma[0] = 0.06;
    parSigma[1] = 0.0;
    parSigma[2] = 2.0;
    break;
  case 9://works for 0-5 %
    nSigmaPar = 3;
    parSigma[0] = 1.88409e-01;
    parSigma[1] = -3.84073e-03;
    parSigma[2] = 3.03110e-05;
    break;
  case 10://works for 0-5 %
    nSigmaPar = 6;
    parSigma[0]=7.96380e-03;
    parSigma[1]=8.09869e-05;
    parSigma[2]=6.29707e-06;
    parSigma[3]=-3.27295e-08;
    parSigma[4]=-1.20200e+03;
    parSigma[5]=-3.97089e+01;
    break;
  case 11://works for 0-5 %
    nSigmaPar = 6;
    parSigma[0]=7.96380e-03;
    parSigma[1]=8.09869e-05;
    parSigma[2]=6.29707e-06;
    parSigma[3]=-3.27295e-08;
    parSigma[4]=-1.20200e+03;
    parSigma[5]=-3.97089e+01;
    break;
  case 12://works for 0-5 %
    nSigmaPar = 3;
    parSigma[0]=7.96380e-03;
    parSigma[1]=8.09869e-05;
    parSigma[2]=6.29707e-06;
    break;
  case 13://works for 0-5 %
    nSigmaPar = 3;
    parSigma[0]=7.96380e-03;
    parSigma[1]=8.09869e-05;
    parSigma[2]=6.29707e-06;
    break;
  case 14:
    nSigmaPar = 2;
    parSigma[0]=7.96380e-03;
    parSigma[1]=8.09869e-05;
    break;
  case 15:
    nSigmaPar = 3;
    parSigma[0]=7.96380e-03;
    parSigma[1]=8.09869e-05;
    break;




  default:

    cout << "optionSigma does not support option: " << optionSigma << endl;
    return;
    break;
  }

  if(initialFitFileName) {
    
  

    TFile* fresBB = FindFileFresh(initialFitFileName);
    if(!fresBB)
      return;
    TF1 *fres=(TF1 *)fresBB->Get("sigmaRes");
    TF1 *fBB=(TF1 *)fresBB->Get("dedxFunc");

    Int_t nparres=fres->GetNpar();
    Int_t nparBB=fBB->GetNpar();
    for(Int_t ires=0;ires<nparres;++ires){
      parSigma[ires] = fres->GetParameter(ires+1);
      //cout<<fres->GetParameter(ires+1)<<endl;
    }
    for(Int_t iBB=0;iBB<nparBB;++iBB){
      parDeDx[iBB] = fBB->GetParameter(iBB+1);
      cout<<parDeDx[iBB]<<endl;

    }
    //Esto lo puse apenas
    fixPlateau  = parDeDx[4];
    fixMIP = fBB->Eval(3.5);
    cout<<"Hello!!"<<"   fixPlateau="<<fixPlateau<<"  fixMIP="<<fixMIP<<endl;
    
    //fixMIP      = fitPar->MIP;
    //return;
  }
 


  const Int_t nGlobalParams = 2  // binStart, binStop, 
    + 2 + nDeDxPar               // optionDeDx, nDeDxPar, dedxpar0 ....
    + 2 + nSigmaPar;             // nSigmaPar, sigmapar0 .....
  
  const Int_t nParams = nBins*nLocalParams + nGlobalParams;
  
  
  TF2* fitAll = new TF2("fitAll", fitf3G, pStart, pStop, 30, 90, nParams);
  Double_t parametersIn[nParams]; 
  
  parametersIn[0] = binStart;
  fitAll->SetParName(0,"binStart");
  fitAll->FixParameter(0, parametersIn[0]);

  parametersIn[1] = binStop;
  fitAll->SetParName(1,"binStop");
  fitAll->FixParameter(1, parametersIn[1]);

  // dE/dx parameters
  parametersIn[2] = nDeDxPar;
  fitAll->SetParName(2,"nDeDxPar");
  fitAll->FixParameter(2, parametersIn[2]);

  parametersIn[3] = optionDeDx;
  fitAll->SetParName(3,"optionDeDx");
  fitAll->FixParameter(3, parametersIn[3]);

  for(Int_t i = 0; i < nDeDxPar; i++) {

    parametersIn[4+i] = parDeDx[i];
    fitAll->SetParName(4+i,Form("dEdXpar%d", i));
    Int_t bit = TMath::Nint(TMath::Power(2, i));
    if(fixParametersDedx>=0)
      if (fixParametersDedx==0 || fixParametersDedx&bit) {
	
	cout << "Fixing dE/dx parameter " << i << endl;
	fitAll->FixParameter(4+i, parametersIn[4+i]);
      }
    
    if(fixAllParBB) {

      fitAll->FixParameter(4+i, parametersIn[4+i]);
      }
  }

  // sigma parameters
  parametersIn[4+nDeDxPar] = nSigmaPar;
  fitAll->SetParName(4+nDeDxPar,"nSigmaPar");
  fitAll->FixParameter(4+nDeDxPar, parametersIn[4+nDeDxPar]);

  parametersIn[5+nDeDxPar] = optionSigma;
  fitAll->SetParName(5+nDeDxPar,"optionSigma");
  fitAll->FixParameter(5+nDeDxPar, parametersIn[5+nDeDxPar]);

  for(Int_t i = 0; i < nSigmaPar; i++) {

    parametersIn[6+nDeDxPar+i] = parSigma[i];
    fitAll->SetParName(6+nDeDxPar+i,Form("sigmapar%d", i));
    Int_t bit = TMath::Nint(TMath::Power(2, i));
    if(fixParametersSigma>=0)
      if (fixParametersSigma==0 || fixParametersSigma&bit) {
	
	cout << "Fixing sigma parameter " << i << endl;
	fitAll->FixParameter(6+nDeDxPar+i, parametersIn[6+nDeDxPar+i]);
      }
    
    if(fixAllParSigma) {
      
      fitAll->FixParameter(6+nDeDxPar+i, parametersIn[6+nDeDxPar+i]);
      }
  }
  
  // Initial yields

  for(Int_t bin = binStart; bin <= binStop; bin++) {
    
    TH1D* hDeDxVsPProj =(TH1D*)hDeDxVsP->ProjectionY("hDeDxVsPProj", bin, bin);
    
    const Int_t offset = nGlobalParams + (bin - binStart)*nLocalParams;
    const Double_t all = hDeDxVsPProj->Integral();

    fitAll->SetParName(offset + 0, Form("piYield%d", bin));
    fitAll->SetParName(offset + 1, Form("kYield%d", bin));
    fitAll->SetParName(offset + 2, Form("pYield%d", bin));
    fitAll->SetParLimits(offset + 0, 0, 10*all);
    fitAll->SetParLimits(offset + 2, 0, 10*all);

    if(fixKtoZero) {
      parametersIn[offset + 0] = (all)*0.5;
      parametersIn[offset + 1] = (all)*0.0;
      fitAll->FixParameter(offset + 1, 0.0);
      parametersIn[offset + 2] = (all)*0.5;
    } else {
      parametersIn[offset + 0] = (all)*0.6;
      parametersIn[offset + 1] = (all)*0.3;
      fitAll->SetParLimits(offset + 1, 0, 10*all);
      parametersIn[offset + 2] = (all)*0.1;
    }    
    // fitAll->SetParLimits(offset + 0, 0., histdEdxvsPpy->GetEntries());
    // fitAll->SetParLimits(offset + 1, 0., histdEdxvsPpy->GetEntries());
    // fitAll->SetParLimits(offset + 2, 0., histdEdxvsPpy->GetEntries());    
  }
  
  fitAll->SetParameters(parametersIn);
  fitAll->Print();

  Bool_t converged = kFALSE;

  if(performFit) {
    for(Int_t i = 0; i < 4; i++) {
      TFitResultPtr fitResultPtr =  hDeDxVsP->Fit(fitAll, "0S", "", pStart+0.01, pStop-0.01);
      if(!fitResultPtr->Status()) {
	//      if(!fitResultPtr->Status() && !fitResultPtr->CovMatrixStatus()) {
	
	converged = kTRUE;
	break;
      }
    }
  }
  // else we just draw how the results looks with the input parameters

  Double_t parametersOut[nParams];
  fitAll->GetParameters(parametersOut);
  const Double_t* parameterErrorsOut = fitAll->GetParErrors();

  // overlay the fitfunction
  

  TF1* fDeDxPi = new TF1("fDeDxPi", FitFunc, 0, 50, nDeDxPar+1); // +1 because of option! 
  fDeDxPi->SetParameters(&parametersOut[3]);
  fDeDxPi->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxPi->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxPi->Draw("SAME");

  TF1* fDeDxK = new TF1("fDeDxK", FitFunc, 0, 50, nDeDxPar+1); 
  fDeDxK->SetParameters(&parametersOut[3]);
  fDeDxK->SetParameter(0, fDeDxK->GetParameter(0)+10);
  fDeDxK->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxK->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxK->Draw("SAME");

  TF1* fDeDxP = new TF1("fDeDxP", FitFunc, 0, 50, nDeDxPar+1); 
  fDeDxP->SetParameters(&parametersOut[3]);
  fDeDxP->SetParameter(0, fDeDxP->GetParameter(0)+20);
  fDeDxP->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxP->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxP->Draw("SAME");

  TF1* fDeDxE = new TF1("fDeDxE", FitFunc, 0, 50, nDeDxPar+1); 
  fDeDxE->SetParameters(&parametersOut[3]);
  fDeDxE->SetParameter(0, fDeDxE->GetParameter(0)+30);
  fDeDxE->SetLineWidth(2);
  cDeDxVsP->cd();
  fDeDxE->Draw("SAME");
  cDeDxVsPLogX->cd();
  fDeDxE->Draw("SAME");

  TF1* fSigma = new TF1("fSigma", SigmaFunc, 0, 50, nSigmaPar+1); 
  fSigma->SetParameters(&parametersOut[5 + nDeDxPar]);

  //  fitAll->Draw("same cont3"); 

  CreateDir(Form("results/%s/calibdedx_%d%d",endingCent[i_cent],etaLow,etaHigh));

  //CreateDir(Form("results/calibdedx_%d%d",etaLow,etaHigh));

  cDeDxVsP->cd();
  cDeDxVsP->Modified();
  cDeDxVsP->Update();
  gROOT->ProcessLine(".x drawText.C");
  cDeDxVsP->SaveAs(Form("results/%s/calibdedx_%d%d/dedx_vs_p.gif",endingCent[i_cent],etaLow,etaHigh));
  cDeDxVsP->SaveAs(Form("results/%s/calibdedx_%d%d/dedx_vs_p.pdf",endingCent[i_cent],etaLow,etaHigh));

  cDeDxVsPLogX->cd();
  cDeDxVsPLogX->Modified();
  cDeDxVsPLogX->Update();
  gROOT->ProcessLine(".x drawText.C");
  cDeDxVsPLogX->SaveAs(Form("results/%s/calibdedx_%d%d/dedx_vs_p_logx.gif",endingCent[i_cent],etaLow,etaHigh));
  cDeDxVsPLogX->SaveAs(Form("results/%s/calibdedx_%d%d/dedx_vs_p_logx.pdf",endingCent[i_cent],etaLow,etaHigh));

  //cross check
  TCanvas* cFits = new TCanvas("cFits", "Fit comparison to data", 1200, 800);

  cFits->Clear();
  cFits->Divide(7, 5);

  TF1* pion = new TF1("pion", "gausn", 30, 90);
  pion->SetLineWidth(2);
  pion->SetLineColor(kRed);
  TF1* kaon = new TF1("kaon", "gausn", 30, 90);
  kaon->SetLineWidth(2);
  kaon->SetLineColor(kGreen);
  TF1* proton = new TF1("proton", "gausn", 30, 90);
  proton->SetLineWidth(2);
  proton->SetLineColor(kBlue);
  TF1* total = new TF1("total", "gausn(0)+gausn(3)+gausn(6)", 30, 90);
  total->SetLineColor(kBlack);
  total->SetLineWidth(2);
  total->SetLineStyle(2);

  TLegend* legend = new TLegend(0.11, 0.6, 0.35, 0.85);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(total, "3-Gauss fit", "L");
  legend->AddEntry(pion, "#pi", "L");
  legend->AddEntry(kaon, "K", "L");
  legend->AddEntry(proton, "p", "L");


  TCanvas* cSingleFit = new TCanvas("cSingleFit", "single fit");
  cSingleFit->SetLogy(1);

  TH1D* hPionRatio =(TH1D*)hDeDxVsP->ProjectionX("hPionRatio", 1, 1);
  hPionRatio->Reset();
  hPionRatio->GetXaxis()->SetRangeUser(pStart+0.001, pStop-0.001);
  hPionRatio->GetYaxis()->SetRangeUser(0.0, 1.0);
  hPionRatio->SetTitle("particle fractions; p [GeV/c]; particle fraction");
  TH1D* hKaonRatio   = (TH1D*)hPionRatio->Clone("hKaonRatio");
  TH1D* hProtonRatio = (TH1D*)hPionRatio->Clone("hProtonRatio");
  
  for(Int_t bin = binStart; bin <= binStop; bin++){
    
    cout << "Making projection for bin: " << bin << endl;
    
    const Int_t j = bin-binStart;
    cFits->cd();
    cFits->cd(j + 1);
    
    TH1D* hDeDxVsPProj =(TH1D*)hDeDxVsP->ProjectionY(Form("hDeDxVsPProj%d", bin), bin, bin);
    //    hDeDxVsPProj->GetXaxis()->SetRangeUser(40, 90);
    hDeDxVsPProj->SetTitle(Form("%.1f<p<%.1f GeV/c", 
				hDeDxVsP->GetXaxis()->GetBinLowEdge(bin),
				hDeDxVsP->GetXaxis()->GetBinUpEdge(bin)));
    hDeDxVsPProj->Draw();
    
    const Int_t offset = nGlobalParams + j*nLocalParams; 
    const Double_t p = hDeDxVsP->GetXaxis()->GetBinCenter(bin);
    const Double_t pKeff = p*piMass/kMass; // corresponding p of a pion with same dE/dx
    const Double_t pPeff = p*piMass/pMass; // corresponding p of a pion with same dE/dx
    const Double_t meanDeDxPi = fDeDxPi->Eval(p);
    const Double_t meanDeDxK  = fDeDxPi->Eval(pKeff);
    const Double_t meanDeDxP  = fDeDxPi->Eval(pPeff);
    Double_t gausParams[9] = { 
      parametersOut[offset + 0], 
      meanDeDxPi, 
      fSigma->Eval(meanDeDxPi) ,
      parametersOut[offset + 1], 
      meanDeDxK, 
      fSigma->Eval(meanDeDxK) ,
      parametersOut[offset + 2], 
      meanDeDxP, 
      fSigma->Eval(meanDeDxP) ,
    };

    for(Int_t i = 0; i < 9; i++) 
      cout << gausParams[i] << ", ";

    cout << endl;
    
    pion->SetParameters(&gausParams[0]);
    pion->DrawCopy("same");
    Double_t all =  hDeDxVsPProj->Integral();
    hPionRatio->SetBinContent(bin, parametersOut[offset + 0]/all);
    hPionRatio->SetBinError(bin, parameterErrorsOut[offset + 0]/all);

    kaon->SetParameters(&gausParams[3]);
    kaon->DrawCopy("same");
    hKaonRatio->SetBinContent(bin, parametersOut[offset + 1]/all);
    hKaonRatio->SetBinError(bin, parameterErrorsOut[offset + 1]/all);
    
    proton->SetParameters(&gausParams[6]);
    proton->DrawCopy("same");
    hProtonRatio->SetBinContent(bin, parametersOut[offset + 2]/all);
    hProtonRatio->SetBinError(bin, parameterErrorsOut[offset + 2]/all);
    
    total->SetParameters(gausParams);
    total->DrawCopy("same");

    cSingleFit->cd();
    cSingleFit->Clear();
    //    cSingleFit->SetLogy();
    hDeDxVsPProj->Draw();
    pion->DrawCopy("same");
    kaon->DrawCopy("same");
    proton->DrawCopy("same");
    total->DrawCopy("same");
    
    gROOT->ProcessLine(".x drawText.C(2)");


    cSingleFit->SaveAs(Form("results/%s/calibdedx_%d%d/ptspectrum_bin%d.gif",endingCent[i_cent], etaLow,etaHigh, bin));
    cSingleFit->SaveAs(Form("results/%s/calibdedx_%d%d/ptspectrum_bin%d.pdf",endingCent[i_cent],etaLow, etaHigh, bin));
    //    legend->Draw();


  }

  TCanvas* cRatio = new TCanvas("cRatio", "ratios/all vs p", 600, 400);
  cRatio->Clear();
  hPionRatio->SetMaximum(0.8);
  hPionRatio->SetMarkerStyle(20);
  hPionRatio->SetMarkerColor(2);
  hPionRatio->Draw("P E");

  hKaonRatio->SetMarkerStyle(20);
  hKaonRatio->SetMarkerColor(3);
  hKaonRatio->Draw("SAME P E");

  hProtonRatio->SetMarkerStyle(20);
  hProtonRatio->SetMarkerColor(4);
  hProtonRatio->Draw("SAME P E");
  gROOT->ProcessLine(".x drawText.C(2)");
  cRatio->SaveAs(Form("results/%s/calibdedx_%d%d/particle_ratios.gif",endingCent[i_cent],etaLow,etaHigh));
  cRatio->SaveAs(Form("results/%s/calibdedx_%d%d/particle_ratios.pdf",endingCent[i_cent],etaLow,etaHigh));



  //
  // Store the <dE/dx> parameters in a file that we can get them back to use for the Delta-pi!
  //
  DeDxFitInfo* fitInfo = new DeDxFitInfo();
  fitInfo->MIP     = fixMIP;
  //fitInfo->plateau = 80; 
  fitInfo->optionDeDx = optionDeDx; 
  fitInfo->nDeDxPar = nDeDxPar; 
  for(Int_t i = 0; i < nDeDxPar; i++) {
    fitInfo->parDeDx[i] = fDeDxPi->GetParameter(i+1); // 1st parameter is option
  }
  //fitInfo->plateau = fDeDxPi->GetParameter(5);//lo puse 16/09/13
  fitInfo->optionSigma = optionSigma; 
  fitInfo->nSigmaPar = nSigmaPar; 
  for(Int_t i = 0; i < nSigmaPar; i++) {
    cout<<"my sugma param="<<fSigma->GetParameter(i+1)<<endl;
    fitInfo->parSigma[i] = fSigma->GetParameter(i+1); // 1st parameter is option 
  }


  fitInfo->Print();


  if(converged) {

    cout << "Fit converged and error matrix was ok" << endl;
  } else {

    cout << "WARNING!!!!!!!!!!!!!!! Fit did not converge, or error matrix was  not ok!" << endl;
  }


  CreateDir(Form("fitparameters/%s",endingCent[i_cent]));
  //cout<<"flag 0"<<"    gSystem->BaseName(calibFileName))="<< gSystem->BaseName(calibFileName)<<endl;
  TFile* outFile = new TFile(Form("fitparameters/%s/%d%d_dataPbPb.root", endingCent[i_cent], etaLow, etaHigh), "RECREATE");
  outFile->cd();
  fitInfo->Write("fitInfo");
  cout<<"flag 1"<<endl;

  outFile->Close();

  if(converged) {

    cout << "Fit converged and error matrix was ok" << endl;
  } else {

    cout << "WARNING!!!!!!!!!!!!!!! Fit did not converge, or error matrix was  not ok!" << endl;
  }
}





void MakeFitsV0s(const Char_t * fileInputName, const Char_t * fileE, const Char_t *outDir, Int_t index_eta){
  /*
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  */
  TFile *f12=TFile::Open(fileE);
  TF1 *fPlateau = (TF1 *)f12->Get(Form("fGauss%d",index_eta));
  /*
  TH1D * hEV0s=(TH1D *)f12->Get("hdedx_vs_p_e");
  TGraphErrors	*gResE=(TGraphErrors	*)f12->Get("gRes_Electrons");
  */

  const Char_t * endNamePos[4]={"02","24","46","68"};
  const Char_t * endNameNeg[4]={"20","42","64","86"};



  TFile *f1=TFile::Open(fileInputName);
  TList * list = (TList *)f1->Get("outputdedx");


  //positive eta
  TH2D *hPi2DPos=(TH2D *)list->FindObject(Form("histPiTof%s",endNamePos[index_eta]));
  TH2D *hPi2DV0sPos=(TH2D *)list->FindObject(Form("histPiV0%s",endNamePos[index_eta]));
  TH2D *hP2DPos=(TH2D *)list->FindObject(Form("histPV0%s",endNamePos[index_eta]));
  TH1D *hPi2DpPos=(TH1D *)list->FindObject(Form("histPiTof1D%s",endNamePos[index_eta]));
  TH1D *hPi2DV0spPos=(TH1D *)list->FindObject(Form("histPiV01D%s",endNamePos[index_eta]));
  TH1D *hP2DpPos=(TH1D *)list->FindObject(Form("histPV01D%s",endNamePos[index_eta]));


  gSystem->Exec("rm *.root");

  TFile *fpos = new TFile("histtmppos.root","recreate");
  fpos->cd();
  hPi2DPos->SetName(Form("histPiTof%s",endNamePos[index_eta]));
  hPi2DV0sPos->SetName(Form("histPiV0%s",endNamePos[index_eta]));
  hP2DPos->SetName(Form("histPV0%s",endNamePos[index_eta]));
  hPi2DpPos->SetName(Form("histPiTof1D%s",endNamePos[index_eta]));
  hPi2DV0spPos->SetName(Form("histPiV01D%s",endNamePos[index_eta]));
  hP2DpPos->SetName(Form("histPV01D%s",endNamePos[index_eta]));
  hPi2DPos->Write();
  hPi2DV0sPos->Write();
  hP2DPos->Write();
  hPi2DpPos->Write();
  hPi2DV0spPos->Write();
  hP2DpPos->Write();
  fpos->Close();


  //negative eta
  TH2D *hPi2DNeg=(TH2D *)list->FindObject(Form("histPiTof%s",endNameNeg[index_eta]));
  TH2D *hPi2DV0sNeg=(TH2D *)list->FindObject(Form("histPiV0%s",endNameNeg[index_eta]));
  TH2D *hP2DNeg=(TH2D *)list->FindObject(Form("histPV0%s",endNameNeg[index_eta]));
  TH1D *hPi2DpNeg=(TH1D *)list->FindObject(Form("histPiTof1D%s",endNameNeg[index_eta]));
  TH1D *hPi2DV0spNeg=(TH1D *)list->FindObject(Form("histPiV01D%s",endNameNeg[index_eta]));
  TH1D *hP2DpNeg=(TH1D *)list->FindObject(Form("histPV01D%s",endNameNeg[index_eta]));



  TFile *fneg = new TFile("histtmpneg.root","recreate");
  fneg->cd();
  hPi2DNeg->SetName(Form("histPiTof%s",endNamePos[index_eta]));
  hPi2DV0sNeg->SetName(Form("histPiV0%s",endNamePos[index_eta]));
  hP2DNeg->SetName(Form("histPV0%s",endNamePos[index_eta]));
  hPi2DpNeg->SetName(Form("histPiTof1D%s",endNamePos[index_eta]));
  hPi2DV0spNeg->SetName(Form("histPiV01D%s",endNamePos[index_eta]));
  hP2DpNeg->SetName(Form("histPV01D%s",endNamePos[index_eta]));
  hPi2DNeg->Write();
  hPi2DV0sNeg->Write();
  hP2DNeg->Write();
  hPi2DpNeg->Write();
  hPi2DV0spNeg->Write();
  hP2DpNeg->Write();
  fneg->Close();

  //pos+neg eta
  TH2D *hPi2D=0;
  TH2D *hPi2DV0s=0;
  TH2D *hP2D=0;
  TH1D *hPi2Dp=0;
  TH1D *hPi2DV0sp=0;
  TH1D *hP2Dp=0;


  gSystem->Exec(Form("hadd tmp%d.root histtmppos.root histtmpneg.root",index_eta));

  TFile *ffinalv0 = TFile::Open(Form("tmp%d.root",index_eta));
  hPi2D=(TH2D *)ffinalv0->Get(Form("histPiTof%s",endNamePos[index_eta]));
  hPi2DV0s=(TH2D *)ffinalv0->Get(Form("histPiV0%s",endNamePos[index_eta]));
  hP2D=(TH2D *)ffinalv0->Get(Form("histPV0%s",endNamePos[index_eta]));
  hPi2Dp=(TH1D *)ffinalv0->Get(Form("histPiTof1D%s",endNamePos[index_eta]));
  hPi2DV0sp=(TH1D *)ffinalv0->Get(Form("histPiV01D%s",endNamePos[index_eta]));
  hP2Dp=(TH1D *)ffinalv0->Get(Form("histPV01D%s",endNamePos[index_eta]));



  TList *l1=new TList();
  l1->SetOwner();


  TGraphErrors* graphRes = new TGraphErrors();
  Int_t nERes = 0;

  TGraphErrors* graphBB = new TGraphErrors();
  Int_t nBB = 0;

  TGraphErrors* graphBBpitof = new TGraphErrors();
  Int_t nBBpitof = 0;
  TGraphErrors* graphBBpiv0s = new TGraphErrors();
  Int_t nBBpiv0s = 0;
  TGraphErrors* graphBBpv0s = new TGraphErrors();
  Int_t nBBpv0s = 0;


 

  TF1* fFit1D = new TF1("fFit1D", "gausn(0)", 40, 90);
 


  //TF1* pionsv0s = new TF1("pionsv0s", "gausn", 56, 80);
  TF1* pionsv0s = new TF1("pionsv0s", "gausn", 40, 80);
  pionsv0s->SetLineColor(kRed);
  pionsv0s->SetLineWidth(2);
  pionsv0s->SetLineStyle(2);

  //TF1* protonsv0s = new TF1("protonsv0s", "gausn", 40, 55);
  TF1* protonsv0s = new TF1("protonsv0s", "gausn", 40, 80);
  protonsv0s->SetLineColor(kBlue);
  protonsv0s->SetLineWidth(2);
  protonsv0s->SetLineStyle(2);


  TF1* total = new TF1("total", "gausn(0)+gausn(3)", 40, 80);
  total->SetLineColor(1);
  total->SetLineWidth(2);
  total->SetLineStyle(2);



  Int_t nbins=hPi2D->GetXaxis()->GetNbins();



  for(Int_t bin=1;bin<=nbins;++bin){
    
 

    TH1D *hproy_pion=0;
    TH1D *hproy_proton=0;

    Double_t p=hPi2D->GetXaxis()->GetBinCenter(bin);
    Double_t lowp = hPi2D->GetXaxis()->GetBinLowEdge(bin);
    Double_t widthp = hPi2D->GetXaxis()->GetBinWidth(bin);
    hPi2Dp->GetXaxis()->SetRangeUser(lowp,lowp+widthp);
    Double_t meanp= hPi2Dp->GetMean();
    Double_t meanpE= hPi2Dp->GetMeanError();


    cout<<"p="<<p<<endl;

    if(p>2.0&&p<9.0){
      
      cout<<bin<<endl;

      hproy_pion=(TH1D *)hPi2D->ProjectionY(Form("hproy_pionTof_%d",bin),bin,bin);

      cout<<bin<<endl;

      hproy_pion->SetTitle(Form("pionsv0s: %.1f < p < %.1f GeV/c; dE/dx (GeV/c); Entries",
				hPi2D->GetXaxis()->GetBinLowEdge(bin),
				hPi2D->GetXaxis()->GetBinUpEdge(bin)));
      cout<<bin<<endl;
      hproy_proton=(TH1D *)hP2D->ProjectionY(Form("hproy_proton_%d",bin),bin,bin);

      cout<<"nbinsproj="<<hproy_proton->GetNbinsX()<<endl;
      //hproy_pion->Draw();


      TSpectrum *s = new TSpectrum(2);
      Int_t nfound = s->Search(hproy_pion,2,"",0.15);
      cout<<"!!!!!!!!!!!!!!         nfound="<<nfound<<endl;
      Int_t npeaks = 0;
      Float_t *xpeaks = s->GetPositionX();


      for (Int_t p1=0;p1<nfound;p1++) {
	//Float_t xp = xpeaks[p1];
	npeaks++;
      }

      if(npeaks==0)
	continue;
      if(npeaks==1){

	fFit1D->SetParameters(hproy_pion->GetEntries(), xpeaks[0], hproy_pion->GetRMS());
	hproy_pion->Fit(fFit1D, "Q");
	hproy_pion->Fit(fFit1D, "Q", "", xpeaks[0]-hproy_pion->GetRMS(),
			xpeaks[0]+hproy_pion->GetRMS());
	
	hproy_pion->Fit(fFit1D, "R", "", fFit1D->GetParameter(1)-1.0*fFit1D->GetParameter(2),
			fFit1D->GetParameter(1)+2*fFit1D->GetParameter(2));




	  graphBB->SetPoint(nBB,meanp/0.14, fFit1D->GetParameter(1));
	  graphBB->SetPointError(nBB,meanpE, fFit1D->GetParError(1));
	  nBB++;
	  
	  graphBBpitof->SetPoint(nBBpitof,meanp/0.14, fFit1D->GetParameter(1));
	  graphBBpitof->SetPointError(nBBpitof,meanpE, fFit1D->GetParError(1));
	  nBBpitof++;
	  
	  if(p>3&&p<9){
	    graphRes->SetPoint(nERes,fFit1D->GetParameter(1),fFit1D->GetParameter(2)/fFit1D->GetParameter(1));
	    graphRes->SetPointError(nERes,fFit1D->GetParError(1),fFit1D->GetParError(2)/fFit1D->GetParameter(1));
	    nERes++;
	  }
	  



      }
      else if(npeaks==2){
	Int_t bin1 = hproy_pion->GetXaxis()->FindBin(xpeaks[0]);
	Int_t bin2 = hproy_pion->GetXaxis()->FindBin(xpeaks[1]);

	Int_t bin1C = hproy_pion->GetBinContent(bin1);
	Int_t bin2C = hproy_pion->GetBinContent(bin2);

	cout<<"bin1="<<bin1<<"  content="<<bin1C<<endl;
	cout<<"bin2="<<bin2<<"  content="<<bin2C<<endl;

	
	if(bin==24){
	  
	  hproy_pion->Fit(fFit1D, "R", "", 65, 85);
	  graphBB->SetPoint(nBB,meanp/0.14, fFit1D->GetParameter(1));
	  graphBB->SetPointError(nBB,meanpE, fFit1D->GetParError(1));
	  nBB++;
	  
	  graphBBpitof->SetPoint(nBBpitof,meanp/0.14, fFit1D->GetParameter(1));
	  graphBBpitof->SetPointError(nBBpitof,meanpE, fFit1D->GetParError(1));
	  nBBpitof++;
	    
	  
	  
	}
	

	else if(bin1C>bin2C){//pion signal larger
	  
	  
	  fFit1D->SetParameters(hproy_pion->GetEntries(), xpeaks[0], hproy_pion->GetRMS());
	  hproy_pion->Fit(fFit1D, "Q");
	  hproy_pion->Fit(fFit1D, "Q", "", xpeaks[0]-hproy_pion->GetRMS(),
			  xpeaks[0]+hproy_pion->GetRMS());
	  
	  
	  
	  
	  hproy_pion->Fit(fFit1D, "R", "", fFit1D->GetParameter(1)-1.0*fFit1D->GetParameter(2),
	  		  fFit1D->GetParameter(1)+2.5*fFit1D->GetParameter(2));

	  
	  
	  
	  cout<<"hello"<<endl;
	  //protonsv0s->SetParameters(hproy_proton->GetEntries(), xpeaks[1], hproy_proton->GetRMS());
	  
	  
	  hproy_proton->Fit(protonsv0s, "Q");

	  
	  hproy_proton->Fit(protonsv0s, "Q", "", xpeaks[1]-hproy_proton->GetRMS(),
			    xpeaks[1]+hproy_proton->GetRMS());
	  
	  
	  hproy_proton->Fit(protonsv0s, "R", "", xpeaks[1]-2*protonsv0s->GetParameter(2),
			  xpeaks[1]+2*protonsv0s->GetParameter(2));
	  
	  protonsv0s->SetParameters(hproy_proton->GetEntries(), xpeaks[1], protonsv0s->GetParameter(2));
	  hproy_proton->Fit(protonsv0s, "R", "", xpeaks[1]-2.0*protonsv0s->GetParameter(2),
			    xpeaks[1]+2.0*protonsv0s->GetParameter(2));
	  
	  total->SetParameter(1,protonsv0s->GetParameter(1));
	  total->SetParameter(2,protonsv0s->GetParameter(2));
	  total->SetParameter(4,fFit1D->GetParameter(1));
	  total->SetParameter(5,fFit1D->GetParameter(2));
	  hproy_pion->Fit(total,"R");

	  Float_t reducedCh2=total->GetChisquare()/total->GetNDF();
	  if(reducedCh2>1.0){
	    hproy_pion->Fit(fFit1D, "R", "", fFit1D->GetParameter(1)-2*fFit1D->GetParameter(2),
			    fFit1D->GetParameter(1)+2*fFit1D->GetParameter(2));
	    //hproy_pion->Fit(fFit1D, "R", "", fFit1D->GetParameter(1)-0.5*fFit1D->GetParameter(2),
	    hproy_pion->Fit(fFit1D, "R", "", fFit1D->GetParameter(1)-1.0*fFit1D->GetParameter(2),
			    fFit1D->GetParameter(1)+2*fFit1D->GetParameter(2));


	    graphBB->SetPoint(nBB,meanp/0.14, fFit1D->GetParameter(1));
	    graphBB->SetPointError(nBB,meanpE, fFit1D->GetParError(1));
	    nBB++;
	    

	    graphBBpitof->SetPoint(nBBpitof,meanp/0.14, fFit1D->GetParameter(1));
	    graphBBpitof->SetPointError(nBBpitof,meanpE, fFit1D->GetParError(1));
	    nBBpitof++;
	    
	    if(p>3.0&&p<9){
	      graphRes->SetPoint(nERes,fFit1D->GetParameter(1),fFit1D->GetParameter(2)/fFit1D->GetParameter(1));
	      graphRes->SetPointError(nERes,fFit1D->GetParError(1),fFit1D->GetParError(2)/fFit1D->GetParameter(1));
	      nERes++;
	    }
	    
	    

	  }else{
	    
	    graphBB->SetPoint(nBB,meanp/0.14, total->GetParameter(4));
	    graphBB->SetPointError(nBB,meanpE, total->GetParError(4));
	    nBB++;
	    
	    graphBBpitof->SetPoint(nBBpitof,meanp/0.14, total->GetParameter(4));
	    graphBBpitof->SetPointError(nBBpitof,meanpE, total->GetParError(4));
	    nBBpitof++;

	    if(p>3.0&&p<9){
	      graphRes->SetPoint(nERes,total->GetParameter(4),total->GetParameter(5)/total->GetParameter(4));
	      
	      graphRes->SetPointError(nERes,total->GetParError(4),total->GetParError(5)/total->GetParameter(4));
	      nERes++;
	    }
	    
	    
	    
	    
	  }
	  
	  
	  
	  
	}else{
	  
	  fFit1D->SetParameters(hproy_pion->GetEntries(), xpeaks[1], hproy_pion->GetRMS());
	  hproy_pion->Fit(fFit1D, "Q");
	  hproy_pion->Fit(fFit1D, "Q", "", xpeaks[1]-hproy_pion->GetRMS(),
			  xpeaks[0]+hproy_pion->GetRMS());
	  
	  hproy_pion->Fit(fFit1D, "R", "", fFit1D->GetParameter(1)-2*fFit1D->GetParameter(2),
			  fFit1D->GetParameter(1)+2*fFit1D->GetParameter(2));
	  

	  protonsv0s->SetParameters(hproy_proton->GetEntries(), xpeaks[0], hproy_proton->GetRMS());
	  hproy_proton->Fit(protonsv0s, "Q");
	  hproy_proton->Fit(protonsv0s, "Q", "", xpeaks[0]-hproy_proton->GetRMS(),
			  xpeaks[0]+hproy_proton->GetRMS());
	  
	  hproy_proton->Fit(protonsv0s, "R", "", xpeaks[0]-2*protonsv0s->GetParameter(2),
			  xpeaks[0]+2*protonsv0s->GetParameter(2));

	  protonsv0s->SetParameters(hproy_proton->GetEntries(), xpeaks[0], protonsv0s->GetParameter(2));
	  hproy_proton->Fit(protonsv0s, "R", "", xpeaks[0]-2*protonsv0s->GetParameter(2),
			  xpeaks[0]+2*protonsv0s->GetParameter(2));


	  //graphBB->SetPoint(nBB,p/0.14, xpeaks[1]);
	  //graphBB->SetPointError(nBB,0, fFit1D->GetParError(1));
	  //nBB++;
	  total->SetParameter(1,protonsv0s->GetParameter(1));
	  total->SetParameter(2,protonsv0s->GetParameter(2));
	  total->SetParameter(4,fFit1D->GetParameter(1));
	  total->SetParameter(5,fFit1D->GetParameter(2));
	  hproy_pion->Fit(total,"R");

	  //graphBB->SetPoint(nBB,p/0.14, xpeaks[1]);
	  //graphBB->SetPointError(nBB,0, total->GetParError(4));
	  graphBB->SetPoint(nBB,meanp/0.14, fFit1D->GetParameter(1));
	  graphBB->SetPointError(nBB,meanpE, fFit1D->GetParError(1));
	  nBB++;

	  graphBBpitof->SetPoint(nBBpitof,meanp/0.14, fFit1D->GetParameter(1));
	  graphBBpitof->SetPointError(nBBpitof,meanpE, fFit1D->GetParError(1));
	  nBBpitof++;
	  

	  if(p>3.0&&p<9){
	    graphRes->SetPoint(nERes,fFit1D->GetParameter(1),fFit1D->GetParameter(2)/fFit1D->GetParameter(1));
	    graphRes->SetPointError(nERes,fFit1D->GetParError(1),fFit1D->GetParError(2)/fFit1D->GetParameter(1));
	    nERes++;
	  }



	}



    







      }
      
      l1->Add(hproy_pion);
      
      
    }
    else
      continue;
    
  }



  graphRes->SetTitle(";#LT dE/dx  #GT;#sigma/#LT dE/dx  #GT");
  graphRes->SetName("gRes_PionsTof");
  
  graphBBpitof->SetTitle(";#beta*#gamma; #LT dE/dx  #GT");
  graphBBpitof->SetName("BB_piTOF");
  graphBBpitof->SetMarkerColor(kGreen+2);
  graphBBpitof->SetLineColor(kGreen+2);
  graphBBpitof->SetMarkerStyle(21);

  l1->Add(graphRes);
  l1->Add(graphBBpitof);

  

  //pi by V0s

  TGraphErrors* graphResV0s = new TGraphErrors();
  nERes = 0;

  nbins=hPi2DV0s->GetXaxis()->GetNbins();



  for(Int_t bin=1;bin<=nbins;++bin){
    
  

    TH1D *hproy_pionV0s=0;
    TH1D *hproy_proton=0;
 
    Double_t p=hPi2DV0s->GetXaxis()->GetBinCenter(bin);

    Double_t lowp = hPi2DV0s->GetXaxis()->GetBinLowEdge(bin);
    Double_t widthp = hPi2DV0s->GetXaxis()->GetBinWidth(bin);
    hPi2DV0sp->GetXaxis()->SetRangeUser(lowp,lowp+widthp);
    Double_t meanp= hPi2DV0sp->GetMean();
    Double_t meanpE= hPi2DV0sp->GetMeanError();



    if(p>2.0&&p<9){

      hproy_pionV0s=(TH1D *)hPi2DV0s->ProjectionY(Form("hproy_pionV0s_%d",bin),bin,bin);
      hproy_pionV0s->SetTitle(Form("pionsv0s: %.1f < p < %.1f GeV/c; dE/dx (GeV/c); Entries",
				hPi2DV0s->GetXaxis()->GetBinLowEdge(bin),
				hPi2DV0s->GetXaxis()->GetBinUpEdge(bin)));

     fFit1D->SetParameters(hproy_pionV0s->GetEntries(), hproy_pionV0s->GetMean(), hproy_pionV0s->GetRMS());
      hproy_pionV0s->Fit(fFit1D, "Q");
      hproy_pionV0s->Fit(fFit1D, "Q", "", hproy_pionV0s->GetMean()-hproy_pionV0s->GetRMS(),
			 hproy_pionV0s->GetMean()+hproy_pionV0s->GetRMS());
      
      
      //      hproy_pionV0s->Fit(fFit1D, "Q", "", fFit1D->GetParameter(1)-1.5*fFit1D->GetParameter(2),
      //	      fFit1D->GetParameter(1)+3.0*fFit1D->GetParameter(2));

      hproy_pionV0s->Fit(fFit1D, "R", "", fFit1D->GetParameter(1)-1.0*fFit1D->GetParameter(2),
			 fFit1D->GetParameter(1)+2.0*fFit1D->GetParameter(2));
      
      
      hproy_proton=(TH1D *)hP2D->ProjectionY(Form("hproy_proton_%d",bin),bin,bin);
      hproy_proton->SetTitle(Form("protonsv0s: %.1f < p < %.1f GeV/c; dE/dx; Entries",
				  hPi2D->GetXaxis()->GetBinLowEdge(bin),
				  hPi2D->GetXaxis()->GetBinUpEdge(bin)));
      
      //protonsv0s->SetParameters(hproy_proton->GetEntries(), hproy_proton->GetMean(), hproy_proton->GetRMS());
      hproy_proton->Fit(protonsv0s,"R");
      //hproy_proton->Fit(protonsv0s,"R", "",protonsv0s->GetParameter(1)-3.0*protonsv0s->GetParameter(2),protonsv0s->GetParameter(1)+1.0*protonsv0s->GetParameter(2));
      hproy_proton->Fit(protonsv0s,"R","",protonsv0s->GetParameter(1)-2.0*protonsv0s->GetParameter(2),protonsv0s->GetParameter(1)+0.8*protonsv0s->GetParameter(2));

      total->SetParameter(0,0.1*hproy_pionV0s->GetEntries());
      total->SetParameter(1,protonsv0s->GetParameter(1));
      total->SetParameter(2,protonsv0s->GetParameter(2));
      total->SetParameter(3,0.9*hproy_pionV0s->GetEntries());
      total->SetParameter(4,fFit1D->GetParameter(1));
      total->SetParameter(5,fFit1D->GetParameter(2));
      
      hproy_pionV0s->Fit(total,"0R");
      
      cout<<"p="<<p<<"   yield0="<<total->GetParameter(0)<<"  yield1="<<total->GetParameter(3)<<endl;
      cout<<"p="<<p<<"   mean0="<<total->GetParameter(1)<<"  mean1="<<total->GetParameter(4)<<endl;
      cout<<"p="<<p<<"   sigma0="<<total->GetParameter(1)<<"  sigma1="<<total->GetParameter(5)<<endl;
      
      total->SetParameter(1,total->GetParameter(1));
      total->SetParameter(2,total->GetParameter(2));
      total->SetParameter(4,total->GetParameter(4));
      total->SetParameter(5,total->GetParameter(5));
      
      hproy_pionV0s->Fit(total,"R");

 
      
      if(p>4.0){


	//return;

	graphResV0s->SetPoint(nERes,total->GetParameter(4),total->GetParameter(5)/total->GetParameter(4));
	graphResV0s->SetPointError(nERes,total->GetParError(4),total->GetParError(5)/total->GetParameter(4));
	nERes++;
	
	//graphBB->SetPoint(nBB,p/0.14, total->GetParameter(4));
	//graphBB->SetPointError(nBB,0, total->GetParError(4));
	//nBB++;
	
      
	graphBBpiv0s->SetPoint(nBBpiv0s,meanp/0.14, total->GetParameter(4));
	graphBBpiv0s->SetPointError(nBBpiv0s,meanpE,total->GetParError(4));
	nBBpiv0s++;
	
	l1->Add(hproy_pionV0s);
	
      }
      
      
      
    }
    
  }
  
  
  graphResV0s->SetTitle(";#LT dE/dx  #GT;#sigma/#LT dE/dx  #GT");
  graphResV0s->SetName("gRes_PionsV0sArmenteros");
  l1->Add(graphResV0s);

  graphBBpiv0s->SetTitle(";#beta*#gamma; #LT dE/dx  #GT");
  graphBBpiv0s->SetName("BB_piV0s");
  graphBBpiv0s->SetMarkerColor(2);
  graphBBpiv0s->SetLineColor(2);
  graphBBpiv0s->SetMarkerStyle(24);
  l1->Add(graphBBpiv0s);
  
  //protons
  
  

  
  TGraphErrors* graphPRes = new TGraphErrors();
  nERes = 0;
  TGraphErrors* graphdEdxvsP = new TGraphErrors();
  TGraphErrors* graphSigmavsP = new TGraphErrors();
  graphdEdxvsP->SetName("graphdEdxvsP_p");
  graphSigmavsP->SetName("graphSigmavsP_p");

  nbins=hPi2D->GetXaxis()->GetNbins();
  for(Int_t bin=1;bin<=nbins;++bin){
    
    TH1D *hproy_pion=0;
    TH1D *hproy_proton=0;

    
    Double_t p=hPi2D->GetXaxis()->GetBinCenter(bin);
    Double_t lowp = hP2D->GetXaxis()->GetBinLowEdge(bin);
    Double_t widthp = hP2D->GetXaxis()->GetBinWidth(bin);
    hP2Dp->GetXaxis()->SetRangeUser(lowp,lowp+widthp);
    Double_t meanp= hP2Dp->GetMean();
    Double_t meanpE= hP2Dp->GetMeanError();


    if(p<2.0||p>5.0)
      continue;


    if(p>2){

      cout<<"p="<<hPi2D->GetXaxis()->GetBinCenter(bin)<<endl;

      hproy_pion=(TH1D *)hPi2D->ProjectionY(Form("hproy_pion_%d",bin),bin,bin);
      hproy_pion->SetTitle(Form("pionsv0s: %.1f < p < %.1f GeV/c; dE/dx (GeV/c); Entries",
				hPi2D->GetXaxis()->GetBinLowEdge(bin),
				hPi2D->GetXaxis()->GetBinUpEdge(bin)));
      
      hproy_pion->Fit(pionsv0s,"R");
      
      hproy_proton=(TH1D *)hP2D->ProjectionY(Form("hproy_proton_%d",bin),bin,bin);
      hproy_proton->SetTitle(Form("protonsv0s: %.1f < p < %.1f GeV/c; dE/dx (GeV/c); Entries",
				  hPi2D->GetXaxis()->GetBinLowEdge(bin),
				  hPi2D->GetXaxis()->GetBinUpEdge(bin)));
      


      hproy_proton->Fit(protonsv0s,"R");

      hproy_proton->Fit(protonsv0s,"R","",protonsv0s->GetParameter(1)-2.0*protonsv0s->GetParameter(2),protonsv0s->GetParameter(1)+0.8*protonsv0s->GetParameter(2));


      Float_t reducedCh2=protonsv0s->GetChisquare()/protonsv0s->GetNDF();
      //if(reducedCh2<3)
      if(reducedCh2>3){
	total->SetParameter(1,protonsv0s->GetParameter(1));
	total->SetParameter(2,protonsv0s->GetParameter(2));
	total->SetParameter(4,pionsv0s->GetParameter(1));
	total->SetParameter(5,pionsv0s->GetParameter(2));
	
	hproy_proton->Fit(total,"R");
      }
      //if(total->GetParameter(1)<=0)continue;	

      if(p>2&&p<7){
	if(reducedCh2<5){



	  if(p>3.0&&p<5){
	    graphPRes->SetPoint(nERes,protonsv0s->GetParameter(1),protonsv0s->GetParameter(2)/protonsv0s->GetParameter(1));
	    graphPRes->SetPointError(nERes,protonsv0s->GetParError(1),protonsv0s->GetParError(2)/protonsv0s->GetParameter(1));
	    nERes++;
	  }
	  graphBB->SetPoint(nBB,meanp/0.943, protonsv0s->GetParameter(1));
	  graphBB->SetPointError(nBB,meanpE,protonsv0s->GetParError(1));

	  graphBBpv0s->SetPoint(nBBpv0s,meanp/0.943, protonsv0s->GetParameter(1));
	  graphBBpv0s->SetPointError(nBBpv0s,meanpE,protonsv0s->GetParError(1));
	  nBBpv0s++;

	  graphdEdxvsP->SetPoint(nBB,p, protonsv0s->GetParameter(1));
	  graphdEdxvsP->SetPointError(nBB,0.001, protonsv0s->GetParError(1));
	  graphSigmavsP->SetPoint(nBB,p, protonsv0s->GetParameter(2));
	  graphSigmavsP->SetPointError(nBB,0.0001, protonsv0s->GetParError(2));
	  nBB++;

	}
	else{

	  if(p>3.0&&p<5){
	    graphPRes->SetPoint(nERes,total->GetParameter(1),total->GetParameter(2)/total->GetParameter(1));
	    graphPRes->SetPointError(nERes,total->GetParError(1),total->GetParError(2)/total->GetParameter(1));
	    nERes++;
	  }

	  //graphBB->SetPoint(nBB,p/0.943, total->GetParameter(1));
	  //graphBB->SetPointError(nBB,0,total->GetParError(1));

	  graphBBpv0s->SetPoint(nBBpv0s,meanp/0.943, protonsv0s->GetParameter(1));
	  graphBBpv0s->SetPointError(nBBpv0s,meanpE,protonsv0s->GetParError(1));
	  nBBpv0s++;

	  
	  graphdEdxvsP->SetPoint(nBB,p, total->GetParameter(1));
	  graphdEdxvsP->SetPointError(nBB,0.001, total->GetParError(1));
	  graphSigmavsP->SetPoint(nBB,p, total->GetParameter(2));
	  graphSigmavsP->SetPointError(nBB,0.0001, total->GetParError(2));

	  graphBB->SetPoint(nBB,meanp/0.943, protonsv0s->GetParameter(1));
	  graphBB->SetPointError(nBB,meanpE,protonsv0s->GetParError(1));

	  nBB++;
	}

	
      }
      /*else{
	
	graphPRes->SetPoint(nERes,protonsv0s->GetParameter(1),protonsv0s->GetParameter(2)/protonsv0s->GetParameter(1));
	graphPRes->SetPointError(nERes,protonsv0s->GetParError(1),protonsv0s->GetParError(2)/protonsv0s->GetParameter(1));
	
	
	nERes++;
	
	}*/

      l1->Add(hproy_proton);


    }
    
 
 
  }

  graphPRes->SetTitle(";#LT dE/dx  #GT;#sigma/#LT dE/dx  #GT");
  graphPRes->SetName("gRes_Protons");

  l1->Add(graphPRes);

  graphBBpv0s->SetTitle(";#beta*#gamma; #LT dE/dx  #GT");
  graphBBpv0s->SetName("BB_pv0s");
  graphBBpv0s->SetMarkerColor(4);
  graphBBpv0s->SetLineColor(4);
  graphBBpv0s->SetMarkerStyle(24);
  l1->Add(graphBBpv0s);
  
  

  graphBB->SetTitle(";#beta*#gamma; #LT dE/dx  #GT");
  graphBB->SetName("gBB");
 

  Double_t dedxPar[10]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  dedxPar[0] =  0;
  dedxPar[1] =  31.384;
  dedxPar[2] =  9.5253;
  dedxPar[3] =  7.3542;
  dedxPar[4] =  1.5;
  
  TF1 *dedxFunc = new TF1("dedxFunc", FitFunc, 2.5, 60, 6);
  dedxFunc->SetParameters(dedxPar);
  // Fit betagamma (option+40)
  dedxFunc->SetParameter(0, 6+40);
  dedxFunc->FixParameter(0, 6+40);
  //parameters for MB
  dedxFunc->SetParameter(1, 33.7538);
  dedxFunc->SetParameter(2, 8.05926);
  dedxFunc->SetParameter(3, 1.69954);
  dedxFunc->SetParameter(4, 1.37793);
  dedxFunc->SetParameter(5, 78.0465);

 
  cout << "!!!!!!!!!!!!!!!!!                    Setting new plateau=" << fPlateau->GetParameter(1) <<  endl;
  dedxFunc->SetParameter(5, fPlateau->GetParameter(1));
  dedxFunc->FixParameter(5, fPlateau->GetParameter(1));
 

  graphBB->Fit(dedxFunc, "0R", "SAME");
  graphBB->Fit(dedxFunc, "0R", "SAME");
  graphBB->Fit(dedxFunc, "0R", "SAME");

  graphBB->Draw();

  for(Int_t i=0;i<dedxFunc->GetNpar();++i)
    dedxFunc->SetParameter(i,dedxFunc->GetParameter(i));
  
  l1->Add(dedxFunc);
  l1->Add(graphBB);



  //Fit to get resolution
  TGraphErrors *g1=0;
  g1=new TGraphErrors();
  g1->SetName("res_allpions");

  Int_t npionts=0;

 
  for(Int_t ipoint = 0; ipoint < graphPRes->GetN(); ipoint++){
    g1->SetPoint(npionts,graphPRes->GetX()[ipoint],graphPRes->GetY()[ipoint]);
    g1->SetPointError(npionts,graphPRes->GetEX()[ipoint],graphPRes->GetEY()[ipoint]);
    npionts++;
  }
  
  //  pions TOF
  
  for(Int_t ipoint = 0; ipoint < graphRes->GetN(); ipoint++){
 
    g1->SetPoint(npionts,graphRes->GetX()[ipoint],graphRes->GetY()[ipoint]);
    g1->SetPointError(npionts,graphRes->GetEX()[ipoint],graphRes->GetEY()[ipoint]);
    npionts++;
    }
  //electrons
  g1->SetPoint(npionts,fPlateau->GetParameter(1),fPlateau->GetParameter(2)/fPlateau->GetParameter(1));
  g1->SetPointError(npionts,fPlateau->GetParError(1),fPlateau->GetParError(2)/fPlateau->GetParameter(1));

 
  Float_t sigmaPar[7];
  sigmaPar[0]=12;
  sigmaPar[1]=7.96380e-03;
  sigmaPar[2]=8.09869e-05;
  sigmaPar[3]=6.29707e-06;
  sigmaPar[4]=-3.27295e-08;
  sigmaPar[5]=-1.20200e+03;
  sigmaPar[6]=-3.97089e+01;

 

  TF1 *sigmarrrrr=new TF1("sigmaRes",SigmaFunc,47.5, 85,4);
  sigmarrrrr->FixParameter(0,sigmaPar[0]);
  g1->Fit(sigmarrrrr, "0R", "SAME");
  g1->Fit(sigmarrrrr, "0R", "SAME");
  g1->Fit(sigmarrrrr, "0R", "SAME");

  TF1 *fdedxvspP=new TF1("fdedxvspP","pol0",3.0,5.5);
  graphdEdxvsP->Fit(fdedxvspP, "R", "SAME");
  TF1 *fsigmavspP=new TF1("fsigmavspP","pol0",3.0,5.5);
  graphSigmavsP->Fit(fsigmavspP, "R", "SAME");



  for(Int_t i=0;i<6;++i)
    sigmarrrrr->SetParameter(i,sigmarrrrr->GetParameter(i));

  l1->Add(graphdEdxvsP);
  l1->Add(graphSigmavsP);

  l1->Add(fdedxvspP);
  l1->Add(fsigmavspP);

  l1->Add(sigmarrrrr);
  l1->Add(g1);

  TFile *fout=new TFile(Form("%s/hres_0_5_%s.root",outDir,endNamePos[index_eta]),"recreate");
  fout->cd();
  l1->Write();
  fout->Close();




}



void MakeFitsExternalData(const Char_t * inFile, const Char_t * outDir){


  TFile *fIn = TFile::Open(inFile);  
  TList * list = (TList *)fIn->Get("outputdedx");
  //plateau, electrons 0.4<p<0.6 GeV/c
  TH2D *hEVsEta = (TH2D *)list->FindObject("hPlateauVsEta");
  hEVsEta->Sumw2();
  //four th1, projections
  TH1D *hdEdxE[4]={0,0,0,0};// |eta|
  TH1D *hdEdxEpos[4]={0,0,0,0};// eta>0
  TF1 *fGaussE[4]={0,0,0,0};
  Int_t index_max = 1;
  for(Int_t i=0; i<4; ++i){

    Int_t index_min_negeta = i+index_max;
    Int_t index_max_negeta = i+index_max+1;
    Int_t index_min_poseta = 17-(i+index_max+1);
    Int_t index_max_poseta = 17-(i+index_max);

    hdEdxE[i] = (TH1D *)hEVsEta->ProjectionY(Form("hdEdxNegEta%d",i),index_min_negeta,index_max_negeta);
    hdEdxEpos[i] = (TH1D *)hEVsEta->ProjectionY(Form("hdEdxPosEta%d",i),index_min_poseta,index_max_poseta);
    hdEdxE[i]->Add(hdEdxEpos[i]);

    //Fitting the signal
    fGaussE[i]=new TF1(Form("fGauss%d",i),"gaus");
    hdEdxE[i]->Fit(fGaussE[i],"0","");
    hdEdxE[i]->Fit(fGaussE[i],"","", fGaussE[i]->GetParameter(1)-1.5*fGaussE[i]->GetParameter(2), fGaussE[i]->GetParameter(1)+1.5*fGaussE[i]->GetParameter(2));

    index_max++;

  }

  CreateDir(outDir);

  TFile *fout = 0;

  if(outDir){

    fout = new TFile(Form("%s/PrimaryElectrons.root",outDir),"recreate");
    fout->cd();
    for(Int_t i=0; i<4; ++i){
      hdEdxE[i]->Write();
      fGaussE[i]->Write();
    }
    fout->Close();


  }



}
//__________________________________________
void PlotQA(const Char_t * inFile){

  //gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  
  const Char_t* ending[8] = { "86", "64", "42", "20", "02", "24", "46", "68" };
  const Char_t* LatexEta[8] = { "-0.8<#eta<-0.6", "-0.6<#eta<-0.4", "-0.4<#eta<-0.2", "-0.2<#eta<0", 
				"0<#eta<0.2", "0.2<#eta<0.4", "0.4<#eta<0.6", "0.6<#eta<0.8" };
  TFile *fIn = TFile::Open(inFile);
  if(!fIn)
    return;

  TList * list = (TList *)fIn->Get("outputdedx");


  TH2D *hMIP[8];
  TH2D *hPlateau[8];
  TProfile *pMIP[8];
  TProfile *pPlateau[8];

  TH2D *hMIPVsNch[8];
  TProfile *pMIPVsNch[8];


  //Secondary pions
  TH2D *hPiV0s2D[8];
  TH1D *hPion3[8];
  TH1D *hPion2GeV = new TH1D("hPion2GeV","hPion2GeV",8,-0.8,0.8);


  for(Int_t index_eta = 0; index_eta < 8; ++index_eta){
    hMIP[index_eta] = (TH2D *)list->FindObject(Form("hMIPVsPhi%s",ending[index_eta]));
    hPlateau[index_eta] = (TH2D *)list->FindObject(Form("hPlateauVsPhi%s",ending[index_eta]));
    pMIP[index_eta] = (TProfile *)list->FindObject(Form("pMIPVsPhi%s",ending[index_eta]));
    pPlateau[index_eta] = (TProfile *)list->FindObject(Form("pPlateauVsPhi%s",ending[index_eta]));
    hPiV0s2D[index_eta] = (TH2D *)list->FindObject(Form("histPiV0%s",ending[index_eta]));
    hPion3[index_eta] = (TH1D *)hPiV0s2D[index_eta]->ProjectionY(Form("Pion3%s",ending[index_eta]),16,16);

    hMIPVsNch[index_eta] = (TH2D *)list->FindObject(Form("hMIPVsNch%s",ending[index_eta]));
    pMIPVsNch[index_eta] = (TProfile *)list->FindObject(Form("pMIPVsNch%s",ending[index_eta]));



    TF1 *fgaus = 0;
    fgaus = new TF1("fgaus","gausn",40,80);
    hPion3[index_eta]->Fit(fgaus,"0","",40,80);

    hPion2GeV->SetBinContent(index_eta+1,fgaus->GetParameter(1));
    hPion2GeV->SetBinError(index_eta+1,fgaus->GetParError(1));
  }



  TH2D *hMIPVsEta = (TH2D *)list->FindObject("hMIPVsEta");
  TProfile *pMIPVsEta = (TProfile *)list->FindObject("pMIPVsEta");
  TProfile *pMIPVsEtaV0s = (TProfile *)list->FindObject("pMIPVsEtaV0s");
  TProfile *pPlateauVsEta = (TProfile *)list->FindObject("pPlateauVsEta");
  TH2D *hPlateauVsEta = (TH2D *)list->FindObject("hPlateauVsEta");



  const Int_t nPadX = 4;
  const Int_t nPadY = 2;
  Float_t noMargin    = 0.000;
  Float_t topMargin   = 0.01;
  Float_t botMargin   = 0.1;
  Float_t leftMargin  = 0.05;
  Float_t rightMargin = 0.01;
  Float_t width       = (1-rightMargin-leftMargin)/nPadX;
  Float_t height      = (1-botMargin-topMargin)/nPadY;
  Float_t shift = 0.05;
    
  
  TCanvas* cMIP = new TCanvas("cMIP2", "Raa Pions", 900, 500);
  TPad* padMIP[nPadX*nPadY];
  //Float_t shift = 0.1;
  //Float_t shift = 0.05;
  const char* yTitleMIP = "d#it{E}/d#it{x}_{MIP}";
  const char* xTitleMIP = "#phi (rad)";
  cMIP->cd();
  
  
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(22);
  latex->SetTextAngle(0);
  latex->SetTextFont(62);
  
  latex->DrawLatex(0.53,0.04,xTitleMIP);
  latex->SetTextAngle(90);
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.025,0.5,yTitleMIP);
  latex->SetTextSize(0.05);

  
  for (Int_t i = 1; i < 9; i++) {
    
    Int_t iX = (i-1)%nPadX;
    Int_t iY = (i-1)/nPadX;
    Float_t x1 = leftMargin + iX*width;
    if(iX==2)
      x1-=0.01;
    else if(iX==1)
      x1+=0.01;
    Float_t x2 = leftMargin + (iX + 1)*width;
     if(iX==0)
       x2+=0.01;
     else if(iX==1)
       x2-=0.01;
     Float_t y1 = 1 - topMargin - (iY +1)*height;
     Float_t y2 = 1 - topMargin - iY*height;
     padMIP[i-1] = new TPad(Form("padRaa%d",i),"", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMargin; 
     Float_t mTop = noMargin; 
     Float_t mLeft = noMargin; 
     Float_t mRight = noMargin; 
     if(iY==0)       mTop   = shift;
     if(iY==nPadY-1) mBot   = shift;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadX-1) mRight = 0.08;
     padMIP[i-1]->SetBottomMargin(mBot);
     padMIP[i-1]->SetTopMargin(mTop);
     padMIP[i-1]->SetLeftMargin(mLeft);
     padMIP[i-1]->SetRightMargin(mRight);
     
     cMIP->cd(); 
     
     padMIP[i-1]->Draw();
     
  }

  latex->SetTextAngle(0);
  latex->SetTextSize(0.085);
  
  Int_t index_dec = 7;

  for(Int_t index_eta = 0; index_eta < 8; ++index_eta){
    padMIP[index_eta]->cd();
    padMIP[index_eta]->SetTickx(1);
    padMIP[index_eta]->SetTicky(1);
    padMIP[index_eta]->SetLogy(0);
    padMIP[index_eta]->SetLogx(0);
    hMIP[index_eta]->GetXaxis()->SetLabelSize(0.06);
    hMIP[index_eta]->GetYaxis()->SetLabelSize(0.06);
    hMIP[index_eta]->GetXaxis()->SetTitleSize(0);
    hMIP[index_eta]->GetYaxis()->SetTitleSize(0);

    if(index_eta<4){
      hMIP[index_eta]->Draw("colz"); 
      pMIP[index_eta]->Draw("samep"); 
      latex->DrawLatex(0.5,0.1,LatexEta[index_eta]);
    }
    else{
      hMIP[index_dec]->Draw("colz");
      pMIP[index_dec]->Draw("samep"); 
      latex->DrawLatex(0.5,0.1,LatexEta[index_dec]);
      index_dec--;
    }

  }



  cMIP->SaveAs("MIPvsPhi.png");

  //Plateau vs phi, different eta intervals 
  TCanvas* cPlateau = new TCanvas("cPlateau2", "Raa Pions", 900, 500);
  TPad* padPlateau[nPadX*nPadY];
 
  const char* yTitlePlateau = "d#it{E}/d#it{x}_{Plateau}";
  const char* xTitlePlateau = "#phi (rad)";
  cPlateau->cd();
  
  
  //TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(22);
  latex->SetTextAngle(0);
  latex->SetTextFont(62);
  
  latex->DrawLatex(0.53,0.04,xTitlePlateau);
  latex->SetTextAngle(90);
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.025,0.5,yTitlePlateau);
  latex->SetTextSize(0.05);

  
  for (Int_t i = 1; i < 9; i++) {
    
    Int_t iX = (i-1)%nPadX;
    Int_t iY = (i-1)/nPadX;
    Float_t x1 = leftMargin + iX*width;
    if(iX==2)
      x1-=0.01;
    else if(iX==1)
      x1+=0.01;
    Float_t x2 = leftMargin + (iX + 1)*width;
     if(iX==0)
       x2+=0.01;
     else if(iX==1)
       x2-=0.01;
     Float_t y1 = 1 - topMargin - (iY +1)*height;
     Float_t y2 = 1 - topMargin - iY*height;
     padPlateau[i-1] = new TPad(Form("padPlateau%d",i),"", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMargin; 
     Float_t mTop = noMargin; 
     Float_t mLeft = noMargin; 
     Float_t mRight = noMargin; 
     if(iY==0)       mTop   = shift;
     if(iY==nPadY-1) mBot   = shift;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadX-1) mRight = 0.08;
     padPlateau[i-1]->SetBottomMargin(mBot);
     padPlateau[i-1]->SetTopMargin(mTop);
     padPlateau[i-1]->SetLeftMargin(mLeft);
     padPlateau[i-1]->SetRightMargin(mRight);
     
     cPlateau->cd(); 
     
     padPlateau[i-1]->Draw();
     
  }

  latex->SetTextAngle(0);
  latex->SetTextSize(0.085);
  
  index_dec = 7;

  for(Int_t index_eta = 0; index_eta < 8; ++index_eta){
    padPlateau[index_eta]->cd();
    padPlateau[index_eta]->SetTickx(1);
    padPlateau[index_eta]->SetTicky(1);
    padPlateau[index_eta]->SetLogy(0);
    padPlateau[index_eta]->SetLogx(0);
    hPlateau[index_eta]->GetXaxis()->SetLabelSize(0.06);

    hPlateau[index_eta]->GetYaxis()->SetLabelSize(0.06);
    hPlateau[index_eta]->GetXaxis()->SetTitleSize(0);
    hPlateau[index_eta]->GetYaxis()->SetTitleSize(0);

    if(index_eta<4){
      hPlateau[index_eta]->Draw("colz"); 
      pPlateau[index_eta]->Draw("samep"); 
      latex->DrawLatex(0.5,0.1,LatexEta[index_eta]);
    }
    else{
      hPlateau[index_dec]->Draw("colz");
      pPlateau[index_dec]->Draw("samep"); 
      latex->DrawLatex(0.5,0.1,LatexEta[index_dec]);
      index_dec--;
    }

  }



  cPlateau->SaveAs("PlateauVsPhi.png");


  //MIP vs Nch
  TCanvas* cMIPVsNch = new TCanvas("cMIPVsNch2", "Raa Pions", 900, 500);
  TPad* padMIPVsNch[nPadX*nPadY];
 
  const char* yTitleMIPVsNch = "d#it{E}/d#it{x}_{MIP}";
  const char* xTitleMIPVsNch = "TPC-track multiplicity";
  cMIPVsNch->cd();
  
  
  latex->SetNDC();
  latex->SetTextAlign(22);
  latex->SetTextAngle(0);
  latex->SetTextFont(62);
  
  latex->DrawLatex(0.53,0.04,xTitleMIPVsNch);
  latex->SetTextAngle(90);
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.025,0.5,yTitleMIPVsNch);
  latex->SetTextSize(0.05);

  
  for (Int_t i = 1; i < 9; i++) {
    
    Int_t iX = (i-1)%nPadX;
    Int_t iY = (i-1)/nPadX;
    Float_t x1 = leftMargin + iX*width;
    if(iX==2)
      x1-=0.01;
    else if(iX==1)
      x1+=0.01;
    Float_t x2 = leftMargin + (iX + 1)*width;
     if(iX==0)
       x2+=0.01;
     else if(iX==1)
       x2-=0.01;
     Float_t y1 = 1 - topMargin - (iY +1)*height;
     Float_t y2 = 1 - topMargin - iY*height;
     padMIPVsNch[i-1] = new TPad(Form("padMIPVsNch%d",i),"", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMargin; 
     Float_t mTop = noMargin; 
     Float_t mLeft = noMargin; 
     Float_t mRight = noMargin; 
     if(iY==0)       mTop   = shift;
     if(iY==nPadY-1) mBot   = shift;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadX-1) mRight = 0.08;
     padMIPVsNch[i-1]->SetBottomMargin(mBot);
     padMIPVsNch[i-1]->SetTopMargin(mTop);
     padMIPVsNch[i-1]->SetLeftMargin(mLeft);
     padMIPVsNch[i-1]->SetRightMargin(mRight);
     
     if(i>=5&&i<9)
       padMIPVsNch[i-1]->SetBottomMargin(mBot+0.05);

     cMIPVsNch->cd(); 
     
     padMIPVsNch[i-1]->Draw();
     
  }

  latex->SetTextAngle(0);
  latex->SetTextSize(0.085);
  
  index_dec = 7;

  for(Int_t index_eta = 0; index_eta < 8; ++index_eta){
    padMIPVsNch[index_eta]->cd();

    padMIPVsNch[index_eta]->SetTickx(1);
    padMIPVsNch[index_eta]->SetTicky(1);
    padMIPVsNch[index_eta]->SetLogy(0);
    padMIPVsNch[index_eta]->SetLogx(1);
    hMIPVsNch[index_eta]->GetXaxis()->SetLabelSize(0.08);

    hMIPVsNch[index_eta]->GetXaxis()->SetRangeUser(10,2000);
    hMIPVsNch[index_eta]->GetYaxis()->SetRangeUser(43.5,54.5);
    hMIPVsNch[index_eta]->GetYaxis()->SetLabelSize(0.06);
    hMIPVsNch[index_eta]->GetXaxis()->SetTitleSize(0);
    hMIPVsNch[index_eta]->GetYaxis()->SetTitleSize(0);

    if(index_eta<4){
      hMIPVsNch[index_eta]->Draw("colz"); 
      pMIPVsNch[index_eta]->Draw("samep"); 
      latex->DrawLatex(0.5,0.1,LatexEta[index_eta]);
    }
    else{
      hMIPVsNch[index_dec]->Draw("colz");
      pMIPVsNch[index_dec]->Draw("samep"); 
      latex->DrawLatex(0.5,0.2,LatexEta[index_dec]);
      index_dec--;
    }

  }








  cMIPVsNch->SaveAs("MIPvsNch.png");









  //eta dependence, scaling
  const Int_t nPadXS = 3;
  const Int_t nPadYS = 1;
  Float_t noMarginS    = 0.0;
  Float_t topMarginS   = 0.01;
  Float_t botMarginS   = 0.1;
  Float_t leftMarginS  = 0.05;
  Float_t rightMarginS = 0.01;
  Float_t widthS       = (1-rightMarginS-leftMarginS)/nPadXS;
  Float_t heightS      = (1-botMarginS-topMarginS)/nPadYS;
  Float_t shiftS = 0.05;
    
  
  TCanvas* cMIPEta = new TCanvas("cMIPEta", "MIP eta", 900, 500);
  TPad* padMIPEta[nPadXS*nPadYS];
  //Float_t shift = 0.1;
  //Float_t shift = 0.05;
  const char* yTitleMIPS = "d#it{E}/d#it{x}";
  const char* xTitleMIPS = "#eta";
  cMIPEta->cd();
  
  
  //TLatexlatex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(22);
  latex->SetTextAngle(0);
  latex->SetTextFont(42);
  
  latex->DrawLatex(0.53,0.04,xTitleMIPS);
  latex->SetTextAngle(90);
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.025,0.5,yTitleMIPS);
  latex->SetTextSize(0.05);

  
  for (Int_t i = 1; i < 4; i++) {
    
    Int_t iX = (i-1)%nPadXS;
    Int_t iY = (i-1)/nPadXS;
    Float_t x1 = leftMarginS + iX*widthS;
    if(iX==2)
      x1-=0.01;
    else if(iX==1)
      x1+=0.01;
    Float_t x2 = leftMarginS + (iX + 1)*widthS;
     if(iX==0)
       x2+=0.01;
     else if(iX==1)
       x2-=0.01;
     Float_t y1 = 1 - topMarginS - (iY +1)*heightS;
     Float_t y2 = 1 - topMarginS - iY*heightS;
     padMIPEta[i-1] = new TPad(Form("padMIPEta%d",i),"", x1, y1, x2, y2, 0, 0, 0);
     Float_t mBot = noMarginS; 
     Float_t mTop = noMarginS; 
     Float_t mLeft = noMarginS; 
     Float_t mRight = noMarginS; 
     if(iY==0)       mTop   = shiftS;
     if(iY==nPadYS-1) mBot   = shiftS;
     if(iX==0)       mLeft  = 0.08;
     if(iX==nPadXS-1) mRight = 0.08;
     padMIPEta[i-1]->SetBottomMargin(mBot);
     padMIPEta[i-1]->SetTopMargin(mTop);
     padMIPEta[i-1]->SetLeftMargin(mLeft);
     padMIPEta[i-1]->SetRightMargin(mRight);
     
     cMIPEta->cd(); 
     
     padMIPEta[i-1]->Draw();
     
  }

  latex->SetTextAngle(0);
  latex->SetTextSize(0.085);
  

  TF1 *fetapos = new TF1("fetapos","pol3",0.0,0.8);
  TF1 *fetaneg = new TF1("fetaneg","pol3",-0.8,0.0);

  pMIPVsEta->Fit(fetaneg,"0","",-0.8,0.0);
  pMIPVsEta->Fit(fetapos,"0","",0.0,0.8);
  TLegend *leg1=new TLegend(0.15,0.5,0.4,0.7);

  fetapos->SetLineColor(1);
  fetapos->SetLineStyle(1);
  fetapos->SetLineWidth(2);

  fetaneg->SetLineColor(1);
  fetaneg->SetLineStyle(1);
  fetaneg->SetLineWidth(2);

  TH1D *hframe20 = new TH1D("hframe20","hframe20",8,-0.8,0.8);

  TH1D *hMIPPiPrim = new TH1D("hMIPPiPrim","",16,-0.8,0.8);
  for(Int_t bin =1; bin <= hMIPPiPrim->GetNbinsX(); ++bin){
    Double_t eta = pMIPVsEta->GetBinCenter(bin);
    Double_t dedx = pMIPVsEta->GetBinContent(bin);
    Double_t e_dedx = pMIPVsEta->GetBinError(bin);
    if(eta<0){
      dedx *= 50/fetaneg->Eval(eta);
      e_dedx *= 50/fetaneg->Eval(eta);
    }else{
      dedx *= 50/fetapos->Eval(eta);
      e_dedx *= 50/fetapos->Eval(eta);
    }
    hMIPPiPrim->SetBinContent(bin,dedx);
    hMIPPiPrim->SetBinError(bin,e_dedx);

  }

  TH1D *hPlateauPiPrim = new TH1D("hPlateauPiPrim","",16,-0.8,0.8);
  for(Int_t bin =1; bin <= hPlateauPiPrim->GetNbinsX(); ++bin){
    Double_t eta = pPlateauVsEta->GetBinCenter(bin);
    Double_t dedx = pPlateauVsEta->GetBinContent(bin);
    Double_t e_dedx = pPlateauVsEta->GetBinError(bin);
    if(eta<0){
      dedx *= 50/fetaneg->Eval(eta);
      e_dedx *= 50/fetaneg->Eval(eta);
    }else{
      dedx *= 50/fetapos->Eval(eta);
      e_dedx *= 50/fetapos->Eval(eta);
    }
    hPlateauPiPrim->SetBinContent(bin,dedx);
    hPlateauPiPrim->SetBinError(bin,e_dedx);

  }


  TH1D *h2GeVPiSec = new TH1D("h2GeVPiSec","",8,-0.8,0.8);
  for(Int_t bin =1; bin <= h2GeVPiSec->GetNbinsX(); ++bin){
    Double_t eta = hPion2GeV->GetBinCenter(bin);
    Double_t dedx = hPion2GeV->GetBinContent(bin);
    Double_t e_dedx = hPion2GeV->GetBinError(bin);
    if(eta<0){
      dedx *= 50/fetaneg->Eval(eta);
      e_dedx *= 50/fetaneg->Eval(eta);
    }else{
      dedx *= 50/fetapos->Eval(eta);
      e_dedx *= 50/fetapos->Eval(eta);
    }
    h2GeVPiSec->SetBinContent(bin,dedx);
    h2GeVPiSec->SetBinError(bin,e_dedx);

  }









  
  for(Int_t index_eta = 0; index_eta < 3; ++index_eta){
    padMIPEta[index_eta]->cd();
    padMIPEta[index_eta]->SetTickx(1);
    padMIPEta[index_eta]->SetTicky(1);
    padMIPEta[index_eta]->SetLogy(0);
    padMIPEta[index_eta]->SetLogx(0);
 
    if(index_eta==0){
      hframe20->GetYaxis()->SetRangeUser(30,95);
      hframe20->GetXaxis()->SetLabelSize(0.06);
      hframe20->GetXaxis()->SetLabelOffset(0.004);
      hframe20->GetYaxis()->SetLabelSize(0.06);
      hframe20->GetXaxis()->SetTitleSize(0);
      hframe20->GetYaxis()->SetTitleSize(0);
      hframe20->Draw();
      hMIPVsEta->Draw("colz same");
      hPlateauVsEta->Draw("colz same");

      latex->SetTextSize(0.05);
      latex->DrawLatex(0.55,0.45,"MIP Pions, 0.4<#it{p}<0.6 GeV/c");
      latex->DrawLatex(0.55,0.85,"Electrons, 0.4<#it{p}<0.6 GeV/c");

      }
 

    if(index_eta==1){
      pMIPVsEtaV0s->GetXaxis()->SetLabelSize(0.06);
      pMIPVsEtaV0s->GetYaxis()->SetLabelSize(0.06);
      pMIPVsEtaV0s->GetXaxis()->SetLabelOffset(0.004);
      pMIPVsEtaV0s->GetXaxis()->SetTitleSize(0);
      pMIPVsEtaV0s->GetYaxis()->SetTitleSize(0);
      pMIPVsEtaV0s->GetYaxis()->SetRangeUser(30,95);
      pMIPVsEtaV0s->SetMarkerStyle(29);
      pMIPVsEtaV0s->SetMarkerColor(4);
      pMIPVsEtaV0s->SetLineColor(4);
      pMIPVsEtaV0s->Draw();

      pMIPVsEta->SetMarkerStyle(25);
      pMIPVsEta->SetMarkerColor(1);
      pMIPVsEta->SetLineColor(1);
      pMIPVsEta->Draw("samep");
      pPlateauVsEta->SetMarkerStyle(21);
      pPlateauVsEta->SetMarkerColor(kBlue);
      pPlateauVsEta->SetLineColor(4);
      pPlateauVsEta->Draw("samep");

      //fetapos->Draw("same");
      //fetaneg->Draw("same");


      pMIPVsEta->Draw("samep");
      hPion2GeV->SetMarkerStyle(20);
      hPion2GeV->SetMarkerColor(4);
      hPion2GeV->SetLineColor(4);
      hPion2GeV->Draw("samep");
      pPlateauVsEta->Draw("samep");

    
      leg1->SetTextFont(42);
      leg1->SetTextSize(0.045);
      leg1->SetLineColor(kWhite);
      leg1->SetLineStyle(3);
      leg1->SetShadowColor(kWhite);
      leg1->SetFillColor(kWhite);
      leg1->SetFillStyle(0);
      leg1->AddEntry(pMIPVsEta,"Primary pions at MIP","P");
      leg1->AddEntry(pMIPVsEtaV0s,"Secondary pions at MIP","P");
      leg1->AddEntry(hPion2GeV,"Secondary pions at #beta#gamma ~ 14","P");
      leg1->AddEntry(pPlateauVsEta,"Electrons,  at #beta#gamma ~ 800","P");
      leg1->Draw();



    }



    if(index_eta==2){
      hMIPPiPrim->GetXaxis()->SetLabelSize(0.06);
      hMIPPiPrim->GetYaxis()->SetLabelSize(0.06);
      hMIPPiPrim->GetXaxis()->SetLabelOffset(0.004);
      hMIPPiPrim->GetXaxis()->SetTitleSize(0);
      hMIPPiPrim->GetYaxis()->SetTitleSize(0);
      hMIPPiPrim->GetYaxis()->SetRangeUser(30,95);
      hMIPPiPrim->SetMarkerStyle(29);
      hMIPPiPrim->SetMarkerColor(4);
      hMIPPiPrim->SetLineColor(4);
      hMIPPiPrim->Draw();


      hPlateauPiPrim->SetMarkerStyle(21);
      hPlateauPiPrim->SetLineColor(4);
      hPlateauPiPrim->SetMarkerColor(4);
      hPlateauPiPrim->Draw("samep");
      h2GeVPiSec->SetMarkerStyle(20);
      h2GeVPiSec->SetLineColor(4);
      h2GeVPiSec->SetMarkerColor(4);
      h2GeVPiSec->Draw("samep");

      latex->DrawLatex(0.45,0.85,"#LT d#it{E}/d#it{x} #GT #times 50 / #LT d#it{E}/d#it{x} #GT_{primary #pi at MIP}");

    }











  }


  cMIPEta->SaveAs("MIPvsEta.png");


  
  
}
TH2D * AddTwoSameBinningTH2D(TH2D *hPos, TH2D *hNeg, const Char_t *nameHist){
  
  TH2D *hout = (TH2D *)hPos->Clone(nameHist);
  hout->Reset();


  for(Int_t binx = 1; binx <= hPos->GetNbinsX(); ++binx){
    
    for(Int_t biny = 1; biny <= hPos->GetNbinsY(); ++biny){
      
      Double_t yield1 =  hPos->GetBinContent(binx,biny);
      Double_t yield2 =  hNeg->GetBinContent(binx,biny);
      Double_t e_yield1 =  hPos->GetBinError(binx,biny);
      Double_t e_yield2 =  hNeg->GetBinError(binx,biny);

      Double_t yield = yield1 + yield2;
      Double_t e_yield = TMath::Sqrt( e_yield1*e_yield1 + e_yield2*e_yield2 );

      hout->SetBinContent(binx,biny,yield);
      hout->SetBinError(binx,biny,e_yield);

    }
    
  }

  hout->SetEntries( hPos->GetEntries()+hNeg->GetEntries() );

  return hout;

}
