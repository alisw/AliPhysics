#include "AliAnalysisMultPbTrackHistoManager.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include "TDatabasePDG.h"
#include "AliPhysicsSelection.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisMultPbCentralitySelector.h"
#include "TLegend.h"
#include "AliBWTools.h"

using namespace std;

AliAnalysisMultPbTrackHistoManager * hManData = 0;
AliAnalysisMultPbTrackHistoManager * hManCorr = 0;
TH2D * hEvStatData = 0;
TH2D * hEvStatCorr = 0;

const Int_t kHistoFitCompoments = 2;
TH1D * gHistoCompoments[kHistoFitCompoments];

void LoadLibs(  Bool_t debug=0);
void LoadData(TString dataFolder, TString correctionFolder);
void SetStyle();
void CheckSecondaries(Double_t & fracWeak, Double_t &fracMaterial);
void CheckSanity(); 
void CheckCorrections();
void ShowAcceptanceInVzSlices() ;
TH1D * GetRatioIntegratedFractions (TH1 * hNum, TH1 * hDenum) ;
TH1D * GetCumulativeHisto (TH1 * h) ;
static Double_t HistoSum(const double * x, const double* p);
TF1 * GetFunctionHistoSum() ;
TF1 * GetMTExp(Float_t norm=68, Float_t t=25) ;
TF1 * GetHagedorn(Float_t norm=68, Float_t pt0=25, Float_t n=13) ;
TF1 * GetLevy(Double_t temp=0.1, Double_t n=7, Double_t norm=10, const char * name="fLevy") ;
void PrintCanvas(TCanvas* c,const TString formats) ;

// global switches
Bool_t doPrint=kFALSE;// disable PrintCanvas
Float_t zmin = -10;
Float_t zmax = 10;
Float_t etaMin = -0.5;
Float_t etaMax = 0.5;
Int_t gCentralityBin = -1;

#define CORRECT_2D

void correct(TString dataFolder = "./output/LHC10g2d_130844_V0M_bin_10/", TString correctionFolder = "./output/LHC10g2a_130844_V0M_bin_10/", 
	     Float_t zminl=-10, Float_t zmaxl=10, Float_t etaMinl = -0.5, Float_t etaMaxl=0.5, Float_t npart=381.3, Float_t scaleWeakByHand = -1, Int_t centrBin = -1) {

  // Set globals 
  zmin = zminl;
  zmax = zmaxl;
  etaMin = etaMinl;
  etaMax = etaMaxl;
  gCentralityBin = centrBin;

  // Load stuff and set some styles
  LoadLibs();
  LoadData(dataFolder,correctionFolder);
  SetStyle();
  // ShowAcceptanceInVzSlices();
  // return;

  // TODO add some cool printout for cuts and centrality selection
  
  CheckSanity();

  Double_t fractionWeak = 1, fractionMaterial=1; 
  CheckSecondaries(fractionWeak, fractionMaterial);
  cout << "Rescaling secondaries correction, weak: " << fractionWeak << ", material: " << fractionMaterial <<endl;
  if(scaleWeakByHand >= 0) {
    fractionWeak = scaleWeakByHand;
    cout << "Scaling Strangeness by hand (factor " << fractionWeak <<")" << endl;
    
  }
  // Some shorthands

#if defined CORRECT_1D
  TH1D * hDataPt   = (TH1D*) hManData->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRec, etaMin,etaMax,zmin,zmax)->Clone("hDataPt");
  TH1D * hMCPtGen  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,         etaMin,etaMax,zmin,zmax); 
  //zTH1D * hMCPtGen  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,         etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtRec  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRec,         etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtPri  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim,     etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtSeM  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat,   etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtSeW  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak,  etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtFak  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecFake,     etaMin,etaMax,zmin,zmax);
#elif defined CORRECT_2D
  TH1  * hDataPt   = (TH2D*) hManData->GetHistoPtVz(AliAnalysisMultPbTrackHistoManager::kHistoRec,         etaMin,etaMax)->Clone("hDataPt");
  TH1  * hMCPtGen  = (TH2D*) hManCorr->GetHistoPtVz(AliAnalysisMultPbTrackHistoManager::kHistoGen,         etaMin,etaMax);
  TH1  * hMCPtRec  = (TH2D*) hManCorr->GetHistoPtVz(AliAnalysisMultPbTrackHistoManager::kHistoRec,         etaMin,etaMax);
  TH1  * hMCPtPri  = (TH2D*) hManCorr->GetHistoPtVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim,     etaMin,etaMax);
  TH1  * hMCPtSeM  = (TH2D*) hManCorr->GetHistoPtVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat,   etaMin,etaMax);
  TH1  * hMCPtSeW  = (TH2D*) hManCorr->GetHistoPtVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak,  etaMin,etaMax);
  TH1  * hMCPtFak  = (TH2D*) hManCorr->GetHistoPtVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake,     etaMin,etaMax);
#endif
 
  //  hMCPtRec->Draw("same");

  TCanvas * cMC = new TCanvas ("cMC", "Monte Carlo");
  cMC->SetLogy();
  cMC->cd();
  // hMCPtSeW->Add(hMCPtSeM);
  // hMCPtSeW->Divide(hMCPtRec);


  hMCPtGen ->Draw();
  hMCPtRec ->Draw("same");
  hMCPtPri ->Draw("same");
  hMCPtSeM ->Draw("same");
  hMCPtSeW ->Draw("same");
  hMCPtFak ->Draw("same");
#ifdef CORRECT_1D  
  hMCPtGen ->GetXaxis()->SetRangeUser(0,4.5);
  hMCPtGen ->GetYaxis()->SetRangeUser(0.1,1e4);
#endif
  TLegend * lMC = new TLegend(0.505034, 0.59965, 0.877517, 0.926573,"Monte Carlo");
  lMC->AddEntry( hMCPtGen, "Generated");
  lMC->AddEntry( hMCPtRec, "All Rec");
  lMC->AddEntry( hMCPtPri, "Rec Primaries");
  lMC->AddEntry( hMCPtSeM, "Rec Sec. Material");
  lMC->AddEntry( hMCPtSeW, "Rec Sec. Weak");
  lMC->AddEntry( hMCPtFak, "Rec Fakes");
  lMC->Draw();


  cout << "Fake/All Rec  = " << hMCPtFak->Integral()/hMCPtRec->Integral()  << endl;
  cout << "SM/All   Rec  = " << hMCPtSeM->Integral()/hMCPtRec->Integral()  << endl;
  cout << "SW/All   Rec  = " << hMCPtSeW->Integral()/hMCPtRec->Integral()  << endl;


  // Compute efficiency and subtract secondaries and fakes after rescaling to data
  // PRIM_DATA = ALL_DATA - SEC_MC/ALL_MC*ALL_DATA - FAK_MC/ALL_MC*ALL_DATA
  // TRUE_DATA = PRIM_DATA * GEN_MC/PRIM_MC

  TH1D * hEffPt = (TH1D*) hMCPtPri->Clone("hEffPt");
  hEffPt->Divide(hMCPtPri,hMCPtGen,1,1,"B");

  TH1D * hCorSeM = (TH1D*) hMCPtSeM->Clone("hCorSeM");
  hCorSeM->Divide(hMCPtSeM,hMCPtRec,1,1,"B");
  hCorSeM->Scale(fractionMaterial);// rescale material correction
  hCorSeM->Multiply(hDataPt); 

  TH1D * hCorSeW = (TH1D*) hMCPtSeW->Clone("hCorSeW");
  hCorSeW->Divide(hMCPtSeW,hMCPtRec,1,1,"B");
  hCorSeW->Scale(fractionWeak);// rescale weak correction
  hCorSeW->Multiply(hDataPt); 

  TH1D * hCorFak = (TH1D*) hMCPtFak->Clone("hCorFak");
  hCorFak->Divide(hMCPtFak,hMCPtRec,1,1,"B");
  hCorFak->Multiply(hDataPt);

  TH1D * hCorrected = (TH1D*) hDataPt->Clone("hCorrected");
  hCorrected->Add(hCorSeM,-1);
  hCorrected->Add(hCorSeW,-1);
  hCorrected->Add(hCorFak,-1);
  hCorrected->Divide(hEffPt);
  hCorrected->SetMarkerStyle(kOpenStar);


  TCanvas * cCorrections = new TCanvas("cCorrections", "cCorrections");

  hEffPt->Draw();
  // hCorSeM->Draw();
  // hCorSeM->SetLineColor(kRed);
  // hCorSeM->SetMarkerColor(kRed);
  // hMCPtSeM->Draw("same");  
  // hCorSeW->Draw("same");
  // hCorSeW->SetLineColor(kRed);
  // hCorSeW->SetMarkerColor(kRed);
  // hMCPtSeW->Draw("same");


  // hDataPt->Draw();
  // return;
  TCanvas * cdata = new TCanvas ("cData", "Data");  
  cdata->SetLogy();
  cdata->cd();  
#ifdef CORRECT_2D
  Int_t minProj = hDataPt->GetYaxis()->FindBin(zmin);
  Int_t maxProj = hDataPt->GetYaxis()->FindBin(zmax-0.0001);

  cout << minProj << "-" << maxProj << endl;
  
  // This correction accounts for the vertex cut;

  hDataPt    = ((TH2D*)hDataPt   )->ProjectionX("_px",minProj,maxProj,"E");
  hCorrected = ((TH2D*)hCorrected)->ProjectionX("_px",minProj,maxProj,"E");
  hMCPtGen   = ((TH2D*)hMCPtGen  )->ProjectionX("_px",minProj,maxProj,"E");

  Float_t vertexCutCorrection = 
    hManCorr->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoGen)->Integral(-1,-1) /
    hManCorr->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoGen)->Integral(minProj,maxProj) ;
  // cout << vertexCutCorrection << " " << hMCPtGen->Integral(-1,-1)  << " " << hMCPtPri->Integral() << endl;
  // vertexCutCorrection /= hMCPtGen->Integral(-1,-1);
  vertexCutCorrection = 1; // FIXME
  cout << "Vertex cut correction " << vertexCutCorrection << " (Efficiency " << 1./vertexCutCorrection << ")" << endl;

  hDataPt    ->Scale(1.,"width");
  hCorrected ->Scale(vertexCutCorrection,"width");
  hMCPtGen   ->Scale(1.,"width");
#endif

  hDataPt->Draw();
  hDataPt ->GetXaxis()->SetRangeUser(0,4.5);
  hDataPt ->GetYaxis()->SetRangeUser(0.1,1e4);
  hCorrected->SetLineColor(kBlack);
  hCorrected->SetMarkerColor(kBlack);
  hCorrected->Draw("same");
  hMCPtGen->DrawClone("same");
  //  TF1 * f = GetLevy();
  TF1 * f = GetHagedorn();
  hCorrected->Fit(f,"", "same");
  hCorrected->Fit(f,"IME", "same",0,2);
  cout << "dN/deta (function)  = " << f->Integral(0,100) << " +- " << f->IntegralError(0,100) << endl;
  cout << "dN/deta (data>0.15) = " <<  hCorrected->Integral(4,-1,"width") << endl;//
  cout << "dN/deta (func+data) = " << f->Integral(0,0.1) + hCorrected->Integral(3,-1,"width") << endl;//
  cout << "dN/deta (func 0.15+data) = " << f->Integral(0,0.15) + hCorrected->Integral(4,-1,"width") << endl;//
  cout << "Generated dN/deta (correction) = " << hMCPtGen->Integral("width") << endl;
  cout << "Generated dN/deta (correction, <0.13) = " << hMCPtGen->Integral(1,2,"width") << endl;
  cout << "-------" << endl;
  cout << "Rescaling secondaries correction, weak: " << fractionWeak << ", material: " << fractionMaterial <<endl;
  Double_t errorData = 0;
  Float_t dNdeta  =  (f->Integral(0,0.1) + hCorrected->IntegralAndError(1,-1,errorData, "width")) / (etaMax-etaMin);
  Float_t dNdetaE = TMath::Sqrt(f->IntegralError(0,0.1)*f->IntegralError(0,0.1) + errorData*errorData) / (etaMax-etaMin);
  cout << "dN/deta              " << dNdeta << " +- " << dNdetaE << endl;
  cout << "(dN/deta)/0.5 Npart  " << dNdeta/npart*2 << " +- " << dNdetaE/npart*2 << endl;
  Double_t mean, meanerr;
  AliBWTools::GetMeanDataAndExtrapolation(hCorrected, f, mean, meanerr,0,1000);
  cout << "Mean Pt " << mean << "+-" << meanerr << endl;
  

  TH1D * hDataGen  = hManData->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,        etaMin,etaMax,zmin,zmax);
  cout << "Generated dN/deta (data) =       " << hDataGen->Integral("width") << endl;
  hDataGen->Draw("same");  
  TLegend * l = new TLegend(0.520134, 0.676573, 0.885906, 0.923077,"137161, p1+++");
  l->AddEntry(hDataPt, "Raw data");
  l->AddEntry(hCorrected, "Corrected data");
  l->AddEntry(hMCPtGen, "Monte Carlo (generated)");
  l->AddEntry(f, "Hagedorn Fit");
  l->Draw();


  Double_t errXCheck = 0, yieldXCheck=0;
  AliBWTools::GetYield(hCorrected,f,yieldXCheck,errXCheck);
  cout << "BWTools crosscheck: " << yieldXCheck << "+-" << errXCheck << endl;
  

}

void CheckSecondaries(Double_t &fracWeak, Double_t &fracMaterial) {
  // Returns the fraction you need to rescale the secondaries from weak decays for.

  // Some shorthands
  TH2D * hDataDCAPt   = hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  //  TH1D * hMCDCAGen  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen       );
  TH2D * hMCDCARecPt  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  TH2D * hMCDCAPriPt  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim   );
  TH2D * hMCDCASWPt   = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak    );
  TH2D * hMCDCASMPt   = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat    );
  TH2D * hMCDCAFakPt  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake   );
 

  TCanvas * cCumulative = new TCanvas("cDCAculumative","DCA cumulative distributions");
  cCumulative->cd();
  GetCumulativeHisto(hMCDCAPriPt )->Draw();
  GetRatioIntegratedFractions(hMCDCASWPt, hMCDCARecPt  )->Draw("same");
  GetRatioIntegratedFractions(hMCDCASMPt, hMCDCARecPt  )->Draw("same");
  GetRatioIntegratedFractions(hMCDCAPriPt,hMCDCARecPt  )->Draw("same");


  TCanvas * c1 = new TCanvas("cDCAFit", "Fit to the DCA distributions");  
  //  c1->Divide(4,2);
  c1->SetLogy();
  // Draw all
  //  hDataDCA->Draw();
  //  hMCDCARec ->Draw("same");
  // hMCDCAPri ->Draw("same");
  // hMCDCASW ->Draw("same");
  // hMCDCASM ->Draw("same");
  // hMCDCAFak ->Draw("same");
  // return;
  
  // Fit the DCA distribution. Uses a TF1 made by summing histograms
  TH2D * hMCPrimSMFakPt = (TH2D*) hMCDCAPriPt->Clone("hMCPrimSMFak");
  hMCPrimSMFakPt->Add(hMCDCASMPt);
  hMCPrimSMFakPt->Add(hMCDCAFakPt);

  // Set the components which are used in HistoSum, the static
  // function for GetFunctionHistoSum
  // Project onti DCA axis
  const Int_t ptBinsFit[] = {3,5,7,9,11,15,19,23,31,-1};
  //const Int_t ptBinsFit[] = {3,20,-1};
  Int_t ibinPt = -1;
  while(ptBinsFit[(++ibinPt)+1]!=-1){
    c1->cd(ibinPt+1);
    gPad->SetLogy();
    Int_t minPtBin=ptBinsFit[ibinPt];
    Int_t maxPtBin=ptBinsFit[ibinPt+1]-1;
  
    TH1D * hDataDCA     = (TH1D*) hDataDCAPt    ->ProjectionX(TString(hDataDCAPt    ->GetName())+long(ibinPt)+"_px",minPtBin,maxPtBin,"E")->Clone();
    TH1D * hMCPrimSMFak = (TH1D*) hMCPrimSMFakPt->ProjectionX(TString(hMCPrimSMFakPt->GetName())+long(ibinPt)+"_px",minPtBin,maxPtBin,"E")->Clone();
    TH1D * hMCDCASW     = (TH1D*) hMCDCASWPt    ->ProjectionX(TString(hMCDCASWPt    ->GetName())+long(ibinPt)+"_px",minPtBin,maxPtBin,"E")->Clone();
    TH1D * hMCDCASM     = (TH1D*) hMCDCASMPt    ->ProjectionX(TString(hMCDCASMPt    ->GetName())+long(ibinPt)+"_px",minPtBin,maxPtBin,"E")->Clone();
    gHistoCompoments[0] = (TH1D*) hMCPrimSMFakPt->ProjectionX(TString(hMCPrimSMFakPt->GetName())+long(ibinPt)+"_px",minPtBin,maxPtBin,"E")->Clone();
    gHistoCompoments[1] = (TH1D*) hMCDCASWPt    ->ProjectionX(TString(hMCDCASWPt    ->GetName())+long(ibinPt)+"_px",minPtBin,maxPtBin,"E")->Clone();
    //  gHistoCompoments[2] = (TH1D*) hMCDCASMPt    ->ProjectionX("_px",minPtBin,maxPtBin,"E")->Clone();
    TF1 * fHistos = GetFunctionHistoSum();
    fHistos->SetParameters(1,1);
    //  fHistos->FixParameter(2,1);
    fHistos->SetLineColor(kRed);
    hDataDCA= (TH1D*) hDataDCAPt->ProjectionX(TString(hDataDCAPt->GetName())+long(ibinPt)+"_px",minPtBin,maxPtBin,"E")->Clone();
    hDataDCA->Draw();
    // Fit!
    hDataDCA->Fit(fHistos,"0Q","",-2,2);
    // Rescale the components and draw to see how it looks like
    hMCPrimSMFak->Scale(fHistos->GetParameter(0));
    hMCDCASW    ->Scale(fHistos->GetParameter(1));
    hMCDCASM    ->Scale(fHistos->GetParameter(2));
    hMCPrimSMFak->Draw("same");
    hMCDCASW    ->Draw("same");
    hMCDCASM    ->Draw("same");
    TH1D* hMCSum = (TH1D*) hMCPrimSMFak->Clone();
    hMCSum->Add(hMCDCASW);
    //    hMCSum->Add(hMCDCASM);
    hMCSum->SetLineColor(kRed);
    hMCSum->SetMarkerStyle(1);
    hMCSum->Draw("lhist,same");
    //    fHistos->Draw("same");
    // compute scaling factors
    //  fracWeak     = fHistos->GetParameter(1)/(fHistos->GetParameter(0)+fHistos->GetParameter(1));
    cout << hDataDCAPt->GetYaxis()->GetBinLowEdge(minPtBin) << " - ";
    cout << hDataDCAPt->GetYaxis()->GetBinLowEdge(maxPtBin+1) << " = ";
    cout << 2*fHistos->GetParameter(1)/(fHistos->GetParameter(0)+fHistos->GetParameter(1));
    cout << " ("<< fHistos->GetChisquare() <<")";
    cout << " " << fHistos->GetParameter(0) << "," << fHistos->GetParameter(1) << endl;
    
    
   fracWeak     = 2*fHistos->GetParameter(1)/(fHistos->GetParameter(0)+fHistos->GetParameter(1));
     //  fracMaterial = fHistos->GetParameter(2)/(fHistos->GetParameter(0);
  }

}

void CheckSanity() {
  // compares various distributions in data and in MC
  TCanvas * c1 = new TCanvas("cVz", "Vertex Z distribution");
  c1->cd();
  TH1 * hData  = hManData->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  TH1 * hCorr  = hManCorr->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  hData->SetLineColor  (kBlack);
  hData->SetMarkerColor(kBlack);
  hData->SetMarkerStyle(kFullCircle);
  hCorr->SetLineColor  (kRed);
  hCorr->SetMarkerColor(kRed);
  hCorr->SetMarkerStyle(kOpenSquare);
  TLegend * l = new TLegend(0.575503, 0.70979, 0.86745, 0.917832);
  l->AddEntry(hData, "Data");
  l->AddEntry(hCorr, "MC");
  hCorr->Draw("");
  hData->Draw("same");
  l->Draw();


  TCanvas * c2 = new TCanvas("cEta", "Eta distribution");
  c2->cd();
  hData  = hManData->GetHistoEta(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  hCorr  = hManCorr->GetHistoEta(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  hData->SetLineColor  (kBlack);
  hData->SetMarkerColor(kBlack);
  hData->SetMarkerStyle(kFullCircle);
  hCorr->SetLineColor  (kRed);
  hCorr->SetMarkerColor(kRed);
  hCorr->SetMarkerStyle(kOpenSquare);
  l = new TLegend(0.575503, 0.70979, 0.86745, 0.917832);
  l->AddEntry(hData, "Data");
  l->AddEntry(hCorr, "MC");
  hData->Draw("");
  hCorr->Draw("same");
  l->Draw();

  TCanvas * c3 = new TCanvas("cMult", "Multiplicity distribution");
  c3->cd();
  hData  = hManData->GetHistoMult(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  hCorr  = hManCorr->GetHistoMult(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  hData->SetLineColor  (kBlack);
  hData->SetMarkerColor(kBlack);
  hData->SetMarkerStyle(kFullCircle);
  hCorr->SetLineColor  (kRed);
  hCorr->SetMarkerColor(kRed);
  hCorr->SetMarkerStyle(kOpenSquare);
  l = new TLegend(0.575503, 0.70979, 0.86745, 0.917832);
  l->AddEntry(hData, "Data");
  l->AddEntry(hCorr, "MC");
  hData->Draw("");
  hCorr->Draw("same");
  l->Draw();

  TCanvas * c4 = new TCanvas("cStatistics", "cStats");
  c4->cd();
  hData  = hManData->GetHistoStats();
  hCorr  = hManCorr->GetHistoStats();
  hData->SetLineColor  (kBlack);
  hData->SetMarkerColor(kBlack);
  //  hData->SetMarkerStyle(kFullCircle);
  hCorr->SetLineColor  (kRed);
  hCorr->SetMarkerColor(kRed);
  //  hCorr->SetMarkerStyle(kOpenSquare);
  l = new TLegend(0.575503, 0.70979, 0.86745, 0.917832);
  l->AddEntry(hData, "Data");
  l->AddEntry(hCorr, "MC");
  hData->Draw("");
  hCorr->Draw("same");
  l->Draw();


}

void CheckCorrections() {

  //  TH1D * hDataPt   = (TH1D*) hManData->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRec, etaMin,etaMax,zmin,zmax)->Clone("hDataPt");
  TH1D * hMCPtGen  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,         etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtRec  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRec,         etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtPri  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim,     etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtSeM  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat,   etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtSeW  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak,  etaMin,etaMax,zmin,zmax);
  TH1D * hMCPtFak  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecFake,     etaMin,etaMax,zmin,zmax);

  hMCPtSeM->Divide(hMCPtRec);
  hMCPtSeW->Divide(hMCPtRec);
  hMCPtPri->Divide(hMCPtGen);
  TCanvas * c1 = new TCanvas("cSecondaries", "Secondaries");
  hMCPtSeM->Draw();
  hMCPtSeW->Draw("same");

  TCanvas * c2 = new TCanvas("cEfficiency", "Efficiency");
  hMCPtPri->Draw();

}

void LoadLibs(  Bool_t debug) {

  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWG2Spectra");
  gSystem->Load("libPWG0base"); 
   
  gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPbPb"));
  gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
  gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG2/SPECTRA/Fit"));

  // Load helper classes
  // TODO: replace this by a list of TOBJStrings
  TString taskName("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisTaskMultPbTracks.cxx+");
  TString histoManName("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisMultPbTrackHistoManager.cxx+");
  TString centrName("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisMultPbCentralitySelector.cxx+");
  TString listName("$ALICE_ROOT/PWG1/background/AliHistoListWrapper.cxx+");

  gSystem->ExpandPathName(taskName);
  gSystem->ExpandPathName(histoManName);
  gSystem->ExpandPathName(centrName);
  gSystem->ExpandPathName(listName);


  gROOT->LoadMacro(listName    +(debug?"+g":""));   
  gROOT->LoadMacro(histoManName+(debug?"+g":""));
  gROOT->LoadMacro(centrName   +(debug?"+g":""));   
  gROOT->LoadMacro(taskName    +(debug?"+g":""));   

  // Histo fitter
  gROOT->LoadMacro("/Users/mfloris/Work/ALICE/ANALYSIS/HistoFitter/fcn.cxx+g");
  gROOT->LoadMacro("/Users/mfloris/Work/ALICE/ANALYSIS/HistoFitter/AliHistoFitter.cxx+g");


}


void SetStyle() {

  hManData->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kBlack);    
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetLineColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetLineColor(kGreen);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetLineColor(kBlue );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetLineColor(kBlue );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetLineColor(kCyan );

  hManData->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kBlack);    
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerColor(kGreen);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerColor(kBlue );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerColor(kCyan );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerColor(kYellow );

  hManData->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kFullCircle);    
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerStyle(kFullSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerStyle(kOpenCircle);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerStyle(kOpenSquare);

 hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kBlack);    
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetLineColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetLineColor(kGreen);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetLineColor(kBlue );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetLineColor(kCyan);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetLineColor(kYellow );

  hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kBlack);    
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerColor(kGreen);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerColor(kBlue );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerColor(kCyan );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerColor(kYellow );

  hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kFullCircle);    
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerStyle(kFullSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerStyle(kOpenSquare);


}

void LoadData(TString dataFolder, TString correctionFolder){

  TString file = "/multPbPbtracks.root";
  if (gCentralityBin != -1 ) file.ReplaceAll(".root", Form("_%2.2d.root", gCentralityBin));
  // Get histo manager for data and MC + stat histos
  TFile * fData = new TFile(dataFolder+file);
  TFile * fCorr = new TFile(correctionFolder+file);
  TFile * fStatData = new TFile(dataFolder+"/event_stat.root");
  TFile * fStatCorr = new TFile(correctionFolder+"/event_stat.root");

  hManData = (AliAnalysisMultPbTrackHistoManager*) fData->Get("histoManager");
  hManCorr = (AliAnalysisMultPbTrackHistoManager*) fCorr->Get("histoManager");
  AliESDtrackCuts * cutsData = (AliESDtrackCuts*) fData->Get("AliESDtrackCuts");
  AliESDtrackCuts * cutsCorr = (AliESDtrackCuts*) fCorr->Get("AliESDtrackCuts");

  // FIXME: implement operator= and Print in AliESDtrackCuts
  // if (cutsData != cutsCorr) {
  //   cout << "ERROR: MC and data do not have the same cuts" << endl;
  //   // }
  cutsData->Print();
  hEvStatData = (TH2D*) fStatData->Get("fHistStatistics");
  hEvStatCorr = (TH2D*) fStatCorr->Get("fHistStatistics");

  AliAnalysisMultPbCentralitySelector * centrData = (AliAnalysisMultPbCentralitySelector*) fData->Get("Cuts");
  if(!centrData) {
    cout << "ERROR:  cannot read centrality data" << endl;
  }
  centrData->Print();
  // Normalize
  Int_t irowGoodTrigger = 1;
  if (hEvStatCorr && hEvStatData) {
    //  hManData->ScaleHistos(75351.36/1.015);// Nint for run 104892 estimated correcting for the trigger efficiency, multiplied for the physics selection efficiency which I'm not correcting for the time being
    // hManData->ScaleHistos(hEvStatData->GetBinContent(AliPhysicsSelection::kStatAccepted,irowGoodTrigger));
    // hManCorr->ScaleHistos(hEvStatCorr->GetBinContent(AliPhysicsSelection::kStatAccepted,irowGoodTrigger));
    // hManData->ScaleHistos(hManData->GetHistoStats()->GetBinContent(AliAnalysisMultPbTrackHistoManager::kStatVtx));
    // hManCorr->ScaleHistos(hManCorr->GetHistoStats()->GetBinContent(AliAnalysisMultPbTrackHistoManager::kStatVtx));
    TH1D* hvzData = hManData->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec);
    TH1D* hvzCorr = hManCorr->GetHistoVzEvent(AliAnalysisMultPbTrackHistoManager::kHistoRec);
    hManData->ScaleHistos(hvzData->Integral(hvzData->FindBin(zmin),hvzData->FindBin(zmax-0.001)));
    hManCorr->ScaleHistos(hvzCorr->Integral(hvzCorr->FindBin(zmin),hvzCorr->FindBin(zmax-0.001)));
  } else {
    cout << "WARNING!!! ARBITRARY SCALING" << endl;
    hManData->ScaleHistos(1000);
    hManCorr->ScaleHistos(1000);    
  }
}

TF1 * GetHagedorn(Float_t norm, Float_t pt0, Float_t n) {

  TF1 * f =0;
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  f=new TF1("fHagedorn",Form("(x/sqrt(x*x+%f*%f))*x*[0]*( 1 + x/[1] )^(-[2])",mass,mass),0,10);
  f->SetParameters(norm, pt0, n);
  f->SetParLimits(1, 0.01, 10);
  f->SetParNames("norm", "pt0", "n");
  f->SetLineWidth(1);
  return f;


}

TF1 * GetMTExp(Float_t norm, Float_t t) {

  TF1 * f =0;
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  f=new TF1("fMTExp",Form("(x/sqrt(x*x+%f*%f))*x*[0]*exp(-sqrt(x*x+%f*%f)/[1])",mass,mass,mass,mass),0,10);
  f->SetParameters(norm, t);
  //  f->SetParLimits(1, 0.01);
  f->SetParNames("norm", "T");
  f->SetLineWidth(1);
  return f;


}

TF1 * GetLevy(Double_t temp, Double_t n, Double_t norm, const char * name) {

  char formula[500];
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  

  sprintf(formula,"(x/sqrt(x*x+[3]*[3]))*x*( [0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (sqrt([3]*[3]+x*x) -[3])/([1]*[2])  )^(-[1])");
  TF1* f=new TF1(name,formula,0,10);
  f->SetParameters(norm, n, temp,mass);
  f->SetParLimits(2, 0.01, 10);
  f->SetParNames("norm (dN/dy)", "n", "T", "mass");
  f->FixParameter(3,mass);
  f->SetLineWidth(1);
  return f;


}


TH1D * GetCumulativeHisto (TH1 * h) { 
  // Returns a cumulative histogram
  // TODO: put this method in a tools class
  TH1D * hInt = (TH1D*) h->Clone(TString(h->GetName())+"_Int");
  hInt->Reset();
  Float_t integral = h->Integral(-1,-1); // Consider under/over flow!
  Int_t nbin = h->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t content = h->Integral(-1,ibin) / integral;
    hInt->SetBinContent(ibin, content);
  }
  return hInt;
}

TH1D * GetRatioIntegratedFractions (TH1 * hNum, TH1 * hDenum) { 
  // Returns the the ratio of integrated histograms 
  // TODO: put this method in a tools class
  TH1D * hRatio = (TH1D*) hNum->Clone(TString(hNum->GetName())+hDenum->GetName()+"_RatioIntegrated");
  hRatio->Reset();
  Int_t nbin = hNum->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t content = hNum->Integral(-1,ibin) / hDenum->Integral(-1,ibin);// consider underflow
    hRatio->SetBinContent(ibin, content);
  }
  return hRatio;
}

void ShowAcceptanceInVzSlices() {
  TCanvas * cvz = new TCanvas("cvz","acc #times eff vs vz");
  for(Int_t ivz = -10; ivz < -6; ivz+=2){
    ivz=0;//FIXME  
    Bool_t first = kTRUE;
    TH1D * hMCPtPri  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim ,   etaMin,etaMax,ivz,ivz+2);
    TH1D * hMCPtGen  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,        etaMin,etaMax,ivz,ivz+2);
    TH1D * hMCPtPriD  = hManData->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim ,   etaMin,etaMax,ivz,ivz+2);
    TH1D * hMCPtGenD  = hManData->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,        etaMin,etaMax,ivz,ivz+2);
    //    hEff= hMCPtGen;
    TH1D * hEff = (TH1D*)hMCPtPri->Clone(Form("hEff_vz_%d_%d",ivz,ivz+2));
    hEff->Divide(hMCPtPri,hMCPtGen,1,1,"B");
    TH1D * hEffD = (TH1D*)hMCPtPriD->Clone(Form("hEffD_vz_%d_%d",ivz,ivz+2));
    hEffD->Divide(hMCPtPriD,hMCPtGenD,1,1,"B");
    hEffD->SetLineColor(kRed);
    cout << "ivz " << ivz << endl;
    if(first) {
      first=kFALSE;
      cout << "First" << endl;
      hEff->Draw();
      hEffD->Draw("same");
      // hMCPtGen->Draw();
      // hMCPtPri->Draw("same");
    }
    else {
      cout << "Same" << endl;
      hEff->Draw("same");
      hEffD->Draw("same");
      // hMCPtGen->Draw("");
      // hMCPtPri->Draw("same");
    }
    cvz->Update();
    //    cvz->WaitPrimitive();
  }
  
  //hManCorr->GetHistoVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim )->Draw();
  //  hManCorr->GetHistoVz(AliAnalysisMultPbTrackHistoManager::kHistoRec )->Draw("");
  // hManCorr->GetHistoVz(AliAnalysisMultPbTrackHistoManager::kHistoGen     )->Draw("same");      


}

TF1 * GetFunctionHistoSum() {

  TF1 * f = new TF1 ("fHistoSum",HistoSum,-1,1,kHistoFitCompoments);
  //  f->SetParNames("Primaries+Fakes", "Sec. Weak decays","Sec. Material");
  f->SetParNames("Primaries+Sec. Material", "Sec. Weak decays");
  f->SetNpx(1000);
  return f;

}

static Double_t HistoSum(const double * x, const double* p){
  // This function uses a global array of histograms, rescaled by some
  // parameters, to define a function
  // The array is called gHistoCompoments
  // The size of the arrays is given by the global constant kHistoFitCompoments
  
  Double_t xx = x[0];
  Double_t value = 0;
  for(Int_t icomp = 0; icomp < kHistoFitCompoments; icomp++){
    //    value += gHistoCompoments[icomp]-Interpolate(xx) * p[icomp];
    Int_t ibin = gHistoCompoments[icomp]->FindBin(xx);
    value += gHistoCompoments[icomp]->GetBinContent(ibin) * p[icomp];
  }
  
  return value;

}

void PrintCanvas(TCanvas* c,const TString formats) {
  // print a canvas in every of the given comma-separated formats
  // ensure the canvas is updated
  if(!doPrint) return;
  c->Update();
  gSystem->ProcessEvents();
  TObjArray * arr = formats.Tokenize(",");
  TIterator * iter = arr->MakeIterator();
  TObjString * element = 0;
  TString name  =c ->GetName();
  name.ReplaceAll(" ","_");
  name.ReplaceAll("+","Plus");
  name.ReplaceAll("-","");
  while ((element = (TObjString*) iter->Next())) {
    c->Print(name+ "."+element->GetString());
  }
}

